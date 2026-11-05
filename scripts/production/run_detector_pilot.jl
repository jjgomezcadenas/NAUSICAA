#!/usr/bin/env julia
# Run shared-pair ion/photon detector pilot from five TOML cards.
using HDF5, Random, SHA, TOML, Dates, Printf
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis
include(joinpath(@__DIR__,"specs.jl"))
using .ProductionSpecs

"""Return the source-energy-weighted transverse reference in original mm."""
function transverse_reference(hits::EventHits)::NTuple{2,Float32}
    weights=Float64.(hits.energy_MeV)
    total=sum(weights)
    total>0 || throw(ArgumentError("Event deposited energy must be positive"))
    return (Float32(sum(Float64.(hits.x_mm).*weights)/total),
            Float32(sum(Float64.(hits.y_mm).*weights)/total))
end

"""Populate top-level source, card, run and completion metadata in one file."""
function write_header(f::HDF5.File,schema::String,sample::String,
        source_path::String,source_sha::String,card_hash::String,
        code_hash::String,seed::Int,branch::String,
        realisation_card::String,branch_card::String)
    attributes(f)["schema"]=schema
    m=create_group(f,"metadata")
    m["complete"]=false
    m["sample"]=sample
    m["source_path"]=source_path
    m["source_sha256"]=source_sha
    m["card_hash"]=card_hash
    m["code_hash"]=code_hash
    m["branch"]=branch
    m["realisation_card_toml"]=realisation_card
    m["branch_card_toml"]=branch_card
    m["signal_unit"]=branch=="ions" ? "ions" : branch=="realisation" ? "pairs" : "incident_photons"
    m["time_unit"]=branch=="ions" ? "ms" : branch=="realisation" ? "" : "us"
    m["seed"]=Int64(seed)
    m["created_utc"]=string(now(UTC))
    create_group(f,"events")
    return nothing
end

"""Write all four event files to temporary paths, then publish complete files."""
function run_sample(spec,sample::String,code_hash::String)
    r=spec.realisation.data;i=spec.ions.data;e=spec.electrons.data
    source_path=spec.sample_paths[sample]
    source_sha=content_hash(source_path;file=true)
    identities=realisation_identities(ROOT,sample,source_sha,r)
    stream_identity=identities.stream_identity;realisation_id=identities.realisation_id
    ion_id=content_hash(join((realisation_id,canonical(i),code_hash),"|"))
    electron_id=content_hash(join((realisation_id,canonical(e),code_hash),"|"))
    base=joinpath(spec.results_dir,sample,realisation_id)
    paths=(realisation=joinpath(base,"realisation.h5"),
           ions=joinpath(base,"detector","ions",ion_id,"ions.h5"),
           pure_xe=joinpath(base,"detector","photons",electron_id,"photons_pure_xe.h5"),
           xe_he15=joinpath(base,"detector","photons",electron_id,"photons_xe_he15.h5"))
    for path in values(paths)
        mkpath(dirname(path))
    end
    temporary=map(path->path*".tmp",paths)
    for path in values(temporary)
        isfile(path) && rm(path)
    end
    seed=Int(r["seed"])
    started=time_ns()
    totals=Dict("events"=>Int64(0),"pairs"=>Int64(0),"ions"=>Int64(0),
                "survivors"=>Int64(0),"pure_xe_emitted"=>Int64(0),
                "pure_xe_incident"=>Int64(0),"xe_he15_emitted"=>Int64(0),
                "xe_he15_incident"=>Int64(0))
    h5open(source_path,"r") do source
        ids=read(source["events/event_id"])
        ids=ids[1:min(length(ids),Int(r["max_events_per_sample"]))]
        h5open(temporary.realisation,"w") do fr
            h5open(temporary.ions,"w") do fi
                h5open(temporary.pure_xe,"w") do fx
                    h5open(temporary.xe_he15,"w") do fh
                        write_header(fr,"nausicaa.event_realisation.pilot.v1",sample,
                            source_path,source_sha,realisation_id,code_hash,seed,
                            "realisation",read(spec.realisation.path,String),"")
                        write_header(fi,"nausicaa.ion_pad_file.pilot.v1",sample,
                            source_path,source_sha,ion_id,code_hash,seed,
                            "ions",read(spec.realisation.path,String),read(spec.ions.path,String))
                        write_header(fx,"nausicaa.photon_sipm_file.pilot.v1",sample,
                            source_path,source_sha,electron_id,code_hash,seed,
                            "pure_xe",read(spec.realisation.path,String),read(spec.electrons.path,String))
                        write_header(fh,"nausicaa.photon_sipm_file.pilot.v1",sample,
                            source_path,source_sha,electron_id,code_hash,seed,
                            "xe_he15",read(spec.realisation.path,String),read(spec.electrons.path,String))
                        for file in (fr,fi,fx,fh)
                            file["metadata/stream_identity"]=stream_identity
                        end
                        for raw_id in ids
                            id=Int64(raw_id)
                            hits,_=read_event(source,id)
                            length_cm=Float32(r["chamber_length_cm"])
                            ion_cm=sample_drift_cm(Xoshiro(event_seed(stream_identity,id,"placement")),
                                Float32(r["min_ion_drift_cm"]),Float32(get(r,"max_ion_drift_cm",length_cm)))
                            placement=EventPlacement(id,ion_cm,length_cm-ion_cm)
                            pairs=sample_fano_pairs(Xoshiro(event_seed(stream_identity,id,"pairs")),
                                hits,Float32(r["w_eV"]),Float32(r["fano_factor"]))
                            n_pairs=foldl(Base.Checked.checked_add,pairs;init=Int64(0))
                            checksum=pair_checksum(pairs)
                            rx,ry=transverse_reference(hits)
                            ion_phase=sample_grid_offset_mm(Xoshiro(event_seed(stream_identity,id,"ion_phase")),
                                Float32(i["pad_pitch_mm"]))
                            electron_phase=sample_grid_offset_mm(Xoshiro(event_seed(stream_identity,id,"electron_phase")),
                                Float32(e["hole_pitch_mm"]))
                            ion_bucket_mm=Float32(i["drift_speed_m_s"]*i["bucket_width_ms"])
                            electron_bucket_mm=Float32(e["drift_speed_mm_us"]*e["bucket_width_us"])
                            ion_grid=GridSpec(Float32(i["pad_pitch_mm"]),ion_bucket_mm,
                                (rx,ry,minimum(hits.z_mm)),ion_phase,Int8(1))
                            electron_grid=GridSpec(Float32(e["hole_pitch_mm"]),electron_bucket_mm,
                                (rx,ry,maximum(hits.z_mm)),electron_phase,Int8(-1))
                            sigma_ion=sigma_t_ion_mm(Float32(i["temperature_K"]),ion_cm,
                                Float32(i["drift_field_V_cm"]))
                            ion=collect_ion_pad_history(Xoshiro(event_seed(stream_identity,id,"ion_diffusion")),
                                hits,pairs,placement,ion_grid;
                                sigma_t_mm=sigma_ion,sigma_l_mm=sigma_ion)
                            drift_ms=Float32(10*placement.electron_drift_cm/
                                (1000*e["drift_speed_mm_us"]))
                            survivors=sample_attached_electrons(
                                Xoshiro(event_seed(stream_identity,id,"attachment")),pairs,
                                drift_ms,Float32(e["lifetime_ms"]))
                            n_survivors=foldl(Base.Checked.checked_add,survivors;init=Int64(0))
                            gas_results=Dict{String,PhotonSiPMHistory}()
                            for gas in ("pure_xe","xe_he15")
                                g=e[gas]
                                scale=sqrt(placement.electron_drift_cm/Float32(e["pressure_bar"]))
                                sigma_t=Float32(g["d_T_mm_sqrtbar_per_sqrtcm"])*scale
                                sigma_l=Float32(g["d_L_mm_sqrtbar_per_sqrtcm"])*scale
                                yield=el_photon_yield(Float32(g["yield_gain_photons_per_electron_kV"]),
                                    Float32(g["yield_threshold_kV_cm_bar"]),
                                    Float32(e["reduced_field_kV_cm_bar"]),Float32(e["pressure_bar"]),
                                    Float32(e["hole_depth_mm"]))
                                acceptance=axial_light_acceptance(Float32(e["hole_diameter_mm"])/2f0,
                                    Float32(e["hole_depth_mm"]))
                                gas_results[gas]=collect_photon_sipm_history(
                                    Xoshiro(event_seed(stream_identity,id,"$(gas)_diffusion_el")),
                                    hits,survivors,placement,electron_grid;
                                    sigma_t_mm=sigma_t,sigma_l_mm=sigma_l,
                                    yield_per_e=yield,J=Float32(e["j_variance_ratio"]),
                                    acceptance=acceptance,gas_model=Symbol(gas),
                                    initial_pairs=n_pairs)
                            end
                            gr=create_group(fr["events"],string(id))
                            gr["event_id"]=id;gr["ion_drift_cm"]=ion_cm
                            gr["electron_drift_cm"]=placement.electron_drift_cm
                            gr["initial_pairs"]=n_pairs;gr["pair_checksum"]=checksum
                            for (file,history) in ((fi,ion),(fx,gas_results["pure_xe"]),
                                                   (fh,gas_results["xe_he15"]))
                                group=create_group(file["events"],string(id))
                                if history isa IonPadHistory
                                    write_ion_history(group,history)
                                else
                                    write_photon_history(group,history)
                                end
                                group["pair_checksum"]=checksum
                            end
                            totals["events"]+=1
                            totals["pairs"]+=n_pairs
                            totals["ions"]+=sum(ion.cells.count)
                            totals["survivors"]+=n_survivors
                            for gas in ("pure_xe","xe_he15")
                                h=gas_results[gas]
                                totals[gas*"_emitted"]+=h.emitted_photons
                                totals[gas*"_incident"]+=sum(h.cells.count)
                            end
                            xe_incident=sum(gas_results["pure_xe"].cells.count)
                            he_incident=sum(gas_results["xe_he15"].cells.count)
                            println("$sample event $id: pairs=$n_pairs, ions=$(sum(ion.cells.count)), "*
                                "incident Xe=$xe_incident, XeHe=$he_incident")
                            flush(stdout)
                        end
                        for file in (fr,fi,fx,fh)
                            write(file["metadata/complete"],true)
                            file["metadata/event_count"]=totals["events"]
                        end
                    end
                end
            end
        end
    end
    for key in keys(paths)
        mv(getproperty(temporary,key),getproperty(paths,key);force=true)
    end
    elapsed=(time_ns()-started)/1e9
    n_events=totals["events"]
    println("$sample: $n_events events in $(round(elapsed;digits=2)) s")
    return (;sample,realisation_id,stream_identity,ion_id,electron_id,
            paths,totals,elapsed,source_sha)
end

"""Run both sample branches and save numerical pilot summary for presentation."""
function main(card::String)
    spec=load_specs(ROOT,card)
    code_hash=content_hash(join((content_hash(joinpath(ROOT,p);file=true) for p in
        ("src/event_analysis/EventAnalysis.jl","src/event_analysis/types.jl",
         "src/event_analysis/io.jl","src/event_analysis/production.jl",
         "src/event_analysis/production_readout.jl","src/event_analysis/grid_readout.jl",
         "src/ion_signal/event_io.jl",
         "scripts/production/specs.jl","scripts/production/run_detector_pilot.jl")),"|"))
    run_id=content_hash(canonical(spec.run_spec)*canonical(spec.configuration)*
        canonical(spec.realisation.data)*canonical(spec.ions.data)*
        canonical(spec.electrons.data)*code_hash)[1:16]
    results=[run_sample(spec,name,code_hash) for name in
             get(spec.run_spec,"selected_samples",["bb0nu","single_electrons"])]
    summary_path=joinpath(spec.results_dir,"pilot_summary_$(run_id).h5")
    mkpath(dirname(summary_path))
    h5open(summary_path*".tmp","w") do f
        attributes(f)["schema"]="nausicaa.detector_pilot_summary.v1"
        f["metadata/complete"]=false
        f["metadata/run_id"]=run_id
        f["metadata/code_hash"]=code_hash
        f["metadata/run_card"]=read(spec.run_path,String)
        f["metadata/realisation_card"]=read(spec.realisation.path,String)
        f["metadata/ion_card"]=read(spec.ions.path,String)
        f["metadata/electron_card"]=read(spec.electrons.path,String)
        f["metadata/peak_rss_bytes"]=Int64(Sys.maxrss())
        r=spec.realisation.data;i=spec.ions.data;e=spec.electrons.data
        d=create_group(f,"derived")
        chamber_cm=Float32(r["chamber_length_cm"])
        pressure=Float32(e["pressure_bar"])
        reduced_field=Float32(e["reduced_field_kV_cm_bar"])
        hole_depth=Float32(e["hole_depth_mm"])
        d["ion_sigma_max_mm"]=sigma_t_ion_mm(Float32(i["temperature_K"]),
            chamber_cm,Float32(i["drift_field_V_cm"]))
        d["ion_bucket_mm"]=Float32(i["drift_speed_m_s"]*i["bucket_width_ms"])
        d["electron_bucket_mm"]=Float32(e["drift_speed_mm_us"]*e["bucket_width_us"])
        d["electron_max_drift_ms"]=Float32(10*chamber_cm/(1000*e["drift_speed_mm_us"]))
        d["electron_min_survival"]=exp(-read(d["electron_max_drift_ms"])/Float32(e["lifetime_ms"]))
        d["optical_acceptance"]=axial_light_acceptance(Float32(e["hole_diameter_mm"])/2f0,hole_depth)
        gate_V=Float32(i["drift_field_V_cm"])*chamber_cm
        d["gate_voltage_V"]=gate_V
        d["anode_voltage_V"]=gate_V+reduced_field*pressure*(hole_depth/10f0)*1000f0
        for gas in ("pure_xe","xe_he15")
            model=e[gas]
            gd=create_group(d,gas)
            gd["sigma_t_max_mm"]=Float32(model["d_T_mm_sqrtbar_per_sqrtcm"])*sqrt(chamber_cm/pressure)
            gd["sigma_l_max_mm"]=Float32(model["d_L_mm_sqrtbar_per_sqrtcm"])*sqrt(chamber_cm/pressure)
            gd["yield_photons_per_e"]=el_photon_yield(
                Float32(model["yield_gain_photons_per_electron_kV"]),
                Float32(model["yield_threshold_kV_cm_bar"]),reduced_field,pressure,hole_depth)
        end
        for result in results
            g=create_group(f,result.sample)
            for (key,value) in result.totals
                g[key]=value
            end
            g["ion_charge_conserved"]=result.totals["ions"]==result.totals["pairs"]
            g["photon_ledgers_valid"]=all(gas->result.totals[gas*"_incident"]<=
                result.totals[gas*"_emitted"],("pure_xe","xe_he15"))
            g["elapsed_s"]=result.elapsed
            g["realisation_id"]=result.realisation_id
            g["stream_identity"]=result.stream_identity
            g["ion_spec_id"]=result.ion_id
            g["electron_spec_id"]=result.electron_id
            g["source_sha256"]=result.source_sha
            for key in keys(result.paths)
                g[string(key)*"_path"]=getproperty(result.paths,key)
            end
        end
        write(f["metadata/complete"],true)
    end
    mv(summary_path*".tmp",summary_path;force=true)
    report_dir=joinpath(spec.report_dir,run_id)
    mkpath(report_dir)
    python=get(ENV,"PYTHON","python3.13")
    run(`$python $(joinpath(@__DIR__,"present_detector_pilot.py")) $summary_path $report_dir`)
    println("Summary: $summary_path")
    println("Report: $(joinpath(report_dir,"formulas.pdf"))")
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)<=1 || error("Usage: julia --project=. scripts/production/run_detector_pilot.jl [run.toml]")
    main(isempty(ARGS) ? "config/specs/runs/production_pilot.toml" : ARGS[1])
end
