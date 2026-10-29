#!/usr/bin/env julia
# Sample diffusion and collect one carrier/sample stream in a single pass.
using HDF5, Random, SHA, TOML, DataFrames
const ROOT=normpath(joinpath(@__DIR__,".."))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis

"""Derive a stable event seed from source identity, event ID and stream label."""
function readout_seed(source::String,id::Int64,base::Int,label::String)::UInt64
    digest=bytes2hex(sha256(join((realpath(source),string(id),string(base),label),"|")))
    return parse(UInt64,digest[1:16];base=16)
end

"""Return the source-hit energy-weighted transverse centre and minimum z."""
function source_reference_mm(hits::EventHits)::NTuple{3,Float32}
    total=sum(hits.energy_MeV)
    isfinite(total) && total>0f0 || throw(ArgumentError("Event energy must be positive"))
    return (sum(hits.x_mm.*hits.energy_MeV)/total,
            sum(hits.y_mm.*hits.energy_MeV)/total,minimum(hits.z_mm))
end

"""Write an Int64 vector with compression when nonempty."""
function write_count_vector(g::HDF5.Group,key::String,values::Vector{Int64})
    if isempty(values)
        g[key]=values
    else
        ds=create_dataset(g,key,Int64,(length(values),);
                          chunk=(min(length(values),4096),),shuffle=true,deflate=3)
        write(ds,values)
    end
end

"""Evaluate one stream's physical and optical settings from repository cards."""
function readout_settings(carrier::String,ion_card,detector,electron_card)
    ion=carrier=="ion"
    pressure=Float32(electron_card["gas"]["pressure_bar"])
    w=Float32(ion_card["transport"]["ionisation_energy_eV"])
    pitch=ion ? Float32(detector["geometry"]["pitch_m"]*1000) :
                Float32(electron_card["elcc"]["pitch_mm"])
    bucket=ion ? Float32(ion_card["pads"]["z_bin_width_mm"]) :
                 Float32(electron_card["readout"]["bucket_width_mm"])
    speed=ion ? Float32(detector["gas"]["bulk_speed_m_s"]) :
                Float32(electron_card["readout"]["electron_drift_speed_m_s"])
    direction=ion ? Int(ion_card["pads"]["drift_direction"]) : -1
    minimum_span=Float32(ion_card["pads"]["minimum_xy_span_mm"])
    acceptance=0f0;pde=0f0;yield=0f0
    if !ion
        el=electron_card["elcc"]
        depth=Float32(el["hole_depth_mm"])
        expected_drop=Float32(el["reduced_field_kV_cm_bar"])*pressure*(depth/10f0)*1000f0
        isapprox(Float32(el["exit_anode_voltage_V"]-el["entry_gate_voltage_V"]),
                 expected_drop;atol=0.01f0) ||
            error("EL electrode voltages disagree with E/P × pressure × depth")
        Float32(el["sipm_diameter_mm"])>=Float32(el["hole_diameter_mm"]) ||
            error("SiPM must cover the EL hole")
        fit=electron_card["el_yield"][carrier]
        yield=el_photon_yield(Float32(fit["gain_photons_per_electron_kV"]),
            Float32(fit["threshold_kV_cm_bar"]),
            Float32(el["reduced_field_kV_cm_bar"]),pressure,depth)
        acceptance=axial_light_acceptance(Float32(el["hole_diameter_mm"])/2f0,depth)
        pde=Float32(electron_card["readout"]["sipm_pde"])
    end
    return (;ion,pressure,w,pitch,bucket,speed,direction,minimum_span,
            acceptance,pde,yield)
end

"""Write one event's complete sparse ion or PE history and exact carrier ledger."""
function write_readout_event(events::HDF5.Group,id::Int64,history,signal,
        carrier::String,ion_cm::Float32,electron_cm::Float32,
        sigma_t::Float32,sigma_l::Float32,location_seed::UInt64,
        diffusion_seed::UInt64,phase::NTuple{2,Float32},settings)
    g=create_group(events,string(id))
    rows=signal.counts
    for field in (:pad_ix,:pad_iy,:slice_index)
        write_count_vector(g,string(field),copy(rows[!,field]))
    end
    write_count_vector(g,"signal_count",copy(rows[!,carrier=="ion" ? :ions : :pe]))
    g["x_edges_mm"]=history.x_edges_mm;g["y_edges_mm"]=history.y_edges_mm
    g["z_edges_mm"]=history.z_edges_mm
    if carrier=="ion"
        g["time_edges_ms"]=history.time_edges_ms
    else
        g["time_edges_us"]=history.time_edges_ms.*1000f0
    end
    g["ion_drift_cm"]=ion_cm;g["electron_drift_cm"]=electron_cm
    g["sigma_t_mm"]=sigma_t;g["sigma_l_mm"]=sigma_l
    g["sampled_carriers"]=history.sampled_ions
    g["input_carriers"]=history.retained_ions
    g["outside_carriers"]=history.outside_ions
    g["total_signal"]=sum(rows[!,carrier=="ion" ? :ions : :pe])
    g["grid_offset_x_mm"]=phase[1];g["grid_offset_y_mm"]=phase[2]
    g["reference_x_mm"]=history.reference_mm[1]
    g["reference_y_mm"]=history.reference_mm[2]
    g["reference_z_mm"]=history.reference_mm[3]
    g["location_seed"]=string(location_seed)
    g["diffusion_seed"]=string(diffusion_seed)
end

"""Sample Gaussian carriers directly into physical grid/time buckets.

Process the first `limit` events of one sample (or the optional explicit
`event_id`) using the same six-stream seed, placement, diffusion and PE models
as the earlier two-stage run. The source hit table is read in Float32 mm/MeV.
Every retained carrier is assigned directly to one pad/hole and 2 mm bucket;
individual carrier coordinates and display histograms are not written. Ion
signals are exact carrier counts. Electron counts per hole/bucket are reduced
to PE with the fitted mean EL yield and binomial optical detection. Each PE
stays in its parent arrival bucket. Output is a versioned compressed HDF5
history with complete axes, exact accounting, saved cards, seeds and source
hash. The file is moved into place only after a completion marker is written.
"""
function sample_collect_stream(sample::String,carrier::String,limit::Int,
                               output::String;event_id::Union{Nothing,Int64}=nothing)
    sample in ("bb0nu","electrons") || error("Unknown sample")
    carrier in ("ion","pure_xe","xe_he") || error("Unknown carrier")
    limit>0 || error("Event limit must be positive")
    ion_card_path=joinpath(ROOT,"config/pixel_charge_1mm.toml")
    detector_path=joinpath(ROOT,"config/detector.toml")
    electron_card_path=joinpath(ROOT,"config/electron_transport.toml")
    ion_card=TOML.parsefile(ion_card_path)
    detector=TOML.parsefile(detector_path)
    electron_card=TOML.parsefile(electron_card_path)
    settings=readout_settings(carrier,ion_card,detector,electron_card)
    source=realpath(joinpath(ROOT,ion_card["samples"][sample]))
    base=Int(ion_card["run"]["seed"])
    maximum_cm=Float32(electron_card["placement"]["l_max_cm"])
    minimum_cm=Float32(electron_card["placement"]["min_ion_drift_cm"])
    maximum_cm==Float32(detector["geometry"]["bulk_reference_length_m"]*100) ||
        error("Placement length disagrees with detector card")
    mkpath(dirname(output));temporary=output*".tmp"
    started=time_ns()
    h5open(source,"r") do input
        ids=Int64.(read(input["events/event_id"]))
        ids=event_id===nothing ? ids[1:min(limit,length(ids))] : Int64[event_id]
        event_id===nothing || event_id in read(input["events/event_id"]) ||
            error("Event ID is absent from source")
        h5open(temporary,"w") do target
            attributes(target)["schema"]="nausicaa.direct_grid_readout.v2"
            target["metadata/sample"]=sample;target["metadata/carrier"]=carrier
            target["metadata/signal_unit"]=settings.ion ? "ions" : "pe"
            target["metadata/time_unit"]=settings.ion ? "ms" : "us"
            target["metadata/source_path"]=source
            target["metadata/source_sha256"]=open(io->bytes2hex(sha256(io)),source)
            target["metadata/ion_card_toml"]=read(ion_card_path,String)
            target["metadata/electron_card_toml"]=read(electron_card_path,String)
            target["metadata/detector_card_toml"]=read(detector_path,String)
            target["metadata/pitch_mm"]=settings.pitch
            target["metadata/bucket_mm"]=settings.bucket
            target["metadata/drift_speed_m_s"]=settings.speed
            target["metadata/drift_direction"]=Int64(settings.direction)
            target["metadata/finite_window_nsigma"]=Float32(ion_card["transport"]["nsigma"])
            target["metadata/el_photons_per_electron"]=settings.yield
            target["metadata/pe_probability_per_photon"]=settings.acceptance*settings.pde
            target["metadata/base_seed"]=Int64(base)
            target["metadata/julia_version"]=string(VERSION)
            target["metadata/implementation_sha256"]=bytes2hex(sha256(read(joinpath(ROOT,"src/event_analysis/grid_readout.jl"))))
            code_paths=sort(vcat(["scripts/sample_collect_stream.jl","Project.toml","Manifest.toml"],
                [relpath(joinpath(dir,name),ROOT) for (dir,_,names) in
                 walkdir(joinpath(ROOT,"src/event_analysis")) for name in names if endswith(name,".jl")]))
            target["metadata/code_paths"]=code_paths
            target["metadata/code_sha256"]=[bytes2hex(sha256(read(joinpath(ROOT,path)))) for path in code_paths]
            events=create_group(target,"events")
            for (index,id) in enumerate(ids)
                hits,_=read_event(input,id)
                counts=ion_counts(hits,settings.w).sampled
                reference=source_reference_mm(hits)
                location_seed=readout_seed(source,id,base,"position")
                ion_cm=sample_drift_cm(Xoshiro(location_seed),minimum_cm,maximum_cm)
                electron_cm=maximum_cm-ion_cm
                drift_cm=settings.ion ? ion_cm : electron_cm
                if settings.ion
                    sigma_t=sigma_l=sigma_t_ion_mm(Float32(detector["gas"]["temperature_K"]),
                        drift_cm,Float32(detector["bias"]["bulk_field_V_m"])/100f0)
                else
                    coeff=electron_card["diffusion"][carrier]
                    scale=sqrt(drift_cm/settings.pressure)
                    sigma_t=Float32(coeff["d_T_mm_sqrtbar_per_sqrtcm"])*scale
                    sigma_l=Float32(coeff["d_L_mm_sqrtbar_per_sqrtcm"])*scale
                end
                diffusion_seed=readout_seed(source,id,base,carrier)
                phase_label=settings.ion ? "ion_grid" : "electron_grid"
                phase=sample_grid_offset_mm(Xoshiro(readout_seed(source,id,base,phase_label)),settings.pitch)
                history=diffuse_to_grid(Xoshiro(diffusion_seed),hits,counts;
                    sigma_t_mm=sigma_t,sigma_l_mm=sigma_l,
                    nsigma=Float32(ion_card["transport"]["nsigma"]),
                    reference_mm=reference,grid_offset_mm=phase,
                    pitch_mm=settings.pitch,bucket_width_mm=settings.bucket,
                    drift_speed_mm_ms=settings.speed,drift_direction=settings.direction,
                    minimum_xy_span_mm=settings.minimum_span)
                signal=if settings.ion
                    history
                else
                    pe_seed=readout_seed(source,id,base,"pe_"*carrier)
                    pe=detect_pe_history(Xoshiro(pe_seed),history;
                        photons_per_electron=settings.yield,
                        acceptance=settings.acceptance,pde=settings.pde)
                    (;counts=pe)
                end
                write_readout_event(events,id,history,signal,carrier,ion_cm,electron_cm,
                    sigma_t,sigma_l,location_seed,diffusion_seed,phase,settings)
                index%100==0 && (println("$sample/$carrier $index/$(length(ids))");flush(stdout))
            end
            target["metadata/event_count"]=Int64(length(ids))
            target["metadata/complete"]=true
        end
    end
    mv(temporary,output;force=true)
    println("completed $sample/$carrier in ",round((time_ns()-started)/1e9;digits=2)," s")
    println("peak_rss_bytes=",Sys.maxrss())
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (4,5) || error("Expected SAMPLE CARRIER LIMIT OUTPUT [EVENT_ID]")
    sample_collect_stream(ARGS[1],ARGS[2],parse(Int,ARGS[3]),abspath(ARGS[4]);
                          event_id=length(ARGS)==5 ? parse(Int64,ARGS[5]) : nothing)
end
