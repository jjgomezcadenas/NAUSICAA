#!/usr/bin/env julia
# Map one saved diffusion stream to ion pad or electron SiPM time histories.
using HDF5, Random, SHA, TOML, DataFrames
const ROOT = normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis

"""Deterministic event seed, independent of process scheduling and event order."""
function grid_seed(source::String,id::Int64,base::Int,label::String)::UInt64
    digest=bytes2hex(sha256(join((realpath(source),string(id),string(base),label),"|")))
    return parse(UInt64,digest[1:16];base=16)
end

"""Read one lossless occupied cloud, retaining Float32 mm and Int64 carriers."""
function read_cloud_event(g::HDF5.Group,id::Int64,w::Float32)::DataFrame
    counts=read(g["ions"])
    x,y,z=(read(g[axis]) for axis in ("x","y","z"))
    sum(counts)==read(g["retained_ions"]) || error("Cloud carrier ledger mismatch: $id")
    return DataFrame(event_id=fill(id,length(counts)),x=x,y=y,z=z,
                     energy=Float32.(counts).*w./1f6,ions=counts)
end

"""Use the undiffused event to place both grids in one common source frame."""
function source_reference(hits::EventHits)::NTuple{3,Float32}
    total=sum(hits.energy_MeV)
    isfinite(total) && total>0f0 || error("Event energy must be positive")
    return (sum(hits.x_mm.*hits.energy_MeV)/total,
            sum(hits.y_mm.*hits.energy_MeV)/total,minimum(hits.z_mm))
end

"""Write a compressed nonempty Int64 vector as one HDF5 dataset."""
function write_vector(g::HDF5.Group,key::String,values::Vector{Int64})
    ds=create_dataset(g,key,Int64,(length(values),);
                      chunk=(max(1,min(length(values),4096)),),shuffle=true,deflate=3)
    write(ds,values)
end

"""Collect one carrier stream with independent physical-grid phases.

For ions use the detector's 2 mm pad pitch and 2 mm (20 ms) bucket. For
electrons assign every retained cloud cell to a 10 mm-pitch EL hole and 2 mm
(2 microsecond) bucket, with reversed arrival order. The hole radius enters
the axial optical acceptance; geometrical focusing sends the full retained
electron cloud to holes. For each occupied hole/bucket sample PE from the
matched pure-Xe or Xe--15% He fitted yield and combined acceptance × PDE.
The fitted yield is fixed at its mean; photon and electron intermediate counts
are not written. Output consists of sparse pad/time signal counts and complete
axes, phase, drift and provenance. Each process owns one output HDF5 file.
"""
function collect_stream(sample::String,carrier::String,limit::Int,input::String,output::String)
    sample in ("bb0nu","electrons") || error("Unknown sample")
    carrier in ("ion","pure_xe","xe_he") || error("Unknown carrier")
    limit>0 || error("Event limit must be positive")
    card=TOML.parsefile(joinpath(ROOT,"config/pixel_charge_1mm.toml"))
    detector=TOML.parsefile(joinpath(ROOT,"config/detector.toml"))
    ec=TOML.parsefile(joinpath(ROOT,"config/electron_transport.toml"))
    w=Float32(card["transport"]["ionisation_energy_eV"])
    base=Int(card["run"]["seed"])
    ion=carrier=="ion"
    pitch=ion ? Float32(detector["geometry"]["pitch_m"]*1000) : Float32(ec["elcc"]["pitch_mm"])
    bucket=ion ? Float32(card["pads"]["z_bin_width_mm"]) :
          Float32(ec["readout"]["photon_bucket_us"]*ec["readout"]["electron_drift_speed_m_s"]/1000)
    speed=ion ? Float32(detector["gas"]["bulk_speed_m_s"]) : Float32(ec["readout"]["electron_drift_speed_m_s"])
    direction=ion ? Int(card["pads"]["drift_direction"]) : -1
    acceptance=0f0;pde=0f0;probability=0f0;yield=0f0
    if !ion
        el=ec["elcc"]
        depth=Float32(el["hole_depth_mm"]);pressure=Float32(ec["gas"]["pressure_bar"])
        expected_drop=Float32(el["reduced_field_kV_cm_bar"])*pressure*(depth/10f0)*1000f0
        isapprox(Float32(el["exit_anode_voltage_V"]-el["entry_gate_voltage_V"]),expected_drop;atol=0.01f0) ||
            error("EL electrode voltages disagree with E/P × pressure × depth")
        Float32(el["sipm_diameter_mm"])>=Float32(el["hole_diameter_mm"]) || error("SiPM must cover the hole")
        fit=ec["el_yield"][carrier]
        yield=el_photon_yield(Float32(fit["gain_photons_per_electron_kV"]),
            Float32(fit["threshold_kV_cm_bar"]),Float32(el["reduced_field_kV_cm_bar"]),pressure,depth)
        acceptance=axial_light_acceptance(Float32(el["hole_diameter_mm"])/2f0,depth)
        pde=Float32(ec["readout"]["sipm_pde"])
        probability=acceptance*pde
    end
    mkpath(dirname(output))
    started=time_ns();temporary=output*".tmp"
    h5open(input,"r") do cloud
        read(cloud["metadata/complete"])==true || error("Input cloud is incomplete")
        read(cloud["metadata/sample"])==sample && read(cloud["metadata/carrier"])==carrier ||
            error("Input cloud identity mismatch")
        source=read(cloud["metadata/source_path"])
        ids=h5open(f->read(f["events/event_id"]),source,"r")
        ids=ids[1:min(limit,length(ids))]
        h5open(source,"r") do mc
            h5open(temporary,"w") do target
                attributes(target)["schema"]="nausicaa.grid_readout.v1"
                target["metadata/sample"]=sample;target["metadata/carrier"]=carrier
                target["metadata/signal_unit"]=ion ? "ions" : "pe"
                target["metadata/time_unit"]=ion ? "ms" : "us"
                target["metadata/source_cloud_path"]=realpath(input)
                target["metadata/source_cloud_sha256"]=open(io->bytes2hex(sha256(io)),input)
                target["metadata/source_hit_path"]=source
                target["metadata/pitch_mm"]=pitch;target["metadata/bucket_mm"]=bucket
                target["metadata/drift_speed_m_s"]=speed
                target["metadata/drift_direction"]=Int64(direction)
                target["metadata/el_photons_per_electron"]=yield
                target["metadata/pe_probability_per_photon"]=probability
                target["metadata/electron_card_toml"]=read(joinpath(ROOT,"config/electron_transport.toml"),String)
                target["metadata/ion_card_toml"]=read(joinpath(ROOT,"config/pixel_charge_1mm.toml"),String)
                events=create_group(target,"events")
                for (index,raw_id) in enumerate(ids)
                    id=Int64(raw_id)
                    key="events/"*string(id)
                    haskey(cloud,key) || error("Missing source cloud for event $id")
                    gcloud=cloud[key]
                    d=read_cloud_event(gcloud,id,w)
                    hits,_=read_event(mc,id)
                    reference=source_reference(hits)
                    phase_label=ion ? "ion_grid" : "electron_grid"
                    phase=sample_grid_offset_mm(Xoshiro(grid_seed(source,id,base,phase_label)),pitch)
                    result=collect_grid_history(d;reference_mm=reference,pitch_mm=pitch,
                        z_bin_width_mm=bucket,grid_offset_mm=phase,
                        drift_speed_mm_ms=speed,drift_direction=direction,
                        minimum_xy_span_mm=160f0)
                    raw=result.history.counts
                    counts=if ion
                        copy(raw.ions)
                    else
                        rng=Xoshiro(grid_seed(source,id,base,"pe_"*carrier))
                        sampled=detect_pe_history(rng,result.history;
                            photons_per_electron=yield,
                            acceptance=acceptance,pde=pde)
                        # Build the sparse result below from sampled row keys.
                        raw=sampled
                        copy(sampled.pe)
                    end
                    g=create_group(events,string(id))
                    for field in ("pad_ix","pad_iy","slice_index")
                        write_vector(g,field,copy(raw[!,Symbol(field)]))
                    end
                    write_vector(g,"signal_count",counts)
                    g["x_edges_mm"]=result.axes.x_edges_mm
                    g["y_edges_mm"]=result.axes.y_edges_mm
                    if ion
                        g["time_edges_ms"]=result.history.time_edges_ms
                    else
                        g["time_edges_us"]=result.history.time_edges_ms.*1000f0
                    end
                    g["ion_drift_cm"]=read(gcloud["ion_drift_cm"])
                    g["electron_drift_cm"]=read(gcloud["electron_drift_cm"])
                    g["input_carriers"]=read(gcloud["retained_ions"])
                    g["total_signal"]=sum(counts)
                    g["grid_offset_x_mm"]=phase[1];g["grid_offset_y_mm"]=phase[2]
                    g["reference_x_mm"]=reference[1];g["reference_y_mm"]=reference[2]
                    g["reference_z_mm"]=reference[3]
                    index%100==0 && (println("$sample/$carrier $index/$(length(ids))");flush(stdout))
                end
                target["metadata/event_count"]=Int64(length(ids))
                target["metadata/complete"]=true
            end
        end
    end
    mv(temporary,output;force=true)
    println("completed $sample/$carrier in ",round((time_ns()-started)/1e9;digits=2)," s")
    println("peak_rss_bytes=",Sys.maxrss())
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==5 || error("Expected SAMPLE CARRIER LIMIT INPUT_CLOUD OUTPUT")
    collect_stream(ARGS[1],ARGS[2],parse(Int,ARGS[3]),abspath(ARGS[4]),abspath(ARGS[5]))
end
