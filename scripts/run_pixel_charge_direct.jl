#!/usr/bin/env julia
# Recompute the earlier ion-only pixel studies with direct carrier assignment.
using HDF5, TOML, SHA, Random, DataFrames
const ROOT=normpath(joinpath(@__DIR__,".."))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis
include(joinpath(@__DIR__,"reduce_pixel_charge.jl"))

"""Convert one pixel-charge card to a versioned direct-readout output path."""
direct_output(card)=joinpath(ROOT,card["run"]["output_directory"]*"_direct")

"""Run one ion sample using exact sampled-carrier grid assignments.

The first RNG replay reduces the retained positions to their ion-weighted x/y
centre and minimum z, preserving the earlier viewer's grid-reference rule.
The second replay assigns each carrier to the physical pad and time bucket.
No per-carrier positions or diffusion histogram are stored. Both passes use
the same deterministic diffusion seed, finite window and Float32 draw order.
The output keeps the earlier pixel-charge event/cell interface with a new HDF5
schema and source/configuration/code provenance, allowing the numerical
histogram reducer and Python presentation to operate on corrected histories.
"""
function run_sample(card_path::String,sample::String)
    sample in ("bb0nu","electrons") || error("Unknown sample")
    card=TOML.parsefile(card_path);tc=card["transport"];pc=card["pads"]
    detector_path=joinpath(ROOT,tc["detector_card"])
    detector=TOML.parsefile(detector_path)
    input=realpath(joinpath(ROOT,card["samples"][sample]))
    output=direct_output(card);mkpath(output)
    target=joinpath(output,sample*".h5");temporary=target*".tmp"
    pitch=Float32(get(pc,"pitch_mm",detector["geometry"]["pitch_m"]*1000))
    bucket=Float32(pc["z_bin_width_mm"])
    speed=Float32(detector["gas"]["bulk_speed_m_s"])
    w=Float32(tc["ionisation_energy_eV"])
    nsigma=Float32(tc["nsigma"])
    n_min=Int64(card["statistics"]["n_min"])
    seed_base=Int(card["run"]["seed"])
    events=DataFrame(event_id=Int64[],drift_cm=Float32[],sigma_mm=Float32[],seed=String[],
        produced_ions=Int64[],retained_ions=Int64[],outside_ions=Int64[],collected_ions=Int64[],
        cells_with_charge=Int64[],pulse_count=Int64[],has_pulse=Bool[],
        maximum_ions=Int64[],minimum_ions=Int64[],
        grid_offset_x_mm=Float32[],grid_offset_y_mm=Float32[])
    cells=DataFrame(event_id=Int64[],pixel_x=Int64[],pixel_y=Int64[],
                    time_slice=Int64[],n_ions=Int64[])
    started=time_ns()
    h5open(input,"r") do source
        ids=Int64.(read(source["events/event_id"]))
        ids=ids[1:min(length(ids),Int(card["run"]["n_events"]))]
        h5open(temporary,"w") do file
            attributes(file)["schema"]="nausicaa.pixel_charge.direct.v2"
            for (i,id) in enumerate(ids)
                hits,_=read_event(source,id)
                token=bytes2hex(sha256(join((input,string(id),string(seed_base),"0"),"|")))
                seed=parse(UInt64,token[1:16];base=16)
                drift=sample_drift_cm(Xoshiro(seed),Float32(tc["min_cm"]),Float32(tc["max_cm"]))
                sigma=sigma_t_ion_mm(Float32(detector["gas"]["temperature_K"]),drift,
                    Float32(detector["bias"]["bulk_field_V_m"])/100f0)
                counts=ion_counts(hits,w).sampled
                diffusion_seed=seed⊻0x9e3779b97f4a7c15
                reference=diffused_carrier_reference(Xoshiro(diffusion_seed),hits,counts;
                    sigma_t_mm=sigma,sigma_l_mm=sigma,nsigma).reference_mm
                phase=sample_grid_offset_mm(Xoshiro(seed⊻0xd1b54a32d192ed03),pitch)
                history=diffuse_to_grid(Xoshiro(diffusion_seed),hits,counts;
                    sigma_t_mm=sigma,sigma_l_mm=sigma,nsigma,reference_mm=reference,
                    grid_offset_mm=phase,pitch_mm=pitch,bucket_width_mm=bucket,
                    drift_speed_mm_ms=speed,drift_direction=Int(pc["drift_direction"]),
                    minimum_xy_span_mm=Float32(pc["minimum_xy_span_mm"]))
                aligned=align_pad_history(history,reference)
                stats=pixel_frame_statistics(history;n_min)
                history.sampled_ions==history.retained_ions+history.outside_ions ||
                    error("Direct diffusion ledger mismatch at event $id")
                for row in eachrow(aligned)
                    push!(cells,(id,row.pad_dx,row.pad_dy,row.arrival_bin,row.ions))
                end
                push!(events,(id,drift,sigma,string(seed),history.sampled_ions,
                    history.retained_ions,history.outside_ions,history.total_ions,
                    nrow(aligned),length(stats.counts),!isempty(stats.counts),
                    something(stats.maximum_ions,0),something(stats.minimum_ions,0),phase...))
                base="event_axes/$(id)/"
                for (key,value) in (("x_edges_mm",history.x_edges_mm),
                    ("y_edges_mm",history.y_edges_mm),("z_edges_mm",history.z_edges_mm),
                    ("time_edges_ms",history.time_edges_ms),
                    ("reference_mm",collect(reference)),("grid_offset_mm",collect(phase)))
                    file[base*key]=value
                end
                i%100==0 && (println("$sample: $i/$(length(ids)) events");flush(stdout))
            end
            for (name,table) in (("events",events),("cells",cells)),
                field in propertynames(table)
                file[name*"/"*string(field)]=table[!,field]
            end
            attributes(file["cells"])["index_convention"]="pixel_x/y relative to reference pixel; time_slice zero-based relative arrival frame"
            file["summary/processed_events"]=length(ids)
            file["summary/produced_ions"]=sum(events.produced_ions)
            file["summary/outside_ions"]=sum(events.outside_ions)
            file["summary/collected_ions"]=sum(events.collected_ions)
            file["settings/n_min"]=n_min
            file["settings/pitch_mm"]=pitch
            file["settings/frame_ms"]=bucket/speed
            file["settings/drift_direction"]=pc["drift_direction"]
            file["provenance/card_toml"]=read(card_path,String)
            file["provenance/detector_toml"]=read(detector_path,String)
            file["provenance/source_path"]=input
            file["provenance/source_sha256"]=open(io->bytes2hex(sha256(io)),input)
            file["provenance/julia_version"]=string(VERSION)
            paths=["scripts/run_pixel_charge_direct.jl","src/event_analysis/grid_readout.jl",
                   "src/event_analysis/production.jl","src/event_analysis/diffusion.jl"]
            file["provenance/code_paths"]=paths
            file["provenance/code_sha256"]=[bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in paths]
            file["timing/elapsed_seconds"]=(time_ns()-started)/1e9
            file["summary/complete"]=true
        end
    end
    mv(temporary,target;force=true)
    println("$sample direct pixel histories: ",round((time_ns()-started)/1e9;digits=2),
            " s, ",nrow(cells)," occupied cells, ",filesize(target)," bytes")
    return output
end

"""Run both samples, then reduce saved counts and write presentation inputs."""
function run_pixel_charge_direct(card_path::String)
    for sample in ("bb0nu","electrons")
        run_sample(card_path,sample)
    end
    output=direct_output(TOML.parsefile(card_path))
    reduce_pixel_charge(output)
    summarize_pixel_charge(output)
    return output
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("Expected PIXEL_CHARGE_CARD")
    run_pixel_charge_direct(abspath(ARGS[1]))
end
