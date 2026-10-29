# Deterministic integration checkpoint: hits -> diffusion -> selection -> pads.
using HDF5, DataFrames, TOML, SHA, Dates, Random, Statistics
include(joinpath(@__DIR__,"../../src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis
const ROOT = normpath(joinpath(@__DIR__,"../.."))
const CARD = isempty(ARGS) ? joinpath(ROOT,"config/event_analysis.toml") : abspath(ARGS[1])
const OUT = joinpath(ROOT,"results/event_analysis/step3_validation.h5")

"""Save typed table columns and their order in the validation HDF5."""
function save_table(parent,name,table)
    g = create_group(parent,name)
    attributes(g)["columns"] = names(table)
    for column in names(table);g[column] = table[!,column];end
    return g
end

"""
Derive a UInt64 seed from length-prefixed (version, run seed, sample ID, event ID,
stage) strings using the first eight SHA-256 bytes in little-endian order. This
keeps drift, diffusion and grid streams independent of processing order and
selection outcome. The scheme and each realised seed are saved in the output.
"""
function stage_seed(seed::UInt64,sample::String,event_id::Int64,stage::String)::UInt64
    fields = ("nausicaa.event.rng.v1",string(seed),sample,string(event_id),stage)
    encoded = join(string(ncodeunits(s))*":"*s for s in fields)
    digest = sha256(encoded)
    return foldl((v,i)->v|(UInt64(digest[i]) << (8(i-1))),1:8;init=UInt64(0))
end

"""Run the reviewed computations for one event with explicit, saved RNG streams."""
function process_event(hits,config,detector,sample)
    run_seed = UInt64(config["random"]["seed"])
    seeds = (;drift=stage_seed(run_seed,sample,hits.event_id,"drift"),
              diffusion=stage_seed(run_seed,sample,hits.event_id,"diffusion"),
              grid=stage_seed(run_seed,sample,hits.event_id,"grid"))
    w = Float32(config["ionization"]["ionisation_energy_eV"])
    counts = ion_counts(hits,w)
    L = sample_drift_cm(Xoshiro(seeds.drift),Float32(config["drift"]["min_cm"]),Float32(config["drift"]["max_cm"]))
    sigma = sigma_t_ion_mm(Float32(detector["gas"]["temperature_K"]),L,Float32(detector["bias"]["bulk_field_V_m"])/100f0)
    model = config["diffusion"]["longitudinal"]
    model in ("thermal","zero") || error("Unknown longitudinal model: $model")
    sigma_z = model == "thermal" ? sigma : 0f0
    b,_,distance = get_voxel_size_and_distance(sigma;
        bin_sigma_scale=Float32(config["graph"]["bin_sigma_scale"]),
        edge_sigma_scale=Float32(config["graph"]["edge_sigma_scale"]))
    nbins = Int(config["diffusion"]["nbins"])
    nsigma = Float32(config["diffusion"]["nsigma"])
    threshold = Int64(config["graph"]["min_ions_per_bin"])
    params = DiffusionParams(L,sigma,sigma_z,b,distance,Float32(threshold)*w/1000f0,nbins,nsigma)
    diffused = diffuse_xyz_image_mc(Xoshiro(seeds.diffusion),hits,counts.sampled;
        sigma_t_mm=sigma,sigma_l_mm=sigma_z,nbins,nsigma,ionisation_energy_eV=w)
    selection = select_diffused_event(diffused.bins;bin_width_mm=b,max_distance_mm=distance,
        min_ions_per_bin=threshold,method=config["graph"]["method"],k=Int(config["graph"]["k"]),diffusion=params)
    pitch = Float32(detector["geometry"]["pitch_m"])*1000f0
    offset = config["pads"]["random_grid_offset"] ? sample_grid_offset_mm(Xoshiro(seeds.grid),pitch) : (0f0,0f0)
    axes,history,backbone = nothing,nothing,nothing
    if selection.accepted
        track = only(selection.tracks)
        backbone = walk_track_from_extremes(track)
        axes = pad_axes(diffused.bins,track;pitch_mm=pitch,
            z_bin_width_mm=Float32(config["pads"]["z_bin_width_mm"]),grid_offset_mm=offset,
            minimum_xy_span_mm=Float32(config["pads"]["minimum_xy_span_mm"]))
        # 1 m/s = 1 mm/ms; the numerical conversion factor is one.
        speed = Float32(detector["gas"]["bulk_speed_m_s"])
        history = compute_occupancy_analysis(diffused.bins;x_edges_mm=axes.x_edges_mm,
            y_edges_mm=axes.y_edges_mm,z_edges_mm=axes.z_edges_mm,drift_speed_mm_ms=speed,
            drift_direction=Int(config["pads"]["drift_direction"]))
        @assert sum(history.counts.ions) == diffused.retained_ions
    end
    @assert diffused.sampled_ions == diffused.retained_ions+diffused.outside_ions
    @assert diffused.retained_ions == selection.graph_retained_ions+selection.graph_below_threshold_ions
    return (;counts,seeds,L,sigma,sigma_z,b,distance,offset,diffused,selection,axes,history,backbone)
end

"""Persist complete numerical event histories; all plots read these values."""
function save_event(parent,id,r)
    g = create_group(parent,string(id))
    for (name,value) in pairs(r.seeds);g["seeds/"*string(name)] = value;end
    dg = save_table(g,"diffusion",r.diffused.bins)
    for name in (:x_edges_mm,:y_edges_mm,:z_edges_mm)
        dg[string(name)] = getproperty(r.diffused,name)
    end
    if r.selection.accepted
        a,h = r.axes,r.history
        hg = save_table(g,"pad_history",h.counts)
        for name in (:x_edges_mm,:y_edges_mm,:z_edges_mm,:time_edges_ms)
            hg[string(name)] = getproperty(h,name)
        end
        hg["drift_direction"] = h.drift_direction
        hg["time_origin_z_mm"] = h.time_origin_z_mm
        hg["reference_mm"] = collect(a.reference_mm)
        hg["grid_offset_mm"] = collect(a.grid_offset_mm)
        hg["first_edge_indices"] = collect(a.first_edge_indices)
        save_table(g,"backbone",r.backbone.path_voxels)
        g["backbone/path_indices"] = r.backbone.path_indices
        g["backbone/length_mm"] = r.backbone.total_length
    end
    return g
end

function main()
    config = TOML.parsefile(CARD)
    detector_path = joinpath(ROOT,config["inputs"]["detector_card"])
    detector = TOML.parsefile(detector_path)
    record = TOML.parsefile(joinpath(ROOT,"work/event_analysis/step3_validation.toml"))
    record["status"] == "passed" || error("Run step-three tests first")
    for (path,hash) in record["files"]
        bytes2hex(sha256(read(joinpath(ROOT,path)))) == hash || error("Changed after tests: $path")
    end
    mkpath(dirname(OUT))
    h5open(OUT,"w") do f
        attributes(f)["schema"] = "nausicaa.event.step3.validation.v1"
        attributes(f)["created_utc"] = string(now(UTC))
        attributes(f)["julia_version"] = string(VERSION)
        attributes(f)["scope"] = "Validation subset: source-frame hits, common event L, per-ion diffusion, exact ion threshold, exactly one component, full retained charge in pads."
        attributes(f)["rng_scheme"] = "nausicaa.event.rng.v1; length-prefixed SHA-256 keys; first eight bytes little endian; Xoshiro"
        f["configuration/card"] = read(CARD,String)
        f["configuration/detector"] = read(detector_path,String)
        f["tests/passed"] = record["counts"]["passed"]
        paths = sort(collect(keys(record["files"])))
        append!(paths,[relpath(CARD,ROOT),relpath(detector_path,ROOT),"scripts/validation/validate_event_step3.jl"])
        f["provenance/paths"] = paths
        f["provenance/sha256"] = [bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in paths]
        for relative in config["inputs"]["samples"]
            input = joinpath(ROOT,relative)
            sample = replace(basename(input),"_hitdf.h5"=>"")
            h5open(input,"r") do source
                sg = create_group(f,"samples/"*sample)
                sg["input_sha256"] = bytes2hex(open(sha256,input))
                ids = read(source["events/event_id"])
                n = min(Int(config["validation"]["selection_events_per_sample"]),length(ids))
                n > 0 || error("selection_events_per_sample must be positive")
                selected = ids[unique(round.(Int,range(1,length(ids);length=n)))]
                event_group = create_group(sg,"events")
                rows = NamedTuple[]
                example_saved = false
                # Compile every branch using a small deterministic accepted event.
                warm = EventHits(0,source_frame,Float32[0],Float32[0],Float32[0],Float32[2.5],Float32[0],Int32[7],Int32[0])
                process_event(warm,config,detector,sample)
                for id in selected
                    hits,_ = read_event(source,id)
                    timed = @timed process_event(hits,config,detector,sample)
                    r = timed.value
                    save_event(event_group,id,r)
                    accepted = r.selection.accepted
                    pad_total = accepted ? r.history.total_ions : Int64(0)
                    push!(rows,(event_id=id,accepted=Int8(accepted),n_tracks=length(r.selection.tracks),
                        energy_MeV=r.counts.event_energy_MeV,expected_ions=r.counts.total_expected,
                        sampled_ions=r.diffused.sampled_ions,retained_ions=r.diffused.retained_ions,
                        outside_ions=r.diffused.outside_ions,graph_retained_ions=r.selection.graph_retained_ions,
                        graph_below_threshold_ions=r.selection.graph_below_threshold_ions,pad_ions=pad_total,
                        drift_cm=r.L,sigma_t_mm=r.sigma,sigma_l_mm=r.sigma_z,
                        bin_width_mm=r.b,max_distance_mm=r.distance,
                        grid_offset_x_mm=r.offset[1],grid_offset_y_mm=r.offset[2],
                        n_pad_entries=accepted ? nrow(r.history.counts) : 0,
                        path_length_mm=accepted ? r.backbone.total_length : 0f0,
                        seconds=timed.time,allocated_bytes=timed.bytes))
                    if accepted && !example_saved
                        replay = process_event(hits,config,detector,sample)
                        @assert replay.history.counts == r.history.counts
                        @assert replay.seeds == r.seeds && replay.L == r.L && replay.offset == r.offset
                        proj = pad_projections(r.history)
                        ex = create_group(sg,"example")
                        ex["event_id"] = id
                        ex["xy_ions"] = proj.xy_ions
                        ex["slice_ions"] = proj.slice_ions
                        # Strongest pad's complete sequence including zero slices.
                        totals = combine(groupby(r.history.counts,[:pad_ix,:pad_iy]),:ions=>sum=>:ions)
                        best = totals[argmax(totals.ions),:]
                        pulse = zeros(Int64,length(proj.slice_ions))
                        for row in eachrow(r.history.counts)
                            row.pad_ix == best.pad_ix && row.pad_iy == best.pad_iy && (pulse[row.slice_index] = row.ions)
                        end
                        ex["strongest_pad_indices"] = Int64[best.pad_ix,best.pad_iy]
                        ex["strongest_pad_ions"] = pulse
                        example_saved = true
                    end
                end
                table = DataFrame(rows)
                save_table(sg,"event_summary",table)
                accepted = table.accepted .== 1
                sg["summary/tested_events"] = nrow(table)
                sg["summary/accepted_events"] = count(accepted)
                sg["summary/acceptance_fraction"] = Float32(count(accepted))/Float32(nrow(table))
                sg["summary/accepted_retained_ions"] = sum(table.retained_ions[accepted])
                sg["summary/pad_ions"] = sum(table.pad_ions)
                sg["summary/accepted_below_graph_threshold_ions"] = sum(table.graph_below_threshold_ions[accepted])
                sg["summary/mean_seconds"] = mean(table.seconds)
                @assert sum(table.pad_ions) == sum(table.retained_ions[accepted])
                println(sample,": ",count(accepted),"/",nrow(table)," accepted; pad ledger exact; mean ",round(mean(table.seconds)*1000;digits=1)," ms")
            end
        end
    end
    # Read-back validation of persisted histories and every event charge ledger.
    h5open(OUT,"r") do f
        for sample in keys(f["samples"])
            sg = f["samples/"*sample];summary = sg["event_summary"]
            for (id,accepted,retained,pad) in zip(read(summary["event_id"]),read(summary["accepted"]),read(summary["retained_ions"]),read(summary["pad_ions"]))
                eg = sg["events/"*string(id)]
                @assert sum(read(eg["diffusion/ions"])) == retained
                if accepted == 1
                    @assert sum(read(eg["pad_history/ions"])) == pad == retained
                else
                    @assert !haskey(eg,"pad_history") && pad == 0
                end
            end
        end
    end
    println(OUT)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
