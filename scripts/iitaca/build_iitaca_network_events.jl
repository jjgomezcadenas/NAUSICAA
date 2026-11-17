#!/usr/bin/env julia
#=
iItaca: training events of one run for the classification network.

  julia --project=. scripts/iitaca/build_iitaca_network_events.jl <run folder> [OUTPUT.h5] [--window=LO,HI]

For each batch it reads features.h5 and hits.h5 (run_iitaca_batch.jl,
prepare_iitaca_features.jl), keeps the events whose event energy is in the ROI
and that pass the track selection (in_roi and selected), and places their
itrack hits on the readout lattice inside the fixed frame, one hit per cell.
Writes results/network/<variant>/events/<sample>/events.h5 by default, in the
format of network_events_builder.jl (schema itaca.network_events.v1) with the
itrack only; network_trainer.py reads it with --input=itrack. The event table
holds the columns of the feature tables; the track features also appear with
the itrack_ prefix used by the trainer (itrack_E1_keV, itrack_E2_keV,
itrack_length_mm).

With --window=LO,HI (keV) the energy cut is the box energy energy_true_keV in
[LO, HI] instead of in_roi, with a selection text of its own; the default
output is then results/network/<variant>/events/<run label>_window/events.h5,
and roi/selected_events counts the events in the ROI that pass the track
selection.
=#
include(joinpath(@__DIR__, "../network/network_events_builder.jl"))     # lattice_cells, lattice_indices, one_hit_per_cell; GuardedFiles

const IITACA_FEATURES = "iitaca.features.v1"
const TRACK_COLUMNS = ("hits", "components", "main_keV", "second_keV", "main_voxels", "length_mm", "extent_x_mm", "extent_y_mm", "extent_z_mm", "E1_keV", "E2_keV")

"""Pitches (x, y, z in mm) of the ion lattice from the simulation inputs of a hit file: pad pitch, and ion velocity x frame."""
function iitaca_pitch(f)
    s = f["configuration/simulation_inputs"]; pad = Float64(read(s["ions.pad_pitch_mm"]))
    (pad, pad, read(s["drift.ion_drift_velocity_cm_s"]) / 100 * read(s["ions.frame_ms"]))
end

function build_iitaca_network_events(run_dir::AbstractString, target::AbstractString = ""; window = nothing)
    window === nothing || (length(window) == 2 && window[1] < window[2]) || error("The window needs LO < HI, found $window")
    run_dir = abspath(run_dir); root = joinpath(run_dir, "batches")
    tags = sort(filter(t -> occursin(r"^b\d{5}$", t) && isfile(joinpath(root, t, "features.h5")), readdir(root))); isempty(tags) && error("No feature tables under $root")
    started = time(); columns = Dict{String,Vector}(); accepted = 0; in_roi_events = 0; roi_selected = 0; shas = String[]
    first = nothing; pitch = nothing; side = 0.0; n = (0, 0, 0)
    a = (offsets = Int64[0], ix = UInt8[], iy = UInt8[], iz = UInt8[], energy = Float32[]); dropped = 0
    for t in tags
        dir = joinpath(root, t); hits = joinpath(dir, "hits.h5")
        h5open(joinpath(dir, "features.h5"), "r") do ft
            read(ft["schema"]) == IITACA_FEATURES || error("Unexpected schema in $dir/features.h5")
            ident = Dict(k => read(ft[k]) for k in ("sample", "variant", "frame_track", "settings", "roi_lower_keV", "roi_upper_keV", "run_label"))
            ident["run_label"] = dirname(dirname(ident["run_label"]))
            first === nothing ? (first = ident) : for k in keys(ident); ident[k] == first[k] || error("Setting $k differs in batch $t"); end
            bytes2hex(open(sha256, hits)) == read(ft["hit_file_sha256"]) || error("Feature table and hit file of $t do not match")
            push!(shas, read(ft["hit_file_sha256"]))
            e = Dict(k => read(ft["events"][k]) for k in keys(ft["events"]))
            accepted += length(e["event_number"]); in_roi_events += count(==(1), e["in_roi"])
            roi_selected += count((e["in_roi"] .== 1) .& (e["selected"] .== 1))
            energy_cut = window === nothing ? e["in_roi"] .== 1 : window[1] .<= e["energy_true_keV"] .<= window[2]
            keep = findall(energy_cut .& (e["selected"] .== 1))
            h5open(hits, "r") do f
                read(f["events/event_number"]) == e["event_number"] || error("Event order differs between features and hits in $t")
                p = iitaca_pitch(f); lower = [read(f["events/frame_$(c)_min_mm"]) for c in "xyz"]; upper = [read(f["events/frame_$(c)_max_mm"]) for c in "xyz"]
                s = Float64(upper[1][1] - lower[1][1]); all(i -> all(isapprox.(upper[i] .- lower[i], s; atol = 1f-2)), 1:3) || error("Frame side varies in $t")
                pitch === nothing ? (pitch = p; side = s; n = ntuple(i -> lattice_cells(side, pitch[i]), 3)) : (p == pitch && s == side || error("Lattice or frame differs in $t"))
                all(<=(256), n) || error("Grid too large for 8-bit indices: $n")
                off = read(f["events/hit_offset"]); cnt = read(f["events/hit_count"]); xyz = [read(f["itrack/$(c)_mm"]) for c in "xyz"]; energy = read(f["itrack/energy_keV"])
                for row in keep
                    r = off[row]+1:off[row]+cnt[row]
                    idx = [lattice_indices(view(xyz[i], r), lower[i][row], pitch[i], n[i]) for i in 1:3]
                    en = view(energy, r); kept = one_hit_per_cell(idx[1], idx[2], idx[3], en); dropped += length(en) - length(kept)
                    append!(a.ix, idx[1][kept]); append!(a.iy, idx[2][kept]); append!(a.iz, idx[3][kept]); append!(a.energy, en[kept]); push!(a.offsets, length(a.energy))
                end
            end
            for (k, v) in e; append!(get!(columns, k, similar(v, 0)), v[keep]); end
        end
        println("batch $t: $(length(columns["event_number"])) of $accepted events kept, $(round(time() - started; digits = 1)) s"); flush(stdout)
    end
    allunique(columns["event_number"]) || error("Duplicate event numbers across batches")
    for k in TRACK_COLUMNS; columns["itrack_$k"] = columns[k]; end
    isempty(target) && (target = joinpath(ROOT, "results/network", first["variant"], "events", window === nothing ? first["sample"] : first["run_label"] * "_window", "events.h5"))
    target = abspath(target)
    isfile(target) && error("Preserve existing events: $target"); mkpath(dirname(target))
    budget = 16 * 1024^2 + 2 * (Base.summarysize(columns) + Base.summarysize(a)); GuardedFiles.space_check(dirname(target), budget, 1024^3); tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = NETWORK_EVENTS_SCHEMA; f["species"] = first["sample"]; f["run_label"] = first["run_label"]; f["batches"] = tags; f["hit_file_sha256"] = shas
            f["selection"] = window === nothing ? "iItaca $(first["variant"]): ROI cut on the event energy, $(first["roi_lower_keV"]) to $(first["roi_upper_keV"]) keV; track selection of the itrack on the readout lattice (selected)" :
                "iItaca $(first["variant"]): window cut on the box energy (energy_true_keV), $(Float64(window[1])) to $(Float64(window[2])) keV; track selection of the itrack on the readout lattice (selected)"
            window === nothing || (f["window/lower_keV"] = Float64(window[1]); f["window/upper_keV"] = Float64(window[2]); f["window/source"] = "features.h5 of the batches (energy_true_keV)")
            f["accepted_events"] = accepted; f["frame_side_mm"] = side
            f["roi/lower_keV"] = first["roi_lower_keV"]; f["roi/upper_keV"] = first["roi_upper_keV"]; f["roi/events"] = in_roi_events; f["roi/source"] = "features.h5 of the batches (in_roi)"; f["roi/selected_events"] = roi_selected
            f["settings"] = first["settings"]; f["variant"] = first["variant"]
            for (k, v) in columns; f["events/$k"] = v; end
            f["itrack/pitch_mm"] = collect(pitch); f["itrack/cells"] = collect(n); f["itrack/offsets"] = a.offsets; f["itrack/hits_dropped_from_shared_cells"] = dropped
            for (k, v) in (("ix", a.ix), ("iy", a.iy), ("iz", a.iz), ("energy_keV", a.energy)); f["itrack/$k", chunk = (min(length(v), 65536),), shuffle = true, deflate = 4] = v; end
            f["time_s"] = time() - started; flush(f); GuardedFiles.check_size(tmp, budget)
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("NETWORK EVENTS $target events=$(length(columns["event_number"])) of $accepted ($in_roi_events in the ROI) hits dropped from shared cells: $dropped bytes=$(filesize(target)) time_s=$(round(time() - started; digits = 1))")
    target
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = filter(a -> startswith(a, "--"), ARGS); args = filter(a -> !startswith(a, "--"), ARGS)
    all(o -> startswith(o, "--window="), opts) && length(opts) <= 1 && 1 <= length(args) <= 2 ||
        error("Usage: build_iitaca_network_events.jl <run folder> [OUTPUT.h5] [--window=LO,HI]")
    window = isempty(opts) ? nothing : parse.(Float64, split(opts[1][length("--window=")+1:end], ","))
    build_iitaca_network_events(args...; window)
end
