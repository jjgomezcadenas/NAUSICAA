#!/usr/bin/env julia
#=
iItaca: selection and track features of the accepted events of one batch.

  julia --project=. scripts/iitaca/prepare_iitaca_features.jl <folder of a batch, with hits.h5>
  julia --project=. scripts/iitaca/prepare_iitaca_features.jl <folder of a run>          # every batch without features
  julia --project=. scripts/iitaca/prepare_iitaca_features.jl <folder> --tracks-variant=<variant card>

Reads hits.h5 (run_iitaca_batch.jl) and writes features.h5 beside it (schema
iitaca.features.v1): one row per accepted event with the columns of the hit
file, the ROI cut on the event energy, the track selection and the features of
the main track. The settings are those of the resolved configuration stored
in the hit file (tracks and energy cards). The table is enough to select
events and to build the statistics of a run without its hit files.
=#
using HDF5, TOML, SHA
const ROOT = normpath(joinpath(@__DIR__, "../.."))
if !isdefined(@__MODULE__, :IItaca)
    include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl")); using .EventAnalysis
    include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
    include(joinpath(ROOT, "src/iitaca/IItaca.jl")); using .IItaca
end
isdefined(@__MODULE__, :GuardedFiles) || include(joinpath(ROOT, "src/ion_signal/guarded_files.jl"))
isdefined(@__MODULE__, :TrackSelection) || include(joinpath(ROOT, "src/track_selection/TrackSelection.jl"))
using .TrackSelection

const FEATURE_SCHEMA = "iitaca.features.v1"
const HITS_SCHEMA = "iitaca.hits.v1"
"Text entries of the configuration copied from the hit file to the feature table."
const COPIED = ("run_label", "sample", "variant", "batch", "frame_track", "resolved")

"""
    compute_features(hits_path) -> Dict

Feature table of one batch. `events`: column name => vector, one entry per
accepted event: the event columns of the hit file, `batch` (number), `in_roi`,
`segments`, `deposits`, `deposit_keV`, `segments_joined`, `selected` (track
selection) and the track features `hits`, `components`, `main_keV`,
`second_keV`, `main_voxels`, `length_mm`, `extent_x_mm`, `extent_y_mm`,
`extent_z_mm`, `E1_keV`, `E2_keV`. `counts`: events per trigger decision.
Also the copied configuration, the generated events, the ROI bounds, the
settings used and the SHA-256 of the hit file.

With `tracks_variant` (a variant card), the track selection and the features
use the tracks card of that variant, which must give the detector of the hit
file but for `tracks.*`: it must differ from the variant of the hit file only
there, with the same saved hit cut; the overrides of the run are kept. `variant`
and `run_label` are then those of that variant (the variant name of the hit
file replaced in the run label), and `hits_run_label` keeps the run of the
hit file.
"""
function compute_features(path::AbstractString; tracks_variant = nothing)
    out = Dict{String,Any}("hit_file_sha256" => bytes2hex(open(sha256, path)))
    h5open(path, "r") do f
        read(f["schema"]) == HITS_SCHEMA || error("Expected schema $HITS_SCHEMA in $path")
        g = f["configuration"]
        for k in COPIED; out[k] = read(g[k]); end
        out["generated_events"] = read(g["generated_events"])
        detector = TOML.parse(out["resolved"])["detector"]; en = detector["energy"]
        out["hits_run_label"] = out["run_label"]
        if tracks_variant !== nothing
            vr = TOML.parse(out["resolved"])["variant"]; card = vr["card"]; isempty(card) && error("The hit file $path names no variant card")
            bytes2hex(open(sha256, joinpath(ROOT, card))) == vr["sha256"] || error("$card changed since the hit file was written")
            index = joinpath(ROOT, "config/iitaca/iitaca.toml")
            v = resolve_detector(index; variant = tracks_variant, root = ROOT); w = resolve_detector(index; variant = joinpath(ROOT, card), root = ROOT)
            all(v.detector[s] == w.detector[s] for s in SECTIONS if s != "tracks") || error("$(tracks_variant) differs from $card outside tracks")
            v.detector["tracks"]["saved_hit_cut_keV"] == detector["tracks"]["saved_hit_cut_keV"] || error("$(tracks_variant) changes the saved hit cut of the hit file")
            occursin(out["variant"], out["run_label"]) || error("Run label $(out["run_label"]) does not hold the variant name $(out["variant"])")
            out["run_label"] = replace(out["run_label"], out["variant"] => v.variant; count = 1); out["variant"] = v.variant
            detector = merge(detector, Dict("tracks" => v.detector["tracks"]))
        end
        selection = selection_settings(detector); features = feature_settings(detector); roi = roi_keV(en)
        out["roi_lower_keV"] = roi.lower; out["roi_upper_keV"] = roi.upper
        out["settings"] = sprint(io -> TOML.print(io, Dict("tracks" => detector["tracks"], "energy" => en); sorted = true))
        out["counts"] = Dict{String,Int64}(k => read(f["selection/counts"][k]) for k in keys(f["selection/counts"]))
        events = Dict{String,Any}(k => read(f["events"][k]) for k in keys(f["events"]))
        n = length(events["event_number"]); n > 0 || error("No accepted events in $path")
        x = read(f["itrack/x_mm"]); y = read(f["itrack/y_mm"]); z = read(f["itrack/z_mm"]); e = read(f["itrack/energy_keV"])
        rows = NamedTuple[]
        for j in 1:n
            k = events["hit_offset"][j]+1:events["hit_offset"][j]+events["hit_count"][j]
            s = track_selection(view(x, k), view(y, k), view(z, k), view(e, k); selection...)
            t = TrackSelection.track_features(x[k], y[k], z[k], e[k]; features...)
            push!(rows, merge((in_roi = Int8(in_roi(en, events["energy_keV"][j])), segments = s.segments, deposits = s.deposits, deposit_keV = s.deposit_keV,
                segments_joined = Int8(s.segments_joined), selected = Int8(s.selected)), t))
        end
        for k in keys(first(rows)); events[string(k)] = [getfield(row, k) for row in rows]; end
        events["batch"] = fill(parse(Int, out["batch"][2:end]), n)
        out["events"] = events
    end
    out
end

"""
    write_features(batch_dir; tracks_variant=nothing)

Write the feature table of a batch beside its hit file; an existing table is
preserved. With `tracks_variant`, the table goes to the batch folder of the run
of that variant, beside the hit file of its run (`<run>/batches/<tag>`, next to
the run of the hit file), with a relative link `hits.h5` to the shared hit file.
"""
function write_features(batch_dir::AbstractString; tracks_variant = nothing)
    started = time(); t = compute_features(joinpath(batch_dir, "hits.h5"); tracks_variant)
    if tracks_variant !== nothing
        tag = basename(batch_dir); run = dirname(dirname(t["run_label"])); hit_run = dirname(dirname(t["hits_run_label"]))
        out = joinpath(dirname(dirname(dirname(batch_dir))), run, "batches", tag); mkpath(out)
        link = joinpath(out, "hits.h5")
        islink(link) || symlink(joinpath("..", "..", "..", hit_run, "batches", tag, "hits.h5"), link)
        bytes2hex(open(sha256, link)) == t["hit_file_sha256"] || error("Link $link does not reach the hit file")
        batch_dir = out
    end
    target = joinpath(batch_dir, "features.h5"); isfile(target) && error("Preserve existing features: $target")
    budget = 16 * 1024^2 + 2 * Base.summarysize(t["events"]); GuardedFiles.space_check(batch_dir, budget, 1024^3); tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = FEATURE_SCHEMA
            for k in (COPIED..., "hits_run_label", "generated_events", "hit_file_sha256", "roi_lower_keV", "roi_upper_keV", "settings"); f[k] = t[k]; end
            for (k, v) in t["events"]; f["events/$k"] = v; end
            for (k, v) in t["counts"]; f["selection/counts/$k"] = v; end
            f["time_s"] = time() - started; flush(f); GuardedFiles.check_size(tmp, budget)
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("FEATURES $target events=$(length(t["events"]["event_number"])) time_s=$(round(time() - started; digits = 1))")
    target
end

"""Read a feature table written by `write_features`, in the form of `compute_features`."""
function read_features(path::AbstractString)
    h5open(path, "r") do f
        read(f["schema"]) == FEATURE_SCHEMA || error("Unexpected schema in $path")
        t = Dict{String,Any}(k => read(f[k]) for k in (COPIED..., "generated_events", "hit_file_sha256", "roi_lower_keV", "roi_upper_keV", "settings"))
        t["hits_run_label"] = haskey(f, "hits_run_label") ? read(f["hits_run_label"]) : t["run_label"]
        t["events"] = Dict{String,Any}(k => read(f["events"][k]) for k in keys(f["events"]))
        t["counts"] = Dict{String,Int64}(k => read(f["selection/counts"][k]) for k in keys(f["selection/counts"]))
        t
    end
end

"""Folders of the batches of a run that hold `file` (the hit file by default), in batch order."""
function batch_folders(run_dir::AbstractString; file::AbstractString = "hits.h5")
    base = joinpath(run_dir, "batches"); isdir(base) || error("No batches folder in $run_dir")
    [joinpath(base, d) for d in sort(readdir(base)) if isfile(joinpath(base, d, file))]
end

"""Feature table of a batch, beside its hit file or, with `tracks_variant`, in the run of that variant."""
function features_target(batch_dir::AbstractString; tracks_variant = nothing)
    tracks_variant === nothing && return joinpath(batch_dir, "features.h5")
    label = h5open(f -> (read(f["configuration/run_label"]), read(f["configuration/variant"])), joinpath(batch_dir, "hits.h5"), "r")
    run = replace(dirname(dirname(label[1])), label[2] => load_variant(tracks_variant).name; count = 1)
    joinpath(dirname(dirname(dirname(batch_dir))), run, "batches", basename(batch_dir), "features.h5")
end

if abspath(PROGRAM_FILE) == @__FILE__
    k = findall(startswith("--tracks-variant="), ARGS); tv = isempty(k) ? nothing : abspath(split(ARGS[k[1]], "=")[2]); args = ARGS[setdiff(eachindex(ARGS), k)]
    length(args) == 1 || error("Usage: prepare_iitaca_features.jl <batch folder or run folder> [--tracks-variant=<variant card>]")
    dir = abspath(args[1])
    if isfile(joinpath(dir, "hits.h5"))
        target = features_target(dir; tracks_variant = tv)
        isfile(target) ? println("Kept $target") : write_features(dir; tracks_variant = tv)
    else
        for b in batch_folders(dir); isfile(features_target(b; tracks_variant = tv)) || write_features(b; tracks_variant = tv); end
    end
end
