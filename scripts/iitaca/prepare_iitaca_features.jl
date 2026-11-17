#!/usr/bin/env julia
#=
iItaca: selection and track features of the accepted events of one batch.

  julia --project=. scripts/iitaca/prepare_iitaca_features.jl <folder of a batch, with hits.h5>
  julia --project=. scripts/iitaca/prepare_iitaca_features.jl <folder of a run>          # every batch without features

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
"""
function compute_features(path::AbstractString)
    out = Dict{String,Any}("hit_file_sha256" => bytes2hex(open(sha256, path)))
    h5open(path, "r") do f
        read(f["schema"]) == HITS_SCHEMA || error("Expected schema $HITS_SCHEMA in $path")
        g = f["configuration"]
        for k in COPIED; out[k] = read(g[k]); end
        out["generated_events"] = read(g["generated_events"])
        detector = TOML.parse(out["resolved"])["detector"]; en = detector["energy"]
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

"""Write the feature table of a batch beside its hit file; an existing table is preserved."""
function write_features(batch_dir::AbstractString)
    target = joinpath(batch_dir, "features.h5"); isfile(target) && error("Preserve existing features: $target")
    started = time(); t = compute_features(joinpath(batch_dir, "hits.h5"))
    budget = 16 * 1024^2 + 2 * Base.summarysize(t["events"]); GuardedFiles.space_check(batch_dir, budget, 1024^3); tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = FEATURE_SCHEMA
            for k in (COPIED..., "generated_events", "hit_file_sha256", "roi_lower_keV", "roi_upper_keV", "settings"); f[k] = t[k]; end
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

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: prepare_iitaca_features.jl <batch folder or run folder>")
    dir = abspath(ARGS[1])
    if isfile(joinpath(dir, "hits.h5"))
        write_features(dir)
    else
        for b in batch_folders(dir); isfile(joinpath(b, "features.h5")) || write_features(b); end
    end
end
