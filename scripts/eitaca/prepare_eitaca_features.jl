#!/usr/bin/env julia
#=
eITACA: selection and track features of the accepted events of one batch, and
the statistics of a run.

  julia --project=. scripts/eitaca/prepare_eitaca_features.jl <folder of a batch, with electrons.h5> [--tracks=CARD]
  julia --project=. scripts/eitaca/prepare_eitaca_features.jl <folder of a run> [--tracks=CARD]   # every batch without features, then the statistics

Reads electrons.h5 (run_eitaca_batch.jl) and the tracks card (default
config/eitaca/tracks.toml) and writes features.h5 beside it (schema
eitaca.features.v1), one row per event accepted by the fast trigger:

 1. The event columns of the electron file and the trigger counts.
 2. energy_keV: the energy in the simulated box (energy_true_keV) smeared with
    the resolution of the tracks card (IItaca.emulated_energy_keV, stream
    "emulated_energy" of the run seed, as iItaca). The reconstructed electron
    energy stays in electron_energy_keV.
 3. in_roi: the ROI cut on energy_keV.
 4. Track selection of the etrack (TrackSelection.track_selection): segments,
    deposits, deposit_keV, segments_joined, selected.
 5. Main track (TrackSelection.track_features): hits, components, main_keV,
    second_keV, main_voxels, length_mm, extents, E1_keV >= E2_keV.

For a run folder it then writes statistics/events.h5 (schema
eitaca.statistics.v1, the layout of iitaca.statistics.v1): the feature tables
joined, the counts summed and the cut flow generated, read, trigger, roi, track.
An existing feature table is preserved; the statistics are rewritten, they hold
nothing that is not in the feature tables.
=#
using HDF5, TOML, SHA, Random
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl")); using .EventAnalysis
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
include(joinpath(ROOT, "src/iitaca/IItaca.jl")); using .IItaca
include(joinpath(ROOT, "src/ion_signal/guarded_files.jl"))
include(joinpath(ROOT, "src/track_selection/TrackSelection.jl")); using .TrackSelection

const FEATURE_SCHEMA = "eitaca.features.v1"
const STATISTICS_SCHEMA = "eitaca.statistics.v1"
const ELECTRONS_SCHEMA = "eitaca.electrons.v1"
const DEFAULT_TRACKS = joinpath(ROOT, "config/eitaca/tracks.toml")
"Event columns of the electron file copied to the feature table."
const COPIED_COLUMNS = ("event_number", "energy_true_keV", "electron_energy_keV", "electron_drift_cm", "sigma_t_mm", "sigma_l_mm",
    "track_distance_to_el_mm", "track_z_length_mm", "dark_pulses", "channels", "buckets", "clipped_buckets", "invalid_electron_buckets", "hit_count")

"""Energy card of the analysis: Qbb of the resolved eITACA energy card and the resolution and ROI of the tracks card."""
analysis_energy(resolved_energy, tracks) = merge(Dict{String,Any}(resolved_energy), tracks["energy"])

"""
    compute_features(electrons_path, tracks_card) -> Dict

Feature table of one batch: `events` (column => vector), `counts` (events per
trigger decision), the configuration, the ROI and the settings used.
"""
function compute_features(path::AbstractString, tracks_card::AbstractString)
    tracks = TOML.parsefile(tracks_card); tracks["schema"] == "nausicaa.eitaca_tracks.v1" || error("Unexpected tracks card $tracks_card")
    out = Dict{String,Any}("electron_file_sha256" => bytes2hex(open(sha256, path)), "tracks_card" => relpath(tracks_card, ROOT),
                           "tracks_card_sha256" => bytes2hex(open(sha256, tracks_card)))
    h5open(path, "r") do f
        read(f["schema"]) == ELECTRONS_SCHEMA || error("Expected schema $ELECTRONS_SCHEMA in $path")
        g = f["configuration"]
        for k in ("run_label", "sample", "variant", "batch", "resolved"); out[k] = read(g[k]); end
        out["generated_events"] = read(g["generated_events"]); seed = read(g["seed"])
        resolved = TOML.parse(out["resolved"]); out["frame_track"] = resolved["detector"]["trigger"]["frame_track"]
        en = analysis_energy(resolved["detector"]["energy"], tracks); roi = roi_keV(en)
        d = Dict("tracks" => tracks); selection = selection_settings(d); features = feature_settings(d)
        out["roi_lower_keV"] = roi.lower; out["roi_upper_keV"] = roi.upper
        out["settings"] = sprint(io -> TOML.print(io, Dict("tracks" => tracks, "energy" => en); sorted = true))
        out["counts"] = Dict{String,Int64}(k => read(f["selection/counts"][k]) for k in keys(f["selection/counts"]))
        ev = f["events"]; events = Dict{String,Any}(k => read(ev[k]) for k in COPIED_COLUMNS)
        offset = read(ev["hit_offset"]); n = length(events["event_number"]); n > 0 || error("No accepted events in $path")
        x = read(f["etrack/x_mm"]); y = read(f["etrack/y_mm"]); z = read(f["etrack/z_mm"]); e = read(f["etrack/energy_keV"])
        rows = NamedTuple[]
        for j in 1:n
            id = events["event_number"][j]
            E = emulated_energy_keV(Xoshiro(IItaca.event_seed(seed, id, "emulated_energy")), events["energy_true_keV"][j], en)
            k = offset[j]+1:offset[j]+events["hit_count"][j]
            s = track_selection(view(x, k), view(y, k), view(z, k), view(e, k); selection...)
            t = TrackSelection.track_features(x[k], y[k], z[k], e[k]; features...)
            push!(rows, merge((energy_keV = E, in_roi = Int8(in_roi(en, E)), segments = s.segments, deposits = s.deposits, deposit_keV = s.deposit_keV,
                segments_joined = Int8(s.segments_joined), selected = Int8(s.selected)), t))
        end
        for k in keys(first(rows)); events[string(k)] = [getfield(row, k) for row in rows]; end
        events["batch"] = fill(parse(Int, out["batch"][2:end]), n)
        out["events"] = events
    end
    out
end

const SAVED = ("run_label", "sample", "variant", "batch", "frame_track", "resolved", "generated_events", "electron_file_sha256",
               "tracks_card", "tracks_card_sha256", "roi_lower_keV", "roi_upper_keV", "settings")

"""Write the feature table of a batch beside its electron file; an existing table is preserved."""
function write_features(batch_dir::AbstractString, tracks_card::AbstractString)
    started = time(); target = joinpath(batch_dir, "features.h5"); isfile(target) && error("Preserve existing features: $target")
    t = compute_features(joinpath(batch_dir, "electrons.h5"), tracks_card)
    budget = 16 * 1024^2 + 2 * Base.summarysize(t["events"]); GuardedFiles.space_check(batch_dir, budget, 1024^3); tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = FEATURE_SCHEMA
            for k in SAVED; f[k] = t[k]; end
            for (k, v) in t["events"]; f["events/$k"] = v; end
            for (k, v) in t["counts"]; f["selection/counts/$k"] = v; end
            f["code_sha256"] = bytes2hex(open(sha256, @__FILE__)); f["time_s"] = time() - started
            flush(f); GuardedFiles.check_size(tmp, budget)
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("FEATURES $target events=$(length(t["events"]["event_number"])) selected=$(count(==(1), t["events"]["selected"])) time_s=$(round(time() - started; digits = 1))")
    target
end

"""Read a feature table written by `write_features`."""
function read_features(path::AbstractString)
    h5open(path, "r") do f
        read(f["schema"]) == FEATURE_SCHEMA || error("Unexpected schema in $path")
        t = Dict{String,Any}(k => read(f[k]) for k in SAVED)
        t["events"] = Dict{String,Any}(k => read(f["events"][k]) for k in keys(f["events"]))
        t["counts"] = Dict{String,Int64}(k => read(f["selection/counts"][k]) for k in keys(f["selection/counts"]))
        t
    end
end

"""Folders of the batches of a run that hold `file`, in batch order."""
function batch_folders(run_dir::AbstractString; file::AbstractString = "electrons.h5")
    base = joinpath(run_dir, "batches"); isdir(base) || error("No batches folder in $run_dir")
    [joinpath(base, d) for d in sort(readdir(base)) if isfile(joinpath(base, d, file))]
end

"""Write statistics/events.h5 of a run from the feature tables of its batches and return its path."""
function write_statistics(run_dir::AbstractString)
    tables = [read_features(joinpath(b, "features.h5")) for b in batch_folders(run_dir; file = "features.h5")]
    isempty(tables) && error("No feature tables in $run_dir"); t1 = first(tables)
    for t in tables, k in ("sample", "variant", "frame_track", "settings", "roi_lower_keV", "roi_upper_keV", "tracks_card_sha256")
        t[k] == t1[k] || error("Batches $(t1["batch"]) and $(t["batch"]) differ in $k")
    end
    events = Dict{String,Any}(k => vcat((t["events"][k] for t in tables)...) for k in keys(t1["events"]))
    allunique(events["event_number"]) || error("Repeated event numbers in $run_dir")
    counts = Dict{String,Int64}(k => sum(t["counts"][k] for t in tables) for k in keys(t1["counts"]))
    counts["accepted_events"] == length(events["event_number"]) || error("Accepted events and table rows differ")
    generated = all(t["generated_events"] > 0 for t in tables) ? sum(t["generated_events"] for t in tables) : 0
    roi = events["in_roi"] .== 1; track = roi .& (events["selected"] .== 1)
    flow = [("generated", generated), ("read", counts["input_events"]), ("trigger", counts["accepted_events"]), ("roi", count(roi)), ("track", count(track))]
    dir = joinpath(run_dir, "statistics"); mkpath(dir); target = joinpath(dir, "events.h5")
    budget = 16 * 1024^2 + 2 * Base.summarysize(events); GuardedFiles.space_check(dir, budget, 1024^3); tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = STATISTICS_SCHEMA; f["run_label"] = dirname(dirname(t1["run_label"]))
            for k in ("sample", "variant", "frame_track", "settings", "resolved", "tracks_card", "tracks_card_sha256", "roi_lower_keV", "roi_upper_keV"); f[k] = t1[k]; end
            f["batches"] = [t["batch"] for t in tables]; f["electron_file_sha256"] = [t["electron_file_sha256"] for t in tables]
            for (k, v) in events; f["events/$k"] = v; end
            for (k, v) in counts; f["selection/counts/$k"] = v; end
            f["cut_flow/step"] = first.(flow); f["cut_flow/events"] = last.(flow)
            flush(f); GuardedFiles.check_size(tmp, budget)
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("STATISTICS $target batches=$(length(tables)) " * join(("$name=$n" for (name, n) in flow), " "))
    target
end

if abspath(PROGRAM_FILE) == @__FILE__
    k = findall(startswith("--tracks="), ARGS); card = isempty(k) ? DEFAULT_TRACKS : abspath(split(ARGS[k[1]], "=")[2])
    args = ARGS[setdiff(eachindex(ARGS), k)]
    length(args) == 1 || error("Usage: prepare_eitaca_features.jl <batch folder or run folder> [--tracks=CARD]")
    dir = abspath(args[1])
    if isfile(joinpath(dir, "electrons.h5"))
        isfile(joinpath(dir, "features.h5")) ? println("Kept $(joinpath(dir, "features.h5"))") : write_features(dir, card)
    else
        for b in batch_folders(dir); isfile(joinpath(b, "features.h5")) || write_features(b, card); end
        write_statistics(dir)
    end
end
