#!/usr/bin/env julia
#=
iItaca: statistics of a run from the feature tables of its batches.

  julia --project=. scripts/iitaca/prepare_iitaca_statistics.jl <folder of a run>

Reads batches/<tag>/features.h5 of every batch and writes statistics/events.h5
(schema iitaca.statistics.v1): the feature tables joined, the selection counts
summed, and the cut flow. An existing file is replaced: it holds nothing that
is not in the feature tables. Python makes the page
(scripts/iitaca/present_iitaca_statistics.py).

Cut flow, each step applied to the events that passed the previous ones:

 1. generated: events generated in Geant4 (run card), when stated.
 2. read: events in the input files.
 3. trigger: accepted by the fast trigger; the rejected ones are counted per decision.
 4. roi: event energy inside the ROI.
 5. track: selected by the track selection.
=#
include("prepare_iitaca_features.jl")

const STATISTICS_SCHEMA = "iitaca.statistics.v1"

"""Join the feature tables of a run: columns, summed counts, cut flow and identity."""
function run_statistics(run_dir::AbstractString)
    tables = [read_features(joinpath(b, "features.h5")) for b in batch_folders(run_dir; file = "features.h5")]      # the hit files may stay on the cluster
    isempty(tables) && error("No feature tables in $run_dir")
    first_table = first(tables)
    for t in tables, k in ("sample", "variant", "frame_track", "settings", "roi_lower_keV", "roi_upper_keV")
        t[k] == first_table[k] || error("Batches $(first_table["batch"]) and $(t["batch"]) differ in $k")
    end
    events = Dict{String,Any}(k => vcat((t["events"][k] for t in tables)...) for k in keys(first_table["events"]))
    allunique(events["event_number"]) || error("Repeated event numbers in $run_dir")
    counts = Dict{String,Int64}(k => sum(t["counts"][k] for t in tables) for k in keys(first_table["counts"]))
    generated = all(t["generated_events"] > 0 for t in tables) ? sum(t["generated_events"] for t in tables) : 0
    roi = events["in_roi"] .== 1; track = roi .& (events["selected"] .== 1)
    flow = [("generated", generated), ("read", counts["input_events"]), ("trigger", counts["accepted_events"]), ("roi", count(roi)), ("track", count(track))]
    counts["accepted_events"] == length(events["event_number"]) || error("Accepted events and table rows differ")
    (events = events, counts = counts, flow = flow, batches = [t["batch"] for t in tables], run_label = dirname(dirname(first_table["run_label"])),
     sample = first_table["sample"], variant = first_table["variant"], settings = first_table["settings"], resolved = first_table["resolved"],
     roi = (first_table["roi_lower_keV"], first_table["roi_upper_keV"]))
end

"""Write statistics/events.h5 of a run and return its path."""
function write_statistics(run_dir::AbstractString)
    s = run_statistics(run_dir); dir = joinpath(run_dir, "statistics"); mkpath(dir); target = joinpath(dir, "events.h5")
    budget = 16 * 1024^2 + 2 * Base.summarysize(s.events); GuardedFiles.space_check(dir, budget, 1024^3); tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = STATISTICS_SCHEMA
            f["run_label"] = s.run_label; f["sample"] = s.sample; f["variant"] = s.variant; f["settings"] = s.settings; f["resolved"] = s.resolved
            f["batches"] = s.batches; f["roi_lower_keV"] = s.roi[1]; f["roi_upper_keV"] = s.roi[2]
            for (k, v) in s.events; f["events/$k"] = v; end
            for (k, v) in s.counts; f["selection/counts/$k"] = v; end
            f["cut_flow/step"] = [name for (name, _) in s.flow]; f["cut_flow/events"] = [n for (_, n) in s.flow]
            flush(f); GuardedFiles.check_size(tmp, budget)
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("STATISTICS $target batches=$(length(s.batches)) " * join(("$name=$n" for (name, n) in s.flow), " "))
    target
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: prepare_iitaca_statistics.jl <run folder>")
    write_statistics(abspath(ARGS[1]))
end
