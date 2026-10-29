#!/usr/bin/env julia
# Run the six independent Julia diffusion streams with bounded process concurrency.
using HDF5, TOML

const ROOT = normpath(joinpath(@__DIR__,"../.."))
const STREAM_SCRIPT = joinpath(@__DIR__,"propagate_cloud_stream.jl")
const STREAMS = [(sample, carrier) for sample in ("bb0nu", "electrons")
                 for carrier in ("ion", "pure_xe", "xe_he")]

"""Parse `--events N --workers N --output DIRECTORY` for a six-stream run."""
function parse_options(args)
    values = Dict("--events" => "100", "--workers" => "6",
                  "--output" => "results/three_clouds")
    length(args) % 2 == 0 || error("Options require values")
    for i in 1:2:length(args)
        haskey(values, args[i]) || error("Unknown option: $(args[i])")
        values[args[i]] = args[i+1]
    end
    events, workers = parse(Int, values["--events"]), parse(Int, values["--workers"])
    events > 0 || error("--events must be positive")
    1 <= workers <= length(STREAMS) || error("--workers must be between 1 and 6")
    output = abspath(joinpath(ROOT, values["--output"]))
    return (; events, workers, output)
end

"""Run one cloud in its own Julia process and return elapsed time and provenance."""
function run_stream(sample::String, carrier::String, events::Int, output::String)
    target = joinpath(output, "$(sample)_$(carrier).h5")
    log = joinpath(output, "$(sample)_$(carrier).log")
    command = addenv(`$(Base.julia_cmd()) --startup-file=no --project=$ROOT $STREAM_SCRIPT $sample $carrier $events $target`, "JULIA_NUM_THREADS" => "1")
    println("Starting $sample/$carrier"); flush(stdout)
    started = time_ns()
    open(log, "w") do io
        run(pipeline(command, stdout=io, stderr=io))
    end
    elapsed = (time_ns()-started)/1e9
    log_text = read(log, String)
    matched = match(r"peak_rss_bytes=(\d+)", log_text)
    matched === nothing && error("Missing process memory measurement in $log")
    println("Finished $sample/$carrier: $(round(elapsed;digits=2)) s"); flush(stdout)
    return (; sample, carrier, elapsed, bytes=filesize(target),
            peak_rss_bytes=parse(Int64, matched.captures[1]))
end

"""Run six HDF5 cloud streams in up to `workers` independent processes.

Each source event uses a deterministic shared location across its three
carrier models. Each process writes one HDF5 file; the launcher writes an
HDF5 run summary after all streams complete. The wall clock includes Julia
startup, source reads, diffusion and HDF5 output.
"""
function propagate_clouds_parallel(events::Int, workers::Int, output::String)
    mkpath(output)
    started = time_ns()
    rows = Channel{NamedTuple}(length(STREAMS))
    @sync for worker in 1:workers
        @async begin
            for index in worker:workers:length(STREAMS)
                sample, carrier = STREAMS[index]
                put!(rows, run_stream(sample, carrier, events, output))
            end
        end
    end
    observations = [take!(rows) for _ in STREAMS]
    sort!(observations; by=r -> findfirst(==((r.sample, r.carrier)), STREAMS))
    wall = (time_ns()-started)/1e9
    summary = joinpath(output, "propagation_summary.h5")
    h5open(summary, "w") do file
        attributes(file)["schema"] = "nausicaa.cloud_propagation_summary.v1"
        file["workers"] = Int64(workers)
        file["events_per_stream"] = Int64(events)
        file["wall_seconds"] = wall
        file["sample"] = [r.sample for r in observations]
        file["carrier"] = [r.carrier for r in observations]
        file["stream_seconds"] = [r.elapsed for r in observations]
        file["output_bytes"] = Int64[r.bytes for r in observations]
        file["peak_rss_bytes"] = Int64[r.peak_rss_bytes for r in observations]
        file["sum_worker_peak_rss_bytes"] = sum(r.peak_rss_bytes for r in observations)
    end
    println("Six streams: $(round(wall;digits=2)) s, ",
            round(sum(r.bytes for r in observations)/1e6;digits=1), " MB HDF5")
    println("Summary: $summary")
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    options = parse_options(ARGS)
    propagate_clouds_parallel(options.events, options.workers, options.output)
end
