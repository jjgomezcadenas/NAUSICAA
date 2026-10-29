#!/usr/bin/env julia
# Run six independent cloud-to-grid readout streams in Julia processes.
using HDF5

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const WORKER=joinpath(@__DIR__,"collect_cloud_stream.jl")
const STREAMS=[(sample,carrier) for sample in ("bb0nu","electrons")
               for carrier in ("ion","pure_xe","xe_he")]

"""Parse `--events`, `--workers`, `--clouds` and `--output` options."""
function grid_options(args)
    values=Dict("--events"=>"100","--workers"=>"6",
                "--clouds"=>"results/three_clouds_8000",
                "--output"=>"results/grid_readout_8000")
    length(args)%2==0 || error("Options require values")
    for i in 1:2:length(args)
        haskey(values,args[i]) || error("Unknown option: $(args[i])")
        values[args[i]]=args[i+1]
    end
    events=parse(Int,values["--events"]);workers=parse(Int,values["--workers"])
    events>0 || error("--events must be positive")
    1<=workers<=6 || error("--workers must be in 1:6")
    clouds=abspath(joinpath(ROOT,values["--clouds"]))
    output=abspath(joinpath(ROOT,values["--output"]))
    return (;events,workers,clouds,output)
end

"""Run one independent HDF5 stream and capture elapsed time and peak RSS."""
function grid_stream(sample,carrier,events,clouds,output)
    input=joinpath(clouds,"$(sample)_$(carrier).h5")
    target=joinpath(output,"$(sample)_$(carrier).h5")
    log=joinpath(output,"$(sample)_$(carrier).log")
    isfile(input) || error("Missing diffusion cloud: $input")
    command=addenv(`$(Base.julia_cmd()) --startup-file=no --project=$ROOT $WORKER $sample $carrier $events $input $target`,"JULIA_NUM_THREADS"=>"1")
    println("Starting $sample/$carrier");flush(stdout)
    started=time_ns()
    open(log,"w") do io
        run(pipeline(command,stdout=io,stderr=io))
    end
    elapsed=(time_ns()-started)/1e9
    match_rss=match(r"peak_rss_bytes=(\d+)",read(log,String))
    match_rss===nothing && error("Missing memory measurement in $log")
    println("Finished $sample/$carrier: $(round(elapsed;digits=2)) s");flush(stdout)
    return (;sample,carrier,elapsed,bytes=filesize(target),
            peak_rss_bytes=parse(Int64,match_rss.captures[1]))
end

"""Process all six clouds, writing one readout HDF5 file per stream.

The launcher limits concurrent Julia processes and writes a numerical HDF5
summary after all six finish. The source diffusion clouds remain reusable for
changes to pitch, bucket size, EL yield, PDE and optical acceptance.
"""
function propagate_grids_parallel(events::Int,workers::Int,clouds::String,output::String)
    mkpath(output)
    started=time_ns()
    channel=Channel{NamedTuple}(6)
    @sync for worker in 1:workers
        @async for i in worker:workers:6
            sample,carrier=STREAMS[i]
            put!(channel,grid_stream(sample,carrier,events,clouds,output))
        end
    end
    rows=[take!(channel) for _ in 1:6]
    sort!(rows;by=r->findfirst(==((r.sample,r.carrier)),STREAMS))
    wall=(time_ns()-started)/1e9
    summary=joinpath(output,"grid_readout_summary.h5")
    h5open(summary,"w") do f
        attributes(f)["schema"]="nausicaa.grid_readout_summary.v1"
        f["events_per_stream"]=Int64(events);f["workers"]=Int64(workers)
        f["wall_seconds"]=wall;f["cloud_directory"]=realpath(clouds)
        f["sample"]=[r.sample for r in rows];f["carrier"]=[r.carrier for r in rows]
        f["stream_seconds"]=[r.elapsed for r in rows]
        f["output_bytes"]=Int64[r.bytes for r in rows]
        f["peak_rss_bytes"]=Int64[r.peak_rss_bytes for r in rows]
        f["sum_worker_peak_rss_bytes"]=sum(r.peak_rss_bytes for r in rows)
    end
    println("Grid readout: $(round(wall;digits=2)) s, ",
            round(sum(r.bytes for r in rows)/1e6;digits=1)," MB")
    println("Summary: $summary")
    return summary
end

if abspath(PROGRAM_FILE)==@__FILE__
    opts=grid_options(ARGS)
    propagate_grids_parallel(opts.events,opts.workers,opts.clouds,opts.output)
end
