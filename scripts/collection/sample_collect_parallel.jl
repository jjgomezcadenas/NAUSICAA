#!/usr/bin/env julia
# Run six independent direct diffusion-to-grid streams in Julia processes.
using HDF5
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const WORKER=joinpath(@__DIR__,"sample_collect_stream.jl")
const STREAMS=[(sample,carrier) for sample in ("bb0nu","electrons")
               for carrier in ("ion","pure_xe","xe_he")]

"""Parse event count, worker count and versioned output directory."""
function options(args)
    values=Dict("--events"=>"100","--workers"=>"6",
                "--output"=>"results/grid_readout_direct_8000",
                "--ion-card"=>"config/pixel_charge_1mm.toml")
    length(args)%2==0 || error("Options require values")
    for i in 1:2:length(args)
        haskey(values,args[i]) || error("Unknown option: $(args[i])")
        values[args[i]]=args[i+1]
    end
    events=parse(Int,values["--events"]);workers=parse(Int,values["--workers"])
    events>0 || error("--events must be positive")
    1<=workers<=6 || error("--workers must be in 1:6")
    return (;events,workers,output=abspath(joinpath(ROOT,values["--output"])),
            ion_card=realpath(joinpath(ROOT,values["--ion-card"])))
end

"""Run one worker with a private HDF5 file and captured timing/peak RSS."""
function run_stream(sample,carrier,events,output,ion_card)
    target=joinpath(output,"$(sample)_$(carrier).h5")
    log=joinpath(output,"$(sample)_$(carrier).log")
    cmd=addenv(`$(Base.julia_cmd()) --startup-file=no --project=$ROOT $WORKER $sample $carrier $events $target --ion-card $ion_card`,
               "JULIA_NUM_THREADS"=>"1")
    println("Starting $sample/$carrier");flush(stdout)
    started=time_ns()
    open(log,"w") do io
        run(pipeline(cmd,stdout=io,stderr=io))
    end
    elapsed=(time_ns()-started)/1e9
    rss=match(r"peak_rss_bytes=(\d+)",read(log,String))
    rss===nothing && error("Missing peak RSS in $log")
    println("Finished $sample/$carrier: $(round(elapsed;digits=2)) s");flush(stdout)
    return (;sample,carrier,elapsed,bytes=filesize(target),
            peak_rss_bytes=parse(Int64,rss.captures[1]))
end

"""Run all six direct streams and save an HDF5 timing/provenance summary.

The workers independently replay event placement, Gaussian diffusion and
physical grid assignment. Their outputs contain ion or PE pad/time histories,
with no per-carrier coordinate cache or intermediate display histogram.
"""
function sample_collect_parallel(events::Int,workers::Int,output::String,ion_card::String)
    mkpath(output)
    started=time_ns()
    channel=Channel{NamedTuple}(6)
    @sync for worker in 1:workers
        @async for index in worker:workers:6
            sample,carrier=STREAMS[index]
            put!(channel,run_stream(sample,carrier,events,output,ion_card))
        end
    end
    rows=[take!(channel) for _ in 1:6]
    sort!(rows;by=r->findfirst(==((r.sample,r.carrier)),STREAMS))
    wall=(time_ns()-started)/1e9
    summary=joinpath(output,"grid_readout_summary.h5")
    h5open(summary,"w") do file
        attributes(file)["schema"]="nausicaa.direct_grid_readout_summary.v2"
        file["events_per_stream"]=Int64(events)
        file["workers"]=Int64(workers)
        file["ion_card_path"]=ion_card
        file["wall_seconds"]=wall
        file["sample"]=[r.sample for r in rows]
        file["carrier"]=[r.carrier for r in rows]
        file["stream_seconds"]=[r.elapsed for r in rows]
        file["output_bytes"]=Int64[r.bytes for r in rows]
        file["peak_rss_bytes"]=Int64[r.peak_rss_bytes for r in rows]
        file["sum_worker_peak_rss_bytes"]=sum(r.peak_rss_bytes for r in rows)
    end
    println("Direct grid readout: $(round(wall;digits=2)) s, ",
            round(sum(r.bytes for r in rows)/1e6;digits=1)," MB")
    println("Summary: $summary")
    return summary
end

if abspath(PROGRAM_FILE)==@__FILE__
    cfg=options(ARGS)
    sample_collect_parallel(cfg.events,cfg.workers,cfg.output,cfg.ion_card)
end
