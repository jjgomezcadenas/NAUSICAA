#!/usr/bin/env julia
# One independent process owns every nth response shard.
include("sipm_production.jl")

"""Process disjoint shard numbers while preserving event/channel RNG identities.

`worker` is in 1:`workers`. Each worker reads the shared photon histories and
writes only its assigned sensor/ASIC files. Completed immutable shards are
verified and reused by the existing production functions.
"""
function run_worker(card,worker,workers)
    1<=worker<=workers || throw(ArgumentError("Invalid worker index"))
    P=SiPMProduction;s=P.settings(card);started=time();processed=0
    println("Worker $worker/$workers PID $(getpid()) started");flush(stdout)
    for branch in P.sources(s)
        for (number,chunk) in enumerate(Iterators.partition(branch.ids,Int(s.run["shard_events"])))
            mod(number-1,workers)==worker-1 || continue
            ids=collect(chunk);start=time()
            sensor=P.sensor_shard(s,branch,ids,number)
            P.asic_shard(s,branch,sensor,ids,number)
            processed+=length(ids)
            println("worker=$worker shard=$number gas=$(branch.src.gas) events=$(length(ids)) last_id=$(last(ids)) processed=$processed shard_s=$(round(time()-start;digits=3)) elapsed_s=$(round(time()-started;digits=3))")
            flush(stdout)
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==3 || error("Usage: RUN.toml WORKER WORKERS")
    run_worker(ARGS[1],parse(Int,ARGS[2]),parse(Int,ARGS[3]))
end
