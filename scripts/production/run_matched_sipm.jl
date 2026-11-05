#!/usr/bin/env julia
module MatchedSiPM
using TOML
include("sipm_production.jl")
const ROOT=SiPMProduction.ROOT
"""Run independent response shards, then assemble the standard SiPM manifest."""
function main(kind)
    kind in ("smoke","full") || error("Choose smoke or full")
    out=joinpath(ROOT,"results",kind=="smoke" ? "elm_matched_smoke" : "elm_matched_production")
    card="config/specs/runs/elm_matched_sipm_$kind.toml";started=time();mkpath(out)
    if kind=="full"
        jobs=Tuple{Base.Process,IO}[]
        try
            for i in 1:4
                io=open(joinpath(out,"sipm_worker_$i.log"),"w")
                script=joinpath(ROOT,"scripts/production/run_sipm_shard_worker.jl")
                process=run(pipeline(`$(Base.julia_cmd()) --startup-file=no --project=$ROOT $script $card $i 4`;stdout=io,stderr=io);wait=false)
                push!(jobs,(process,io));println("SiPM worker $i PID $(getpid(process))");flush(stdout)
            end
            for (process,_) in jobs;wait(process);success(process) || error("SiPM worker failed");end
        finally
            for (_,io) in jobs;close(io);end
        end
    end
    manifest=SiPMProduction.run(card)
    open(joinpath(out,"sipm_files.toml"),"w") do io
        TOML.print(io,Dict("complete"=>true,"manifest"=>manifest,"card"=>card,"response_wall_s"=>time()-started))
    end
    println("Matched SiPM complete in $(round(time()-started;digits=2)) s");manifest
end
end
if abspath(PROGRAM_FILE)==@__FILE__;MatchedSiPM.main(only(ARGS));end
