#!/usr/bin/env julia
"""Run full Xe–He photon replay and four independent response shard workers."""
module ELMFullProduction
using HDF5,TOML,SHA
module Replay
include("run_electron_250ns.jl")
end
include("sipm_production.jl")
const ROOT=Replay.ROOT
const OUT=joinpath(ROOT,"results/elm_sipm_pitch4mm_production")

"""Archive blue photon histories beside the shared primary realisation.

The new index preserves the replay's source/cards and records immutable archived
paths. Existing archives must match the source contents exactly. Ion files are
reused through the detector configuration and are never opened for writing.
"""
function archive_photons(source_index)
    mkpath(OUT)
    target=joinpath(OUT,"photons_ep2.h5")
    cp(source_index,target*".tmp";force=true)
    h5open(target*".tmp","r+") do f
        g=f["bb0nu"]
        base=dirname(String(read(g["realisation_path"])))
        dir=joinpath(base,"detector","photons",String(read(g["electron_id"])));mkpath(dir)
        for gas in ("xe_he15",)
            key="$(gas)_path";source=String(read(g[key]));dest=joinpath(dir,basename(source))
            if isfile(dest)
                sha256(read(dest))==sha256(read(source)) || error("Archived photon contents differ")
            else
                cp(source,dest)
            end
            delete_object(g,key);g[key]=dest
        end
        f["metadata/replay_index"]=String(source_index)
        f["metadata/photon_band"]="blue_tpb"
    end
    mv(target*".tmp",target;force=true)
    target
end

"""Run disjoint restartable shards in four independent Julia processes.

Workers share read-only photon input and write distinct shard files. Event and
channel seeds are independent of worker assignment; the standard production
runner validates completed shards and produces the aggregate manifest afterward.
"""
function run_workers(card;workers=4)
    jobs=Tuple{Base.Process,IO}[]
    try
        for i in 1:workers
            logfile=joinpath(OUT,"worker_$(i).log")
            io=open(logfile,"w")
            script=joinpath(ROOT,"scripts/production/run_sipm_shard_worker.jl")
            cmd=`$(Base.julia_cmd()) --startup-file=no --project=$ROOT $script $card $i $workers`
            process=run(pipeline(cmd;stdout=io,stderr=io);wait=false)
            push!(jobs,(process,io))
            println("Response worker $i PID $(getpid(process)); $logfile");flush(stdout)
        end
        for (process,_) in jobs
            wait(process);success(process) || error("Response worker failed; inspect worker logs")
        end
    finally
        for (_,io) in jobs;close(io);end
    end
end

"""Generate the full selected photon sample and baseline response."""
function main()
    mkpath(OUT)
    started=time()
    println("ELM full production PID $(getpid()), start $(time())");flush(stdout)
    index=Replay.main_electrons("config/specs/runs/elm_photons_pitch4mm_ep2_full.toml",8000)
    archive_photons(index)
    photon_wall=time()-started
    manifests=Dict{String,Any}()
    response_started=time()
    for tag in ("12p5",)
        start=time()
        card="config/specs/runs/elm_pitch4mm_bb0nu_full.toml"
        run_workers(card)
        manifest=SiPMProduction.run(card)
        manifests[tag]=Dict("manifest"=>manifest,"card"=>card,"elapsed_s"=>time()-start)
        println("ELM recovery $tag ns complete: $(round(time()-start;digits=2)) s")
    end
    response_wall=time()-response_started
    open(joinpath(OUT,"files.toml"),"w") do io
        TOML.print(io,Dict("complete"=>true,"photon_index"=>joinpath(OUT,"photons_ep2.h5"),
            "responses"=>manifests,"photon_wall_s"=>photon_wall,
            "response_wall_s"=>response_wall,"total_wall_s"=>time()-started,
            "events_per_gas"=>8000,"sample"=>"bb0nu"))
    end
    println("ELM photons: $(round(photon_wall;digits=3)) s; response: $(round(response_wall;digits=3)) s")
end
end
if abspath(PROGRAM_FILE)==@__FILE__;ELMFullProduction.main();end
