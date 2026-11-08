#!/usr/bin/env julia
"""Run the ELM 100-event photon, finite-cell SiPM and ASIC pilot."""
module ELMSiPMPilot
using HDF5,TOML,SHA
module Replay
include("run_electron_250ns.jl")
end
include("sipm_production.jl")
include(joinpath(@__DIR__,"../../src/ion_signal/guarded_files.jl"))
const ROOT=Replay.ROOT
const OUT=joinpath(ROOT,"results/elm_sipm_pilot")

"""Archive blue photon histories beside the shared primary realisation.

The new index preserves the replay's source/cards and records immutable archived
paths. Existing archives must match the source contents exactly. Ion files are
reused through the detector configuration and are never opened for writing.
"""
function archive_photons(source_index)
    mkpath(OUT)
    target=joinpath(OUT,"photons_ep2.h5")
    GuardedFiles.checked_copy(source_index,target*".tmp";force=true)
    h5open(target*".tmp","r+") do f
        g=f["bb0nu"]
        base=dirname(String(read(g["realisation_path"])))
        dir=joinpath(base,"detector","photons",String(read(g["electron_id"])));mkpath(dir)
        for gas in ("pure_xe","xe_he15")
            key="$(gas)_path";source=String(read(g[key]));dest=joinpath(dir,basename(source))
            if isfile(dest)
                sha256(read(dest))==sha256(read(source)) || error("Archived photon contents differ")
            else
                GuardedFiles.checked_copy(source,dest)
            end
            delete_object(g,key);g[key]=dest
        end
        f["metadata/replay_index"]=String(source_index)
        f["metadata/photon_band"]="blue_tpb"
    end
    mv(target*".tmp",target;force=true)
    target
end

"""Generate/reuse matched first-100 blue histories and independently saved response shards."""
function main()
    mkpath(OUT)
    started=time()
    index=Replay.main_electrons("config/specs/runs/elm_photons_ep2_drift10cm.toml",100)
    archive_photons(index)
    photon_wall=time()-started
    manifests=Dict{String,Any}()
    response_started=time()
    for tag in ("12p5","25")
        start=time()
        card="config/specs/runs/elm_bb0nu_tau$(tag)ns_pilot100.toml"
        manifest=SiPMProduction.run(card)
        manifests[tag]=Dict("manifest"=>manifest,"card"=>card,"elapsed_s"=>time()-start)
        println("ELM recovery $tag ns complete: $(round(time()-start;digits=2)) s")
    end
    response_wall=time()-response_started
    open(joinpath(OUT,"files.toml"),"w") do io
        TOML.print(io,Dict("complete"=>true,"photon_index"=>joinpath(OUT,"photons_ep2.h5"),
            "responses"=>manifests,"photon_wall_s"=>photon_wall,
            "response_wall_s"=>response_wall,"total_wall_s"=>time()-started,
            "events_per_gas"=>100,"sample"=>"bb0nu"))
    end
    println("ELM photons: $(round(photon_wall;digits=3)) s; response: $(round(response_wall;digits=3)) s")
end
end
if abspath(PROGRAM_FILE)==@__FILE__;ELMSiPMPilot.main();end
