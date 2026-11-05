#!/usr/bin/env julia
"""Run the ELM 100-event photon, linear-SiPM and ASIC clipping pilot."""
module ELMPilot
using HDF5,TOML,SHA
module Replay
include("run_electron_250ns.jl")
end
include("sipm_production.jl")
const ROOT=Replay.ROOT
const OUT=joinpath(ROOT,"results/elm_pilot")

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
        for gas in ("pure_xe","xe_he15")
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

"""Generate/reuse matched first-100 blue histories and independently saved response shards."""
function main()
    mkpath(OUT)
    started=time()
    index=Replay.main_electrons("config/specs/runs/elm_photons_ep2.toml",100)
    archive_photons(index)
    photon_wall=time()-started
    response_started=time()
    manifest=SiPMProduction.run("config/specs/runs/elm_bb0nu_ep2_pilot100.toml")
    response_wall=time()-response_started
    open(joinpath(OUT,"files.toml"),"w") do io
        TOML.print(io,Dict("complete"=>true,"photon_index"=>joinpath(OUT,"photons_ep2.h5"),
            "response_manifest"=>manifest,"photon_wall_s"=>photon_wall,
            "response_wall_s"=>response_wall,"total_wall_s"=>time()-started,
            "events_per_gas"=>100,"sample"=>"bb0nu"))
    end
    println("ELM photons: $(round(photon_wall;digits=3)) s; response: $(round(response_wall;digits=3)) s")
end
end
if abspath(PROGRAM_FILE)==@__FILE__;ELMPilot.main();end
