#!/usr/bin/env julia
include("run_electron_250ns.jl")
include("sipm_production.jl")
"""Replay photons and produce both gas responses for one 250 ns operating point.

The detector card selects E/P and J. The saved primary-pair realisation and
ion response are reused. Generated response cards bind to the completed photon
index, and the resulting manifest identifies the reusable SiPM/ASIC histories.
"""
function run_operating_point(tag)
    tag in ("ep2p5","ep3") || error("Choose ep2p5 or ep3")
    index=main_electrons("config/specs/runs/electron_250ns_$tag.toml")
    pilot=TOML.parsefile(joinpath(ROOT,"config/specs/runs/sipm_response_pilot.toml"))
    pilot["photon_index"]=relpath(index,ROOT);pilot["samples"]=["bb0nu"]
    pilotpath="config/specs/runs/sipm_response_$(tag).toml"
    open(joinpath(ROOT,pilotpath),"w") do io;TOML.print(io,pilot);end
    production=TOML.parsefile(joinpath(ROOT,"config/specs/runs/sipm_bb0nu_production.toml"))
    production["photon_index"]=relpath(index,ROOT);production["pilot_card"]=pilotpath
    production["results_directory"]="results/sipm_response_$tag"
    prodpath="config/specs/runs/sipm_bb0nu_$(tag).toml"
    open(joinpath(ROOT,prodpath),"w") do io;TOML.print(io,production);end
    manifest=SiPMProduction.run(prodpath)
    out=joinpath(ROOT,"results/electron_operating_points");mkpath(out)
    open(joinpath(out,"$tag.toml.tmp"),"w") do io
        TOML.print(io,Dict("photon_index"=>relpath(index,ROOT),"sipm_manifest"=>relpath(manifest,ROOT),"complete"=>true))
    end
    mv(joinpath(out,"$tag.toml.tmp"),joinpath(out,"$tag.toml");force=true)
    println("Operating point complete: $tag")
end
if abspath(PROGRAM_FILE)==@__FILE__
    run_operating_point(only(ARGS))
end
