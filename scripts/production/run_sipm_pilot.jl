include("sipm_pipeline.jl")
root=normpath(joinpath(@__DIR__,"../.."))
SiPMPipeline.run(root,isempty(ARGS) ? "config/specs/runs/sipm_response_pilot.toml" : only(ARGS))
