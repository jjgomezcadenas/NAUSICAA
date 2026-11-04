"""Run or publish the ion CSA stage using a resolved production run card."""
root=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(root,"src/ion_signal/ion_csa_pipeline.jl"))
card=isempty(ARGS) ? "config/specs/runs/ion_csa_v3_pilot_all.toml" : only(ARGS)
IonCSAProduction.run(root,card)
