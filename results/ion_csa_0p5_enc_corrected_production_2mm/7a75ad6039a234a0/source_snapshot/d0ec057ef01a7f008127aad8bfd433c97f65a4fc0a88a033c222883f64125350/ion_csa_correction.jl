root=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(root,"src/ion_signal/ion_csa_correction.jl"))
const C=IonCSACorrection
card=isempty(ARGS) ? "config/specs/runs/ion_csa_v3_corrected_pilot.toml" : ARGS[1]
C.correct(C.setup(root,card))
