# Calibrated final circuit; preserve the historical results/ltspice/averaging_scan.h5.
include(joinpath(@__DIR__,"../src/ion_signal/paper_readout.jl"))
PaperReadout.analyze(normpath(joinpath(@__DIR__,"..")))
