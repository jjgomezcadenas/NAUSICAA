include(joinpath(@__DIR__,"../src/ion_signal/pad_history_spice.jl"))
PadHistorySpice.analyze(normpath(joinpath(@__DIR__,"..")))
PadHistorySpice.summarize(normpath(joinpath(@__DIR__,"..")))
