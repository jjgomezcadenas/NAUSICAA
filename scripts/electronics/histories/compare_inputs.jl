root=normpath(joinpath(@__DIR__,"../../.."))
include(joinpath(root,"src/ion_signal/direct_history_comparison.jl"))
DirectHistoryComparison.run(root)
