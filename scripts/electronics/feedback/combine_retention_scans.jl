root=normpath(joinpath(@__DIR__,"../../.."))
include(joinpath(root,"src/ion_signal/retention_comparison.jl"))
card=isempty(ARGS) ? "config/retention_comparison.toml" : only(ARGS)
RetentionComparison.combine(root,card)
