include(joinpath(@__DIR__,"../src/ion_signal/feedback_optimization.jl"))
root=normpath(joinpath(@__DIR__,".."))
FeedbackOptimization.analyze(root;names=isempty(ARGS) ? nothing : ARGS)
FeedbackOptimization.checks(root)
