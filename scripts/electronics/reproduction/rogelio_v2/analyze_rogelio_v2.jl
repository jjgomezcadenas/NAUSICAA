include(joinpath(@__DIR__,"../../../../src/ion_signal/rogelio_v2.jl"))
root=normpath(joinpath(@__DIR__,"../../../.."))
if "--diagnostics" in ARGS
    RogelioV2.diagnostics(root)
else
    RogelioV2.analyze(root;config=isempty(ARGS) ? "config/rogelio_v2.toml" : first(ARGS))
end
