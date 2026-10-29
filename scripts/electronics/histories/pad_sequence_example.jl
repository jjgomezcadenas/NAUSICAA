include(joinpath(@__DIR__,"../../../src/ion_signal/pad_sequence_example.jl"))
root=normpath(joinpath(@__DIR__,"../../.."))
if "--build" in ARGS
 PadSequenceExample.build(root)
else
 PadSequenceExample.analyze(root)
end
