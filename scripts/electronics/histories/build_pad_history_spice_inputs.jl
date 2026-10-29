get!(ENV,"NAUSICAA_PAD_HISTORY_CARD",normpath(joinpath(@__DIR__,"../../../config/pad_history_readout_direct.toml")))
include(joinpath(@__DIR__,"../../../src/ion_signal/pad_history_readout.jl"))
PadHistoryReadout.spice_inputs(normpath(joinpath(@__DIR__,"../../..")))
