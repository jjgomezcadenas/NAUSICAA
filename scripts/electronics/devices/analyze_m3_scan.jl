include(joinpath(@__DIR__,"../../../src/ion_signal/m4_scan.jl"))
M4Scan.analyze(normpath(joinpath(@__DIR__,"../../..")); study="m3_scan",
    reference="m4_scan",reference_group="scale_4_gate_0.5568034172",gate_node="V(vnload)")
