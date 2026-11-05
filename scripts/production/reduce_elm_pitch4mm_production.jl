#!/usr/bin/env julia
# Reduce the full Xe–He production using streaming quiet-noise statistics.
include("reduce_elm_sipm_pilot.jl")
if abspath(PROGRAM_FILE)==@__FILE__
    ELMSiPMReduction.main(output_dir=joinpath(ELMSiPMReduction.ROOT,"results/elm_sipm_pitch4mm_production"),
        report_name="elm_sipm_pitch4mm_production",store_inclusive=false)
end
