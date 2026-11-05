#!/usr/bin/env julia
# Reduce the Xe–He 4 mm pitch pilot and compare matched saved 5 mm events.
include("reduce_elm_sipm_pilot.jl")
if abspath(PROGRAM_FILE)==@__FILE__
    ELMSiPMReduction.main(output_dir=joinpath(ELMSiPMReduction.ROOT,"results/elm_sipm_pitch4mm_pilot"),
        report_name="elm_sipm_pitch4mm_pilot",
        reference_results=joinpath(ELMSiPMReduction.ROOT,"results/elm_sipm_pilot/elm_sipm_pilot.h5"))
end
