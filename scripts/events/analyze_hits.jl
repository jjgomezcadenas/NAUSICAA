#!/usr/bin/env julia
include(joinpath(@__DIR__,"../../src/ion_signal/event_analysis.jl"))
using .EventAnalysis

if "--help" in ARGS || "-h" in ARGS
    println("Usage: julia --project=. scripts/events/analyze_hits.jl [config/ion_tracks.toml]")
    println("Analyze saved hit DataFrames in Julia; write results/ion_tracks/<sample>/analysis.h5.")
elseif length(ARGS)>1
    error("Expected at most one configuration path")
else
    isempty(ARGS) ? run_analysis() : run_analysis(abspath(only(ARGS)))
end
