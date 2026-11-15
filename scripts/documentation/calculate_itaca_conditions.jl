#!/usr/bin/env julia
#=
Resolve the ITACA detector set and save its derived quantities to HDF5.

  julia --project=. scripts/documentation/calculate_itaca_conditions.jl            # default detector
  julia --project=. scripts/documentation/calculate_itaca_conditions.jl <run card> # with the run's overrides

Default output: results/itaca_detector/itaca_conditions.h5. With a run card the
output is results/readout/<run_label>/configuration/itaca_conditions.h5. The file
holds the card texts and SHA-256, every detector value under /inputs (after
overrides), the overrides under /overrides, and every derived value under /derived
(Float64 scalar with unit and note attributes). Python makes the page from it
(scripts/documentation/present_itaca_detector.py).
=#
using HDF5, TOML
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
using .EventAnalysis, .ItacaConditions

index = joinpath(ROOT, "config/detector/itaca.toml")
if isempty(ARGS)
    set = load_detector_set(index); detector = set.cards; overrides = Dict{String,Any}(); run_path = ""
    out = joinpath(ROOT, "results/itaca_detector/itaca_conditions.h5")
else
    r = resolve_run(abspath(ARGS[1]); root = ROOT)
    set = r.set; detector = r.detector; overrides = r.overrides; run_path = relpath(r.run_path, ROOT)
    out = joinpath(ROOT, "results/readout", r.run["run_label"], "configuration/itaca_conditions.h5")
end
derived = derived_conditions(detector)

mkpath(dirname(out))
h5open(out, "w") do f
    attrs(f)["detector_index"] = relpath(set.index_path, ROOT)
    attrs(f)["detector_name"] = set.name
    attrs(f)["detector_date"] = set.date
    attrs(f)["run_card"] = run_path
    attrs(f)["generator"] = "scripts/documentation/calculate_itaca_conditions.jl"
    cards = create_group(f, "cards")
    for (name, text) in set.texts
        cards[name] = text
        attrs(cards[name])["sha256"] = set.sha256[name]
    end
    inputs = create_group(f, "inputs")
    for (key, value) in ItacaConditions.flatten(detector)
        if value isa Dict
            inputs[key] = sprint(io -> TOML.print(io, value))
        elseif value isa Real && !(value isa Bool)
            inputs[key] = Float64(value)
        else
            inputs[key] = string(value)
        end
    end
    ov = create_group(f, "overrides")
    for (key, value) in overrides
        ov[key] = value isa Real && !(value isa Bool) ? Float64(value) : string(value)
    end
    g = create_group(f, "derived")
    for (group, name, value, unit, note) in derived
        ds = "$group.$name"
        g[ds] = value
        attrs(g[ds])["unit"] = unit
        attrs(g[ds])["note"] = note
    end
end
println("Saved $(length(derived)) derived values to $(relpath(out, ROOT)); overrides: $(length(overrides))")
for (group, name, value, unit, note) in derived
    println(rpad("$group.$name", 52), " ", rpad(string(round(value, sigdigits=6)), 14), " ", unit)
end
