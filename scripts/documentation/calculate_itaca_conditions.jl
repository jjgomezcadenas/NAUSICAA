#!/usr/bin/env julia
#=
Derive the ITACA operating quantities from config/itaca.toml and save them to HDF5.

Run from the repository with
  JULIA_DEPOT_PATH=work/julia-depot julia --project=. scripts/documentation/calculate_itaca_conditions.jl
or `make itaca-detector`. Output: results/itaca_detector/itaca_conditions.h5 with
the card text and SHA-256, every card value under /inputs and every derived value
under /derived (Float64 scalar with unit and note attributes). Python makes the
page from this file (scripts/documentation/present_itaca_detector.py).
=#
using HDF5
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
using .EventAnalysis, .ItacaConditions

card_path = joinpath(ROOT, "config/itaca.toml")
card = load_itaca_card(card_path)
derived = derived_conditions(card)

out = joinpath(ROOT, "results/itaca_detector/itaca_conditions.h5")
mkpath(dirname(out))
h5open(out, "w") do f
    attrs(f)["card_path"] = relpath(card_path, ROOT)
    attrs(f)["card_sha256"] = card.sha256
    attrs(f)["card_date"] = string(card.data["date"])
    attrs(f)["schema"] = string(card.data["schema"])
    attrs(f)["generator"] = "scripts/documentation/calculate_itaca_conditions.jl"
    f["card_toml"] = card.text
    inputs = create_group(f, "inputs")
    function store!(group, prefix, table)
        for (k, v) in table
            name = prefix == "" ? k : "$prefix.$k"
            if v isa Dict
                store!(group, name, v)
            elseif v isa Real
                group[name] = Float64(v)
            else
                group[name] = string(v)
            end
        end
    end
    for (section, table) in card.data
        table isa Dict || continue
        store!(inputs, section, table)
    end
    g = create_group(f, "derived")
    for (group, name, value, unit, note) in derived
        ds = "$group.$name"
        g[ds] = value
        attrs(g[ds])["unit"] = unit
        attrs(g[ds])["note"] = note
    end
end
println("Saved $(length(derived)) derived values to $(relpath(out, ROOT)) (card sha256 $(card.sha256[1:12]))")
for (group, name, value, unit, note) in derived
    println(rpad("$group.$name", 52), " ", rpad(string(round(value, sigdigits=6)), 14), " ", unit)
end
