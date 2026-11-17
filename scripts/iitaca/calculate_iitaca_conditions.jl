#!/usr/bin/env julia
#=
Resolve the iItaca detector set with each of its variants and save the derived quantities to HDF5.

  julia --project=. scripts/iitaca/calculate_iitaca_conditions.jl

Output: results/iitaca/detector/conditions.h5. It holds the card texts and
SHA-256 under /cards and, for every variant card of config/iitaca/variants/
(in name order), a group /variants/<name> with its card text, every detector
value after the overrides under inputs, the overrides, and every derived value
under derived (Float64 scalar with unit and note attributes). Python makes the
page from it (scripts/iitaca/present_iitaca_detector.py).
=#
using HDF5, TOML
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
include(joinpath(ROOT, "src/iitaca/IItaca.jl"))
using .IItaca

scalar(value) = value isa Dict ? sprint(io -> TOML.print(io, value)) : value isa Real && !(value isa Bool) ? Float64(value) : string(value)

index = joinpath(ROOT, "config/iitaca/iitaca.toml")
paths = sort(filter(p -> endswith(p, ".toml"), readdir(joinpath(ROOT, "config/iitaca/variants"); join = true)))
resolved = [resolve_detector(index; variant = p, root = ROOT) for p in paths]
out = joinpath(ROOT, "results/iitaca/detector/conditions.h5")
mkpath(dirname(out))
h5open(out, "w") do f
    set = first(resolved).set
    attrs(f)["schema"] = "nausicaa.iitaca_conditions.v1"
    attrs(f)["detector_index"] = relpath(set.index_path, ROOT)
    attrs(f)["detector_name"] = set.name
    attrs(f)["detector_date"] = set.date
    attrs(f)["generator"] = "scripts/iitaca/calculate_iitaca_conditions.jl"
    cards = create_group(f, "cards")
    for (name, text) in set.texts
        cards[name] = text
        attrs(cards[name])["sha256"] = set.sha256[name]
    end
    variants = create_group(f, "variants")
    for (order, r) in enumerate(resolved)
        g = create_group(variants, r.variant)
        attrs(g)["order"] = order
        attrs(g)["card"] = relpath(r.variant_path, ROOT)
        g["card_text"] = r.variant_text
        inputs = create_group(g, "inputs")
        for (key, value) in ItacaConditions.flatten(r.detector)
            inputs[key] = scalar(value)
        end
        ov = create_group(g, "overrides")
        for (key, value) in r.overrides
            ov[key] = scalar(value)
        end
        derived = create_group(g, "derived")
        for (group, name, value, unit, note) in derived_conditions(r.detector)
            derived["$group.$name"] = value
            attrs(derived["$group.$name"])["unit"] = unit
            attrs(derived["$group.$name"])["note"] = note
        end
    end
end
println("Saved $(length(resolved)) variants to $(relpath(out, ROOT))")
for r in resolved
    println("\n", r.variant, "  overrides: ", length(r.overrides))
    for (group, name, value, unit, note) in derived_conditions(r.detector)
        println(rpad("$group.$name", 46), " ", rpad(string(round(value, sigdigits = 6)), 14), " ", unit)
    end
end
