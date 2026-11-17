#!/usr/bin/env julia
#=
eITACA: NEXUS macros of a variant, written from the detector cards.

  julia --project=. scripts/eitaca/write_eitaca_macros.jl                    # every variant of config/eitaca/variants
  julia --project=. scripts/eitaca/write_eitaca_macros.jl <variant card>     # one variant

Writes simulation/nexus/macros/eitaca/<variant>/ with the same macros as the
iItaca writer (geometry; bb0nu, electron, bi214, tl208 templates), from the
resolved eITACA detector: gas (xenon with its helium fraction), pressure,
temperature, radius and length, and the generation card. The ItacaChamber
geometry places z = 0 at one end and z = length at the other; for eITACA the
cathode is at z = 0 and the gate at z = length.
=#
include(joinpath(@__DIR__, "../iitaca/write_iitaca_macros.jl"))      # geometry_macro, config_macro, init_macro, SAMPLES
isdefined(@__MODULE__, :EItaca) || include(joinpath(ROOT, "src/eitaca/EItaca.jl"))

const EITACA_MACRO_ROOT = "simulation/nexus/macros/eitaca"

"""All macros of an eITACA variant: file name => text."""
function eitaca_macros(variant_card::AbstractString)
    r = EItaca.resolve_detector(joinpath(ROOT, "config/eitaca/eitaca.toml"); variant = variant_card, root = ROOT)
    sim = ItacaConditions.flatten(r.detector)
    files = Dict{String,String}("ItacaChamber.mac" => geometry_macro(sim, r.variant; family = "eitaca"))
    for s in SAMPLES
        files["ItacaChamber.$s.config.mac"] = config_macro(sim, r.variant, s; family = "eitaca")
        files["ItacaChamber.$s.init.mac"] = init_macro(sim, r.variant, s; family = "eitaca")
    end
    r.variant, files
end

eitaca_variant_cards() = sort(filter(p -> endswith(p, ".toml"), readdir(joinpath(ROOT, "config/eitaca/variants"); join = true)))

if abspath(PROGRAM_FILE) == @__FILE__
    for card in (isempty(ARGS) ? eitaca_variant_cards() : abspath.(ARGS))
        variant, files = eitaca_macros(card); dir = joinpath(ROOT, EITACA_MACRO_ROOT, variant); mkpath(dir)
        for (name, text) in files; write(joinpath(dir, name), text); end
        println("$(length(files)) macros in $(relpath(dir, ROOT))")
    end
end
