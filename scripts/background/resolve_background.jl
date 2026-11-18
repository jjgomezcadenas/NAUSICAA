#=
Shared by the background scripts: load a shielding card, resolve the detector
variant it names (iItaca or eITACA) and build the shield geometry.

  include("resolve_background.jl"); r = resolved_shield("config/background/ics_iitaca.toml")

`resolved_shield` returns a named tuple: `card` (ShieldCard), `geometry`
(ShieldGeometry), `sim` (flat simulation inputs of the detector: gas,
geometry, generation), `variant` (name) and `family`.
=#
const ROOT = normpath(joinpath(@__DIR__, "../.."))
if !isdefined(@__MODULE__, :IItaca)
    include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
    include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
    include(joinpath(ROOT, "src/iitaca/IItaca.jl"))
end
isdefined(@__MODULE__, :EItaca) || include(joinpath(ROOT, "src/eitaca/EItaca.jl"))
isdefined(@__MODULE__, :Background) || include(joinpath(ROOT, "src/background/Background.jl"))
using .Background

"""Shielding card, resolved detector inputs and shield geometry of a card path."""
function resolved_shield(card_path::AbstractString)
    card = load_shielding(card_path); c = card.card; variant = joinpath(ROOT, c["variant"])
    isfile(variant) || error("Variant card $variant of $card_path not found")
    if c["family"] == "iitaca"
        r = IItaca.resolve_detector(joinpath(ROOT, "config/iitaca/iitaca.toml"); variant, root = ROOT)
        sim = IItaca.simulation_inputs(r.detector)
    else
        r = EItaca.resolve_detector(joinpath(ROOT, "config/eitaca/eitaca.toml"); variant, root = ROOT)
        sim = ItacaConditions.flatten(r.detector)
    end
    geometry = shield_geometry(card, sim["geometry.fiducial_diameter_cm"] / 2, sim["geometry.fiducial_length_cm"])
    (card = card, geometry = geometry, sim = sim, variant = r.variant, family = c["family"])
end

shield_cards() = sort(filter(p -> endswith(p, ".toml"), readdir(joinpath(ROOT, "config/background"); join = true)))
