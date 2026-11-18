#!/usr/bin/env julia
#=
Radiogenic background: NEXUS macros of the ICS runs, written from the shielding cards.

  julia --project=. scripts/background/write_background_macros.jl                 # every card of config/background
  julia --project=. scripts/background/write_background_macros.jl <card> ...     # some cards

Writes simulation/nexus/macros/background/<name>/: the geometry macro
ItacaChamber.mac (fiducial gas of the variant, cavity, layers and copper) and,
for each isotope (bi214, tl208), the templates ItacaChamber.<isotope>.config.mac
and ItacaChamber.<isotope>.init.mac, in which the job script replaces @SEED@,
@START_ID@, @OUTPUT@ and @CONFIG@ (simulation/nexus/slurm/run_iitaca_nexus_batch.sh
with MACRO_SET=background), and the delayed macro ItacaChamber.<isotope>.decay.mac
that limits the radioactive decay to the generated nucleus. The generator is the
IonGenerator in region ICS: full decays, every daughter followed; an event is
saved when the energy deposited in the fiducial gas lies in the saved window of
the generation card. The macros are committed; tests/background/test_background.jl
requires them to equal what the cards give.
=#
include(joinpath(@__DIR__, "resolve_background.jl"))
include(joinpath(ROOT, "scripts/iitaca/write_iitaca_macros.jl"))     # header, number

const BACKGROUND_MACRO_ROOT = "simulation/nexus/macros/background"
const ISOTOPES = Dict("bi214" => (Z = 83, A = 214, name = "214Bi"), "tl208" => (Z = 81, A = 208, name = "208Tl"))

"""Geometry macro: the gas cylinder of the variant and the shield of the card."""
function shield_geometry_macro(r)
    sim = r.sim; c = r.card.card; ics = c["ics"]
    face(f) = "/Geometry/ItacaChamber/$(f)_gap $(number(c[f]["gap_cm"])) cm\n/Geometry/ItacaChamber/$(f)_layer_thickness $(number(c[f]["layer_thickness_cm"])) cm\n/Geometry/ItacaChamber/$(f)_layer_material $(c[f]["layer_material"])\n"
    header("ItacaChamber.mac", c["name"], ["Fiducial gas of $(r.variant): radius $(number(sim["geometry.fiducial_diameter_cm"] / 2)) cm, length $(number(sim["geometry.fiducial_length_cm"])) cm,",
        "enriched xenon with a helium mole fraction of $(sim["gas.helium_mole_fraction"]) at $(number(sim["gas.pressure_bar"])) bar and $(number(sim["gas.temperature_K"])) K.",
        "Shield: $(number(ics["thickness_cm"])) cm of $(ics["material"]) around a xenon cavity with one passive layer per face;",
        "region ICS samples the decay depth as exp(-depth/$(number(ics["sampling_length_cm"])) cm)."]; family = "background") * """

/Geometry/ItacaChamber/radius $(number(sim["geometry.fiducial_diameter_cm"] / 2)) cm
/Geometry/ItacaChamber/length $(number(sim["geometry.fiducial_length_cm"])) cm
/Geometry/ItacaChamber/pressure $(number(sim["gas.pressure_bar"])) bar
/Geometry/ItacaChamber/temperature $(number(sim["gas.temperature_K"])) kelvin
/Geometry/ItacaChamber/helium_mole_fraction $(sim["gas.helium_mole_fraction"])
/Geometry/ItacaChamber/step_limit $(number(sim["generation.step_limit_mm"])) mm

/Geometry/ItacaChamber/shield_material $(ics["material"])
/Geometry/ItacaChamber/shield_thickness $(number(ics["thickness_cm"])) cm
/Geometry/ItacaChamber/shield_sampling_length $(number(ics["sampling_length_cm"])) cm
$(face("barrel"))$(face("anode"))$(face("cathode"))"""
end

"""Saved energy window of the card, in keV, from the generation card of the detector."""
saved_window(r) = r.card.card["saved_window"] == "gamma" ?
    (min = r.sim["generation.gamma_saved_energy_min_keV"], max = r.sim["generation.gamma_saved_energy_max_keV"]) :
    (min = r.sim["generation.saved_energy_min_keV"], max = r.sim["generation.saved_energy_max_keV"])

function shield_config_macro(r, isotope)
    iso = ISOTOPES[isotope]; name = r.card.card["name"]; w = saved_window(r)
    header("ItacaChamber.$isotope.config.mac", name, ["Template. The job script replaces @SEED@, @START_ID@ and @OUTPUT@.",
        "$(iso.name) decays uniform in the copper shield (region ICS, exponential depth sampling with weights),",
        "every decay product followed. Events are saved when the energy deposited in the fiducial gas lies in the window below."]; family = "background") * """

/run/verbose 1
/event/verbose 0
/tracking/verbose 0
/process/em/verbose 0

/control/execute macros/background/$name/ItacaChamber.mac

/Generator/IonGenerator/atomic_number $(iso.Z)
/Generator/IonGenerator/mass_number $(iso.A)
/Generator/IonGenerator/region ICS
/Generator/IonGenerator/decay_at_time_zero true

/Actions/DefaultEventAction/min_energy $(number(w.min)) keV
/Actions/DefaultEventAction/max_energy $(number(w.max)) keV

/nexus/random_seed @SEED@
/nexus/persistency/start_id @START_ID@
/nexus/persistency/output_file @OUTPUT@
/nexus/persistency/event_type background
"""
end

function shield_init_macro(r, isotope)
    name = r.card.card["name"]
    header("ItacaChamber.$isotope.init.mac", name, ["Template. The job script replaces @CONFIG@ with the configuration macro of the job."]; family = "background") * """

/control/execute macros/physics/DefaultPhysicsList.mac

/nexus/RegisterGeometry ItacaChamber
/nexus/RegisterGenerator IonGenerator

/nexus/RegisterPersistencyManager PersistencyManager

/nexus/RegisterRunAction DefaultRunAction
/nexus/RegisterEventAction DefaultEventAction
/nexus/RegisterTrackingAction DefaultTrackingAction

/nexus/RegisterMacro @CONFIG@
/nexus/RegisterDelayedMacro macros/background/$name/ItacaChamber.$isotope.decay.mac
"""
end

"""Delayed macro: only the generated nucleus decays (its daughters are stable for Geant4)."""
function shield_decay_macro(r, isotope)
    iso = ISOTOPES[isotope]
    header("ItacaChamber.$isotope.decay.mac", r.card.card["name"], ["Radioactive decay limited to $(iso.name): the daughter nucleus does not decay."]; family = "background") * """

/process/had/rdm/nucleusLimits $(iso.A) $(iso.A) $(iso.Z) $(iso.Z)
"""
end

"""All macros of a shielding card: file name => text."""
function shield_macros(card_path::AbstractString)
    r = resolved_shield(card_path)
    files = Dict{String,String}("ItacaChamber.mac" => shield_geometry_macro(r))
    for iso in sort(collect(keys(ISOTOPES)))
        files["ItacaChamber.$iso.config.mac"] = shield_config_macro(r, iso)
        files["ItacaChamber.$iso.init.mac"] = shield_init_macro(r, iso)
        files["ItacaChamber.$iso.decay.mac"] = shield_decay_macro(r, iso)
    end
    r.card.card["name"], files
end

if abspath(PROGRAM_FILE) == @__FILE__
    for card in (isempty(ARGS) ? shield_cards() : abspath.(ARGS))
        name, files = shield_macros(card); dir = joinpath(ROOT, BACKGROUND_MACRO_ROOT, name); mkpath(dir)
        for (file, text) in files; write(joinpath(dir, file), text); end
        println("$(length(files)) macros in $(relpath(dir, ROOT))")
    end
end
