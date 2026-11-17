#!/usr/bin/env julia
#=
iItaca: NEXUS macros of a variant, written from the detector cards.

  julia --project=. scripts/iitaca/write_iitaca_macros.jl                    # every variant of config/iitaca/variants
  julia --project=. scripts/iitaca/write_iitaca_macros.jl <variant card>     # one variant

Writes simulation/nexus/macros/iitaca/<variant>/: the geometry macro
ItacaChamber.mac and, for each sample (bb0nu, electron, bi214, tl208), the
templates ItacaChamber.<sample>.config.mac and ItacaChamber.<sample>.init.mac,
in which the job script replaces @SEED@, @START_ID@, @OUTPUT@ and @CONFIG@
(simulation/nexus/slurm/run_iitaca_nexus_batch.sh). Gas, pressure, temperature,
radius and length come from the resolved detector; step limit, saved energy
window, energies and regions from its generation card. The macros are
committed; tests/iitaca/test_iitaca_macros.jl requires them to equal what the
cards give.
=#
const ROOT = normpath(joinpath(@__DIR__, "../.."))
if !isdefined(@__MODULE__, :IItaca)
    include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
    include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
    include(joinpath(ROOT, "src/iitaca/IItaca.jl")); using .IItaca
end

const SAMPLES = ("bb0nu", "electron", "bi214", "tl208")
const MACRO_ROOT = "simulation/nexus/macros/iitaca"
number(x) = (s = string(Float64(x)); endswith(s, ".0") ? s[1:end-1] : s)       # 200.0 -> "200."

header(name, variant, lines) = join(vcat("## " * "-"^76, "## nexus | $name", "##", ["## " * l for l in lines]...,
    "## Written by scripts/iitaca/write_iitaca_macros.jl from the cards of config/iitaca", "## (variant $variant). Do not edit: change the cards and write it again.", "## " * "-"^76), "\n") * "\n"

"""Geometry macro of the gas cylinder from the simulation inputs."""
geometry_macro(sim, variant) = header("ItacaChamber.mac", variant, ["iItaca gas volume: radius $(number(sim["geometry.fiducial_diameter_cm"] / 2)) cm, length $(number(sim["geometry.fiducial_length_cm"])) cm,",
        "enriched xenon with a helium mole fraction of $(sim["gas.helium_mole_fraction"]) at $(number(sim["gas.pressure_bar"])) bar and $(number(sim["gas.temperature_K"])) K."]) * """

/Geometry/ItacaChamber/radius $(number(sim["geometry.fiducial_diameter_cm"] / 2)) cm
/Geometry/ItacaChamber/length $(number(sim["geometry.fiducial_length_cm"])) cm
/Geometry/ItacaChamber/pressure $(number(sim["gas.pressure_bar"])) bar
/Geometry/ItacaChamber/temperature $(number(sim["gas.temperature_K"])) kelvin
/Geometry/ItacaChamber/helium_mole_fraction $(sim["gas.helium_mole_fraction"])
/Geometry/ItacaChamber/step_limit $(number(sim["generation.step_limit_mm"])) mm
"""

"""Generator lines and description of a sample."""
function generator(sim, sample)
    gamma(E) = ("/Generator/SingleParticle/particle gamma\n/Generator/SingleParticle/min_energy $E keV\n/Generator/SingleParticle/max_energy $E keV\n" *
                "/Generator/SingleParticle/region $(sim["generation.gamma_region"])", "SingleParticleGenerator", "background")
    if sample == "bb0nu"
        ("/Generator/Decay0Interface/inputFile none\n/Generator/Decay0Interface/Xe136DecayMode 1\n/Generator/Decay0Interface/Ba136FinalState 0\n" *
         "/Generator/Decay0Interface/region $(sim["generation.volume_region"])", "Decay0Interface", "signal",
         ["Neutrinoless double beta decays of Xe-136 to the ground state, uniform in the gas."])
    elseif sample == "electron"
        E = sim["generation.single_electron_energy_keV"]
        ("/Generator/SingleParticle/particle e-\n/Generator/SingleParticle/min_energy $E keV\n/Generator/SingleParticle/max_energy $E keV\n" *
         "/Generator/SingleParticle/region $(sim["generation.volume_region"])", "SingleParticleGenerator", "background",
         ["Single electrons of $E keV, isotropic, uniform in the gas."])
    elseif sample == "bi214"
        (gamma(sim["generation.bi214_gamma_keV"])..., ["Gammas of the $(sim["generation.bi214_gamma_keV"]) keV line of 214Bi, emitted isotropically from the inner", "surface of the gas cylinder."])
    elseif sample == "tl208"
        (gamma(sim["generation.tl208_gamma_keV"])..., ["Gammas of the $(sim["generation.tl208_gamma_keV"]) keV line of 208Tl, emitted isotropically from the inner", "surface of the gas cylinder."])
    else
        throw(ArgumentError("Unknown sample $sample"))
    end
end

"""Configuration template of a sample: geometry, generator, saved energy window, seed and output."""
function config_macro(sim, variant, sample)
    lines, _, event_type, description = generator(sim, sample)
    header("ItacaChamber.$sample.config.mac", variant, vcat("Template. The job script replaces @SEED@, @START_ID@ and @OUTPUT@.", description,
        "Events are saved when the energy deposited in the gas lies in the window below.")) * """

/run/verbose 1
/event/verbose 0
/tracking/verbose 0
/process/em/verbose 0

/control/execute macros/iitaca/$variant/ItacaChamber.mac

$lines

/Actions/DefaultEventAction/min_energy $(number(sim["generation.saved_energy_min_keV"])) keV
/Actions/DefaultEventAction/max_energy $(number(sim["generation.saved_energy_max_keV"])) keV

/nexus/random_seed @SEED@
/nexus/persistency/start_id @START_ID@
/nexus/persistency/output_file @OUTPUT@
/nexus/persistency/event_type $event_type
"""
end

"""Initialisation template of a sample: physics list, geometry, generator and actions."""
function init_macro(sim, variant, sample)
    _, name, _, _ = generator(sim, sample)
    header("ItacaChamber.$sample.init.mac", variant, ["Template. The job script replaces @CONFIG@ with the configuration macro of the job."]) * """

/control/execute macros/physics/DefaultPhysicsList.mac

/nexus/RegisterGeometry ItacaChamber
/nexus/RegisterGenerator $name

/nexus/RegisterPersistencyManager PersistencyManager

/nexus/RegisterRunAction DefaultRunAction
/nexus/RegisterEventAction DefaultEventAction
/nexus/RegisterTrackingAction DefaultTrackingAction

/nexus/RegisterMacro @CONFIG@
"""
end

"""All macros of a variant: file name => text."""
function variant_macros(variant_card::AbstractString)
    r = resolve_detector(joinpath(ROOT, "config/iitaca/iitaca.toml"); variant = variant_card, root = ROOT); sim = simulation_inputs(r.detector)
    files = Dict{String,String}("ItacaChamber.mac" => geometry_macro(sim, r.variant))
    for s in SAMPLES
        files["ItacaChamber.$s.config.mac"] = config_macro(sim, r.variant, s); files["ItacaChamber.$s.init.mac"] = init_macro(sim, r.variant, s)
    end
    r.variant, files
end

variant_cards() = sort(filter(p -> endswith(p, ".toml"), readdir(joinpath(ROOT, "config/iitaca/variants"); join = true)))

if abspath(PROGRAM_FILE) == @__FILE__
    for card in (isempty(ARGS) ? variant_cards() : abspath.(ARGS))
        variant, files = variant_macros(card); dir = joinpath(ROOT, MACRO_ROOT, variant); mkpath(dir)
        for (name, text) in files; write(joinpath(dir, name), text); end
        println("$(length(files)) macros in $(relpath(dir, ROOT))")
    end
end
