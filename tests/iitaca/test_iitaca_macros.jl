# The NEXUS macros of the repository are those the detector cards give, and the run cards of the pilots resolve.
using Test
include(joinpath(@__DIR__, "../../scripts/iitaca/write_iitaca_macros.jl"))

@testset "Macros of the variants" begin
    for card in variant_cards()
        variant, files = variant_macros(card); dir = joinpath(ROOT, MACRO_ROOT, variant)
        @test sort(readdir(dir)) == sort(collect(keys(files)))
        for (name, text) in files
            @test read(joinpath(dir, name), String) == text
            @test !occursin("@", replace(text, r"@(SEED|START_ID|OUTPUT|CONFIG)@" => ""))        # no placeholder left over
        end
        r = resolve_detector(joinpath(ROOT, "config/iitaca/iitaca.toml"); variant = card, root = ROOT).detector
        g = files["ItacaChamber.mac"]
        @test occursin("/Geometry/ItacaChamber/radius $(number(r["geometry"]["fiducial_diameter_cm"] / 2)) cm", g)
        @test occursin("/Geometry/ItacaChamber/length $(number(r["geometry"]["fiducial_length_cm"])) cm", g)
        @test occursin("/Geometry/ItacaChamber/pressure 16. bar", g) && occursin("/Geometry/ItacaChamber/helium_mole_fraction 0.0", g)
        for s in SAMPLES
            c = files["ItacaChamber.$s.config.mac"]
            @test occursin("/control/execute macros/iitaca/$variant/ItacaChamber.mac", c)
            @test occursin("/Actions/DefaultEventAction/min_energy 2300. keV", c) && occursin("/Actions/DefaultEventAction/max_energy 2470. keV", c)
        end
        @test occursin("min_energy 2614.511 keV", files["ItacaChamber.tl208.config.mac"]) && occursin("region INNER_SURFACE", files["ItacaChamber.tl208.config.mac"])
        @test occursin("min_energy 2447.86 keV", files["ItacaChamber.bi214.config.mac"]) && occursin("particle gamma", files["ItacaChamber.bi214.config.mac"])
        @test occursin("particle e-", files["ItacaChamber.electron.config.mac"]) && occursin("region GAS", files["ItacaChamber.electron.config.mac"])
        @test occursin("Xe136DecayMode 1", files["ItacaChamber.bb0nu.config.mac"]) && occursin("RegisterGenerator Decay0Interface", files["ItacaChamber.bb0nu.init.mac"])
    end
    @test variant_macros(variant_cards()[1])[2]["ItacaChamber.mac"] != variant_macros(variant_cards()[2])[2]["ItacaChamber.mac"]
end

@testset "Run cards of the pilots" begin
    for (key, sample, track, n) in (("bb0nu", "bb0nu", "primary_and_delta_electrons", 5000), ("single_electrons", "single_electrons", "primary_and_delta_electrons", 5000),
            ("bi214", "bi214", "leading_gamma_electron", 100_000), ("tl208", "tl208", "leading_gamma_electron", 100_000))
        r = resolve_run(joinpath(ROOT, "config/iitaca/runs/II_L1000_D500_$(key)_pilot.toml"); root = ROOT)
        @test r.resolved.variant == "II_L1000_D500" && isempty(r.resolved.overrides)
        @test (r.run["sample"], r.run["frame_track"], r.run["generated_events_per_batch"]) == (sample, track, n)
        @test occursin("{batch}", r.run["hits"]) && occursin("{batch}", r.run["particles"]) && occursin("II_L1000_D500/pilot", r.run["hits"])
        p = resolve_run(joinpath(ROOT, "config/iitaca/runs/II_L1000_D500_$(key).toml"); root = ROOT)          # the production card of the same sample
        @test p.resolved.detector == r.resolved.detector && p.run["frame_track"] == track && p.run["sample"] == sample
        @test p.run["seed"] == r.run["seed"] && occursin("II_L1000_D500/prod", p.run["hits"]) && p.run["run_label"] == "II_L1000_D500_$(key)"
    end
end
