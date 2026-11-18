# Shielding cards, shield geometry, the ICS vertex sampler and its weights, and the NEXUS macros.
using Test, Random, Statistics
include(joinpath(@__DIR__, "../../scripts/background/write_background_macros.jl"))

@testset "Shielding cards and geometry" begin
    @test basename.(shield_cards()) == ["ics_eitaca.toml", "ics_iitaca.toml"]
    ii = resolved_shield(joinpath(ROOT, "config/background/ics_iitaca.toml")); ie = resolved_shield(joinpath(ROOT, "config/background/ics_eitaca.toml"))
    @test (ii.family, ii.variant) == ("iitaca", "II_L1500_D350") && (ie.family, ie.variant) == ("eitaca", "IE_L1500_D400")
    for r in (ii, ie)
        g = r.geometry
        @test (g.radius_cm, g.length_cm, g.thickness_cm, g.sampling_length_cm) == (200.0, 150.0, 20.0, 2.9)
        @test g.inner_cm == 201.5 && outer_cm(g) == 221.5
    end
    @test (ii.geometry.zmin_cm, ii.geometry.zmax_cm) == (-21.5, 161.5)        # IFG and MARS region; transfer gap, boards and 10 cm of HDPE
    @test (ie.geometry.zmin_cm, ie.geometry.zmax_cm) == (-10.0, 151.0)        # HV buffer; EL gap and EL structure
    g = ii.geometry
    @test volume_barrel_cm3(g) ≈ π * (221.5^2 - 201.5^2) * 183.0 && volume_disk_cm3(g) ≈ π * 221.5^2 * 20.0
    @test copper_mass_kg(g) ≈ (volume_barrel_cm3(g) + 2 * volume_disk_cm3(g)) * 8.96 / 1000
    @test 95_000 < copper_mass_kg(g) < 105_000                                  # about 99 t: 4.9 m3 of barrel and 6.2 m3 of disks
    @test decays_per_year(g, 1.2e-6) ≈ copper_mass_kg(g) * 1.2e-6 * 365.25 * 86400
    @test_throws ArgumentError load_shielding(joinpath(ROOT, "config/iitaca/gas.toml"))
end

@testset "ICS sampler and weights" begin
    r = resolved_shield(joinpath(ROOT, "config/background/ics_iitaca.toml")); g = r.geometry
    rng = MersenneTwister(20261118); n = 400_000
    v = [sample_vertex(g, rng) for _ in 1:n]
    w = [weight(g, p...) for p in v]; d = [depth_cm(g, p...) for p in v]; k = [piece(g, p...) for p in v]
    @test all(0 .<= d .<= g.thickness_cm) && all(>(0), w)
    tol(x) = 4 * std(x) / sqrt(n)                                              # four standard errors of a weighted mean
    @test isapprox(mean(w), 1.0; atol = tol(w))                                 # the weights restore the uniform distribution
    frac(sel) = (x = w .* sel; (sum(x) / n, tol(x)))
    for (p, expected) in ((:barrel, volume_barrel_cm3(g)), (:cathode, volume_disk_cm3(g)), (:anode, volume_disk_cm3(g)))
        f, e = frac(k .== p); @test isapprox(f, expected / volume_cm3(g); atol = e)
    end
    # weighted depths follow the copper volume: in the barrel the volume grows with the depth as (inner + depth)
    inner = g.inner_cm; out = outer_cm(g); L = cavity_length_cm(g); t = g.thickness_cm
    volume_below(x) = (π * ((inner + x)^2 - inner^2) * L + 2π * out^2 * x) / volume_cm3(g)
    mean_depth = (2π * L * (inner * t^2 / 2 + t^3 / 3) + 2π * out^2 * t^2 / 2) / volume_cm3(g)
    @test isapprox(sum(w .* d) / n, mean_depth; atol = tol(w .* d)) && 10.0 < mean_depth < 10.3
    for x in (5.0, 10.0, 15.0)
        f, e = frac(d .< x); @test isapprox(f, volume_below(x); atol = e)
    end
    @test mean(d) < 3.0 && count(<(3.0), d) / n > 0.6                          # most decays sampled near the inner face
    # the depth density integrates to one and matches the sampled depths
    λ = g.sampling_length_cm; t = g.thickness_cm
    @test isapprox(sum(depth_density(g, x) for x in 0.005:0.01:t) * 0.01, 1.0; atol = 1e-3)
    @test isapprox(count(<(λ), d) / n, (1 - exp(-1)) / (1 - exp(-t / λ)); atol = 0.005)   # truncated exponential below one sampling length
    # uniform depth when the sampling length is zero: the disks have weight one; the barrel
    # keeps the radial factor 2 (inner + depth) / (inner + outer), which averages to one
    u = ShieldGeometry(g.radius_cm, g.length_cm, g.inner_cm, g.zmin_cm, g.zmax_cm, g.thickness_cm, 0.0, g.density_g_cm3)
    vu = [sample_vertex(u, rng) for _ in 1:20_000]
    for p in vu
        expected = piece(u, p...) == :barrel ? 2 * (u.inner_cm + depth_cm(u, p...)) / (u.inner_cm + outer_cm(u)) : 1.0
        @test isapprox(weight(u, p...), expected; rtol = 1e-9)
    end
    wu = [weight(u, p...) for p in vu]; @test isapprox(mean(wu), 1.0; atol = 4 * std(wu) / sqrt(length(wu))) && all(0.95 .< wu .< 1.05)
    @test_throws ArgumentError piece(g, 0.0, 0.0, 50.0)                        # inside the cavity
    @test_throws ArgumentError piece(g, 0.0, 0.0, g.zmax_cm + g.thickness_cm + 1.0)
end

@testset "Macros of the shielding cards" begin
    for card in shield_cards()
        name, files = shield_macros(card); dir = joinpath(ROOT, BACKGROUND_MACRO_ROOT, name)
        @test sort(readdir(dir)) == sort(collect(keys(files)))
        for (file, text) in files
            @test read(joinpath(dir, file), String) == text
            @test !occursin("@", replace(text, r"@(SEED|START_ID|OUTPUT|CONFIG)@" => ""))
        end
        m = files["ItacaChamber.mac"]
        @test occursin("/Geometry/ItacaChamber/radius 200. cm", m) && occursin("/Geometry/ItacaChamber/length 150. cm", m)
        @test occursin("/Geometry/ItacaChamber/shield_thickness 20. cm", m) && occursin("/Geometry/ItacaChamber/shield_sampling_length 2.9 cm", m)
        @test occursin("/Geometry/ItacaChamber/barrel_layer_material G4_TEFLON", m)
        for iso in ("bi214", "tl208")
            c = files["ItacaChamber.$iso.config.mac"]
            @test occursin("/Generator/IonGenerator/region ICS", c) && occursin("/control/execute macros/background/$name/ItacaChamber.mac", c)
            @test occursin("/Actions/DefaultEventAction/min_energy 2420. keV", c) && occursin("/Actions/DefaultEventAction/max_energy 2490. keV", c)
            @test occursin("RegisterGenerator IonGenerator", files["ItacaChamber.$iso.init.mac"])
            @test occursin("RegisterDelayedMacro macros/background/$name/ItacaChamber.$iso.decay.mac", files["ItacaChamber.$iso.init.mac"])
        end
        @test occursin("atomic_number 83\n/Generator/IonGenerator/mass_number 214", files["ItacaChamber.bi214.config.mac"])
        @test occursin("nucleusLimits 214 214 83 83", files["ItacaChamber.bi214.decay.mac"]) && occursin("nucleusLimits 208 208 81 81", files["ItacaChamber.tl208.decay.mac"])
    end
    ii = shield_macros(joinpath(ROOT, "config/background/ics_iitaca.toml"))[2]["ItacaChamber.mac"]; ie = shield_macros(joinpath(ROOT, "config/background/ics_eitaca.toml"))[2]["ItacaChamber.mac"]
    @test occursin("helium_mole_fraction 0.0", ii) && occursin("anode_layer_thickness 10. cm", ii) && occursin("cathode_gap 20.5 cm", ii)
    @test occursin("helium_mole_fraction 0.158", ie) && occursin("cathode_layer_thickness 0. cm", ie) && occursin("cathode_gap 10. cm", ie)
end
