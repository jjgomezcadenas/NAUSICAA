using Test, TOML, SHA
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl")); using .EventAnalysis: axial_light_acceptance, elm_light_acceptance
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
include(joinpath(ROOT, "src/eitaca/EItaca.jl")); using .EItaca
const IC = ItacaConditions

INDEX = joinpath(ROOT, "config/eitaca/eitaca.toml")
VARIANT = joinpath(ROOT, "config/eitaca/variants/IE_L1500_D400.toml")
table(d) = Dict("$g.$n" => v for (g, n, v, u, note) in derived_conditions(d; root = ROOT))
base = resolve_detector(INDEX; variant = VARIANT, root = ROOT)
a = table(base.detector)

@testset "Optics without guides" begin
    # Equal inputs: the guide capture 1 - n_gap/n_guide = 1/3 is the only difference
    for (r, d) in ((1.25f0, 5f0), (2.5f0, 5f0), (1.0f0, 10f0))
        direct = direct_light_acceptance(r, d; tpb_efficiency = 0.65f0, forward_fraction = 0.5f0)
        guided = elm_light_acceptance(r, d; tpb_efficiency = 0.65f0, forward_fraction = 0.5f0, guide_index = 1.5f0, gap_index = 1.0f0)
        @test direct / guided ≈ 3 rtol = 1e-6
        @test direct ≈ axial_light_acceptance(r, d) * 0.65f0 * 0.5f0
    end
    # Independent quadrature of the axial acceptance: solid angle of the exit disk from points on the axis
    r = 1.25; d = 5.0; n = 200_000
    quad = sum(0.5 * (1 - u / hypot(u, r)) for u in ((i - 0.5) * d / n for i in 1:n)) / n
    @test direct_light_acceptance(Float32(r), Float32(d); tpb_efficiency = 1f0, forward_fraction = 1f0) ≈ quad rtol = 1e-5
    @test_throws ArgumentError direct_light_acceptance(1.25f0, 5f0; tpb_efficiency = 1.2f0, forward_fraction = 0.5f0)
    # The set uses the direct optics; the reference card the guides, with the same factor 3 at equal yield
    ref = resolve_detector(INDEX; variant = VARIANT, overrides = Dict("el.optics" => "light_guides"), root = ROOT)
    @test light_acceptance(base.detector) / light_acceptance(ref.detector) ≈ 3 rtol = 1e-6
end

@testset "Design values and derived quantities" begin
    set = load_eitaca_set(INDEX)
    @test Set(keys(set.cards)) == Set(SECTIONS)
    @test base.variant == "IE_L1500_D400" && isempty(base.overrides)
    @test base.detector == resolve_detector(INDEX; root = ROOT).detector
    for s in SECTIONS, k in keys(IC.flatten(set.cards[s]))
        @test !occursin(r"^(potential_|drift_voltage|reduced_field|photons_per_electron|.*_mass_kg|channel_count|electron_.*_diffusion_mm_sqrtm)", k)
    end
    @test a["gas.xenon_density_kg_m3"] ≈ 82.2 rtol = 1e-3
    @test a["geometry.xe136_mass_kg"] ≈ 1401 rtol = 1e-3
    @test a["voltages.potential_cathode_kV"] ≈ -75 && a["voltages.potential_gate_kV"] ≈ -15 && a["voltages.potential_el_anode_kV"] == 0
    @test a["drift.reduced_field_V_cm_bar"] ≈ 25
    @test a["el.reduced_field_kV_cm_bar"] ≈ 1.875
    @test a["el.photons_per_electron"] ≈ 134 * (1.875 - 1.083) * 16 * 0.5 rtol = 1e-6
    @test a["drift.electron_transverse_diffusion_mm_sqrtm"] ≈ 1094 / 4 * 10 / 1000 rtol = 1e-12     # D*/sqrt(16 bar)
    @test a["drift.electron_longitudinal_diffusion_mm_sqrtm"] ≈ 1489 / 4 * 10 / 1000 rtol = 1e-12
    @test a["sipm.asic_limit_pe_per_bucket"] ≈ 1e-3 * 250e-9 / (3.6e5 * 1.602176634e-19) rtol = 1e-12
    s = electron_spreads_mm(base.detector, 100.0); D = electron_diffusion_mm_sqrtm(base.detector)
    @test s.transverse ≈ D.transverse && electron_spreads_mm(base.detector, 25.0).transverse ≈ D.transverse / 2     # sqrt(L) scaling
end

@testset "Checks" begin
    bad(o) = resolve_detector(INDEX; variant = VARIANT, overrides = o, root = ROOT)
    @test_throws ArgumentError bad(Dict("drift.field_V_cm" => 350.0))          # off the field of the transport coefficients
    @test_throws ArgumentError bad(Dict("el.optics" => "mirrors"))
    @test_throws ArgumentError bad(Dict("el.hole_diameter_mm" => 3.5))        # the SiPM no longer covers the hole
    @test_throws ArgumentError bad(Dict("el.hole_diameter_mm" => 2.0))        # illuminated microcells of the response card differ
    @test_throws ArgumentError bad(Dict("voltages.gate_V" => 1000.0))
end

@testset "Run cards and seeds" begin
    r = resolve_run(joinpath(ROOT, "config/eitaca/runs/IE_pilot_bb0nu_guides.toml"); root = ROOT)
    d = r.resolved.detector
    @test d["geometry"]["fiducial_length_cm"] == 130.0 && d["el"]["optics"] == "light_guides"
    @test el_photons_per_electron(d) ≈ 139 * (1.875 - 0.80) * 16 * 0.5 rtol = 1e-6
    @test occursin("nausicaa.eitaca_resolved_run.v1", resolved_text(r, ROOT))
    # the rule of the hybrid production (seedfor) and of IItaca.event_seed
    @test event_seed(42, 7, "pairs") == parse(UInt64, bytes2hex(sha256("framed-v1|42|7|pairs"))[1:16]; base = 16)
end
