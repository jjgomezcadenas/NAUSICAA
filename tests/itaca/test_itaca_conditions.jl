using Test
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
using .EventAnalysis, .ItacaConditions
const IC = ItacaConditions

card = load_itaca_card(joinpath(ROOT, "config/itaca.toml"))
model = card.data["gas"]["density_model"]
M_xe = enriched_xenon_molar_mass_g_mol(card.data["gas"]["xenon_isotopic_fractions"])
derived = Dict("$g.$n" => v for (g, n, v, u, note) in derived_conditions(card))

@testset "Card" begin
    @test card.data["schema"] == "nausicaa.itaca_conditions.v1"
    @test length(card.sha256) == 64
    @test M_xe ≈ 135.685 atol = 0.002           # docs/handoffs/itaca_detector_page.md
    @test derived["gas.xe136_mass_fraction"] ≈ 0.9041 atol = 0.0002
end

# Reference densities: docs/handoffs/itaca_detector_page.md, "Correction to the gas density"
# (additive partial volumes), Xe-15% He, fiducial volume 12.064 m^3.
@testset "Gas density model" begin
    V = pi * 1.6^2 * 1.5
    @test V ≈ 12.064 atol = 0.001
    for (P, T, rho, m_xe, m_136, m_he, P_eq) in [(15, 290, 77.9, 935, 845, 4.9, 16.195), (16, 290, 83.5, 1002, 906, 5.2, 17.367),
            (17, 290, 89.2, 1070, 968, 5.6, 18.551), (17.5, 290, 92.1, 1105, 999, 5.8, 19.148),
            (18, 290, 95.0, 1140, 1030, 5.9, 19.749), (15, 273.15, 84.0, 1008, 911, 5.2, 16.450)]
        s = gas_state(P, T, 0.15, M_xe, model)
        @test s.density_kg_m3 ≈ rho atol = 0.06
        @test s.xenon_density_kg_m3 * V ≈ m_xe atol = 1
        @test s.xenon_density_kg_m3 * V * 0.9041 ≈ m_136 atol = 1
        @test s.helium_density_kg_m3 * V ≈ m_he atol = 0.06
        @test s.ideal_equivalent_pressure_bar ≈ P_eq atol = 0.001
    end
    s16 = gas_state(16, 290, 0.15, M_xe, model)
    @test s16.xenon_compressibility ≈ 0.9074 atol = 0.0001
    @test s16.number_density_relative_to_STP ≈ 16.14 atol = 0.01
    pure = gas_state(15, 290, 0.0, M_xe, model)
    @test pure.density_kg_m3 * V ≈ 1115 atol = 2               # pure enriched xenon, 15 bar, 290 K
    @test_throws ArgumentError gas_state(16, 290, 1.0, M_xe, model)
    @test_throws ArgumentError enriched_xenon_molar_mass_g_mol(Dict("Xe136" => 0.5))
end

# Final configuration of docs/handoffs/itaca_detector_page.md: 16 bar, 290 K, 360 cm x 130 cm, 350 V/cm.
@testset "Derived values at the card conditions" begin
    @test derived["gas.density_kg_m3"] ≈ 83.50 atol = 0.01
    @test derived["gas.xenon_density_kg_m3"] ≈ 83.07 atol = 0.01
    @test derived["gas.helium_density_kg_m3"] ≈ 0.43 atol = 0.01
    @test derived["gas.ideal_equivalent_pressure_bar"] ≈ 17.367 atol = 0.001
    @test derived["geometry.fiducial_volume_m3"] ≈ 13.23 atol = 0.005
    @test derived["geometry.end_area_m2"] ≈ 10.18 atol = 0.005
    @test derived["geometry.gas_mass_kg"] ≈ 1105 atol = 1
    @test derived["geometry.xenon_mass_kg"] ≈ 1099 atol = 1
    @test derived["geometry.xe136_mass_kg"] ≈ 994 atol = 1
    @test derived["geometry.helium_mass_kg"] ≈ 5.7 atol = 0.05
    @test derived["geometry.track_length_relative_to_pure_xenon_15bar"] ≈ 0.11 atol = 0.005
    @test derived["drift.drift_voltage_kV"] ≈ 45.5 atol = 0.01
    @test derived["drift.reduced_field_V_cm_bar"] ≈ 21.9 atol = 0.05
    @test derived["drift.reduced_field_Td"] ≈ 0.0807 atol = 0.0001
    @test derived["drift.electron_max_drift_time_ms"] ≈ 2.29 atol = 0.005
    @test derived["drift.electron_max_attachment_loss"] ≈ 0.023 atol = 0.0005
    @test derived["drift.electron_transverse_spread_mm"] ≈ 2.29 atol = 0.005
    @test derived["drift.electron_longitudinal_spread_mm"] ≈ 3.90 atol = 0.005
    @test derived["drift.electron_longitudinal_spread_us"] ≈ 6.9 atol = 0.05
    @test derived["drift.ion_max_drift_time_s"] ≈ 7.9 atol = 0.05
    @test derived["drift.ion_thermal_spread_mm"] ≈ 1.36 atol = 0.005
    @test derived["drift.ion_reduced_mobility_cm2_V_s"] ≈ 0.65 * 1.17 atol = 0.01
    @test derived["drift.potential_drift_end_kV"] ≈ 53.5 atol = 0.01
    @test derived["drift.potential_el_entry_kV"] ≈ 53.8 atol = 0.01
    @test derived["drift.potential_el_exit_kV"] ≈ 68.8 atol = 0.01
    @test derived["drift.potential_sipm_plane_kV"] == 0
    @test derived["ion_readout.ifg_field_kV_cm"] ≈ 5.0
    @test derived["ion_readout.pad_gap_field_kV_cm"] ≈ 6.0
    @test derived["el.field_kV_cm"] ≈ 30.0
    @test derived["el.reduced_field_kV_cm_bar"] ≈ 1.875 atol = 0.001
    @test derived["el.photons_per_electron"] ≈ 1195 atol = 1
    @test derived["elm.grading_field_kV_cm"] ≈ 6.88 atol = 0.005
    @test derived["elm.grading_reduced_field_kV_cm_bar"] ≈ 0.43 atol = 0.005
    @test abs(derived["elm.grading_field_card_difference_kV_cm"]) < 0.005
    @test derived["elm.hole_count"] ≈ 10.18 / 0.004^2 rtol = 0.001
    @test derived["energy.electrons_at_qbb"] ≈ 2457.83e3 / 21.9 rtol = 1e-6
    @test derived["energy.fano_fwhm_at_qbb_keV"] ≈ 2.355 * sqrt(0.15 * 21.9 * 2457.83e3) / 1000 rtol = 1e-6
end

# Earlier EL yields quoted in the handoff: 15 bar and 15 kV give 1251 photons per electron;
# 17.5 bar: 17.5 kV -> 1460, 16 kV -> 1251, 15 kV -> 1112.
@testset "EL yield model" begin
    yield(P, V) = el_photon_yield(139f0, 0.80f0, Float32(V / 0.5 / P), Float32(P), 5f0)
    @test yield(15, 15) ≈ 1251 atol = 1
    @test yield(17.5, 17.5) ≈ 1460 atol = 1
    @test yield(17.5, 16) ≈ 1251 atol = 1
    @test yield(17.5, 15) ≈ 1112 atol = 1
end
