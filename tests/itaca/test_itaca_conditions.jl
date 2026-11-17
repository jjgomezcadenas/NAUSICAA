using Test, TOML
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "scripts/production/run_framed_detector_response.jl"))   # loads EventAnalysis, DetectorResponse, ItacaConditions and event_settings
const IC = ItacaConditions

INDEX = joinpath(ROOT, "config/detector/itaca.toml")
set = load_detector_set(INDEX)
model = set.cards["gas"]["density_model"]
M_xe = enriched_xenon_molar_mass_g_mol(set.cards["gas"]["xenon_isotopic_fractions"])
table(d) = Dict("$g.$n" => v for (g, n, v, u, note) in derived_conditions(d))
derived = table(set.cards)

@testset "Detector set" begin
    @test Set(keys(set.cards)) == Set(IC.SECTIONS)
    @test all(length(set.sha256[s]) == 64 for s in IC.SECTIONS)
    @test M_xe ≈ 135.685 atol = 0.002
    @test derived["gas.xe136_mass_fraction"] ≈ 0.9041 atol = 0.0002
    # no derived quantity is typed in the cards
    for s in IC.SECTIONS, k in keys(IC.flatten(set.cards[s]))
        @test !occursin(r"^(potential_|drift_voltage|reduced_field|photons_per_electron|grading_field|chamber_length|z_el_entry|geant4_z_offset)", k)
    end
    # a run card with an unknown section, key, or wrong type is refused
    @test_throws ArgumentError IC.apply_overrides(set, Dict("drift.field_V_cm" => 350.0, "nowhere.x" => 1.0))
    @test_throws ArgumentError IC.apply_overrides(set, Dict("drift.field_Vcm" => 350.0))
    @test_throws ArgumentError IC.apply_overrides(set, Dict("drift.field_V_cm" => "350"))
    @test IC.apply_overrides(set, Dict("gas.density_model.xenon_virial_a" => 0.09))["gas"]["density_model"]["xenon_virial_a"] == 0.09
    @test set.cards["gas"]["density_model"]["xenon_virial_a"] == 0.083   # the set itself is untouched
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
@testset "Derived values of the default detector" begin
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
    @test !haskey(derived, "geometry.geant4_z_offset_cm")   # the Geant4 volume uses detector coordinates
    @test derived["geometry.z_el_entry_cm"] == 130.0
    @test derived["drift.drift_voltage_kV"] ≈ 45.5 atol = 0.01
    @test derived["drift.reduced_field_V_cm_bar"] ≈ 21.9 atol = 0.05
    @test derived["drift.reduced_field_Td"] ≈ 0.0807 atol = 0.0001
    @test derived["drift.electron_max_drift_time_ms"] ≈ 2.29 atol = 0.005
    @test derived["drift.electron_max_attachment_loss"] ≈ 0.023 atol = 0.0005
    @test derived["drift.electron_transverse_diffusion_mm_sqrtbar_sqrtcm"] ≈ 0.8046 atol = 0.0001
    @test derived["drift.electron_longitudinal_diffusion_mm_sqrtbar_sqrtcm"] ≈ 1.3691 atol = 0.0001
    @test derived["drift.electron_transverse_spread_mm"] ≈ 2.29 atol = 0.005
    @test derived["drift.electron_longitudinal_spread_mm"] ≈ 3.90 atol = 0.005
    @test derived["drift.electron_longitudinal_spread_us"] ≈ 6.9 atol = 0.05
    @test derived["drift.ion_max_drift_time_s"] ≈ 130/15 atol = 1e-9
    @test derived["drift.ion_thermal_spread_mm"] ≈ 1.36 atol = 0.005
    @test derived["drift.ion_reduced_mobility_cm2_V_s"] ≈ 0.65 * 1.17 * 15/16.5 atol = 0.01   # adopted velocity of 15.0 cm/s against 16.5 from the mobility
    @test derived["drift.ion_bucket_length_mm"] ≈ 3.0 atol = 1e-9
    @test derived["voltages.potential_drift_end_kV"] ≈ 53.5 atol = 0.01
    @test derived["voltages.potential_el_entry_kV"] ≈ 53.8 atol = 0.01
    @test derived["voltages.potential_el_exit_kV"] ≈ 68.8 atol = 0.01
    @test derived["voltages.potential_sipm_plane_kV"] == 0
    @test derived["voltages.ifg_field_kV_cm"] ≈ 5.0
    @test derived["voltages.pad_gap_field_kV_cm"] ≈ 6.0
    @test derived["voltages.grading_field_kV_cm"] ≈ 6.88 atol = 0.005
    @test derived["voltages.grading_reduced_field_kV_cm_bar"] ≈ 0.43 atol = 0.005
    @test derived["el.field_kV_cm"] ≈ 30.0
    @test derived["el.reduced_field_kV_cm_bar"] ≈ 1.875 atol = 0.001
    @test derived["el.photons_per_electron"] ≈ 1195 atol = 1
    @test derived["el.hole_count"] ≈ 10.18 / 0.004^2 rtol = 0.001
    @test derived["sipm.sipm_count"] == derived["el.hole_count"]
    @test derived["ions.pad_count"] ≈ 10.18 / 0.003^2 rtol = 0.001
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

# The pilot run card overrides the default detector to the conditions of the current
# production. The detector-response settings it produces must equal, field by field, the
# settings the production card gave through the earlier flat keys (Float32, as the code uses).
@testset "Run cards and production settings" begin
    r = resolve_run(joinpath(ROOT, "config/runs/bb0nu_vh1_legacy_conditions_pilot100.toml"); root = ROOT)
    @test r.detector["gas"]["pressure_bar"] == 15.0
    @test r.set.cards["gas"]["pressure_bar"] == 16.0
    old = TOML.parsefile(joinpath(ROOT, "config/readout/bb0nu_detector_response_3mm_mainframe20.toml"))
    @test r.run["seed"] == old["seed"] && r.run["hits"] == old["hits"]
    p = EventPlacement(1, 37.5f0, 112.5f0)
    s = event_settings(r, p)
    F(k) = Float32(old[k]); scale = sqrt(p.electron_drift_cm / F("pressure_bar"))
    @test s.chamber_length_mm == F("chamber_length_mm") && s.W_eV == F("W_eV") && s.fano_factor == F("fano_factor")
    @test s.ion_grid.pitch_mm == F("ion_pitch_mm") && s.electron_grid.pitch_mm == F("EL_pitch_mm")
    @test s.ion_bucket_ms == F("ion_bucket_ms") && s.electron_bucket_us == F("electron_bucket_us")
    @test s.ion_transport.ion_drift_velocity_mm_ms == F("ion_drift_velocity_mm_ms")
    @test s.ion_transport.sigma_t_mm == sigma_t_ion_mm(F("temperature_K"), p.ion_drift_cm, F("drift_field_V_cm"))
    @test s.electron_transport.electron_drift_velocity_mm_us == F("electron_drift_velocity_mm_us")
    @test s.electron_transport.sigma_t_mm == F("electron_dT_mm_sqrtbar_sqrtcm") * scale
    @test s.electron_transport.sigma_l_mm == F("electron_dL_mm_sqrtbar_sqrtcm") * scale
    @test s.electron_lifetime_ms == F("electron_lifetime_ms")
    @test s.el_photons_per_electron == el_photon_yield(F("EL_gain_photons_e_kV"), F("EL_threshold_kV_cm_bar"), F("EL_reduced_field_kV_cm_bar"), F("pressure_bar"), F("hole_depth_mm"))
    @test s.el_variance_ratio == F("EL_variance_ratio")
    @test s.optical_acceptance == elm_light_acceptance(F("hole_diameter_mm") / 2, F("hole_depth_mm"); tpb_efficiency = F("tpb_efficiency"),
        forward_fraction = F("forward_fraction"), guide_index = F("guide_index"), gap_index = F("gap_index"),
        propagation_transmission = F("propagation_transmission"), coupling_transmission = F("coupling_transmission"), electrode_transmission = F("electrode_transmission"))
    @test s.photon_wavelength == Symbol(old["photon_wavelength"])
    # grid offsets come from the run seed and the event number, as before
    s2 = event_settings(r, p)
    @test s2.ion_grid.edge_mm == s.ion_grid.edge_mm
    # every run card resolves, names existing response cards, and defines the stages
    for card in readdir(joinpath(ROOT, "config/runs"); join = true)
        rr = resolve_run(card; root = ROOT)
        @test haskey(rr.run, "hit_production") && haskey(rr.run["hit_production"], "output_subdir")
        for (k, v) in rr.run["response"]
            v isa AbstractString && endswith(v, ".toml") && @test isfile(joinpath(ROOT, v))
        end
        ec = electronics_card(rr)
        @test haskey(ec, "csa_card") && haskey(ec, "adc_bits") && haskey(ec, "seed")
    end
    # a stray physical key at the top level of a run card is refused
    bad = tempname() * ".toml"
    write(bad, replace(read(r.run_path, String), "seed = 42" => "seed = 42\npressure_bar = 15.0"))
    @test_throws ArgumentError resolve_run(bad; root = ROOT)
    # the VH1 frame must equal the ion bucket of the resolved detector
    pilot15 = joinpath(ROOT, "config/runs/bb0nu_itaca_pilot_10x100_ion15.toml")
    mismatch = tempname() * ".toml"
    write(mismatch, replace(read(pilot15, String), "vh1_device_changes_cf20f.toml" => "vh1_device_changes_cf20f_15ms.toml"))
    @test_throws ArgumentError resolve_run(mismatch; root = ROOT)
    f15 = resolve_run(joinpath(ROOT, "config/runs/bb0nu_itaca_pilot_10x100_ion20_15ms.toml"); root = ROOT)
    @test f15.detector["ions"]["ion_bucket_ms"] == 15.0 && f15.detector["drift"]["ion_drift_velocity_cm_s"] == 20.0
    @test f15.set.cards["ions"]["ion_bucket_ms"] == 20.0 && length(f15.overrides) == 2
    @test event_settings(f15, p).ion_bucket_ms == 15f0 && event_settings(f15, p).ion_transport.ion_drift_velocity_mm_ms == 0.2f0
    @test event_settings(f15, p).ion_transport.sigma_t_mm == event_settings(resolve_run(pilot15; root = ROOT), p).ion_transport.sigma_t_mm
    # configuration writing: detector stage is shared by run cards with the same detector response, hits stage is per output
    tmp = mktempdir()
    reldir = write_configuration(r, tmp)
    @test reldir == "results/readout/bb0nu_newconfig_legacy_conditions_pilot100/configuration"
    @test write_configuration(r, tmp) == reldir
    resolved = TOML.parsefile(joinpath(tmp, reldir, "resolved.toml"))
    @test resolved["detector"]["gas"]["pressure_bar"] == 15.0
    @test resolved["overrides"]["drift.field_V_cm"] == 300.0
    @test resolved["detector_set"]["sha256"]["drift"] == r.set.sha256["drift"]
    @test !haskey(resolved["run"], "hit_production") && haskey(resolved["run"], "selection")
    hitsdir = write_configuration(r, tmp; stage = "hits")
    @test endswith(hitsdir, "vh1_hit_output_100_cf20f/configuration")
    @test haskey(TOML.parsefile(joinpath(tmp, hitsdir, "resolved.toml"))["run"], "hit_production")
    write(joinpath(tmp, reldir, "resolved.toml"), "changed = true\n")
    @test_throws ErrorException write_configuration(r, tmp)
    texts = configuration_texts(r, ROOT)
    @test haskey(texts, "run") && haskey(texts, "resolved") && haskey(texts, "vh1_card") && haskey(texts, "detector/drift")
    # the two production run cards share the detector-response configuration
    a = resolve_run(joinpath(ROOT, "config/runs/bb0nu_vh1_cf20f_full_6156.toml"); root = ROOT)
    b = resolve_run(joinpath(ROOT, "config/runs/bb0nu_vh1_cf20f_1000.toml"); root = ROOT)
    @test resolved_text(a, ROOT; full = false) == resolved_text(b, ROOT; full = false)
    @test resolved_text(a, ROOT) != resolved_text(b, ROOT)
end
