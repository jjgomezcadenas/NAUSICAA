using Test, TOML
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
include(joinpath(ROOT, "src/iitaca/IItaca.jl"))
using .IItaca
const IC = ItacaConditions

INDEX = joinpath(ROOT, "config/iitaca/iitaca.toml")
VARIANT(name) = joinpath(ROOT, "config/iitaca/variants/$name.toml")
table(d) = Dict("$g.$n" => v for (g, n, v, u, note) in derived_conditions(d))
short = resolve_detector(INDEX; variant = VARIANT("II_L1000_D500"), root = ROOT)
long = resolve_detector(INDEX; variant = VARIANT("II_L1500_D350"), root = ROOT)
a = table(short.detector); b = table(long.detector)

@testset "Detector set and variants" begin
    set = load_iitaca_set(INDEX)
    @test Set(keys(set.cards)) == Set(SECTIONS)
    @test all(length(set.sha256[s]) == 64 for s in SECTIONS)
    @test short.variant == "II_L1000_D500" && isempty(short.overrides)       # the default set is that variant
    @test short.detector == resolve_detector(INDEX; root = ROOT).detector
    @test long.variant == "II_L1500_D350"
    @test long.detector["geometry"]["fiducial_length_cm"] == 150.0 && long.detector["drift"]["field_V_cm"] == 350.0
    @test set.cards["geometry"]["fiducial_length_cm"] == 100.0                # the set itself is untouched
    # no derived quantity is typed in the cards
    for s in SECTIONS, k in keys(IC.flatten(set.cards[s]))
        @test !occursin(r"^(potential_|drift_voltage|reduced_field|photons_per_electron|roi_lower|roi_upper|.*_mass_kg|pad_count|channel_count)", k)
    end
    # overrides of the run are applied after those of the variant
    r = resolve_detector(INDEX; variant = VARIANT("II_L1500_D350"), overrides = Dict("geometry.fiducial_diameter_cm" => 360.0, "drift.field_V_cm" => 300.0), root = ROOT)
    @test r.detector["drift"]["field_V_cm"] == 300.0 && r.detector["geometry"]["fiducial_diameter_cm"] == 360.0
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("drift.field_Vcm" => 300.0))
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("sipm.pde" => 0.3))
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("drift.field_V_cm" => "300"))
end

@testset "Consistency checks" begin
    # velocity x frame must equal the Z cell
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("drift.ion_drift_velocity_cm_s" => 15.0))
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("ions.frame_ms" => 20.0))
    @test resolve_detector(INDEX; overrides = Dict("drift.ion_drift_velocity_cm_s" => 15.0, "ions.frame_ms" => 20.0)).detector["ions"]["frame_ms"] == 20.0
    @test resolve_detector(INDEX; overrides = Dict("ions.frame_ms" => 20.0, "ions.z_cell_mm" => 4.0)).detector["ions"]["z_cell_mm"] == 4.0
    # the VH1 card must have the frame of the detector (checked with the repository root)
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("drift.ion_drift_velocity_cm_s" => 15.0, "ions.frame_ms" => 20.0), root = ROOT)
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("ions.vh1_card" => "config/specs/response/ion_csa/none.toml"), root = ROOT)
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("geometry.fiducial_length_cm" => -1.0))
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("voltages.documentation.pad_face_ifg_V" => 9000.0))
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("energy.resolution_scaling" => "linear"))
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("ions.pad_side_mm" => 4.0))
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("el.documentation.hole_diameter_mm" => 7.5))   # guide and gas gap exceed the pitch
end

# The simulation reads no documentation value: changing any of them leaves its inputs unchanged.
@testset "Simulation and documentation inputs" begin
    sim = simulation_inputs(short.detector); doc = documentation_inputs(short.detector)
    @test isempty(intersect(keys(sim), keys(doc)))
    @test length(sim) + length(doc) == length(IC.flatten(short.detector))
    @test all(startswith(k, "voltages.documentation.") for k in keys(doc) if startswith(k, "voltages."))
    @test all(startswith(k, "el.documentation.") for k in keys(doc) if startswith(k, "el."))
    @test !any(startswith(k, "voltages.") || startswith(k, "el.") for k in keys(sim))
    for k in ("gas.pressure_bar", "gas.temperature_K", "gas.helium_mole_fraction", "geometry.fiducial_diameter_cm", "geometry.fiducial_length_cm",
              "drift.field_V_cm", "drift.ion_drift_velocity_cm_s", "ions.pad_pitch_mm", "ions.frame_ms", "ions.vh1_card",
              "energy.resolution_fwhm_at_qbb_keV", "energy.roi_width_fwhm", "trigger.frame_side_mm")
        @test haskey(sim, k)
    end
    changed = Dict{String,Any}(k => (v isa Real && !(v isa Bool) ? 1.01 * v : v) for (k, v) in doc if !(k in ("voltages.documentation.pads_V", "voltages.documentation.sipm_plane_V")))
    changed["el.documentation.gap_refractive_index"] = 1.0     # stays below the guide index
    @test simulation_inputs(resolve_detector(INDEX; overrides = changed).detector) == sim
    # the two variants differ in the simulation only by length, field, velocity, frame and VH1 card
    simb = simulation_inputs(long.detector)
    @test sort([k for k in keys(sim) if sim[k] != simb[k]]) == ["drift.field_V_cm", "drift.ion_drift_velocity_cm_s", "geometry.fiducial_length_cm", "ions.frame_ms", "ions.vh1_card"]
end

# Light budget of docs/elm.tex written out for a 5 mm hole 5 mm deep and a 3 mm SiPM.
@testset "EL module and light budget" begin
    hole = 2.5 / (5 + 2.5 + hypot(5, 2.5))                       # axial acceptance of the hole
    guide = hole * (1 - 1 / 1.5) * 0.65 * 0.5
    sensor = 9 / (pi * 2.5^2)
    for t in (a, b)
        @test t["el.hole_open_fraction"] ≈ pi * 2.5^2 / 64
        @test t["el.sipm_coverage"] ≈ 9 / 64
        @test t["el.threshold_voltage_kV"] ≈ 0.8286 * 8
        @test t["el.guide_exit_photons_per_el_photon"] ≈ guide rtol = 1e-5
        @test t["el.sensor_fraction_of_guide_exit"] ≈ sensor
        @test t["el.photoelectrons_per_el_photon"] ≈ guide * sensor * 0.32 rtol = 1e-5
        @test t["el.photoelectrons_at_qbb"] ≈ t["el.photons_at_qbb"] * guide * sensor * 0.32 rtol = 1e-5
        @test t["energy.photoelectron_fwhm_at_qbb_percent"] ≈ 100 * FWHM_PER_SIGMA / sqrt(t["el.photoelectrons_at_qbb"])
        @test t["energy.expected_fwhm_at_qbb_percent"]^2 ≈ t["energy.fano_fwhm_at_qbb_percent"]^2 + t["energy.el_fluctuation_fwhm_at_qbb_percent"]^2 + t["energy.photoelectron_fwhm_at_qbb_percent"]^2
    end
    @test a["voltages.grading_field_kV_cm"] ≈ 7.33 && b["voltages.grading_field_kV_cm"] ≈ 7.58
    big = table(resolve_detector(INDEX; overrides = Dict("el.documentation.sipm_side_mm" => 6.0)).detector)
    @test big["el.sensor_fraction_of_guide_exit"] == 1.0          # a 6 mm sensor covers the 5 mm guide
end

# Reference values: independent calculation of work/iitaca/iitaca_numbers.py (2026-11-17) and
# the formulas written here.
@testset "Derived values, both variants" begin
    for t in (a, b)
        @test t["gas.xenon_compressibility"] ≈ 0.9074 atol = 0.0001
        @test t["gas.density_kg_m3"] ≈ 99.22 atol = 0.01
        @test t["gas.xe136_mass_fraction"] ≈ 0.9041 atol = 0.0001
        @test t["gas.number_density_relative_to_STP"] ≈ 16.39 atol = 0.01
        @test t["geometry.end_area_m2"] ≈ pi * 2.0^2
        @test t["ions.pad_count"] ≈ 1.396e6 rtol = 1e-3
        @test t["el.channel_count"] ≈ 196_350 rtol = 1e-4
        @test t["el.reduced_field_kV_cm_bar"] ≈ 1.875
        @test t["el.photons_per_electron"] ≈ 140 * (15 - 0.8286 * 16 * 0.5) rtol = 1e-5      # 1172
        @test t["energy.electrons_at_qbb"] ≈ 112_230 atol = 1
        @test t["energy.fano_fwhm_at_qbb_keV"] ≈ 6.69 atol = 0.01
        @test t["energy.resolution_fwhm_at_qbb_percent"] ≈ 0.4719 atol = 0.0001
        @test t["energy.roi_lower_keV"] ≈ 2452.0306 atol = 1e-4                               # bounds of the earlier comparison file
        @test t["energy.roi_upper_keV"] ≈ 2463.6294 atol = 1e-4
        @test t["drift.ion_frame_length_mm"] ≈ 3.0
        @test t["voltages.ifg_field_kV_cm"] ≈ 5.0 && t["voltages.pad_gap_field_kV_cm"] ≈ 6.0
    end
    @test a["geometry.fiducial_volume_m3"] ≈ 12.566 atol = 0.001
    @test b["geometry.fiducial_volume_m3"] ≈ 18.850 atol = 0.001
    @test a["geometry.xenon_mass_kg"] ≈ 1247 atol = 1
    @test b["geometry.xenon_mass_kg"] ≈ 1870 atol = 1
    @test a["geometry.xe136_mass_kg"] ≈ 1127.3 atol = 0.5
    @test b["geometry.xe136_mass_kg"] ≈ 1690.9 atol = 0.5
    @test a["drift.drift_voltage_kV"] ≈ 50.0 && b["drift.drift_voltage_kV"] ≈ 52.5
    @test a["drift.reduced_field_V_cm_bar"] ≈ 31.25 && b["drift.reduced_field_V_cm_bar"] ≈ 21.875
    @test [a["voltages.potential_$(k)_kV"] for k in ("tpc_face_ifg", "drift_end", "el_entry", "el_exit")] ≈ [8.0, 58.0, 58.3, 73.3]
    @test [b["voltages.potential_$(k)_kV"] for k in ("tpc_face_ifg", "drift_end", "el_entry", "el_exit")] ≈ [8.0, 60.5, 60.8, 75.8]
    # thermal ion spread sqrt(2 k T L / (e E)), k/e = 8.6173e-5 V/K
    @test a["drift.ion_thermal_spread_mm"] ≈ 1000 * sqrt(2 * 8.6173e-5 * 290 * 1.0 / 5.0e4) rtol = 1e-5      # 1.000 mm
    @test b["drift.ion_thermal_spread_mm"] ≈ 1000 * sqrt(2 * 8.6173e-5 * 290 * 1.5 / 3.5e4) rtol = 1e-5      # 1.464 mm
    @test a["drift.ion_thermal_spread_mm"] ≈ 1.00 atol = 0.005
    @test b["drift.ion_thermal_spread_mm"] ≈ 1.46 atol = 0.005
    @test a["drift.ion_max_drift_time_s"] ≈ 5.0 && b["drift.ion_max_drift_time_s"] ≈ 10.0
    @test a["drift.ion_frames_in_longest_drift"] ≈ 1000 / 3 && b["drift.ion_frames_in_longest_drift"] ≈ 500.0
    @test a["drift.ion_reduced_mobility_cm2_V_s"] ≈ 0.656 atol = 0.001
    @test b["drift.ion_reduced_mobility_cm2_V_s"] ≈ 0.702 atol = 0.001
    @test a["drift.ion_reference_velocity_cm_s"] ≈ 19.8 atol = 0.05
    @test b["drift.ion_reference_velocity_cm_s"] ≈ 13.9 atol = 0.05
    @test a["drift.electron_transverse_spread_mm"] ≈ 9.8 && b["drift.electron_transverse_spread_mm"] ≈ 9.95 * sqrt(1.5)
    @test a["drift.electron_max_drift_time_ms"] ≈ 1000 / 0.9 / 1000 && b["drift.electron_max_drift_time_ms"] ≈ 1500 / 0.9 / 1000
    @test a["drift.electron_max_attachment_loss"] ≈ 1 - exp(-(1 / 0.9) / 100)
    @test b["drift.electron_max_attachment_loss"] ≈ 0.0165 atol = 0.0001
end

# 2 mm pads: the detector of II_L1000_D500 on a cubic 2 mm lattice, with the track selection scaled to it.
@testset "Variant with 2 mm pads" begin
    p2 = resolve_detector(INDEX; variant = VARIANT("II_L1000_D500_P2"), root = ROOT)
    c = table(p2.detector); sim = simulation_inputs(short.detector); sim2 = simulation_inputs(p2.detector)
    @test sort([k for k in keys(sim) if sim[k] != sim2[k]]) == ["ions.frame_ms", "ions.pad_pitch_mm", "ions.vh1_card", "ions.z_cell_mm",
        "tracks.connection_mm", "tracks.hit_cut_keV", "tracks.segment_connection_mm", "tracks.segment_min_hits", "tracks.voxel_width_mm"]
    @test c["drift.ion_frame_length_mm"] ≈ 2.0 && p2.detector["ions"]["pad_pitch_mm"] == 2.0
    @test c["ions.pad_count"] ≈ a["ions.pad_count"] * 9 / 4
    @test c["drift.ion_spread_over_pad_pitch"] ≈ a["drift.ion_spread_over_pad_pitch"] * 3 / 2
    @test c["drift.ion_frames_in_longest_drift"] ≈ 500.0
    t = p2.detector["tracks"]
    @test sqrt(3) * 2.0 < t["segment_connection_mm"] < 4.0        # touching cells, not the next cell but one along an axis
    @test (t["connection_mm"], t["segment_min_hits"], t["voxel_width_mm"]) == (5.3, 9, 2.0)
    @test (t["deposit_min_keV"], t["blob_radius_mm"], t["hit_cut_keV"], t["saved_hit_cut_keV"]) == (10.0, 7.0, 3.0, 3.0)
    # the same readout with the radii of the 3 mm selection
    r3 = resolve_detector(INDEX; variant = VARIANT("II_L1000_D500_P2_R3"), root = ROOT).detector
    @test all(r3[s] == p2.detector[s] for s in SECTIONS if s != "tracks")
    @test sort([k for (k, v) in r3["tracks"] if v != t[k]]) == ["connection_mm", "segment_connection_mm"]
    @test (r3["tracks"]["segment_connection_mm"], r3["tracks"]["connection_mm"]) == (short.detector["tracks"]["segment_connection_mm"], short.detector["tracks"]["connection_mm"])
end

@testset "Energy resolution and ROI" begin
    en = short.detector["energy"]
    @test resolution_fwhm_keV(en, 2457.83) ≈ 11.5988
    @test resolution_fwhm_keV(en, 2457.83 / 4) ≈ 11.5988 / 2          # sqrt(E)
    @test resolution_fwhm_keV(en, 0) == 0
    @test_throws ArgumentError resolution_fwhm_keV(en, -1)
    wide = resolve_detector(INDEX; overrides = Dict("energy.resolution_fwhm_at_qbb_keV" => 8.0, "energy.roi_width_fwhm" => 2.0)).detector["energy"]
    @test collect(roi_keV(wide)) ≈ [2457.83 - 8.0, 2457.83 + 8.0]
    # a ROI reaching beyond the energy window of the saved Geant4 events is refused
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("energy.resolution_fwhm_at_qbb_keV" => 24.5783, "energy.roi_width_fwhm" => 2.0))
    @test resolve_detector(INDEX; overrides = Dict("energy.resolution_fwhm_at_qbb_keV" => 24.5783, "energy.roi_width_fwhm" => 2.0, "generation.saved_energy_max_keV" => 2500.0)) isa ResolvedDetector
end

# Scaling with the free parameters: mass with R^2 L, spread with sqrt(L/E), ladder with E L.
@testset "Scaling" begin
    r = table(resolve_detector(INDEX; overrides = Dict("geometry.fiducial_diameter_cm" => 360.0, "geometry.fiducial_length_cm" => 130.0, "drift.field_V_cm" => 350.0)).detector)
    @test r["geometry.xe136_mass_kg"] / a["geometry.xe136_mass_kg"] ≈ (1.8 / 2.0)^2 * 1.3
    @test r["drift.ion_thermal_spread_mm"] / a["drift.ion_thermal_spread_mm"] ≈ sqrt(1.3 * 500 / 350) rtol = 1e-5
    @test r["voltages.potential_el_exit_kV"] ≈ 8 + 45.5 + 0.3 + 15
    hot = table(resolve_detector(INDEX; overrides = Dict("gas.temperature_K" => 300.0)).detector)
    @test hot["gas.density_kg_m3"] < a["gas.density_kg_m3"] && hot["drift.ion_thermal_spread_mm"] > a["drift.ion_thermal_spread_mm"]
end

@testset "Run cards" begin
    card = joinpath(ROOT, "config/iitaca/runs/hybrid_equivalence_bb0nu_pilot.toml")
    r = resolve_run(card; root = ROOT)
    @test r.run["run_label"] == "hybrid_equivalence_bb0nu_pilot" && !haskey(r.run, "overrides")
    @test r.resolved.detector["geometry"]["fiducial_length_cm"] == 130.0 && r.resolved.detector["ions"]["frame_ms"] == 20.0
    text = resolved_text(r, ROOT); again = TOML.parse(text)
    @test again["detector"]["drift"]["field_V_cm"] == 350.0 && again["overrides"]["drift.field_V_cm"] == 350.0
    @test length(again["run_card"]["sha256"]) == 64 && again["detector_set"]["index"] == "config/iitaca/iitaca.toml"
    @test resolved_text(resolve_run(card; root = ROOT), ROOT) == text
    # the random streams of an event: fixed, distinct per event, part and seed, and those of the ITACA production
    @test event_seed(42, 7, "ions") == event_seed(42, 7, "ions")
    @test length(unique([event_seed(42, 7, "ions"), event_seed(42, 8, "ions"), event_seed(42, 7, "pairs"), event_seed(43, 7, "ions")])) == 4
    @test event_seed(42, 7, "ions") == parse(UInt64, bytes2hex(IItaca.sha256("framed-v1|42|7|ions"))[1:16]; base = 16)
    mktempdir() do dir
        write(joinpath(dir, "a.toml"), replace(read(card, String), "frame_track = \"primary_and_delta_electrons\"" => "frame_track = \"all\""))
        @test_throws ArgumentError resolve_run(joinpath(dir, "a.toml"); root = ROOT)
        write(joinpath(dir, "b.toml"), replace(read(card, String), "seed = 42\n" => "seed = 42\nfield_V_cm = 300.0\n"))
        @test_throws ArgumentError resolve_run(joinpath(dir, "b.toml"); root = ROOT)         # a detector key outside [overrides]
        write(joinpath(dir, "c.toml"), replace(read(card, String), "electronics_seed = 20261113\n" => ""))
        @test_throws ArgumentError resolve_run(joinpath(dir, "c.toml"); root = ROOT)
    end
end
