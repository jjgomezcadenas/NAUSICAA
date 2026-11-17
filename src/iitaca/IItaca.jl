"""iItaca detector configuration: the detector set under `config/iitaca/`, its
variants, the run cards that reference them, and the quantities derived from
the resolved configuration.

Three layers:

1. Detector set: an index (`config/iitaca/iitaca.toml`) naming ten cards,
   gas, geometry, drift, voltages, ions, energy, el, trigger, tracks, generation. Design values
   only; nothing derived is typed in them. Values under a `[documentation]`
   table are not read by the simulation: `simulation_inputs` leaves them out
   and `documentation_inputs` returns them.
2. Variant: a card with `[overrides]`, dotted keys into the set
   (`config/iitaca/variants/`).
3. Run card: `detector` (the index), optionally `variant`, and its own
   `[overrides]`, applied after those of the variant.

`resolve_detector` returns the detector values after the overrides and checks
them; `derived_conditions` computes the physics of a resolved detector.
`frame_energy_keV` and `emulated_energy_keV` give the event energy of the
electron side, which is not simulated; `selection_settings` and
`feature_settings` pass the tracks card to `TrackSelection`.

Load `src/event_analysis/EventAnalysis.jl` and `src/itaca/ItacaConditions.jl`
first: the thermal ion spread, the EL yield, the gas state and the override
rules come from there. Units are in the key names; Float64 throughout.
"""
module IItaca
using TOML, SHA
using Random
using ..EventAnalysis: el_photon_yield, sigma_t_ion_mm, elm_light_acceptance, MCEventHits, EventFrame
using ..ItacaConditions: DetectorSet, apply_overrides, flatten, gas_state, enriched_xenon_molar_mass_g_mol,
                         XENON_ISOTOPE_MASS_g_mol, K_B_J_K, TD_V_m2

export SECTIONS, load_iitaca_set, load_variant, ResolvedDetector, resolve_detector, check_detector,
       simulation_inputs, documentation_inputs, potentials_kV, resolution_fwhm_keV, roi_keV, in_roi, derived_conditions,
       FWHM_PER_SIGMA, frame_energy_keV, emulated_energy_keV, selection_settings, feature_settings,
       ResolvedRun, resolve_run, resolved_text, event_seed

const SET_SCHEMA = "nausicaa.iitaca_set.v1"
const VARIANT_SCHEMA = "nausicaa.iitaca_variant.v1"
const RUN_SCHEMA = "nausicaa.iitaca_run.v1"
const SECTIONS = ("gas", "geometry", "drift", "voltages", "ions", "energy", "el", "trigger", "tracks", "generation")
"FWHM of a Gaussian over its sigma, 2 sqrt(2 ln 2)."
const FWHM_PER_SIGMA = 2 * sqrt(2 * log(2))
const CARD_SCHEMAS = Dict(s => "nausicaa.iitaca_$s.v1" for s in SECTIONS)

# ---------------------------------------------------------------- detector set

"""
    load_iitaca_set(index_path)::DetectorSet

Read the index and its ten cards. Require the set schema, every section of
`SECTIONS` exactly once, and each card's own schema.
"""
function load_iitaca_set(index_path::AbstractString)::DetectorSet
    text = read(index_path, String)
    index = TOML.parse(text)
    get(index, "schema", "") == SET_SCHEMA || throw(ArgumentError("Expected schema $SET_SCHEMA in $index_path"))
    paths = get(index, "cards", Dict{String,Any}())
    sort(collect(keys(paths))) == sort(collect(SECTIONS)) ||
        throw(ArgumentError("The index must name exactly the cards $(join(SECTIONS, ", "))"))
    cards = Dict{String,Dict{String,Any}}(); texts = Dict{String,String}("index" => text)
    shas = Dict{String,String}("index" => bytes2hex(sha256(text)))
    for s in SECTIONS
        p = joinpath(dirname(index_path), paths[s])
        isfile(p) || throw(ArgumentError("Missing card $p"))
        t = read(p, String); d = TOML.parse(t)
        get(d, "schema", "") == CARD_SCHEMAS[s] || throw(ArgumentError("Expected schema $(CARD_SCHEMAS[s]) in $p"))
        delete!(d, "schema")
        cards[s] = d; texts[s] = t; shas[s] = bytes2hex(sha256(t))
    end
    DetectorSet(String(index_path), string(get(index, "name", "")), string(get(index, "date", "")), cards, texts, shas)
end

"""
    load_variant(path) -> (name, overrides, text)

Read a variant card: its name and its `[overrides]` (dotted keys into the set).
"""
function load_variant(path::AbstractString)
    text = read(path, String); d = TOML.parse(text)
    get(d, "schema", "") == VARIANT_SCHEMA || throw(ArgumentError("Expected schema $VARIANT_SCHEMA in $path"))
    haskey(d, "name") || throw(ArgumentError("Variant card without name: $path"))
    (name = string(d["name"]), overrides = Dict{String,Any}(get(d, "overrides", Dict{String,Any}())), text = text)
end

"""A detector resolved against the set: the values after the overrides of the
variant and of the run, those overrides, and the provenance."""
struct ResolvedDetector
    set::DetectorSet
    variant::String                       # name of the variant, "" without one
    variant_path::String
    variant_text::String
    overrides::Dict{String,Any}           # variant overrides, then the further ones
    detector::Dict{String,Dict{String,Any}}
end

"""
    resolve_detector(index_path; variant=nothing, overrides=Dict(), root=nothing)::ResolvedDetector

Load the set, apply the overrides of the variant card and then `overrides`, and
check the result with `check_detector`. `root` is the repository root, needed
to check the VH1 card; without it that check is skipped.
"""
function resolve_detector(index_path::AbstractString; variant = nothing, overrides = Dict{String,Any}(), root = nothing)::ResolvedDetector
    set = load_iitaca_set(index_path)
    name = ""; vpath = ""; vtext = ""; all = Dict{String,Any}()
    if variant !== nothing
        v = load_variant(variant); name = v.name; vpath = String(variant); vtext = v.text
        merge!(all, v.overrides)
    end
    merge!(all, Dict{String,Any}(overrides))
    detector = apply_overrides(set, all)
    check_detector(detector; root)
    ResolvedDetector(set, name, vpath, vtext, all, detector)
end

"""A run card resolved against the detector set and its variant."""
struct ResolvedRun
    run_path::String
    run::Dict{String,Any}                 # the run card without its overrides
    resolved::ResolvedDetector
end

"Keys every run card holds."
const RUN_KEYS = ("detector", "run_label", "sample", "hits", "particles", "seed", "electronics_seed", "frame_track")

"""
    resolve_run(run_path; root)::ResolvedRun

Load a run card (schema `nausicaa.iitaca_run.v1`): the detector index, an
optional variant card and `[overrides]`, applied after those of the variant.
The other keys are the run parameters: `run_label`, `sample`, the input files
`hits` and `particles`, `seed` (transport and emulated energy),
`electronics_seed` (VH1 noise and converter), `frame_track`
(`primary_and_delta_electrons` or `leading_gamma_electron`) and optionally
`max_events` (0: all), `max_vector_bytes` and `generated_events_per_batch`
(events generated in Geant4 per batch, when the input files hold only part of
them). Paths are relative to `root`.
A detector key at the top level of the run card is refused.
"""
function resolve_run(run_path::AbstractString; root::AbstractString)::ResolvedRun
    run = TOML.parsefile(run_path)
    get(run, "schema", "") == RUN_SCHEMA || throw(ArgumentError("Expected schema $RUN_SCHEMA in $run_path"))
    for k in RUN_KEYS
        haskey(run, k) || throw(ArgumentError("Run card without '$k': $run_path"))
    end
    run["frame_track"] in ("primary_and_delta_electrons", "leading_gamma_electron") ||
        throw(ArgumentError("frame_track must be primary_and_delta_electrons or leading_gamma_electron"))
    variant = haskey(run, "variant") ? joinpath(root, run["variant"]) : nothing
    r = resolve_detector(joinpath(root, run["detector"]); variant, overrides = get(run, "overrides", Dict{String,Any}()), root)
    physical = Set(k for s in SECTIONS for k in keys(r.set.cards[s]))
    stray = [k for k in keys(run) if k in physical]
    isempty(stray) || throw(ArgumentError("Detector keys $(join(stray, ", ")) in the run card; use [overrides]"))
    ResolvedRun(String(run_path), Dict{String,Any}(k => v for (k, v) in run if k != "overrides"), r)
end

"""
    resolved_text(r::ResolvedRun, root)::String

The resolved configuration as TOML: identity and SHA-256 of the detector cards,
of the variant card and of the run card, the overrides, the detector values
after the overrides, and the run parameters.
"""
function resolved_text(r::ResolvedRun, root::AbstractString)::String
    d = r.resolved
    out = Dict{String,Any}(
        "schema" => "nausicaa.iitaca_resolved_run.v1",
        "detector_set" => Dict("index" => relpath(d.set.index_path, root), "name" => d.set.name, "date" => d.set.date, "sha256" => d.set.sha256),
        "variant" => Dict("name" => d.variant, "card" => d.variant_path == "" ? "" : relpath(d.variant_path, root),
                          "sha256" => d.variant_text == "" ? "" : bytes2hex(sha256(d.variant_text))),
        "overrides" => d.overrides, "detector" => d.detector, "run" => r.run,
        "run_card" => Dict("path" => relpath(r.run_path, root), "sha256" => bytes2hex(sha256(read(r.run_path, String)))))
    sprint(io -> TOML.print(io, out; sorted = true))
end

"""
    event_seed(seed, event_number, part)::UInt64

Seed of one random stream of one event: the first 64 bits of the SHA-256 of
`"framed-v1|seed|event_number|part"`. Each part (`pairs`, `ions`, `ion_grid`,
`vh1_noise`, `vh1_adc`, `frame_noise`, `emulated_energy`) has its own stream,
so a stage can be repeated or left out without changing the others. The same
rule as the ITACA production (`seedfor`), which makes the two chains draw the
same numbers for the same event.
"""
event_seed(seed, event_number, part) = parse(UInt64, bytes2hex(sha256("framed-v1|$seed|$event_number|$part"))[1:16]; base = 16)

"""
    simulation_inputs(detector)::Dict{String,Any}

The detector values the simulation may read, as dotted keys: everything outside
the `[documentation]` tables. The simulation driver is given this dictionary
and nothing else of the detector.
"""
simulation_inputs(d::AbstractDict) = Dict{String,Any}(k => v for (k, v) in flatten(d) if !occursin(".documentation.", k))

"""
    documentation_inputs(detector)::Dict{String,Any}

The values under the `[documentation]` tables, as dotted keys.
"""
documentation_inputs(d::AbstractDict) = Dict{String,Any}(k => v for (k, v) in flatten(d) if occursin(".documentation.", k))

"""
    check_detector(detector; root=nothing)

Consistency of a resolved detector:

- positive pressure, temperature, lengths, pitch, field, velocity and frame;
- track settings: saved hit cut at most the hit cut, segment radius at most
  the connection radius;
- the length drifted in one frame, velocity x frame, equals `ions.z_cell_mm`;
- the frame of the VH1 card equals `ions.frame_ms` (with `root`);
- potentials rise from the pads to the IFG;
- the guide with its gas gap fits in the EL cell;
- the ROI inside the energy windows of the saved Geant4 events (general and gamma);
- known resolution scaling and frame placement.
"""
function check_detector(d::AbstractDict; root = nothing)
    gas = d["gas"]; geo = d["geometry"]; gdoc = geo["documentation"]; dr = d["drift"]; v = d["voltages"]["documentation"]
    ions = d["ions"]; en = d["energy"]; el = d["el"]["documentation"]; tr = d["trigger"]
    positive = (("gas.pressure_bar", gas["pressure_bar"]), ("gas.temperature_K", gas["temperature_K"]),
                ("geometry.fiducial_diameter_cm", geo["fiducial_diameter_cm"]), ("geometry.fiducial_length_cm", geo["fiducial_length_cm"]),
                ("geometry.documentation.ifg_thickness_mm", gdoc["ifg_thickness_mm"]), ("geometry.documentation.pad_gap_mm", gdoc["pad_gap_mm"]),
                ("geometry.documentation.transfer_gap_mm", gdoc["transfer_gap_mm"]), ("geometry.documentation.el_depth_mm", gdoc["el_depth_mm"]),
                ("geometry.documentation.guide_length_mm", gdoc["guide_length_mm"]),
                ("drift.field_V_cm", dr["field_V_cm"]), ("drift.ion_drift_velocity_cm_s", dr["ion_drift_velocity_cm_s"]),
                ("drift.documentation.electron_drift_velocity_mm_us", dr["documentation"]["electron_drift_velocity_mm_us"]),
                ("drift.documentation.electron_lifetime_ms", dr["documentation"]["electron_lifetime_ms"]),
                ("ions.pad_pitch_mm", ions["pad_pitch_mm"]), ("ions.pad_side_mm", ions["pad_side_mm"]), ("ions.frame_ms", ions["frame_ms"]),
                ("ions.z_cell_mm", ions["z_cell_mm"]), ("el.documentation.hole_pitch_mm", el["hole_pitch_mm"]),
                ("el.documentation.hole_diameter_mm", el["hole_diameter_mm"]), ("el.documentation.sipm_side_mm", el["sipm_side_mm"]),
                ("voltages.documentation.el_voltage_kV", v["el_voltage_kV"]),
                ("voltages.documentation.transfer_gap_field_V_cm", v["transfer_gap_field_V_cm"]), ("energy.W_eV", en["W_eV"]), ("energy.qbb_keV", en["qbb_keV"]),
                ("energy.resolution_fwhm_at_qbb_keV", en["resolution_fwhm_at_qbb_keV"]), ("energy.roi_width_fwhm", en["roi_width_fwhm"]),
                ("trigger.frame_side_mm", tr["frame_side_mm"]))
    for (key, value) in positive
        value > 0 || throw(ArgumentError("$key must be positive, found $value"))
    end
    ions["pad_side_mm"] <= ions["pad_pitch_mm"] || throw(ArgumentError("ions.pad_side_mm exceeds the pad pitch"))
    el["hole_diameter_mm"] + 2 * el["gas_gap_radial_mm"] <= el["hole_pitch_mm"] || throw(ArgumentError("The guide with its gas gap exceeds the EL pitch"))
    cell = dr["ion_drift_velocity_cm_s"] / 100 * ions["frame_ms"]    # cm/s to mm/ms
    isapprox(cell, ions["z_cell_mm"]; rtol = 1e-6) ||
        throw(ArgumentError("One frame of $(ions["frame_ms"]) ms at $(dr["ion_drift_velocity_cm_s"]) cm/s spans $(round(cell, digits = 4)) mm, not ions.z_cell_mm = $(ions["z_cell_mm"]) mm"))
    if root !== nothing
        path = joinpath(root, ions["vh1_card"])
        isfile(path) || throw(ArgumentError("Missing VH1 card $path"))
        frame = Float64(TOML.parsefile(path)["frame_ms"])
        frame == ions["frame_ms"] || throw(ArgumentError("VH1 frame of $frame ms in $(ions["vh1_card"]) differs from ions.frame_ms = $(ions["frame_ms"]) ms"))
    end
    v["pads_V"] < v["pad_face_ifg_V"] < v["tpc_face_ifg_V"] || throw(ArgumentError("Potentials must rise from the pads to the TPC-facing IFG face"))
    en["resolution_scaling"] == "sqrt_energy" || throw(ArgumentError("Unsupported energy.resolution_scaling '$(en["resolution_scaling"])'"))
    tr["frame_placement"] == "inside_detector" || throw(ArgumentError("Unsupported trigger.frame_placement '$(tr["frame_placement"])'"))
    0 <= tr["simulation_margin_mm"] && 0 <= tr["el_exclusion_mm"] && 0 <= tr["gamma_energy_cut_keV"] ||
        throw(ArgumentError("Trigger margins and cuts must not be negative"))
    tk = d["tracks"]
    for key in ("segment_connection_mm", "connection_mm", "voxel_width_mm", "blob_radius_mm", "segment_min_hits", "max_segments")
        tk[key] > 0 || throw(ArgumentError("tracks.$key must be positive, found $(tk[key])"))
    end
    0 <= tk["saved_hit_cut_keV"] <= tk["hit_cut_keV"] || throw(ArgumentError("Require 0 <= tracks.saved_hit_cut_keV <= tracks.hit_cut_keV"))
    0 <= tk["deposit_min_keV"] || throw(ArgumentError("tracks.deposit_min_keV must not be negative"))
    tk["segment_connection_mm"] <= tk["connection_mm"] || throw(ArgumentError("tracks.segment_connection_mm exceeds tracks.connection_mm"))
    gn = d["generation"]
    gn["step_limit_mm"] > 0 && 0 <= gn["saved_energy_min_keV"] < gn["saved_energy_max_keV"] || throw(ArgumentError("Invalid generation step limit or saved energy window"))
    r = roi_keV(en)
    gn["saved_energy_min_keV"] < r.lower && r.upper < gn["saved_energy_max_keV"] || throw(ArgumentError("The ROI is not inside the saved energy window of the generation"))
    gn["gamma_saved_energy_min_keV"] < r.lower && r.upper < gn["gamma_saved_energy_max_keV"] || throw(ArgumentError("The ROI is not inside the saved energy window of the gamma samples"))
    all(gn[k] > 0 for k in ("single_electron_energy_keV", "bi214_gamma_keV", "tl208_gamma_keV")) || throw(ArgumentError("Generation energies must be positive"))
    nothing
end

# --------------------------------------------------------------------- physics

"""Potential ladder from the pads to the SiPM plane, in kV: end of drift = TPC-facing
IFG face + E_d L; EL entry = end of drift + E_g g over the transfer gap; EL exit = entry + V_EL."""
function potentials_kV(d::AbstractDict)
    g = d["geometry"]; gdoc = g["documentation"]; v = d["voltages"]["documentation"]; dr = d["drift"]
    ifg = v["tpc_face_ifg_V"] / 1000
    drift_end = ifg + dr["field_V_cm"] * g["fiducial_length_cm"] / 1000
    el_entry = drift_end + v["transfer_gap_field_V_cm"] * gdoc["transfer_gap_mm"] / 10 / 1000
    (pads = v["pads_V"] / 1000, pad_face_ifg = v["pad_face_ifg_V"] / 1000, tpc_face_ifg = ifg,
     drift_end = drift_end, el_entry = el_entry, el_exit = el_entry + v["el_voltage_kV"], sipm_plane = v["sipm_plane_V"] / 1000)
end

"""
    resolution_fwhm_keV(energy_card, E_keV)::Float64

Energy resolution (FWHM) at `E_keV`: FWHM(Qbb) sqrt(E / Qbb).
"""
function resolution_fwhm_keV(en::AbstractDict, E_keV::Real)::Float64
    E_keV >= 0 || throw(ArgumentError("Energy must not be negative"))
    en["resolution_fwhm_at_qbb_keV"] * sqrt(E_keV / en["qbb_keV"])
end

"""
    roi_keV(energy_card) -> (lower, upper)

Region of interest: Qbb plus and minus half of `roi_width_fwhm` FWHM.
"""
function roi_keV(en::AbstractDict)
    half = en["roi_width_fwhm"] * en["resolution_fwhm_at_qbb_keV"] / 2
    (lower = en["qbb_keV"] - half, upper = en["qbb_keV"] + half)
end

"""
    in_roi(energy_card, E_keV)::Bool

True for an energy strictly inside the region of interest.
"""
in_roi(en::AbstractDict, E_keV::Real) = (r = roi_keV(en); r.lower < E_keV < r.upper)

"""
    frame_energy_keV(hits::MCEventHits, box::EventFrame)::Float64

Energy deposited inside the box: the sum, in Float64, of the Monte Carlo hits
whose position lies in the box, edges included. Hits and box must be in the
same coordinate frame. With the box given to the ion response, this is the
energy of the deposits that the ion readout sees.
"""
function frame_energy_keV(h::MCEventHits, box::EventFrame)::Float64
    h.frame == box.frame || throw(ArgumentError("Coordinate frames differ"))
    total = 0.0
    @inbounds for i in eachindex(h.x_mm)
        p = (h.x_mm[i], h.y_mm[i], h.z_mm[i])
        all(box.lower_mm[k] <= p[k] <= box.upper_mm[k] for k in 1:3) && (total += Float64(h.energy_keV[i]))
    end
    total
end

"""
    emulated_energy_keV(rng, E_keV, energy_card)::Float64

Energy measured by the electron side, emulated: `E_keV` plus a Gaussian
deviate of sigma = FWHM(E) / (2 sqrt(2 ln 2)), with `resolution_fwhm_keV`.
One `randn` is drawn from `rng`. Zero energy is returned unchanged.
"""
function emulated_energy_keV(rng::AbstractRNG, E_keV::Real, en::AbstractDict)::Float64
    E_keV + resolution_fwhm_keV(en, E_keV) / FWHM_PER_SIGMA * randn(rng)
end

"""Keywords of `TrackSelection.track_selection` from the tracks card."""
selection_settings(d::AbstractDict) = (t = d["tracks"]; (hit_cut_keV = t["hit_cut_keV"], segment_connection_mm = t["segment_connection_mm"],
    segment_min_hits = t["segment_min_hits"], deposit_min_keV = t["deposit_min_keV"], max_segments = t["max_segments"], connection_mm = t["connection_mm"]))

"""Keywords of `TrackSelection.track_features` from the tracks card."""
feature_settings(d::AbstractDict) = (t = d["tracks"]; (hit_cut_keV = t["hit_cut_keV"], voxel_width_mm = t["voxel_width_mm"],
    connection_mm = t["connection_mm"], blob_radius_mm = t["blob_radius_mm"]))

"""
    derived_conditions(detector::Dict)::Vector{Tuple{String,String,Float64,String,String}}

Compute the derived quantities of a resolved detector. Each entry is
(group, name, value, unit, note). Conventions:

- Gas state from the density model of the card (`ItacaConditions.gas_state`).
- Reduced drift field in Td uses the real-gas number density; the value in
  V/cm/bar uses the plain pressure.
- Ion spread: thermal, sigma = sqrt(2 k T L / (e E)) (`sigma_t_ion_mm`),
  independent of the ion velocity. The reference velocity is
  K0 E N0 / N with the reference reduced mobility of the card.
- Electron spreads sigma = D sqrt(L) with the documentation coefficients.
- Light budget (docs/elm.tex): `elm_light_acceptance` gives the probability
  that an EL photon produces a blue photon at the guide exit (hole acceptance,
  TPB conversion, forward fraction, guide capture). The sensor sees the
  fraction of the guide exit it covers, min(1, side^2 / (pi r^2)), for uniform
  illumination of the exit; then the PDE.
- Energy: electrons at Qbb N = Q / W; photoelectrons N_pe; expected resolution
  FWHM = 2.355 sqrt((F + J) / N + 1 / N_pe), without attachment, SiPM excess
  noise or calibration terms; sigma of the emulated energy = FWHM / 2.355.
"""
function derived_conditions(d::AbstractDict)
    gas = d["gas"]; model = gas["density_model"]; geo = d["geometry"]; gdoc = geo["documentation"]; dr = d["drift"]; ddoc = dr["documentation"]
    v = d["voltages"]["documentation"]; el = d["el"]["documentation"]; ions = d["ions"]; en = d["energy"]
    out = Tuple{String,String,Float64,String,String}[]
    put!(g, n, x, u, note = "") = push!(out, (g, n, Float64(x), u, note))

    M_xe = enriched_xenon_molar_mass_g_mol(gas["xenon_isotopic_fractions"])
    f136 = gas["xenon_isotopic_fractions"]["Xe136"] * XENON_ISOTOPE_MASS_g_mol["Xe136"] / M_xe
    P = gas["pressure_bar"]; T = gas["temperature_K"]
    s = gas_state(P, T, gas["helium_mole_fraction"], M_xe, model)
    put!("gas", "xenon_molar_mass_g_mol", M_xe, "g/mol", "enriched xenon, isotopic fractions of the card")
    put!("gas", "xe136_mass_fraction", f136, "", "")
    put!("gas", "xenon_compressibility", s.xenon_compressibility, "", "second virial, at the card pressure and temperature")
    put!("gas", "density_kg_m3", s.density_kg_m3, "kg/m^3", "")
    put!("gas", "number_density_m3", s.number_density_m3, "1/m^3", "")
    put!("gas", "number_density_relative_to_STP", s.number_density_relative_to_STP, "", "N/N0, N0 at 273.15 K and 101325 Pa")

    radius_cm = geo["fiducial_diameter_cm"] / 2; L_cm = geo["fiducial_length_cm"]
    area_m2 = pi * (radius_cm / 100)^2; volume_m3 = area_m2 * L_cm / 100
    put!("geometry", "fiducial_radius_cm", radius_cm, "cm", "")
    put!("geometry", "end_area_m2", area_m2, "m^2", "")
    put!("geometry", "fiducial_volume_m3", volume_m3, "m^3", "")
    put!("geometry", "xenon_mass_kg", s.xenon_density_kg_m3 * volume_m3, "kg", "density x volume")
    put!("geometry", "xe136_mass_kg", s.xenon_density_kg_m3 * volume_m3 * f136, "kg", "xenon mass x Xe-136 mass fraction")
    put!("geometry", "pads_to_sipm_plane_mm", gdoc["pad_gap_mm"] + gdoc["ifg_thickness_mm"] + L_cm * 10 + gdoc["transfer_gap_mm"] + gdoc["el_depth_mm"] + gdoc["guide_length_mm"], "mm",
        "pad gap + IFG + drift + transfer gap + EL + guides")

    E = dr["field_V_cm"]; v_i = dr["ion_drift_velocity_cm_s"]
    spread = sigma_t_ion_mm(Float32(T), Float32(L_cm), Float32(E))
    put!("drift", "drift_voltage_kV", E * L_cm / 1000, "kV", "E_d L")
    put!("drift", "reduced_field_V_cm_bar", E / P, "V/cm/bar", "plain pressure")
    put!("drift", "reduced_field_Td", E * 100 / s.number_density_m3 / TD_V_m2, "Td", "real-gas number density")
    put!("drift", "ion_thermal_spread_mm", spread, "mm", "sqrt(2 k T L / (e E)), full drift, arXiv:2604.02443 Eq. (2.1)")
    put!("drift", "ion_spread_over_pad_pitch", spread / ions["pad_pitch_mm"], "", "")
    put!("drift", "ion_max_drift_time_s", L_cm / v_i, "s", "L / v_ion")
    put!("drift", "ion_reduced_mobility_cm2_V_s", v_i * s.number_density_relative_to_STP / E, "cm^2/(V s)", "K0 = v N / (N0 E), implied by the card velocity")
    put!("drift", "ion_reference_velocity_cm_s", ddoc["ion_reference_reduced_mobility_cm2_V_s"] * E / s.number_density_relative_to_STP, "cm/s",
        "K0 E N0 / N with the reference mobility of the card")
    put!("drift", "ion_frame_length_mm", v_i / 100 * ions["frame_ms"], "mm", "ion velocity x frame")
    put!("drift", "ion_frames_in_longest_drift", L_cm * 10 / (v_i / 100 * ions["frame_ms"]), "", "drift length over the frame length")
    put!("drift", "electron_transverse_spread_mm", ddoc["electron_transverse_diffusion_mm_sqrtm"] * sqrt(L_cm / 100), "mm", "D_T sqrt(L), full drift; not simulated")
    put!("drift", "electron_longitudinal_spread_mm", ddoc["electron_longitudinal_diffusion_mm_sqrtm"] * sqrt(L_cm / 100), "mm", "D_L sqrt(L), full drift; not simulated")
    t_e_ms = L_cm * 10 / ddoc["electron_drift_velocity_mm_us"] / 1000
    put!("drift", "electron_max_drift_time_ms", t_e_ms, "ms", "L / v_e; not simulated")
    put!("drift", "electron_max_attachment_loss", 1 - exp(-t_e_ms / ddoc["electron_lifetime_ms"]), "", "1 - exp(-t_max / tau); not simulated")

    pot = potentials_kV(d)
    for (k, note) in (("pads", ""), ("pad_face_ifg", ""), ("tpc_face_ifg", ""), ("drift_end", "TPC-facing IFG face + E_d L"),
            ("el_entry", "+ transfer-gap field x gap"), ("el_exit", "+ EL voltage; the highest potential"), ("sipm_plane", ""))
        put!("voltages", "potential_$(k)_kV", getfield(pot, Symbol(k)), "kV", note)
    end
    put!("voltages", "ifg_field_kV_cm", (v["tpc_face_ifg_V"] - v["pad_face_ifg_V"]) / 1000 / (gdoc["ifg_thickness_mm"] / 10), "kV/cm", "across the IFG")
    put!("voltages", "pad_gap_field_kV_cm", (v["pad_face_ifg_V"] - v["pads_V"]) / 1000 / (gdoc["pad_gap_mm"] / 10), "kV/cm", "between the IFG and the pads")
    grading = (pot.el_exit - pot.sipm_plane) / (gdoc["guide_length_mm"] / 10)
    put!("voltages", "grading_field_kV_cm", grading, "kV/cm", "EL exit potential over the guide length")
    put!("voltages", "grading_reduced_field_kV_cm_bar", grading / P, "kV/cm/bar", "")

    Q = en["qbb_keV"]; N = Q * 1000 / en["W_eV"]
    d_cm = gdoc["el_depth_mm"] / 10; E_el = v["el_voltage_kV"] / d_cm; r_mm = el["hole_diameter_mm"] / 2
    photons = el_photon_yield(Float32(el["yield_gain_photons_per_electron_kV"]), Float32(el["yield_threshold_kV_cm_bar"]),
        Float32(E_el / P), Float32(P), Float32(gdoc["el_depth_mm"]))
    guide = elm_light_acceptance(Float32(r_mm), Float32(gdoc["el_depth_mm"]); tpb_efficiency = Float32(el["tpb_conversion_efficiency"]),
        forward_fraction = Float32(el["forward_fraction"]), guide_index = Float32(el["guide_refractive_index"]), gap_index = Float32(el["gap_refractive_index"]))
    sensor = min(1.0, el["sipm_side_mm"]^2 / (pi * r_mm^2))
    detection = guide * sensor * el["sipm_pde"]
    N_pe = photons * N * detection
    put!("el", "field_kV_cm", E_el, "kV/cm", "V_EL / d")
    put!("el", "reduced_field_kV_cm_bar", E_el / P, "kV/cm/bar", "plain pressure")
    put!("el", "threshold_voltage_kV", el["yield_threshold_kV_cm_bar"] * P * d_cm, "kV", "EL voltage at zero yield")
    put!("el", "photons_per_electron", photons, "", "gain x (V_EL - threshold x P x d)")
    put!("el", "photons_at_qbb", photons * N, "", "photons per electron x electrons at Qbb")
    put!("el", "channel_count", area_m2 / (el["hole_pitch_mm"] / 1000)^2, "", "end area over the pitch squared; one SiPM per hole")
    put!("el", "hole_open_fraction", pi * r_mm^2 / el["hole_pitch_mm"]^2, "", "hole area over cell area")
    put!("el", "sipm_coverage", (el["sipm_side_mm"] / el["hole_pitch_mm"])^2, "", "SiPM area over cell area")
    put!("el", "guide_exit_photons_per_el_photon", guide, "", "hole acceptance x TPB conversion x forward fraction x guide capture, docs/elm.tex")
    put!("el", "sensor_fraction_of_guide_exit", sensor, "", "min(1, side^2 / (pi r^2)), uniform illumination of the guide exit")
    put!("el", "photoelectrons_per_el_photon", detection, "", "guide exit x sensor fraction x PDE")
    put!("el", "photoelectrons_per_electron", photons * detection, "", "")
    put!("el", "photoelectrons_at_qbb", N_pe, "", "before attachment")

    put!("ions", "pad_count", area_m2 / (ions["pad_pitch_mm"] / 1000)^2, "", "end area over the pad pitch squared")
    put!("ions", "pad_fill_fraction", (ions["pad_side_mm"] / ions["pad_pitch_mm"])^2, "", "collecting electrode over cell area")

    roi = roi_keV(en); fwhm = en["resolution_fwhm_at_qbb_keV"]; F = en["fano_factor"]; J = el["variance_ratio_J"]
    put!("energy", "electrons_at_qbb", N, "", "Q / W")
    put!("energy", "fano_fwhm_at_qbb_keV", FWHM_PER_SIGMA * sqrt(F * en["W_eV"] * Q * 1000) / 1000, "keV", "2.355 sqrt(F W Q)")
    put!("energy", "fano_fwhm_at_qbb_percent", 100 * FWHM_PER_SIGMA * sqrt(F / N), "%", "2.355 sqrt(F / N)")
    put!("energy", "el_fluctuation_fwhm_at_qbb_percent", 100 * FWHM_PER_SIGMA * sqrt(J / N), "%", "2.355 sqrt(J / N)")
    put!("energy", "photoelectron_fwhm_at_qbb_percent", 100 * FWHM_PER_SIGMA / sqrt(N_pe), "%", "2.355 / sqrt(N_pe)")
    put!("energy", "expected_fwhm_at_qbb_percent", 100 * FWHM_PER_SIGMA * sqrt((F + J) / N + 1 / N_pe), "%", "2.355 sqrt((F + J) / N + 1 / N_pe)")
    put!("energy", "resolution_fwhm_at_qbb_percent", 100 * fwhm / Q, "%", "FWHM of the card over Qbb; used by the simulation")
    put!("energy", "resolution_sigma_at_qbb_keV", fwhm / FWHM_PER_SIGMA, "keV", "FWHM / 2.355")
    put!("energy", "roi_lower_keV", roi.lower, "keV", "Qbb - roi_width_fwhm x FWHM / 2")
    put!("energy", "roi_upper_keV", roi.upper, "keV", "Qbb + roi_width_fwhm x FWHM / 2")
    out
end

end # module
