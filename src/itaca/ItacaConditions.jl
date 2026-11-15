"""ITACA detector configuration: the detector set under `config/detector/`,
the run cards that reference it with explicit overrides, the quantities
derived from the resolved configuration, and the adapter that presents a
resolved run to the production scripts in the shape of their earlier cards.

Three layers (notes/html/detector_configuration.html):

1. Detector set: an index (`config/detector/itaca.toml`) naming eight cards,
   gas, geometry, drift, voltages, el, sipm, ions, energy. Design values only;
   nothing derived is typed in them.
2. Response models: device-model cards (VH1, CSA, SiPM, ASIC) named by the
   run card and loaded by their own functions.
3. Run card: `detector` (the index), `[overrides]` with dotted keys into the
   set, and the run parameters (events, seeds, selection, stages).

`resolve_run` loads the set, applies and checks the overrides and returns a
`ResolvedRun`. `derived_conditions` computes the physics of a resolved
detector. The production scripts read the resolved detector sections directly;
`write_configuration` saves the resolved configuration with the hashes of every
input beside the results, and `configuration_texts` gives the texts they store
in their outputs.

Load `src/event_analysis/EventAnalysis.jl` first; the EL yield and the thermal
ion spread come from there. Units are in the key names; Float64 throughout.
"""
module ItacaConditions
using TOML, SHA
using ..EventAnalysis: el_photon_yield, sigma_t_ion_mm

export DetectorSet, load_detector_set, ResolvedRun, resolve_run, resolve_detector,
       enriched_xenon_molar_mass_g_mol, xenon_compressibility, GasState, gas_state,
       derived_conditions, electronics_card, resolved_text, write_configuration, configuration_texts, QBB_XE136_keV

"Molar gas constant, J/(mol K) (CODATA 2018)."
const R_J_mol_K = 8.314462618
"Boltzmann constant, J/K (CODATA 2018)."
const K_B_J_K = 1.380649e-23
"Loschmidt number density at 273.15 K and 101325 Pa, 1/m^3 (CODATA 2018)."
const N0_STP_m3 = 2.686780111e25
"One Townsend, V m^2."
const TD_V_m2 = 1e-21
"Q-value of the Xe-136 double beta decay, keV (Redshaw et al., PRL 98 (2007) 053003); the card may override it."
const QBB_XE136_keV = 2457.83
"Molar mass of helium-4, g/mol, as in simulation/nexus/source/geometries/ItacaChamber.cc."
const HELIUM_MOLAR_MASS_g_mol = 4.0026
"Atomic masses of the xenon isotopes, g/mol (AME2020)."
const XENON_ISOTOPE_MASS_g_mol = Dict("Xe124"=>123.905885, "Xe126"=>125.904298, "Xe128"=>127.903531,
    "Xe129"=>128.904781, "Xe130"=>129.903509, "Xe131"=>130.905084, "Xe132"=>131.904155,
    "Xe134"=>133.905395, "Xe136"=>135.907214)

const SET_SCHEMA = "nausicaa.detector_set.v1"
const RUN_SCHEMA = "nausicaa.run.v1"
const SECTIONS = ("gas", "geometry", "drift", "voltages", "el", "sipm", "ions", "energy")
const CARD_SCHEMAS = Dict(s => "nausicaa.detector_$s.v1" for s in SECTIONS)

# ---------------------------------------------------------------- detector set

"""The detector set: the index, the eight cards as parsed, their texts and SHA-256."""
struct DetectorSet
    index_path::String
    name::String
    date::String
    cards::Dict{String,Dict{String,Any}}
    texts::Dict{String,String}
    sha256::Dict{String,String}
end

"""
    load_detector_set(index_path)::DetectorSet

Read the index and its eight cards. Require the set schema, every section of
`SECTIONS` present exactly once, and each card's own schema.
"""
function load_detector_set(index_path::AbstractString)::DetectorSet
    text = read(index_path, String)
    index = TOML.parse(text)
    get(index, "schema", "") == SET_SCHEMA || throw(ArgumentError("Expected schema $SET_SCHEMA in $index_path"))
    paths = get(index, "cards", Dict{String,Any}())
    sort(collect(keys(paths))) == sort(collect(SECTIONS)) ||
        throw(ArgumentError("The index must name exactly the cards $(join(SECTIONS, ", "))"))
    cards = Dict{String,Dict{String,Any}}(); texts = Dict{String,String}(); shas = Dict{String,String}()
    texts["index"] = text; shas["index"] = bytes2hex(sha256(text))
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

"""Flatten nested tables to dotted keys, e.g. `gas.density_model.xenon_virial_a`."""
function flatten(d::AbstractDict, prefix::String = "")
    out = Dict{String,Any}()
    for (k, v) in d
        name = prefix == "" ? string(k) : "$prefix.$k"
        if v isa AbstractDict && !(k == "xenon_isotopic_fractions")
            merge!(out, flatten(v, name))
        else
            out[name] = v
        end
    end
    out
end

# -------------------------------------------------------------------- run card

"""A run resolved against the detector set: the detector values after the
overrides, the overrides themselves, the run card, and the provenance."""
struct ResolvedRun
    run_path::String
    run::Dict{String,Any}
    set::DetectorSet
    detector::Dict{String,Dict{String,Any}}
    overrides::Dict{String,Any}
end

"""
    apply_overrides(set::DetectorSet, overrides)::Dict

Return a deep copy of the set's cards with the overrides applied. Each key
is `section.key` or `section.table.key`; it must exist in the set and the new
value must have the type class of the old one (number, string or boolean).
"""
function apply_overrides(set::DetectorSet, overrides::AbstractDict)
    cards = Dict{String,Dict{String,Any}}(s => deepcopy(d) for (s, d) in set.cards)
    for (key, value) in overrides
        parts = split(key, ".")
        length(parts) >= 2 || throw(ArgumentError("Override key '$key' must be section.key"))
        haskey(cards, parts[1]) || throw(ArgumentError("Override '$key': unknown section '$(parts[1])'"))
        table = cards[parts[1]]
        for p in parts[2:end-1]
            (haskey(table, p) && table[p] isa AbstractDict) || throw(ArgumentError("Override '$key': no table '$p'"))
            table = table[p]
        end
        last = parts[end]
        haskey(table, last) || throw(ArgumentError("Override '$key' does not exist in the detector set"))
        old = table[last]
        same = (old isa Real && !(old isa Bool) && value isa Real && !(value isa Bool)) ||
               (old isa AbstractString && value isa AbstractString) || (old isa Bool && value isa Bool)
        same || throw(ArgumentError("Override '$key': type $(typeof(value)) does not match $(typeof(old))"))
        table[last] = value
    end
    cards
end

"""
    resolve_detector(index_path; overrides=Dict())::Dict

The detector set with overrides applied, without a run card.
"""
resolve_detector(index_path::AbstractString; overrides = Dict{String,Any}()) =
    apply_overrides(load_detector_set(index_path), overrides)

"""
    resolve_run(run_path; root=dirname of the repository)::ResolvedRun

Load a run card (schema `nausicaa.run.v1`), its detector set and apply its
`[overrides]`. The run card may hold no physical key outside `[overrides]`;
the keys of the detector cards are refused at the top level of the run card.
"""
function resolve_run(run_path::AbstractString; root::AbstractString = repository_root(run_path))::ResolvedRun
    run = TOML.parsefile(run_path)
    get(run, "schema", "") == RUN_SCHEMA || throw(ArgumentError("Expected schema $RUN_SCHEMA in $run_path"))
    haskey(run, "detector") || throw(ArgumentError("Run card without 'detector' index: $run_path"))
    set = load_detector_set(joinpath(root, run["detector"]))
    overrides = Dict{String,Any}(get(run, "overrides", Dict{String,Any}()))
    physical = Set(k for s in SECTIONS for k in keys(set.cards[s]))
    stray = [k for k in keys(run) if k in physical]
    isempty(stray) || throw(ArgumentError("Physical keys $(join(stray, ", ")) in the run card; use [overrides]"))
    ResolvedRun(String(run_path), run, set, apply_overrides(set, overrides), overrides)
end

"""Repository root: the parent of the `config/` directory that holds the card."""
function repository_root(path::AbstractString)
    p = abspath(path)
    while basename(p) != "config" && p != dirname(p)
        p = dirname(p)
    end
    basename(p) == "config" ? dirname(p) : throw(ArgumentError("Cannot locate the repository root from $path"))
end

# --------------------------------------------------------------------- physics

"""
    enriched_xenon_molar_mass_g_mol(fractions)::Float64

Molar mass of xenon with the given isotopic fractions (keys `Xe129` ... `Xe136`,
fractions summing to one within 1e-6).
"""
function enriched_xenon_molar_mass_g_mol(fractions::AbstractDict)::Float64
    total = sum(values(fractions))
    abs(total - 1) < 1e-6 || throw(ArgumentError("Isotopic fractions sum to $total"))
    sum(f * XENON_ISOTOPE_MASS_g_mol[k] for (k, f) in fractions)
end

"""
    xenon_compressibility(P_bar, T_K, model)::Float64

Second-virial compressibility z = 1 + (a - b / T_r^1.6) P_r / T_r with
T_r = T / T_c and P_r = P / P_c, the constants from `[density_model]` of the
gas card (NEXT DocDB 1183; simulation/nexus/source/geometries/ItacaChamber.cc).
"""
function xenon_compressibility(P_bar::Real, T_K::Real, model::AbstractDict)::Float64
    T_r = T_K / model["xenon_critical_temperature_K"]
    P_r = P_bar / model["xenon_critical_pressure_bar"]
    1 + (model["xenon_virial_a"] - model["xenon_virial_b"] / T_r^1.6) * P_r / T_r
end

"""Gas state of the Xe-He mixture at one pressure and temperature.

Additive partial volumes for mole fractions x_xe and x_he, xenon with
compressibility z, helium ideal:
rho = P (x_xe M_xe + x_he M_he) / ((x_xe z + x_he) R T) and the number density
n = P / ((x_xe z + x_he) k T). `ideal_equivalent_pressure_bar` is the pressure
of an ideal gas with the same number density at the same temperature.
"""
struct GasState
    pressure_bar::Float64
    temperature_K::Float64
    helium_mole_fraction::Float64
    xenon_molar_mass_g_mol::Float64
    xenon_compressibility::Float64
    mixture_compressibility::Float64
    density_kg_m3::Float64
    xenon_density_kg_m3::Float64
    helium_density_kg_m3::Float64
    number_density_m3::Float64
    xenon_number_density_m3::Float64
    number_density_relative_to_STP::Float64
    ideal_equivalent_pressure_bar::Float64
end

"""
    gas_state(P_bar, T_K, x_he, M_xe_g_mol, model)::GasState

Evaluate the gas state with the density model of the card. `x_he = 0` gives
pure xenon.
"""
function gas_state(P_bar::Real, T_K::Real, x_he::Real, M_xe_g_mol::Real, model::AbstractDict)::GasState
    0 <= x_he < 1 && P_bar > 0 && T_K > 0 || throw(ArgumentError("Require 0 <= x_he < 1, P > 0, T > 0"))
    model["helium"] == "ideal" && model["mixing_rule"] == "additive_partial_volumes" ||
        throw(ArgumentError("Unsupported density model"))
    x_xe = 1 - x_he
    z = xenon_compressibility(P_bar, T_K, model)
    zmix = x_xe * z + x_he
    P = P_bar * 1e5
    n = P / (zmix * K_B_J_K * T_K)
    rho_xe = P * x_xe * M_xe_g_mol * 1e-3 / (zmix * R_J_mol_K * T_K)
    rho_he = P * x_he * HELIUM_MOLAR_MASS_g_mol * 1e-3 / (zmix * R_J_mol_K * T_K)
    GasState(P_bar, T_K, x_he, M_xe_g_mol, z, zmix, rho_xe + rho_he, rho_xe, rho_he, n, x_xe * n,
        n / N0_STP_m3, P_bar / zmix)
end

"""Potential ladder from the pads to the SiPM plane, in kV, from the resolved detector."""
function potentials_kV(d::AbstractDict)
    g = d["geometry"]; v = d["voltages"]; dr = d["drift"]
    ifg = v["tpc_face_ifg_V"] / 1000
    drift_end = ifg + dr["field_V_cm"] * g["fiducial_length_cm"] / 1000
    el_entry = drift_end + v["transfer_gap_field_V_cm"] * g["transfer_gap_mm"] / 10 / 1000
    el_exit = el_entry + v["el_voltage_kV"]
    (pads = v["pads_V"] / 1000, pad_face_ifg = v["pad_face_ifg_V"] / 1000, tpc_face_ifg = ifg,
     drift_end = drift_end, el_entry = el_entry, el_exit = el_exit, sipm_plane = v["sipm_plane_V"] / 1000)
end

"""Diffusion coefficient from mm/sqrt(m) to the mm sqrt(bar)/sqrt(cm) of the production cards."""
diffusion_card_units(D_mm_sqrtm::Real, P_bar::Real) = D_mm_sqrtm * sqrt(P_bar) / 10

"""
    derived_conditions(detector::Dict)::Vector{Tuple{String,String,Float64,String,String}}

Compute the derived quantities of a resolved detector. Each entry is
(group, name, value, unit, note). Conventions:

- Pure-xenon reference for the track length: enriched xenon at 15 bar and
  290 K with the same density model. Track length scales with the inverse of
  the effective atom density, helium counted as 2/54 of a xenon atom (ratio of
  atomic numbers).
- Reduced drift field in Td uses the real-gas number density; the value in
  V/cm/bar uses the plain pressure.
- Potentials: sensors at ground; drift end = TPC-facing IFG + E_d L; EL entry
  = drift end + E_g g over the transfer gap; EL exit = entry + V_EL; graded to
  the SiPM plane along the guides (docs/elm.tex).
- Electron spreads sigma = D sqrt(L) with D in mm/sqrt(m); the longitudinal
  spread in time is sigma_L / v_e. Ion spread from `sigma_t_ion_mm`.
- Energy: electrons at Q_bb = Q / W; Fano-limited FWHM = 2.355 sqrt(F W Q).
"""
function derived_conditions(d::AbstractDict)
    gas = d["gas"]; model = gas["density_model"]; geo = d["geometry"]; dr = d["drift"]; v = d["voltages"]
    el = d["el"]; sipm = d["sipm"]; ions = d["ions"]; en = d["energy"]
    out = Tuple{String,String,Float64,String,String}[]
    put!(g, n, x, u, note = "") = push!(out, (g, n, Float64(x), u, note))

    M_xe = enriched_xenon_molar_mass_g_mol(gas["xenon_isotopic_fractions"])
    f136 = gas["xenon_isotopic_fractions"]["Xe136"] * XENON_ISOTOPE_MASS_g_mol["Xe136"] / M_xe
    P = gas["pressure_bar"]; T = gas["temperature_K"]; x_he = gas["helium_mole_fraction"]
    s = gas_state(P, T, x_he, M_xe, model)
    ref = gas_state(15.0, 290.0, 0.0, M_xe, model)
    put!("gas", "xenon_molar_mass_g_mol", M_xe, "g/mol", "enriched xenon, isotopic fractions of the card")
    put!("gas", "xe136_mass_fraction", f136, "", "")
    put!("gas", "xenon_compressibility", s.xenon_compressibility, "", "second virial, at the card pressure and temperature")
    put!("gas", "density_kg_m3", s.density_kg_m3, "kg/m^3", "additive partial volumes")
    put!("gas", "xenon_density_kg_m3", s.xenon_density_kg_m3, "kg/m^3", "")
    put!("gas", "helium_density_kg_m3", s.helium_density_kg_m3, "kg/m^3", "")
    put!("gas", "number_density_m3", s.number_density_m3, "1/m^3", "")
    put!("gas", "number_density_relative_to_STP", s.number_density_relative_to_STP, "", "N/N0, N0 at 273.15 K and 101325 Pa")
    put!("gas", "ideal_equivalent_pressure_bar", s.ideal_equivalent_pressure_bar, "bar", "ideal gas with the same number density at the card temperature")
    put!("gas", "reference_pure_xenon_density_kg_m3", ref.density_kg_m3, "kg/m^3", "enriched xenon, 15 bar, 290 K")

    radius_cm = geo["fiducial_diameter_cm"] / 2; L_cm = geo["fiducial_length_cm"]
    area_m2 = pi * (radius_cm / 100)^2; volume_m3 = area_m2 * L_cm / 100
    put!("geometry", "fiducial_radius_cm", radius_cm, "cm", "")
    put!("geometry", "end_area_m2", area_m2, "m^2", "")
    put!("geometry", "fiducial_volume_m3", volume_m3, "m^3", "")
    put!("geometry", "gas_mass_kg", s.density_kg_m3 * volume_m3, "kg", "")
    put!("geometry", "xenon_mass_kg", s.xenon_density_kg_m3 * volume_m3, "kg", "")
    put!("geometry", "xe136_mass_kg", s.xenon_density_kg_m3 * volume_m3 * f136, "kg", "")
    put!("geometry", "helium_mass_kg", s.helium_density_kg_m3 * volume_m3, "kg", "")
    put!("geometry", "reference_pure_xenon_mass_kg", ref.density_kg_m3 * volume_m3, "kg", "same volume, enriched xenon at 15 bar, 290 K")
    effective = s.xenon_number_density_m3 + (s.number_density_m3 - s.xenon_number_density_m3) * 2 / 54
    put!("geometry", "track_length_relative_to_pure_xenon_15bar", ref.number_density_m3 / effective - 1, "",
        "inverse effective atom density, helium counted as 2/54 of a xenon atom; reference enriched xenon at 15 bar, 290 K")
    put!("geometry", "z_el_entry_cm", L_cm, "cm", "z of the EL entry, z = 0 at the pad plane")
    put!("geometry", "geant4_z_offset_cm", L_cm / 2, "cm", "z = z_Geant4 + offset; the Geant4 cylinder is centred on its origin")
    stack = geo["pad_gap_mm"] + geo["ifg_thickness_mm"] + L_cm * 10 + geo["transfer_gap_mm"] + geo["el_depth_mm"] + geo["guide_length_mm"]
    put!("geometry", "pads_to_sipm_plane_mm", stack, "mm", "pad gap + IFG + drift + transfer gap + EL + guides")

    E = dr["field_V_cm"]; v_e = dr["electron_drift_velocity_mm_us"]; tau = dr["electron_lifetime_ms"]
    t_e_ms = L_cm * 10 / v_e / 1000
    put!("drift", "drift_voltage_kV", E * L_cm / 1000, "kV", "E_d L")
    put!("drift", "reduced_field_V_cm_bar", E / P, "V/cm/bar", "plain pressure")
    put!("drift", "reduced_field_Td", E * 100 / s.number_density_m3 / TD_V_m2, "Td", "real-gas number density")
    put!("drift", "reduced_field_ideal_gas_Td", E * 100 / (P * 1e5 / (K_B_J_K * T)) / TD_V_m2, "Td", "ideal-gas number density at the card pressure and temperature")
    put!("drift", "electron_max_drift_time_ms", t_e_ms, "ms", "L / v_e")
    put!("drift", "electron_max_attachment_loss", 1 - exp(-t_e_ms / tau), "", "1 - exp(-t_max / tau)")
    DT = dr["electron_transverse_diffusion_mm_sqrtm"]; DL = dr["electron_longitudinal_diffusion_mm_sqrtm"]
    put!("drift", "electron_transverse_diffusion_mm_sqrtbar_sqrtcm", diffusion_card_units(DT, P), "mm sqrt(bar)/sqrt(cm)", "D sqrt(P) / 10, the unit of the production cards")
    put!("drift", "electron_longitudinal_diffusion_mm_sqrtbar_sqrtcm", diffusion_card_units(DL, P), "mm sqrt(bar)/sqrt(cm)", "D sqrt(P) / 10, the unit of the production cards")
    put!("drift", "electron_transverse_spread_mm", DT * sqrt(L_cm / 100), "mm", "D_T sqrt(L), full drift")
    put!("drift", "electron_longitudinal_spread_mm", DL * sqrt(L_cm / 100), "mm", "D_L sqrt(L), full drift")
    put!("drift", "electron_longitudinal_spread_us", DL * sqrt(L_cm / 100) / v_e, "us", "sigma_L / v_e")
    v_i = dr["ion_drift_velocity_cm_s"]
    put!("drift", "ion_drift_velocity_mm_ms", v_i / 100, "mm/ms", "the unit of the production cards")
    put!("drift", "ion_max_drift_time_s", L_cm / v_i, "s", "L / v_ion")
    put!("drift", "ion_thermal_spread_mm", sigma_t_ion_mm(Float32(T), Float32(L_cm), Float32(E)), "mm", "thermal, full drift, arXiv:2604.02443 Eq. (2.1)")
    put!("drift", "ion_reduced_mobility_cm2_V_s", v_i * s.number_density_relative_to_STP / E, "cm^2/(V s)", "K0 = v N / (N0 E), implied by the card velocity")
    put!("drift", "ion_bucket_length_mm", v_i / 100 * ions["ion_bucket_ms"], "mm", "ion velocity x ion bucket")
    put!("drift", "electron_bucket_length_mm", v_e * sipm["electron_bucket_us"], "mm", "electron velocity x electron bucket")

    pot = potentials_kV(d)
    for (k, note) in (("pads", ""), ("pad_face_ifg", ""), ("tpc_face_ifg", ""), ("drift_end", "TPC-facing IFG face + E_d L"),
            ("el_entry", "+ transfer-gap field x gap"), ("el_exit", "+ EL voltage"), ("sipm_plane", ""))
        put!("voltages", "potential_$(k)_kV", getfield(pot, Symbol(k)), "kV", note)
    end
    put!("voltages", "ifg_field_kV_cm", (v["tpc_face_ifg_V"] - v["pad_face_ifg_V"]) / 1000 / (geo["ifg_thickness_mm"] / 10), "kV/cm", "across the IFG")
    put!("voltages", "pad_gap_field_kV_cm", (v["pad_face_ifg_V"] - v["pads_V"]) / 1000 / (geo["pad_gap_mm"] / 10), "kV/cm", "between the IFG and the pads")
    grading = (pot.el_exit - pot.sipm_plane) / (geo["guide_length_mm"] / 10)
    put!("voltages", "grading_field_kV_cm", grading, "kV/cm", "EL exit potential over the guide length, docs/elm.tex")
    put!("voltages", "grading_reduced_field_kV_cm_bar", grading / P, "kV/cm/bar", "")

    d_cm = geo["el_depth_mm"] / 10
    E_el = v["el_voltage_kV"] / d_cm
    put!("el", "field_kV_cm", E_el, "kV/cm", "V_EL / d")
    put!("el", "reduced_field_kV_cm_bar", E_el / P, "kV/cm/bar", "plain pressure")
    put!("el", "photons_per_electron", el_photon_yield(Float32(el["yield_gain_photons_per_electron_kV"]),
        Float32(el["yield_threshold_kV_cm_bar"]), Float32(E_el / P), Float32(P), Float32(geo["el_depth_mm"])), "",
        "gain x (V_EL - threshold x P x d)")
    put!("el", "threshold_voltage_kV", el["yield_threshold_kV_cm_bar"] * P * d_cm, "kV", "EL voltage at zero yield")
    put!("el", "hole_count", area_m2 / (el["hole_pitch_mm"] / 1000)^2, "", "end area over the pitch squared")
    put!("el", "hole_open_fraction", pi * (el["hole_diameter_mm"] / 2)^2 / el["hole_pitch_mm"]^2, "", "hole area over cell area")

    put!("sipm", "sipm_count", area_m2 / (el["hole_pitch_mm"] / 1000)^2, "", "one per hole")
    put!("sipm", "coverage", (sipm["side_mm"] / el["hole_pitch_mm"])^2, "", "SiPM area over cell area")
    put!("ions", "pad_count", area_m2 / (ions["pad_pitch_mm"] / 1000)^2, "", "end area over the pad pitch squared")
    put!("ions", "pad_fill_fraction", (ions["pad_side_mm"] / ions["pad_pitch_mm"])^2, "", "collecting electrode over cell area")

    Q = get(en, "qbb_keV", QBB_XE136_keV)
    N = Q * 1000 / en["W_eV"]
    put!("energy", "electrons_at_qbb", N, "", "Q / W")
    put!("energy", "fano_fwhm_at_qbb_keV", 2.355 * sqrt(en["fano_factor"] * en["W_eV"] * Q * 1000) / 1000, "keV", "2.355 sqrt(F W Q)")
    put!("energy", "fano_fwhm_at_qbb_relative", 2.355 * sqrt(en["fano_factor"] / N), "", "")
    out
end

# ------------------------------------------------------- run configuration

"""Keys of a run card that do not affect the detector response: the later stages
and the response models. Everything else, with the detector and the overrides,
defines the detector-response stage."""
const STAGE_KEYS = ("response", "electronics", "hit_production", "ions", "electrons", "track_analysis",
                    "hit_output", "controls", "zs_controls")

"""
    electronics_card(r::ResolvedRun)::Dict

The response-model paths and the acquisition settings of the run in one flat
dictionary: the `[response]` and `[electronics]` sections merged, as
`prepare_electronics_run` reads them.
"""
electronics_card(r::ResolvedRun) = merge(Dict{String,Any}(get(r.run, "response", Dict{String,Any}())),
                                         Dict{String,Any}(get(r.run, "electronics", Dict{String,Any}())))

"""
    resolved_configuration(r::ResolvedRun, root; full=true)::Dict

The resolved configuration as a table: detector set identity and hashes, run
card path and hash, overrides, detector values after overrides, and the run
card. With `full = false` the stage sections (`STAGE_KEYS`) and the run card
identity are left out: what remains defines the detector-response stage, and
every run card sharing it produces the same text.
"""
function resolved_configuration(r::ResolvedRun, root::AbstractString; full::Bool = true)
    run = Dict{String,Any}(k => v for (k, v) in r.run if k != "overrides" && (full || !(k in STAGE_KEYS)))
    out = Dict{String,Any}(
        "schema" => "nausicaa.resolved_run.v1",
        "detector_set" => Dict("index" => relpath(r.set.index_path, root), "name" => r.set.name, "date" => r.set.date,
                               "sha256" => r.set.sha256),
        "overrides" => r.overrides,
        "detector" => r.detector,
        "run" => run)
    full && (out["run_card"] = Dict("path" => relpath(r.run_path, root), "sha256" => bytes2hex(sha256(read(r.run_path, String)))))
    out
end

resolved_text(r::ResolvedRun, root::AbstractString; full::Bool = true) =
    sprint(io -> TOML.print(io, resolved_configuration(r, root; full); sorted = true))

"""
    write_configuration(r::ResolvedRun, root; stage="detector")::String

Save `resolved.toml` for a stage of the run and return its directory relative
to `root`. The detector stage writes `results/readout/<run_label>/configuration/`
without the later-stage sections, so every run card sharing the same detector
response writes the same file; the hits stage writes the full configuration into
`<run_label>/<output_subdir>/configuration/`, and the electronics stage into its
own output folder. An existing file with the same content is left; a different
content is an error, so one folder holds one configuration.
"""
function write_configuration(r::ResolvedRun, root::AbstractString; stage::AbstractString = "detector")::String
    label = r.run["run_label"]
    if stage == "detector"
        reldir = joinpath("results/readout", label, "configuration"); text = resolved_text(r, root; full = false)
    elseif stage == "hits"
        reldir = joinpath("results/readout", label, r.run["hit_production"]["output_subdir"], "configuration"); text = resolved_text(r, root)
    elseif stage == "electronics"
        reldir = joinpath("results/readout", label, get(r.run["electronics"], "output_subdir", "electronics"), "configuration"); text = resolved_text(r, root)
    else
        throw(ArgumentError("Unknown stage $stage"))
    end
    dir = joinpath(root, reldir); mkpath(dir)
    path = joinpath(dir, "resolved.toml")
    if isfile(path)
        read(path, String) == text || error("$path exists with a different configuration")
    else
        write(path, text)
    end
    reldir
end

"""
    configuration_texts(r::ResolvedRun, root)::Dict{String,String}

The texts to store with an output: the run card, the resolved configuration,
the detector cards and the response-model cards named in `[response]`.
"""
function configuration_texts(r::ResolvedRun, root::AbstractString)
    texts = Dict{String,String}("run" => read(r.run_path, String), "resolved" => resolved_text(r, root))
    for (name, text) in r.set.texts
        texts["detector/" * name] = text
    end
    for (k, v) in get(r.run, "response", Dict{String,Any}())
        v isa AbstractString && endswith(v, ".toml") && (texts[k] = read(joinpath(root, v), String))
    end
    texts
end

end # module
