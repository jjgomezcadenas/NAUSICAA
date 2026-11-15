"""ITACA operating conditions: the card `config/itaca.toml` and the quantities
derived from it (gas density and masses, reduced fields, EL yield, electrode
potentials, drift times and spreads).

Load `src/event_analysis/EventAnalysis.jl` first; the EL yield and the thermal
ion spread come from there. Formulas and sources are given at each function.
Inputs and outputs are Float64 in the units named in the key or field.
"""
module ItacaConditions
using TOML, SHA
using ..EventAnalysis: el_photon_yield, sigma_t_ion_mm

export ItacaCard, load_itaca_card, enriched_xenon_molar_mass_g_mol, xenon_compressibility,
       GasState, gas_state, derived_conditions, QBB_XE136_keV

"Molar gas constant, J/(mol K) (CODATA 2018)."
const R_J_mol_K = 8.314462618
"Boltzmann constant, J/K (CODATA 2018)."
const K_B_J_K = 1.380649e-23
"Loschmidt number density at 273.15 K and 101325 Pa, 1/m^3 (CODATA 2018)."
const N0_STP_m3 = 2.686780111e25
"One Townsend, V m^2."
const TD_V_m2 = 1e-21
"Q-value of the Xe-136 double beta decay, keV (Redshaw et al., PRL 98 (2007) 053003)."
const QBB_XE136_keV = 2457.83
"Molar mass of helium-4, g/mol, as in simulation/nexus/source/geometries/ItacaChamber.cc."
const HELIUM_MOLAR_MASS_g_mol = 4.0026
"Atomic masses of the xenon isotopes, g/mol (AME2020)."
const XENON_ISOTOPE_MASS_g_mol = Dict("Xe124"=>123.905885, "Xe126"=>125.904298, "Xe128"=>127.903531,
    "Xe129"=>128.904781, "Xe130"=>129.903509, "Xe131"=>130.905084, "Xe132"=>131.904155,
    "Xe134"=>133.905395, "Xe136"=>135.907214)

"""The card as parsed, with its path and SHA-256 for provenance."""
struct ItacaCard
    path::String
    sha256::String
    text::String
    data::Dict{String,Any}
end

"""
    load_itaca_card(path)::ItacaCard

Read `config/itaca.toml`. Require the schema `nausicaa.itaca_conditions.v1`
and the sections gas, geometry, drift, el, elm, ion_readout and energy.
"""
function load_itaca_card(path::AbstractString)::ItacaCard
    text = read(path, String)
    data = TOML.parse(text)
    get(data, "schema", "") == "nausicaa.itaca_conditions.v1" || throw(ArgumentError("Unexpected schema in $path"))
    for section in ("gas", "geometry", "drift", "el", "elm", "ion_readout", "energy")
        haskey(data, section) || throw(ArgumentError("Missing [$section] in $path"))
    end
    ItacaCard(String(path), bytes2hex(sha256(text)), text, data)
end

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
T_r = T / T_c and P_r = P / P_c, the constants from `[gas.density_model]`
(NEXT DocDB 1183; simulation/nexus/source/geometries/ItacaChamber.cc).
"""
function xenon_compressibility(P_bar::Real, T_K::Real, model::AbstractDict)::Float64
    T_r = T_K / model["xenon_critical_temperature_K"]
    P_r = P_bar / model["xenon_critical_pressure_bar"]
    1 + (model["xenon_virial_a"] - model["xenon_virial_b"] / T_r^1.6) * P_r / T_r
end

"""Gas state of the Xe-He mixture at one pressure and temperature.

Additive partial volumes (Amagat) for mole fractions x_xe and x_he, xenon
with compressibility z, helium ideal:
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

"""
    derived_conditions(card::ItacaCard)::Vector{Tuple{String,String,Float64,String,String}}

Compute the derived quantities of the detector page. Each entry is
(group, name, value, unit, note). Groups: gas, geometry, drift, el, elm,
ion_readout, energy. Reference conditions and conventions:

- Pure-xenon reference for the track length: enriched xenon at 15 bar and
  290 K with the same density model. Track length scales with the inverse of
  the effective atom density, helium counted as 2/54 of a xenon atom (ratio of
  atomic numbers), the convention of docs/handoffs/itaca_detector_page.md.
- Reduced drift field in Td uses the real-gas number density; the value in
  V/cm/bar uses the plain pressure.
- Potentials: sensors at ground; TPC-facing IFG face at the card potential;
  drift end = IFG + E_d L; EL entry = drift end + E_g g over the transfer gap;
  EL exit = entry + V_EL; graded to the SiPM plane along the guides
  (docs/elm.tex, Section "Geometry and nominal potentials").
- Electron spreads sigma = D sqrt(L) with D in mm/sqrt(m); the longitudinal
  spread in time is sigma_L / v_e. Ion spread from `sigma_t_ion_mm`.
- Energy: electrons at Q_bb = Q / W; Fano-limited FWHM = 2.355 sqrt(F W Q).
"""
function derived_conditions(card::ItacaCard)
    d = card.data
    gas = d["gas"]; model = gas["density_model"]; geo = d["geometry"]; dr = d["drift"]
    el = d["el"]; elm = d["elm"]; ir = d["ion_readout"]; en = d["energy"]
    out = Tuple{String,String,Float64,String,String}[]
    put!(g, n, v, u, note="") = push!(out, (g, n, Float64(v), u, note))

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

    E = dr["field_V_cm"]; v_e = dr["electron_drift_velocity_mm_us"]; tau = dr["electron_lifetime_ms"]
    V_drift_kV = E * L_cm / 1000
    t_e_ms = L_cm * 10 / v_e / 1000
    put!("drift", "drift_voltage_kV", V_drift_kV, "kV", "E_d L")
    put!("drift", "reduced_field_V_cm_bar", E / P, "V/cm/bar", "plain pressure")
    put!("drift", "reduced_field_Td", E * 100 / s.number_density_m3 / TD_V_m2, "Td", "real-gas number density")
    put!("drift", "reduced_field_ideal_gas_Td", E * 100 / (P * 1e5 / (K_B_J_K * T)) / TD_V_m2, "Td", "ideal-gas number density at the card pressure and temperature")
    put!("drift", "electron_max_drift_time_ms", t_e_ms, "ms", "L / v_e")
    put!("drift", "electron_max_attachment_loss", 1 - exp(-t_e_ms / tau), "", "1 - exp(-t_max / tau)")
    sT = dr["electron_transverse_diffusion_mm_sqrtm"] * sqrt(L_cm / 100)
    sL = dr["electron_longitudinal_diffusion_mm_sqrtm"] * sqrt(L_cm / 100)
    put!("drift", "electron_transverse_spread_mm", sT, "mm", "D_T sqrt(L), full drift")
    put!("drift", "electron_longitudinal_spread_mm", sL, "mm", "D_L sqrt(L), full drift")
    put!("drift", "electron_longitudinal_spread_us", sL / v_e, "us", "sigma_L / v_e")
    put!("drift", "electron_transverse_spread_relative_error", dr["electron_transverse_diffusion_relative_error"], "", "statistical error of the Magboltz coefficient")
    put!("drift", "electron_longitudinal_spread_relative_error", dr["electron_longitudinal_diffusion_relative_error"], "", "statistical error of the Magboltz coefficient")
    put!("drift", "electron_drift_time_relative_error", dr["electron_drift_velocity_relative_error"], "", "statistical error of the Magboltz velocity")
    v_i = dr["ion_drift_velocity_cm_s"]
    put!("drift", "ion_max_drift_time_s", L_cm / v_i, "s", "L / v_ion")
    put!("drift", "ion_thermal_spread_mm", sigma_t_ion_mm(Float32(T), Float32(L_cm), Float32(E)), "mm", "thermal, full drift, arXiv:2604.02443 Eq. (2.1)")
    put!("drift", "ion_reduced_mobility_cm2_V_s", v_i * s.number_density_relative_to_STP / E, "cm^2/(V s)", "K0 = v N / (N0 E), implied by the card velocity")

    V_ifg_kV = ir["tpc_face_ifg_V"] / 1000
    V_end_kV = V_ifg_kV + V_drift_kV
    V_in_kV = V_end_kV + elm["transfer_gap_field_V_cm"] * elm["transfer_gap_mm"] / 10 / 1000
    V_out_kV = V_in_kV + el["voltage_kV"]
    put!("drift", "potential_pads_kV", ir["pads_V"] / 1000, "kV", "")
    put!("drift", "potential_pad_face_ifg_kV", ir["pad_face_ifg_V"] / 1000, "kV", "")
    put!("drift", "potential_tpc_face_ifg_kV", V_ifg_kV, "kV", "")
    put!("drift", "potential_drift_end_kV", V_end_kV, "kV", "TPC-facing IFG face + E_d L")
    put!("drift", "potential_el_entry_kV", V_in_kV, "kV", "+ transfer-gap field x gap")
    put!("drift", "potential_el_exit_kV", V_out_kV, "kV", "+ EL voltage")
    put!("drift", "potential_sipm_plane_kV", el["sipm_plane_V"] / 1000, "kV", "")
    put!("ion_readout", "ifg_field_kV_cm", (ir["tpc_face_ifg_V"] - ir["pad_face_ifg_V"]) / 1000 / (ir["ifg_thickness_mm"] / 10), "kV/cm", "across the IFG")
    put!("ion_readout", "pad_gap_field_kV_cm", (ir["pad_face_ifg_V"] - ir["pads_V"]) / 1000 / (ir["pad_gap_mm"] / 10), "kV/cm", "between the IFG and the pads")

    d_cm = el["depth_mm"] / 10
    E_el = el["voltage_kV"] / d_cm
    put!("el", "field_kV_cm", E_el, "kV/cm", "V_EL / d")
    put!("el", "reduced_field_kV_cm_bar", E_el / P, "kV/cm/bar", "plain pressure")
    put!("el", "photons_per_electron", el_photon_yield(Float32(el["yield_gain_photons_per_electron_kV"]),
        Float32(el["yield_threshold_kV_cm_bar"]), Float32(E_el / P), Float32(P), Float32(el["depth_mm"])), "",
        "gain x (V_EL - threshold x P x d)")
    put!("el", "threshold_voltage_kV", el["yield_threshold_kV_cm_bar"] * P * d_cm, "kV", "EL voltage at zero yield")
    grading = (V_out_kV - el["sipm_plane_V"] / 1000) / (elm["guide_length_mm"] / 10)
    put!("elm", "grading_field_kV_cm", grading, "kV/cm", "EL exit potential over the guide length, docs/elm.tex Eq. (grading)")
    put!("elm", "grading_reduced_field_kV_cm_bar", grading / P, "kV/cm/bar", "")
    put!("elm", "grading_field_card_difference_kV_cm", elm["grading_field_kV_cm"] - grading, "kV/cm", "card value minus the value computed from the potentials")
    put!("elm", "hole_count", area_m2 / (elm["hole_pitch_mm"] / 1000)^2, "", "end area over the pitch squared; one SiPM per hole")
    put!("elm", "hole_open_fraction", pi * (elm["hole_diameter_mm"] / 2)^2 / elm["hole_pitch_mm"]^2, "", "hole area over cell area")
    put!("elm", "sipm_coverage", (elm["sipm_side_mm"] / elm["hole_pitch_mm"])^2, "", "SiPM area over cell area")

    N = QBB_XE136_keV * 1000 / en["W_eV"]
    put!("energy", "qbb_keV", QBB_XE136_keV, "keV", "Xe-136, Redshaw et al. (2007)")
    put!("energy", "electrons_at_qbb", N, "", "Q / W")
    put!("energy", "fano_fwhm_at_qbb_keV", 2.355 * sqrt(en["fano_factor"] * en["W_eV"] * QBB_XE136_keV * 1000) / 1000, "keV", "2.355 sqrt(F W Q)")
    put!("energy", "fano_fwhm_at_qbb_relative", 2.355 * sqrt(en["fano_factor"] / N), "", "")
    out
end

end # module
