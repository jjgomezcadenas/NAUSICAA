"""eITACA detector configuration: the detector set under `config/eitaca/`, its
variants, the run cards that reference them, the optics without light guides,
and the quantities derived from the resolved configuration.

Three layers, as in `src/iitaca/IItaca.jl`:

1. Detector set: an index (`config/eitaca/eitaca.toml`) naming nine cards,
   gas, geometry, drift, voltages, el, sipm, energy, trigger, reconstruction.
   Design values only; nothing derived is typed in them.
2. Variant: a card with `[overrides]`, dotted keys into the set
   (`config/eitaca/variants/`).
3. Run card: `detector` (the index), optionally `variant`, and its own
   `[overrides]`, applied after those of the variant.

`resolve_detector` returns the detector values after the overrides and checks
them; `derived_conditions` computes the physics of a resolved detector.
`light_acceptance`, `el_photons_per_electron` and `electron_spreads_mm` give
the numbers the electron response reads.

Load `src/event_analysis/EventAnalysis.jl` and `src/itaca/ItacaConditions.jl`
first: the axial hole acceptance, the guide optics, the EL yield law, the gas
state and the override rules come from there. Units are in the key names.
"""
module EItaca
using TOML, SHA
using ..EventAnalysis: el_photon_yield, axial_light_acceptance, elm_light_acceptance
using ..ItacaConditions: DetectorSet, apply_overrides, flatten, gas_state, enriched_xenon_molar_mass_g_mol,
                         XENON_ISOTOPE_MASS_g_mol, TD_V_m2

export SECTIONS, load_eitaca_set, load_variant, ResolvedDetector, resolve_detector, check_detector,
       ResolvedRun, resolve_run, resolved_text, event_seed, direct_light_acceptance, light_acceptance,
       el_reduced_field_kV_cm_bar, el_photons_per_electron, electron_diffusion_mm_sqrtm, electron_spreads_mm,
       asic_limit_pe, potentials_kV, derived_conditions, response_cards

const SET_SCHEMA = "nausicaa.eitaca_set.v1"
const VARIANT_SCHEMA = "nausicaa.eitaca_variant.v1"
const RUN_SCHEMA = "nausicaa.eitaca_run.v1"
const SECTIONS = ("gas", "geometry", "drift", "voltages", "el", "sipm", "energy", "trigger", "reconstruction", "generation")
const CARD_SCHEMAS = Dict(s => "nausicaa.eitaca_$s.v1" for s in SECTIONS)
const OPTICS = ("direct_tpb_on_sipm", "light_guides")
"Elementary charge, C."
const E_CHARGE_C = 1.602176634e-19

# ---------------------------------------------------------------- detector set

"""
    load_eitaca_set(index_path)::DetectorSet

Read the index and its ten cards. Require the set schema, every section of
`SECTIONS` exactly once, and each card's own schema.
"""
function load_eitaca_set(index_path::AbstractString)::DetectorSet
    text = read(index_path, String); index = TOML.parse(text)
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
to check the SiPM and ASIC cards; without it that check is skipped.
"""
function resolve_detector(index_path::AbstractString; variant = nothing, overrides = Dict{String,Any}(), root = nothing)::ResolvedDetector
    set = load_eitaca_set(index_path)
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
const RUN_KEYS = ("detector", "run_label", "sample", "hits", "particles", "seed", "electronics_seed")

"""
    resolve_run(run_path; root)::ResolvedRun

Load a run card (schema `nausicaa.eitaca_run.v1`): the detector index, an
optional variant card and `[overrides]`, applied after those of the variant.
The other keys are the run parameters: `run_label`, `sample`, the input files
`hits` and `particles`, `seed` (pairs, transport, EL and optics),
`electronics_seed` (SiPM and ASIC), and optionally `max_events` (0: all),
`max_vector_bytes` and `generated_events_per_batch`. Paths are relative to
`root`. A detector key at the top level of the run card is refused.
"""
function resolve_run(run_path::AbstractString; root::AbstractString)::ResolvedRun
    run = TOML.parsefile(run_path)
    get(run, "schema", "") == RUN_SCHEMA || throw(ArgumentError("Expected schema $RUN_SCHEMA in $run_path"))
    for k in RUN_KEYS
        haskey(run, k) || throw(ArgumentError("Run card without '$k': $run_path"))
    end
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
        "schema" => "nausicaa.eitaca_resolved_run.v1",
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
`"framed-v1|seed|event_number|part"`, the rule of `IItaca.event_seed` and of
the hybrid production (`seedfor`). With the seeds of the hybrid run cards, the
streams `pairs`, `electrons` and `electron_grid` repeat the primary pairs,
diffusion draws and SiPM-lattice offset of the hybrid production.
"""
event_seed(seed, event_number, part) = parse(UInt64, bytes2hex(sha256("framed-v1|$seed|$event_number|$part"))[1:16]; base = 16)

"""
    response_cards(detector, root) -> (sipm, asic)

The SiPM response card and the ASIC card named in the sipm card, parsed.
"""
response_cards(d::AbstractDict, root::AbstractString) =
    (TOML.parsefile(joinpath(root, d["sipm"]["response_card"])), TOML.parsefile(joinpath(root, d["sipm"]["asic_card"])))

"""
    check_detector(detector; root=nothing)

Consistency of a resolved detector:

- positive pressure, temperature, lengths, pitch, fields, velocity, lifetime;
- the reduced drift field equals the field of the transport coefficients
  (`drift.transport_reduced_field_V_cm_bar`) within 1 %;
- the anode above the gate (electrons drift to the anode);
- the hole fits the pitch; the SiPM square covers the hole and fits the pitch;
- a known optics model;
- with `root`: the response card has the PDE and SiPM side of the sipm card,
  the ASIC bucket equals the electron bucket, and the illuminated microcells
  are the microcells under the hole disk, round(M pi r^2 / side^2);
- known frame placement and non-negative trigger margins.
"""
function check_detector(d::AbstractDict; root = nothing)
    gas = d["gas"]; geo = d["geometry"]; dr = d["drift"]; v = d["voltages"]; el = d["el"]; si = d["sipm"]; tr = d["trigger"]
    positive = (("gas.pressure_bar", gas["pressure_bar"]), ("gas.temperature_K", gas["temperature_K"]),
        ("geometry.fiducial_diameter_cm", geo["fiducial_diameter_cm"]), ("geometry.fiducial_length_cm", geo["fiducial_length_cm"]),
        ("geometry.el_depth_mm", geo["el_depth_mm"]), ("drift.field_V_cm", dr["field_V_cm"]),
        ("drift.electron_drift_velocity_mm_us", dr["electron_drift_velocity_mm_us"]), ("drift.electron_lifetime_ms", dr["electron_lifetime_ms"]),
        ("el.hole_pitch_mm", el["hole_pitch_mm"]), ("el.hole_diameter_mm", el["hole_diameter_mm"]), ("sipm.side_mm", si["side_mm"]),
        ("sipm.electron_bucket_us", si["electron_bucket_us"]), ("trigger.frame_side_mm", tr["frame_side_mm"]))
    for (key, value) in positive
        value > 0 || throw(ArgumentError("$key must be positive, found $value"))
    end
    isapprox(dr["field_V_cm"] / gas["pressure_bar"], dr["transport_reduced_field_V_cm_bar"]; rtol = 0.01) ||
        throw(ArgumentError("Drift field $(dr["field_V_cm"] / gas["pressure_bar"]) V/cm/bar differs from the field of the transport coefficients, $(dr["transport_reduced_field_V_cm_bar"])"))
    v["el_anode_V"] > v["gate_V"] || throw(ArgumentError("The EL anode must be above the gate"))
    el["hole_diameter_mm"] < el["hole_pitch_mm"] || throw(ArgumentError("The hole exceeds the pitch"))
    si["side_mm"] >= el["hole_diameter_mm"] || throw(ArgumentError("The SiPM square does not cover the hole"))
    si["side_mm"] <= el["hole_pitch_mm"] || throw(ArgumentError("The SiPM exceeds the pitch"))
    el["optics"] in OPTICS || throw(ArgumentError("Unknown el.optics '$(el["optics"])'"))
    el["photon_band"] == "blue_tpb" || throw(ArgumentError("Both optics models give blue TPB photons"))
    if root !== nothing
        sc, ac = response_cards(d, root)
        sc["pde"] == si["pde"] || throw(ArgumentError("PDE of $(si["response_card"]) differs from sipm.pde"))
        isapprox(sc["side_m"] * 1000, si["side_mm"]; rtol = 1e-9) || throw(ArgumentError("SiPM side of the response card differs from sipm.side_mm"))
        isapprox(ac["bucket_width_s"], si["electron_bucket_us"] * 1e-6; rtol = 1e-9) || throw(ArgumentError("ASIC bucket differs from sipm.electron_bucket_us"))
        r = el["hole_diameter_mm"] / 2
        round(Int, sc["microcell_count"] * pi * r^2 / si["side_mm"]^2) == sc["illuminated_microcell_count"] ||
            throw(ArgumentError("Illuminated microcells of the response card differ from the hole disk"))
    end
    tr["frame_placement"] == "inside_detector" || throw(ArgumentError("Unsupported trigger.frame_placement '$(tr["frame_placement"])'"))
    tr["frame_track"] == "primary_and_delta_electrons" || throw(ArgumentError("Unsupported trigger.frame_track '$(tr["frame_track"])'"))
    0 <= tr["simulation_margin_mm"] && 0 <= tr["el_exclusion_mm"] && 0 <= tr["gamma_energy_cut_keV"] ||
        throw(ArgumentError("Trigger margins and cuts must not be negative"))
    rc = d["reconstruction"]
    rc["rebin_factor"] isa Integer && rc["rebin_factor"] > 0 || throw(ArgumentError("reconstruction.rebin_factor must be a positive integer"))
    nothing
end

# --------------------------------------------------------------------- physics

"""
    direct_light_acceptance(radius_mm, depth_mm; tpb_efficiency, forward_fraction)::Float32

Probability that one EL photon emitted in a hole produces a blue photon
entering the SiPM, without light guides. Assumptions: VUV emission isotropic
and uniform along the hole axis; absorbing hole walls; the SiPM, with TPB on
its window, directly against the hole exit, its square covering the hole, so
every VUV photon leaving the exit reaches the TPB; TPB converts with
`tpb_efficiency` and re-emits isotropically, `forward_fraction` towards the
SiPM. Then

    acceptance = axial_light_acceptance(r, d) x tpb_efficiency x forward_fraction,

the factors of `elm_light_acceptance` (docs/elm.tex) without the guide capture
1 - n_gap / n_guide. The PDE is applied by the SiPM response.
"""
function direct_light_acceptance(radius_mm::Float32, depth_mm::Float32; tpb_efficiency::Float32, forward_fraction::Float32)::Float32
    all(x -> isfinite(x) && 0f0 <= x <= 1f0, (tpb_efficiency, forward_fraction)) ||
        throw(ArgumentError("Optical probabilities must lie in [0,1]"))
    axial_light_acceptance(radius_mm, depth_mm) * tpb_efficiency * forward_fraction
end

"""
    light_acceptance(detector)::Float32

Blue photons at the SiPM per EL photon for the optics of `el.optics`:
`direct_light_acceptance` (no guides) or `elm_light_acceptance` (guides, with
the refractive indices of the el card).
"""
function light_acceptance(d::AbstractDict)::Float32
    el = d["el"]; F = Float32; r = F(el["hole_diameter_mm"] / 2); depth = F(d["geometry"]["el_depth_mm"])
    if el["optics"] == "direct_tpb_on_sipm"
        direct_light_acceptance(r, depth; tpb_efficiency = F(el["tpb_conversion_efficiency"]), forward_fraction = F(el["forward_fraction"]))
    elseif el["optics"] == "light_guides"
        elm_light_acceptance(r, depth; tpb_efficiency = F(el["tpb_conversion_efficiency"]), forward_fraction = F(el["forward_fraction"]),
            guide_index = F(el["guide_refractive_index"]), gap_index = F(el["gap_refractive_index"]))
    else
        throw(ArgumentError("Unknown el.optics '$(el["optics"])'"))
    end
end

"EL reduced field, (anode - gate) / d / P, in kV/cm/bar."
el_reduced_field_kV_cm_bar(d::AbstractDict) =
    (d["voltages"]["el_anode_V"] - d["voltages"]["gate_V"]) / 1000 / (d["geometry"]["el_depth_mm"] / 10) / d["gas"]["pressure_bar"]

"""
    el_photons_per_electron(detector)::Float32

Mean EL photons per electron: gain x (E/p - threshold) x P x d (`el_photon_yield`).
"""
el_photons_per_electron(d::AbstractDict)::Float32 = (el = d["el"]; F = Float32;
    el_photon_yield(F(el["yield_gain_photons_per_electron_kV"]), F(el["yield_threshold_kV_cm_bar"]), F(el_reduced_field_kV_cm_bar(d)),
        F(d["gas"]["pressure_bar"]), F(d["geometry"]["el_depth_mm"])))

"""
    electron_diffusion_mm_sqrtm(detector) -> (transverse, longitudinal)

D = D* / sqrt(P) at the card pressure, from um sqrt(bar)/sqrt(cm) to mm/sqrt(m):
D[mm/sqrt(m)] = D*[um sqrt(bar)/sqrt(cm)] / sqrt(P[bar]) x 10 / 1000.
"""
function electron_diffusion_mm_sqrtm(d::AbstractDict)
    dr = d["drift"]; P = d["gas"]["pressure_bar"]
    k = 10 / 1000 / sqrt(P)
    (transverse = dr["electron_reduced_transverse_diffusion_um_sqrtbar_sqrtcm"] * k,
     longitudinal = dr["electron_reduced_longitudinal_diffusion_um_sqrtbar_sqrtcm"] * k)
end

"""
    electron_spreads_mm(detector, drift_cm) -> (transverse, longitudinal)

sigma = D sqrt(L) for a drift of `drift_cm`, in mm; L in m.
"""
electron_spreads_mm(d::AbstractDict, drift_cm::Real) = (D = electron_diffusion_mm_sqrtm(d); L = sqrt(drift_cm / 100);
    (transverse = D.transverse * L, longitudinal = D.longitudinal * L))

"""
    asic_limit_pe(sipm_card, asic_card)::Float64

Upper clipping of the ASIC in fully recovered avalanches per bucket:
current_limit x bucket / (gain x e).
"""
asic_limit_pe(sc::AbstractDict, ac::AbstractDict) = ac["current_limit_A"] * ac["bucket_width_s"] / (sc["gain"] * E_CHARGE_C)

"""Potential ladder in kV: cathode = gate - E_d L; gate; anode; SiPM plane."""
function potentials_kV(d::AbstractDict)
    v = d["voltages"]; L = d["geometry"]["fiducial_length_cm"]
    (cathode = (v["gate_V"] - d["drift"]["field_V_cm"] * L) / 1000, gate = v["gate_V"] / 1000, el_anode = v["el_anode_V"] / 1000,
     sipm_plane = v["sipm_plane_V"] / 1000)
end

"""
    derived_conditions(detector; root=nothing)::Vector{Tuple{String,String,Float64,String,String}}

Compute the derived quantities of a resolved detector. Each entry is
(group, name, value, unit, note). With `root` the SiPM and ASIC cards are read
for the saturation scales:

- ASIC limit in PE per bucket (`asic_limit_pe`).
- Recovery: with recovery time tau and M illuminated microcells, a bucket of
  dt holds about M dt / tau recovered-cell equivalents; the mean inverse of
  the reconstruction, m / (1 - alpha m) with alpha = tau / (M dt)
  (`reco_electron_signal`), diverges at m = 1 / alpha.
"""
function derived_conditions(d::AbstractDict; root = nothing)
    gas = d["gas"]; geo = d["geometry"]; dr = d["drift"]; el = d["el"]; si = d["sipm"]; en = d["energy"]
    out = Tuple{String,String,Float64,String,String}[]
    put!(g, n, x, u, note = "") = push!(out, (g, n, Float64(x), u, note))
    M_xe = enriched_xenon_molar_mass_g_mol(gas["xenon_isotopic_fractions"])
    f136 = gas["xenon_isotopic_fractions"]["Xe136"] * XENON_ISOTOPE_MASS_g_mol["Xe136"] / M_xe
    P = gas["pressure_bar"]; T = gas["temperature_K"]
    s = gas_state(P, T, gas["helium_mole_fraction"], M_xe, gas["density_model"])
    put!("gas", "xenon_molar_mass_g_mol", M_xe, "g/mol", "enriched xenon, isotopic fractions of the card")
    put!("gas", "xe136_mass_fraction", f136, "", "")
    put!("gas", "density_kg_m3", s.density_kg_m3, "kg/m^3", "mixture")
    put!("gas", "xenon_density_kg_m3", s.xenon_density_kg_m3, "kg/m^3", "xenon partial density")
    put!("gas", "number_density_m3", s.number_density_m3, "1/m^3", "")

    radius_cm = geo["fiducial_diameter_cm"] / 2; L_cm = geo["fiducial_length_cm"]
    area_m2 = pi * (radius_cm / 100)^2; volume_m3 = area_m2 * L_cm / 100
    put!("geometry", "fiducial_radius_cm", radius_cm, "cm", "")
    put!("geometry", "fiducial_volume_m3", volume_m3, "m^3", "")
    put!("geometry", "xenon_mass_kg", s.xenon_density_kg_m3 * volume_m3, "kg", "xenon density x volume")
    put!("geometry", "xe136_mass_kg", s.xenon_density_kg_m3 * volume_m3 * f136, "kg", "xenon mass x Xe-136 mass fraction")

    E = dr["field_V_cm"]; v_e = dr["electron_drift_velocity_mm_us"]; D = electron_diffusion_mm_sqrtm(d); sp = electron_spreads_mm(d, L_cm)
    t_ms = L_cm * 10 / v_e / 1000
    put!("drift", "drift_voltage_kV", E * L_cm / 1000, "kV", "E_d L")
    put!("drift", "reduced_field_V_cm_bar", E / P, "V/cm/bar", "plain pressure")
    put!("drift", "reduced_field_Td", E * 100 / s.number_density_m3 / TD_V_m2, "Td", "real-gas number density")
    put!("drift", "electron_transverse_diffusion_mm_sqrtm", D.transverse, "mm/sqrt(m)", "D*_T / sqrt(P)")
    put!("drift", "electron_longitudinal_diffusion_mm_sqrtm", D.longitudinal, "mm/sqrt(m)", "D*_L / sqrt(P)")
    put!("drift", "electron_transverse_spread_mm", sp.transverse, "mm", "D_T sqrt(L), full drift")
    put!("drift", "electron_longitudinal_spread_mm", sp.longitudinal, "mm", "D_L sqrt(L), full drift")
    put!("drift", "electron_max_drift_time_ms", t_ms, "ms", "L / v_e")
    put!("drift", "electron_max_attachment_loss", 1 - exp(-t_ms / dr["electron_lifetime_ms"]), "", "1 - exp(-t_max / tau)")
    put!("drift", "electron_bucket_length_mm", v_e * si["electron_bucket_us"], "mm", "v_e x bucket")

    pot = potentials_kV(d)
    for (k, note) in (("cathode", "gate - E_d L"), ("gate", "EL entry, end of the drift"), ("el_anode", "EL exit"), ("sipm_plane", ""))
        put!("voltages", "potential_$(k)_kV", getfield(pot, Symbol(k)), "kV", note)
    end
    Er = el_reduced_field_kV_cm_bar(d)
    put!("el", "el_voltage_kV", pot.el_anode - pot.gate, "kV", "anode - gate")
    put!("el", "field_kV_cm", Er * P, "kV/cm", "(anode - gate) / d")
    put!("el", "reduced_field_kV_cm_bar", Er, "kV/cm/bar", "plain pressure")
    y = el_photons_per_electron(d); a = light_acceptance(d)
    put!("el", "photons_per_electron", y, "", "gain x (E/p - threshold) x P x d")
    put!("el", "axial_hole_acceptance", axial_light_acceptance(Float32(el["hole_diameter_mm"] / 2), Float32(geo["el_depth_mm"])), "",
        "VUV photons leaving the hole exit per EL photon, absorbing walls")
    put!("el", "blue_photons_per_el_photon", a, "", el["optics"] == "direct_tpb_on_sipm" ? "axial acceptance x TPB x forward fraction" :
        "axial acceptance x TPB x forward fraction x (1 - n_gap/n_guide)")
    put!("el", "photoelectrons_per_electron", y * a * si["pde"], "", "photons per electron x acceptance x PDE")
    put!("el", "channel_count", area_m2 / (el["hole_pitch_mm"] / 1000)^2, "", "end area over the pitch squared; one SiPM per hole")
    put!("el", "hole_open_fraction", pi * (el["hole_diameter_mm"] / 2)^2 / el["hole_pitch_mm"]^2, "", "hole area over cell area")
    N = en["qbb_keV"] * 1000 / en["W_eV"]
    put!("energy", "electrons_at_qbb", N, "", "Q / W")
    put!("energy", "photoelectrons_at_qbb", N * y * a * si["pde"], "", "before attachment")
    if root !== nothing
        sc, ac = response_cards(d, root)
        limit = asic_limit_pe(sc, ac); M = sc["illuminated_microcell_count"]; dt = ac["bucket_width_s"]; tau = sc["recovery_s"]
        put!("sipm", "asic_limit_pe_per_bucket", limit, "PE", "current limit x bucket / (gain x e)")
        put!("sipm", "asic_noise_rms_pe", ac["noise_rms_pe"], "PE", "per channel and bucket")
        put!("sipm", "illuminated_microcells", M, "", "microcells under the 2.5 mm hole disk")
        put!("sipm", "recovered_cell_equivalents_per_bucket", M * dt / tau, "", "M dt / tau")
        put!("sipm", "recovery_occupancy_at_asic_limit", limit * tau / (M * dt), "", "alpha m at the ASIC limit, alpha = tau / (M dt)")
        put!("sipm", "electrons_per_bucket_at_asic_limit", limit / (y * a * si["pde"]), "", "drift electrons in one hole and one bucket that reach the limit")
    end
    out
end

end # module
