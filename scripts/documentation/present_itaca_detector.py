"""Write the ITACA detector page and the configuration architecture page.

Inputs: results/itaca_detector/itaca_conditions.h5 (default detector) or the
itaca_conditions.h5 of a run's configuration folder, produced by
scripts/documentation/calculate_itaca_conditions.jl. This script only formats:
no physical quantity is computed here.

  python scripts/documentation/present_itaca_detector.py                 # default detector and architecture pages
  python scripts/documentation/present_itaca_detector.py <conditions.h5> <page.html>   # one run's page
"""
from __future__ import annotations

import html
import math
import sys
import tomllib
from pathlib import Path

import h5py

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results/itaca_detector/itaca_conditions.h5"
TARGET = ROOT / "notes/html/itaca_detector.html"
ARCHITECTURE = ROOT / "notes/html/detector_configuration.html"
EXAMPLE_RUN = ROOT / "config/runs/bb0nu_vh1_legacy_conditions_pilot100.toml"

STYLE = ("body{font:17px/1.65 system-ui;color:#203343;background:#f3f6f8;margin:0}header{background:#17394b;color:white;"
         "padding:30px max(24px,calc((100vw - 1100px)/2))}header a{color:#d4ecf6}main{max-width:1100px;margin:auto;padding:28px 24px}"
         "section{background:white;border:1px solid #d8e3e8;border-radius:10px;padding:25px;margin-bottom:24px;scroll-margin-top:15px}"
         "h1,h2,h3{line-height:1.25}h2{color:#185b83}a{color:#185b83}pre{background:#edf3f6;padding:18px;overflow:auto;font:14px/1.6 ui-monospace,monospace}"
         "code{font-family:ui-monospace,monospace;font-size:.9em}table{border-collapse:collapse;width:100%;font-size:15px}"
         "th,td{text-align:left;padding:10px;border-bottom:1px solid #d8e3e8;vertical-align:top}th{background:#edf3f6}td:nth-child(2){white-space:nowrap}"
         ".over{background:#fff4d8}.over td:nth-child(2)::after{content:' (overridden)';color:#985a00;font-size:.85em}"
         "nav{display:flex;flex-wrap:wrap;gap:10px 20px;margin-bottom:25px}footer{padding:20px;color:#526776}"
         "svg text{font:14px system-ui;fill:#203343}svg .box{fill:#fff;stroke:#185b83;stroke-width:1.5}svg .store{fill:#edf3f6;stroke:#526776}"
         "svg .arrow{stroke:#526776;stroke-width:1.5;fill:none;marker-end:url(#a)}"
         " @media(max-width:650px){main{padding:15px}section{padding:18px}table{display:block;overflow:auto}}")


def decode(v):
    return v.decode() if isinstance(v, bytes) else v


def read(path: Path):
    with h5py.File(path, "r") as f:
        meta = {k: str(decode(v)) for k, v in f.attrs.items()}
        cards = {k: str(decode(ds[()])) for k, ds in f["cards"].items()}
        inputs = {k: (str(decode(ds[()])) if ds.dtype.kind in "OSU" else float(ds[()])) for k, ds in f["inputs"].items()}
        overrides = {k: (str(decode(ds[()])) if ds.dtype.kind in "OSU" else float(ds[()])) for k, ds in f["overrides"].items()}
        derived, units, notes = {}, {}, {}
        for k, ds in f["derived"].items():
            derived[k] = float(ds[()])
            units[k] = str(decode(ds.attrs["unit"]))
            notes[k] = str(decode(ds.attrs["note"]))
    defaults = {}
    for section, text in cards.items():
        if section == "index":
            continue
        d = tomllib.loads(text)
        d.pop("schema", None)
        defaults.update(flatten(d, section))
    return meta, cards, inputs, overrides, derived, units, notes, defaults


def flatten(d, prefix):
    out = {}
    for k, v in d.items():
        name = f"{prefix}.{k}"
        if isinstance(v, dict) and k != "xenon_isotopic_fractions":
            out.update(flatten(v, name))
        else:
            out[name] = v
    return out


def num(x: float, digits: int = 4) -> str:
    if x == 0:
        return "0"
    if abs(x) >= 1e6 or abs(x) < 1e-3:
        return f"{x:.{digits - 1}e}".replace("e+0", "e").replace("e+", "e")
    magnitude = math.floor(math.log10(abs(x)))
    decimals = max(digits - 1 - magnitude, 0)
    s = f"{x:,.{decimals}f}"
    return s.rstrip("0").rstrip(".") if "." in s else s


def pct(x: float, digits: int = 2) -> str:
    return f"{100 * x:.{digits}g} %"


UNIT_HTML = {"kg/m^3": "kg/m³", "m^3": "m³", "m^2": "m²", "1/m^3": "m⁻³", "cm^2/(V s)": "cm²/(V·s)", "us": "µs",
             "mm/sqrt(m)": "mm/√m", "mm sqrt(bar)/sqrt(cm)": "mm·√bar/√cm"}


def unit(u: str) -> str:
    return UNIT_HTML.get(u, u)


def unit_from_key(key: str) -> str:
    for suffix, u in (("_V_cm", "V/cm"), ("_kV_cm_bar", "kV/cm/bar"), ("_kV_cm", "kV/cm"), ("_mm_us", "mm/µs"), ("_cm_s", "cm/s"),
                      ("_per_electron_kV", "photons/(electron·kV)"), ("_mm_sqrtm", "mm/√m"), ("_mm_sqrtbar_sqrtcm", "mm·√bar/√cm"),
                      ("_Hz", "Hz"), ("_um", "µm"), ("_us", "µs"), ("_ns", "ns"), ("_ms", "ms"), ("_bar", "bar"), ("_K", "K"),
                      ("_cm", "cm"), ("_mm", "mm"), ("_kV", "kV"), ("_V", "V"), ("_eV", "eV"), ("_keV", "keV")):
        if key.endswith(suffix):
            return u
    return ""


def table(rows, head=("Quantity", "Value", "Unit", "Definition")) -> str:
    out = ["<table><tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in head) + "</tr>"]
    for r in rows:
        cls = ' class="over"' if len(r) > 4 and r[4] else ""
        out.append(f"<tr{cls}>" + "".join(f"<td>{c}</td>" for c in r[:4]) + "</tr>")
    out.append("</table>")
    return "".join(out)


def page(title: str, subtitle: str, intro: str, sections, footer: str) -> str:
    nav = "".join(f'<a href="#{sid}">{html.escape(t.split(". ", 1)[-1])}</a>' for sid, t, _ in sections)
    body = "".join(f'<section id="{sid}"><h2>{t}</h2>{c}</section>' for sid, t, c in sections)
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>{html.escape(title)}</title><style>{STYLE}</style></head><body><header><a href="index.html">Documentation index</a>'
            f'<h1>{html.escape(title)}</h1><p>{subtitle}</p></header><main><nav>{nav}</nav>{intro}{body}</main><footer>{footer}</footer></body></html>\n')


def detector_page(source: Path, target: Path) -> None:
    meta, cards, inp, overrides, der, units, notes, defaults = read(source)

    def I(key, label, fmt=num, note=""):
        v = inp[key]
        shown = fmt(v) if isinstance(v, float) else html.escape(str(v))
        over = key in overrides
        if over:
            d = defaults.get(key)
            note = (note + "; " if note else "") + f"default {fmt(d) if isinstance(d, (int, float)) else d}"
        return (label, shown, unit_from_key(key), html.escape(note), over)

    def D(key, label, fmt=num, note=None):
        return (label, fmt(der[key]), unit(units[key]), html.escape(note if note is not None else notes[key]), False)

    rel = lambda k: f"± {pct(inp[k])}"  # noqa: E731
    iso = ", ".join(f"{k[2:]}: {pct(inp[f'gas.xenon_isotopic_fractions.{k}'], 4)}" for k in
                    ("Xe129", "Xe130", "Xe131", "Xe132", "Xe134", "Xe136")) if "gas.xenon_isotopic_fractions.Xe136" in inp else inp.get("gas.xenon_isotopic_fractions", "")
    sections = []
    sections.append(("gas", "1. Gas", f"""
<p>Xenon enriched in Xe-136 with helium. Densities follow the model of the Geant4 geometry: xenon as a real gas with a second-virial compressibility, helium ideal, additive partial volumes.</p>
{table([
    I("gas.mixture", "Mixture"), I("gas.helium_mole_fraction", "Helium mole fraction", pct),
    ("Xenon isotopic fractions", html.escape(str(iso)), "", "NEXUS GXeEnriched", False),
    I("gas.pressure_bar", "Pressure"), I("gas.temperature_K", "Temperature"),
    D("gas.xenon_molar_mass_g_mol", "Molar mass of the enriched xenon"), D("gas.xe136_mass_fraction", "Xe-136 mass fraction", pct),
    D("gas.xenon_compressibility", "Xenon compressibility z"), D("gas.density_kg_m3", "Gas density"),
    D("gas.xenon_density_kg_m3", "Xenon density"), D("gas.helium_density_kg_m3", "Helium density"),
    D("gas.number_density_m3", "Number density", lambda x: f"{x:.3e}"), D("gas.number_density_relative_to_STP", "Number density relative to STP, N/N₀"),
    D("gas.ideal_equivalent_pressure_bar", "Ideal-gas equivalent pressure"),
])}
<p>Density model: z = 1 + (a − b / T<sub>r</sub><sup>1.6</sup>) P<sub>r</sub> / T<sub>r</sub>, with a = {num(inp['gas.density_model.xenon_virial_a'])}, b = {num(inp['gas.density_model.xenon_virial_b'])}, T<sub>r</sub> = T / {num(inp['gas.density_model.xenon_critical_temperature_K'])} K, P<sub>r</sub> = P / {num(inp['gas.density_model.xenon_critical_pressure_bar'])} bar; ρ = P (x<sub>Xe</sub> M<sub>Xe</sub> + x<sub>He</sub> M<sub>He</sub>) / ((x<sub>Xe</sub> z + x<sub>He</sub>) R T). NEXT DocDB 1183, as in the Geant4 geometry <code>ItacaChamber.cc</code>.</p>"""))

    stack = [
        ("SiPM plane", "—", num(der["voltages.potential_sipm_plane_kV"]), "sensors at ground"),
        ("Light guides", num(inp["geometry.guide_length_mm"]), f"{num(der['voltages.potential_el_exit_kV'])} → {num(der['voltages.potential_sipm_plane_kV'])}", f"average grading field {num(der['voltages.grading_field_kV_cm'])} kV/cm"),
        ("EL region", num(inp["geometry.el_depth_mm"]), f"{num(der['voltages.potential_el_entry_kV'])} → {num(der['voltages.potential_el_exit_kV'])}", f"{num(inp['voltages.el_voltage_kV'])} kV across, field {num(der['el.field_kV_cm'])} kV/cm"),
        ("Transfer gap", num(inp["geometry.transfer_gap_mm"]), f"{num(der['voltages.potential_drift_end_kV'])} → {num(der['voltages.potential_el_entry_kV'])}", f"{num(inp['voltages.transfer_gap_field_V_cm'])} V/cm"),
        ("Drift region", num(inp["geometry.fiducial_length_cm"] * 10), f"{num(der['voltages.potential_tpc_face_ifg_kV'])} → {num(der['voltages.potential_drift_end_kV'])}", f"{num(inp['drift.field_V_cm'])} V/cm, {num(der['drift.drift_voltage_kV'])} kV across"),
        ("Ion focusing grid (IFG)", num(inp["geometry.ifg_thickness_mm"]), f"{num(der['voltages.potential_pad_face_ifg_kV'])} → {num(der['voltages.potential_tpc_face_ifg_kV'])}", f"{num(der['voltages.ifg_field_kV_cm'])} kV/cm"),
        ("Pad gap", num(inp["geometry.pad_gap_mm"]), f"{num(der['voltages.potential_pads_kV'])} → {num(der['voltages.potential_pad_face_ifg_kV'])}", f"{num(der['voltages.pad_gap_field_kV_cm'])} kV/cm"),
        ("Pads (z = 0)", "—", num(der["voltages.potential_pads_kV"]), "ground"),
    ]
    sections.append(("geometry", "2. Fiducial volume, masses and potentials", f"""
{table([
    I("geometry.fiducial_diameter_cm", "Fiducial diameter"), I("geometry.fiducial_length_cm", "Fiducial length (drift)"),
    D("geometry.end_area_m2", "Area of each end"), D("geometry.fiducial_volume_m3", "Fiducial volume"),
    D("geometry.gas_mass_kg", "Gas mass"), D("geometry.xenon_mass_kg", "Xenon mass"), D("geometry.xe136_mass_kg", "Xe-136 mass"),
    D("geometry.helium_mass_kg", "Helium mass"), D("geometry.reference_pure_xenon_mass_kg", "Pure enriched xenon in the same volume at 15 bar, 290 K"),
    D("geometry.track_length_relative_to_pure_xenon_15bar", "Track length relative to pure xenon at 15 bar, 290 K", lambda x: f"+{pct(x, 3)}"),
    I("geometry.z_origin", "Origin of z"), D("geometry.z_el_entry_cm", "z of the EL entry"), D("geometry.geant4_z_offset_cm", "Geant4 z offset"),
    D("geometry.pads_to_sipm_plane_mm", "Pads to SiPM plane"),
])}
<h3>Potentials along the axis, anode at top</h3>
<p>The pads and the SiPM plane are at ground. The potential rises from the TPC-facing face of the IFG through the drift region, the transfer gap and the EL region, and is graded back to ground along the light guides. Only the independent potentials are card values; the ladder is derived.</p>
{table(stack, head=("Region", "Length (mm)", "Potential (kV)", "Field"))}
{table([
    I("voltages.pads_V", "Pads"), I("voltages.pad_face_ifg_V", "Pad-facing IFG face"), I("voltages.tpc_face_ifg_V", "TPC-facing IFG face"),
    I("voltages.transfer_gap_field_V_cm", "Transfer-gap field"), I("voltages.el_voltage_kV", "Voltage across the EL"), I("voltages.sipm_plane_V", "SiPM plane"),
    D("voltages.ifg_field_kV_cm", "Field across the IFG"), D("voltages.pad_gap_field_kV_cm", "Field in the pad gap"),
    D("voltages.grading_field_kV_cm", "Average grading field along the guides"), D("voltages.grading_reduced_field_kV_cm_bar", "Reduced grading field"),
], head=("Potential or field", "Value", "Unit", "Definition"))}"""))

    sections.append(("drift", "3. Drift region", f"""
<p>Electron transport from Magboltz for Xe/He 85/15 at the default pressure, temperature and field, with the real-gas density; the errors are the statistical errors of the run. The ion velocity comes from the literature mobility.</p>
{table([
    I("drift.field_V_cm", "Drift field"), D("drift.drift_voltage_kV", "Drift voltage"),
    D("drift.reduced_field_V_cm_bar", "Reduced drift field"), D("drift.reduced_field_Td", "Reduced drift field, field per atom"),
    D("drift.reduced_field_ideal_gas_Td", "Reduced drift field, ideal-gas density"),
    I("drift.electron_drift_velocity_mm_us", "Electron drift velocity", note=f"{rel('drift.electron_drift_velocity_relative_error')}, Magboltz"),
    D("drift.electron_max_drift_time_ms", "Maximum electron drift time"), I("drift.electron_lifetime_ms", "Electron lifetime"),
    D("drift.electron_max_attachment_loss", "Maximum attachment loss", pct),
    I("drift.electron_transverse_diffusion_mm_sqrtm", "Transverse diffusion coefficient", note=f"{rel('drift.electron_transverse_diffusion_relative_error')}, Magboltz"),
    I("drift.electron_longitudinal_diffusion_mm_sqrtm", "Longitudinal diffusion coefficient", note=f"{rel('drift.electron_longitudinal_diffusion_relative_error')}, Magboltz"),
    D("drift.electron_transverse_diffusion_mm_sqrtbar_sqrtcm", "Transverse diffusion, card units"),
    D("drift.electron_longitudinal_diffusion_mm_sqrtbar_sqrtcm", "Longitudinal diffusion, card units"),
    D("drift.electron_transverse_spread_mm", "Electron transverse spread at full drift"), D("drift.electron_longitudinal_spread_mm", "Electron longitudinal spread at full drift"),
    D("drift.electron_longitudinal_spread_us", "Electron longitudinal spread in time"), D("drift.electron_bucket_length_mm", "Length of one electron bucket"),
    I("drift.ion_species", "Ion"), I("drift.ion_drift_velocity_cm_s", "Ion drift velocity", note="reduced mobility 0.65 cm²/(V·s) in xenon, raised 17 % for helium by Blanc's law"),
    D("drift.ion_drift_velocity_mm_ms", "Ion drift velocity, card units"), D("drift.ion_reduced_mobility_cm2_V_s", "Reduced ion mobility implied by the velocity"),
    D("drift.ion_max_drift_time_s", "Maximum ion drift time"), D("drift.ion_thermal_spread_mm", "Ion thermal spread at full drift"),
    D("drift.ion_bucket_length_mm", "Length of one ion bucket"),
])}"""))

    sections.append(("el", "4. EL region", f"""
{table([
    I("geometry.el_depth_mm", "EL depth"), I("voltages.el_voltage_kV", "Voltage across the EL"), D("el.field_kV_cm", "EL field"),
    D("el.reduced_field_kV_cm_bar", "Reduced EL field"), I("el.yield_gain_photons_per_electron_kV", "Yield gain", note="Fernandes et al. (2020), Table 1, 15 % He"),
    I("el.yield_threshold_kV_cm_bar", "Yield threshold", note="Fernandes et al. (2020), Table 1, 15 % He"), D("el.threshold_voltage_kV", "EL voltage at zero yield"),
    D("el.photons_per_electron", "Photons per electron"), I("el.variance_ratio_J", "EL variance ratio J"),
])}
<p>Photons per electron = gain × (V<sub>EL</sub> − threshold × P × d), with V<sub>EL</sub> in kV, P in bar and d in cm. The reduced EL field applies inside the holes and is a different quantity from the reduced drift field of Section 3.</p>"""))

    sections.append(("elm", "5. Electroluminescent module and SiPM plane", f"""
{table([
    I("el.hole_pitch_mm", "Hole pitch"), I("el.hole_diameter_mm", "Hole diameter"), D("el.hole_open_fraction", "Open fraction of the plane", pct),
    D("el.hole_count", "Holes", lambda x: num(x, 3)), I("geometry.transfer_gap_mm", "Transfer gap"), I("geometry.guide_length_mm", "Light guide length", note=f"policy {inp['geometry.guide_length_policy']}"),
    I("el.guide_refractive_index", "Guide refractive index"), I("el.gap_refractive_index", "Gap refractive index"), I("el.gas_gap_radial_mm", "Radial gas gap around the guides"),
    I("el.photon_band", "Photon band at the SiPM"), I("el.tpb_conversion_efficiency", "TPB conversion efficiency"), I("el.forward_fraction", "Forward fraction after TPB"),
    I("el.propagation_transmission", "Propagation transmission"), I("el.coupling_transmission", "Coupling transmission"), I("el.electrode_transmission", "Electrode transmission"),
])}
<h3>SiPMs, one per hole</h3>
{table([
    I("sipm.device", "Device"), I("sipm.side_mm", "Side"), D("sipm.coverage", "SiPM area over cell area", pct), D("sipm.sipm_count", "SiPMs", lambda x: num(x, 3)),
    I("sipm.microcell_count", "Microcells", lambda x: num(x, 5)), I("sipm.microcell_pitch_um", "Microcell pitch"), I("sipm.pde", "PDE"),
    I("sipm.gain", "Gain", lambda x: f"{x:.2e}"), I("sipm.dark_rate_Hz", "Dark rate", lambda x: f"{x:.0e}"), I("sipm.recovery_ns", "Recovery time"),
    I("sipm.electron_bucket_us", "Electron bucket"),
])}
<p>The response model of the SiPM and the ASIC are run cards under <code>config/specs/response/sipm/</code> and <code>config/specs/response/sipm_electronics/</code>.</p>"""))

    sections.append(("ions", "6. Ion readout", f"""
{table([
    I("ions.readout_chip", "Readout chip"), I("ions.pad_pitch_mm", "Pad pitch"), I("ions.pad_side_mm", "Collecting electrode side"),
    D("ions.pad_fill_fraction", "Electrode area over cell area", pct), D("ions.pad_count", "Pads", lambda x: num(x, 3)),
    I("ions.ion_bucket_ms", "Ion bucket"), I("geometry.ifg_thickness_mm", "IFG thickness"), I("geometry.pad_gap_mm", "Pad gap"),
    D("voltages.ifg_field_kV_cm", "Field across the IFG"), D("voltages.pad_gap_field_kV_cm", "Field in the pad gap"),
])}
<p>The charge-sensitive amplifier and the converter are run cards under <code>config/specs/response/ion_csa/</code>.</p>"""))

    sections.append(("energy", "7. Energy", f"""
{table([
    I("energy.W_eV", "Energy per electron–ion pair, W"), I("energy.fano_factor", "Fano factor"), I("energy.qbb_keV", "Q<sub>ββ</sub> of Xe-136"),
    D("energy.electrons_at_qbb", "Electrons at Q<sub>ββ</sub>", lambda x: num(x, 6)), D("energy.fano_fwhm_at_qbb_keV", "Fano-limited FWHM at Q<sub>ββ</sub>"),
    D("energy.fano_fwhm_at_qbb_relative", "Fano-limited FWHM, relative", lambda x: pct(x, 3)),
])}"""))

    sections.append(("references", "8. References", """
<ul>
<li>Gas density: NEXT DocDB 1183 (xenon compressibility), implemented in <code>simulation/nexus/source/geometries/ItacaChamber.cc</code>; enriched xenon isotopic fractions of NEXUS <code>GXeEnriched</code>.</li>
<li>Electron drift velocity and diffusion: Magboltz 11.21, run of 2026-11-15 (<code>magboltz/</code>, <code>docs/handoffs/magboltz_xe_he.md</code>).</li>
<li>Ion mobility: P. N. B. Neves et al., J. Chem. Phys. 133 (2010) 124316, as plotted in A. F. V. Cortez et al., JINST 12 (2017) P06012, Fig. 3; helium correction by Blanc's law.</li>
<li>EL yield: A. F. M. Fernandes et al., JHEP 04 (2020) 034, Table 1 (15 % He). EL variance ratio J: user reading of Oliveira et al., Fig. 6 (2026-11-04).</li>
<li>Potentials and grading field: <code>manuscripts/itaca3/itaca3.tex</code> (IFG at 8 kV, sensors at ground) and <code>docs/elm.tex</code> (regions and grading equation).</li>
<li>Ion thermal spread: ITACA revisited, arXiv:2604.02443, Eq. (2.1). Q<sub>ββ</sub>: M. Redshaw et al., PRL 98 (2007) 053003.</li>
<li>Design values without a literature source are user design inputs of 2026-11-15.</li>
</ul>"""))

    run = meta.get("run_card", "")
    if run:
        ov = "".join(f"<li><code>{html.escape(k)}</code> = {html.escape(num(v) if isinstance(v, float) else str(v))} (default {html.escape(num(defaults[k]) if isinstance(defaults.get(k), (int, float)) else str(defaults.get(k)))})</li>" for k, v in sorted(overrides.items()))
        intro = (f'<p>Configuration of the run <code>{html.escape(run)}</code>: the default detector <code>{html.escape(meta["detector_index"])}</code> with {len(overrides)} overrides, highlighted in the tables.</p><ul>{ov}</ul>')
        title = "ITACA detector · run configuration"
    else:
        intro = (f'<p>Default detector: <code>{html.escape(meta["detector_index"])}</code> ({html.escape(meta["detector_name"])}, {html.escape(meta["detector_date"])}), eight cards described in <a href="detector_configuration.html">the configuration page</a>. '
                 f'Derived values: <code>scripts/documentation/calculate_itaca_conditions.jl</code> → <code>results/itaca_detector/itaca_conditions.h5</code>; this page is written by <code>scripts/documentation/present_itaca_detector.py</code> (<code>make itaca-detector</code>). A run card overrides any value explicitly; its page shows the overrides highlighted.</p>')
        title = "ITACA detector and operating conditions"
    text = page(title, "Gas, fiducial volume, potentials, drift and EL regions, ELM, SiPMs, ion readout and energy constants, with the quantities derived from them.", intro, sections, "ITACA · Detector and operating conditions")
    target.write_text(text)
    shown = target.resolve().relative_to(ROOT) if target.resolve().is_relative_to(ROOT) else target
    print(f"Wrote {shown} ({len(text)} bytes)")
    return cards


DIAGRAM = """<svg viewBox="0 0 1040 420" width="100%" role="img" aria-label="Configuration architecture">
<defs><marker id="a" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#526776"/></marker></defs>
<rect class="box" x="20" y="20" width="300" height="150" rx="8"/><text x="36" y="44" font-weight="600">1 · Detector set  config/detector/</text>
<text x="36" y="70">itaca.toml (index)</text><text x="36" y="92">gas · geometry · drift · voltages</text><text x="36" y="114">el · sipm · ions · energy</text><text x="36" y="150" fill="#526776">design values only, nothing derived</text>
<rect class="box" x="20" y="200" width="300" height="100" rx="8"/><text x="36" y="224" font-weight="600">2 · Response models  config/specs/response/</text>
<text x="36" y="250">VH1 · CSA · SiPM · ASIC cards</text><text x="36" y="272">device models, hashes, calibrations</text>
<rect class="box" x="20" y="320" width="300" height="80" rx="8"/><text x="36" y="344" font-weight="600">3 · Run card  config/runs/</text>
<text x="36" y="368">detector · [overrides] · events, seeds,</text><text x="36" y="388">selection, stages, response cards</text>
<rect class="box" x="400" y="150" width="220" height="120" rx="8"/><text x="416" y="176" font-weight="600">resolve_run</text>
<text x="416" y="200">load the set</text><text x="416" y="220">apply + check overrides</text><text x="416" y="240">event_settings</text><text x="416" y="260">derived_conditions</text>
<rect class="store" x="700" y="20" width="320" height="90" rx="8"/><text x="716" y="44" font-weight="600">results/readout/&lt;run&gt;/configuration/</text>
<text x="716" y="68">resolved.toml per stage: detector values after</text><text x="716" y="90">overrides, run card, SHA-256 of every input</text>
<rect class="store" x="700" y="140" width="320" height="70" rx="8"/><text x="716" y="164" font-weight="600">production scripts</text><text x="716" y="188">run_framed_detector_response.jl · run_hit_production.jl</text>
<rect class="store" x="700" y="240" width="320" height="70" rx="8"/><text x="716" y="264" font-weight="600">pages</text><text x="716" y="288">itaca_detector.html (default and per run)</text>
<rect class="store" x="700" y="340" width="320" height="60" rx="8"/><text x="716" y="364" font-weight="600">Geant4 macro (planned)</text><text x="716" y="386">same resolver, same values</text>
<path class="arrow" d="M320,95 C360,95 360,190 400,190"/><path class="arrow" d="M320,250 C360,250 360,220 400,220"/><path class="arrow" d="M320,360 C360,360 360,240 400,240"/>
<path class="arrow" d="M620,190 C660,190 660,65 700,65"/><path class="arrow" d="M620,200 C660,200 660,175 700,175"/><path class="arrow" d="M620,220 C660,220 660,275 700,275"/><path class="arrow" d="M620,240 C660,240 660,370 700,370"/>
</svg>"""


def architecture_page(cards: dict, target: Path) -> None:
    order = ["index", "gas", "geometry", "drift", "voltages", "el", "sipm", "ions", "energy"]
    card_blocks = "".join(f'<h3><code>config/detector/{"itaca" if s == "index" else s}.toml</code></h3><pre><code>{html.escape(cards[s].rstrip())}</code></pre>' for s in order if s in cards)
    run_text = EXAMPLE_RUN.read_text() if EXAMPLE_RUN.exists() else ""
    sections = [
        ("layers", "1. Three layers", f"""
{DIAGRAM}
<ol>
<li><strong>Detector set.</strong> <code>config/detector/itaca.toml</code> is an index naming eight cards, one per subsystem. They hold design values only. Everything that follows from them (densities, masses, the potential ladder, reduced fields, yields, drift times, spreads, counts) is computed by the resolver and shown on the <a href="itaca_detector.html">detector page</a>, never typed in a card.</li>
<li><strong>Response models.</strong> The VH1 and CSA amplifiers, the SiPM response model and the ASIC are device models with simulation settings, calibrations and hashes. They stay as their own cards under <code>config/specs/response/</code>, named by the run card and loaded by their own functions.</li>
<li><strong>Run card.</strong> Names the detector index, lists its <code>[overrides]</code> as dotted keys into the set, and carries the run parameters: events, seeds, selection, the response cards, the stages and their settings. A run is "the default detector plus what this run changes".</li>
</ol>"""),
        ("rules", "2. Resolution rules", """
<ul>
<li>An override key must exist in the detector set, with the same type class (number, string or boolean). A misspelt key is an error, not a new key.</li>
<li>A physical key at the top level of a run card is refused; changes go through <code>[overrides]</code>.</li>
<li>Each stage writes <code>resolved.toml</code> beside its results: the detector-response stage in <code>results/readout/&lt;run_label&gt;/configuration/</code>, without the later-stage sections, so every run card that shares a detector response writes the same file; the hit-production stage in <code>&lt;run_label&gt;/&lt;output_subdir&gt;/configuration/</code> with the full run. The file holds the detector values after overrides, the overrides, the run card and the SHA-256 of every input. A different content in an existing folder is an error.</li>
<li>The production scripts read the resolved detector sections directly (<code>event_settings</code> in <code>run_framed_detector_response.jl</code>) and the response cards from <code>[response]</code>; they accept run cards only. The cards under <code>config/readout/</code> and <code>config/specs/runs/</code> predate the scheme and stay as the records of their productions, with their copies saved in the HDF5 outputs.</li>
<li>Historical cards under <code>config/readout/</code> and <code>config/specs/</code> are not edited; their productions keep their saved copies.</li>
</ul>"""),
        ("cards", "3. The detector cards", f"<p>Shown as they are in the repository, comments included.</p>{card_blocks}"),
        ("run", "4. A run card", f"""
<p><code>{html.escape(str(EXAMPLE_RUN.relative_to(ROOT)))}</code>: the default detector overridden to the conditions of the current production. Sections: run identity and inputs at the top level; <code>[overrides]</code>; <code>[selection]</code> (fast trigger and frame); <code>[response]</code> (VH1, CSA, SiPM and ASIC cards); <code>[electronics]</code> (present-chain acquisition and converter, and its own pilot stage); <code>[hit_production]</code> (events, output folder); <code>[ions]</code>, <code>[electrons]</code>, <code>[track_analysis]</code>, <code>[hit_output]</code>, <code>[controls]</code>, <code>[zs_controls]</code> for the reconstruction stages.</p>
<pre><code>{html.escape(run_text.rstrip())}</code></pre>
<h3>Making a run</h3>
<ol><li>Copy a run card under <code>config/runs/</code>; change <code>run_label</code>, the event counts and the output folder.</li>
<li>Put every changed detector value in <code>[overrides]</code>; nothing else in the card is physics.</li>
<li><code>make vh1-production RUN=config/runs/&lt;card&gt;.toml</code> runs the detector response and the VH1 hits; <code>julia --project=. scripts/documentation/calculate_itaca_conditions.jl config/runs/&lt;card&gt;.toml</code> then <code>present_itaca_detector.py &lt;h5&gt; &lt;html&gt;</code> make the run's own detector page with the overrides highlighted.</li></ol>
<p>Check of the scheme (2026-11-15): the run card above, the default detector plus nine overrides, reproduced the detector response of the current production on 145 events and the VH1 hits of its first 100 accepted events bit for bit (<code>scripts/validation/compare_detector_controls.jl</code>).</p>"""),
        ("code", "5. Code, commands and tests", """
<table><tr><th>Item</th><th>Path</th></tr>
<tr><td>Resolver and physics</td><td><code>src/itaca/ItacaConditions.jl</code>: <code>load_detector_set</code>, <code>resolve_run</code>, <code>derived_conditions</code>, <code>electronics_card</code>, <code>write_configuration</code>, <code>configuration_texts</code></td></tr>
<tr><td>Production scripts</td><td><code>scripts/production/run_framed_detector_response.jl</code> (<code>event_settings</code> builds the detector-response settings from the resolved set), <code>run_hit_production.jl</code>, <code>run_framed_electronics.jl</code>, <code>run_framed_reconstruction.jl</code>; all take a run card</td></tr>
<tr><td>Run cards</td><td><code>config/runs/</code>: the two bb0nu productions (6156 and 1000 events), the single-electron production, the pilot card</td></tr>
<tr><td>Derived values to HDF5</td><td><code>scripts/documentation/calculate_itaca_conditions.jl [run card]</code></td></tr>
<tr><td>Pages</td><td><code>scripts/documentation/present_itaca_detector.py</code>; <code>make itaca-detector</code> runs both steps for the default detector</td></tr>
<tr><td>Tests</td><td><code>tests/itaca/test_itaca_conditions.jl</code> (<code>make test-itaca</code>): card schemas and sections; no derived key in a card; override checks; the density model against the handoff references; the derived values of the default detector; the detector-response settings of the pilot card against the production card field by field at Float32; every run card resolves and names existing response cards; configuration writing per stage</td></tr>
</table>"""),
    ]
    text = page("Detector configuration", "How the detector is defined: detector cards, response-model cards and run cards with explicit overrides.",
                "", sections, "ITACA · Detector configuration")
    target.write_text(text)
    print(f"Wrote {target.relative_to(ROOT)} ({len(text)} bytes)")


def main() -> None:
    if len(sys.argv) == 3:
        detector_page(Path(sys.argv[1]), Path(sys.argv[2]))
    else:
        cards = detector_page(SOURCE, TARGET)
        architecture_page(cards, ARCHITECTURE)


if __name__ == "__main__":
    main()
