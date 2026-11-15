"""Write notes/html/itaca_detector.html from results/itaca_detector/itaca_conditions.h5.

The HDF5 file is produced by scripts/documentation/calculate_itaca_conditions.jl
from config/itaca.toml. This script only formats: no physical quantity is
computed here. Run `make itaca-detector` for both steps.
"""
from __future__ import annotations

import html
import math
from pathlib import Path

import h5py

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results/itaca_detector/itaca_conditions.h5"
TARGET = ROOT / "notes/html/itaca_detector.html"


def read(path: Path):
    inputs, derived, units, notes = {}, {}, {}, {}
    with h5py.File(path, "r") as f:
        meta = {k: (v.decode() if isinstance(v, bytes) else str(v)) for k, v in f.attrs.items()}
        for k, ds in f["inputs"].items():
            v = ds[()]
            inputs[k] = v.decode() if isinstance(v, bytes) else float(v)
        for k, ds in f["derived"].items():
            derived[k] = float(ds[()])
            units[k] = ds.attrs["unit"].decode() if isinstance(ds.attrs["unit"], bytes) else str(ds.attrs["unit"])
            notes[k] = ds.attrs["note"].decode() if isinstance(ds.attrs["note"], bytes) else str(ds.attrs["note"])
    return meta, inputs, derived, units, notes


def num(x: float, digits: int = 4) -> str:
    """Format with `digits` significant figures, trailing zeros removed, thousands separated."""
    if x == 0:
        return "0"
    if abs(x) >= 1e5 or abs(x) < 1e-3:
        return f"{x:.{digits - 1}e}".replace("e+0", "e").replace("e+", "e")
    magnitude = math.floor(math.log10(abs(x)))
    decimals = max(digits - 1 - magnitude, 0)
    s = f"{x:,.{decimals}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s


def pct(x: float, digits: int = 2) -> str:
    return f"{100 * x:.{digits}g} %"


UNIT_HTML = {"kg/m^3": "kg/m³", "m^3": "m³", "m^2": "m²", "1/m^3": "m⁻³", "cm^2/(V s)": "cm²/(V·s)",
             "us": "µs", "mm/sqrt(m)": "mm/√m", "mm sqrt(bar)/sqrt(cm)": "mm·√bar/√cm"}


def unit(u: str) -> str:
    return UNIT_HTML.get(u, u)


def table(rows, head=("Quantity", "Value", "Unit", "Definition")) -> str:
    out = ["<table><tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in head) + "</tr>"]
    for r in rows:
        out.append("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>")
    out.append("</table>")
    return "".join(out)


def main() -> None:
    meta, inp, der, units, notes = read(SOURCE)

    def I(key, label, fmt=num, unit_override=None, note=""):
        """Row for a card value."""
        v = inp[key]
        shown = fmt(v) if isinstance(v, float) else html.escape(str(v))
        u = unit_override if unit_override is not None else unit_from_key(key)
        return (label, shown, u, html.escape(note))

    def D(key, label, fmt=num, note=None):
        """Row for a derived value."""
        return (label, fmt(der[key]), unit(units[key]), html.escape(note if note is not None else notes[key]))

    def unit_from_key(key: str) -> str:
        tail = key.rsplit("_", 1)[-1]
        table_ = {"bar": "bar", "K": "K", "cm": "cm", "mm": "mm", "us": "µs", "ms": "ms", "s": "s", "V": "V", "kV": "kV",
                  "eV": "eV", "sqrtm": "mm/√m", "sqrtcm": "mm·√bar/√cm", "J": "", "error": "", "fraction": "",
                  "mixture": "", "enrichment": "", "fractions": "", "species": "", "policy": "", "a": "", "b": ""}
        if key.endswith("_V_cm"):
            return "V/cm"
        if key.endswith("_kV_cm"):
            return "kV/cm"
        if key.endswith("_kV_cm_bar"):
            return "kV/cm/bar"
        if key.endswith("_mm_us"):
            return "mm/µs"
        if key.endswith("_cm_s"):
            return "cm/s"
        if key.endswith("_per_electron_kV"):
            return "photons/(electron·kV)"
        return table_.get(tail, "")

    rel = lambda k: f"± {pct(inp[k])}"  # noqa: E731
    sections = []

    # 1. Gas
    iso = ", ".join(f"{k[2:]}: {pct(inp[f'gas.xenon_isotopic_fractions.{k}'], 4)}" for k in
                    ("Xe129", "Xe130", "Xe131", "Xe132", "Xe134", "Xe136"))
    sections.append(("gas", "1. Gas", f"""
<p>Xenon enriched in Xe-136 with 15 % helium by mole. Densities follow the model of the Geant4 geometry: xenon as a real gas with a second-virial compressibility, helium ideal, additive partial volumes.</p>
{table([
    I("gas.mixture", "Mixture"),
    I("gas.helium_mole_fraction", "Helium mole fraction", pct),
    ("Xenon isotopic fractions", iso, "", "NEXUS GXeEnriched"),
    I("gas.pressure_bar", "Pressure"),
    I("gas.temperature_K", "Temperature"),
    D("gas.xenon_molar_mass_g_mol", "Molar mass of the enriched xenon"),
    D("gas.xe136_mass_fraction", "Xe-136 mass fraction", pct),
    D("gas.xenon_compressibility", "Xenon compressibility z"),
    D("gas.density_kg_m3", "Gas density"),
    D("gas.xenon_density_kg_m3", "Xenon density"),
    D("gas.helium_density_kg_m3", "Helium density"),
    D("gas.number_density_m3", "Number density", lambda x: f"{x:.3e}"),
    D("gas.number_density_relative_to_STP", "Number density relative to STP, N/N₀"),
    D("gas.ideal_equivalent_pressure_bar", "Ideal-gas equivalent pressure"),
])}
<p>Density model: z = 1 + (a − b / T<sub>r</sub><sup>1.6</sup>) P<sub>r</sub> / T<sub>r</sub>, with a = {num(inp['gas.density_model.xenon_virial_a'])}, b = {num(inp['gas.density_model.xenon_virial_b'])}, T<sub>r</sub> = T / {num(inp['gas.density_model.xenon_critical_temperature_K'])} K, P<sub>r</sub> = P / {num(inp['gas.density_model.xenon_critical_pressure_bar'])} bar; ρ = P (x<sub>Xe</sub> M<sub>Xe</sub> + x<sub>He</sub> M<sub>He</sub>) / ((x<sub>Xe</sub> z + x<sub>He</sub>) R T). NEXT DocDB 1183, as in the Geant4 geometry <code>ItacaChamber.cc</code>.</p>
"""))

    # 2. Geometry and potentials
    stack = [
        ("SiPM plane", "—", num(der["drift.potential_sipm_plane_kV"]), "sensors at ground"),
        ("Light guides", num(inp["elm.guide_length_mm"]), f"{num(der['drift.potential_el_exit_kV'])} → 0", f"average grading field {num(der['elm.grading_field_kV_cm'])} kV/cm"),
        ("EL region", num(inp["el.depth_mm"]), f"{num(der['drift.potential_el_entry_kV'])} → {num(der['drift.potential_el_exit_kV'])}", f"{num(inp['el.voltage_kV'])} kV across, field {num(der['el.field_kV_cm'])} kV/cm"),
        ("Transfer gap", num(inp["elm.transfer_gap_mm"]), f"{num(der['drift.potential_drift_end_kV'])} → {num(der['drift.potential_el_entry_kV'])}", f"{num(inp['elm.transfer_gap_field_V_cm'])} V/cm"),
        ("Drift region", num(inp["geometry.fiducial_length_cm"] * 10), f"{num(der['drift.potential_tpc_face_ifg_kV'])} → {num(der['drift.potential_drift_end_kV'])}", f"{num(inp['drift.field_V_cm'])} V/cm, {num(der['drift.drift_voltage_kV'])} kV across"),
        ("Ion focusing grid (IFG)", num(inp["ion_readout.ifg_thickness_mm"]), f"{num(der['drift.potential_pad_face_ifg_kV'])} → {num(der['drift.potential_tpc_face_ifg_kV'])}", f"{num(der['ion_readout.ifg_field_kV_cm'])} kV/cm"),
        ("Pad gap", num(inp["ion_readout.pad_gap_mm"]), f"{num(der['drift.potential_pads_kV'])} → {num(der['drift.potential_pad_face_ifg_kV'])}", f"{num(der['ion_readout.pad_gap_field_kV_cm'])} kV/cm"),
        ("Pads", "—", num(der["drift.potential_pads_kV"]), "ground"),
    ]
    sections.append(("geometry", "2. Fiducial volume, masses and potentials", f"""
{table([
    I("geometry.fiducial_diameter_cm", "Fiducial diameter"),
    I("geometry.fiducial_length_cm", "Fiducial length (drift)"),
    D("geometry.end_area_m2", "Area of each end"),
    D("geometry.fiducial_volume_m3", "Fiducial volume"),
    D("geometry.gas_mass_kg", "Gas mass"),
    D("geometry.xenon_mass_kg", "Xenon mass"),
    D("geometry.xe136_mass_kg", "Xe-136 mass"),
    D("geometry.helium_mass_kg", "Helium mass"),
    D("geometry.reference_pure_xenon_mass_kg", "Pure enriched xenon in the same volume at 15 bar, 290 K"),
    D("geometry.track_length_relative_to_pure_xenon_15bar", "Track length relative to pure xenon at 15 bar, 290 K", lambda x: f"+{pct(x, 3)}"),
])}
<h3>Potentials along the axis, anode at top</h3>
<p>The pads and the SiPM plane are at ground. The high voltage is on the anode side: the potential rises from the TPC-facing face of the IFG through the drift region, the transfer gap and the EL region, and is graded back to ground along the light guides.</p>
{table(stack, head=("Region", "Length (mm)", "Potential (kV)", "Field"))}
"""))

    # 3. Drift region
    sections.append(("drift", "3. Drift region", f"""
<p>Electron transport from Magboltz for Xe/He 85/15 at {num(inp['gas.pressure_bar'])} bar, {num(inp['gas.temperature_K'])} K and {num(inp['drift.field_V_cm'])} V/cm, with the real-gas density; the errors are the statistical errors of the run. The ion velocity comes from the literature mobility.</p>
{table([
    I("drift.field_V_cm", "Drift field"),
    D("drift.drift_voltage_kV", "Drift voltage"),
    D("drift.reduced_field_V_cm_bar", "Reduced drift field"),
    D("drift.reduced_field_Td", "Reduced drift field, field per atom"),
    D("drift.reduced_field_ideal_gas_Td", "Reduced drift field, ideal-gas density"),
    I("drift.electron_drift_velocity_mm_us", "Electron drift velocity", note=f"{rel('drift.electron_drift_velocity_relative_error')}, Magboltz"),
    D("drift.electron_max_drift_time_ms", "Maximum electron drift time"),
    I("drift.electron_lifetime_ms", "Electron lifetime"),
    D("drift.electron_max_attachment_loss", "Maximum attachment loss", pct),
    I("drift.electron_transverse_diffusion_mm_sqrtm", "Transverse diffusion coefficient", note=f"{rel('drift.electron_transverse_diffusion_relative_error')}, Magboltz"),
    I("drift.electron_longitudinal_diffusion_mm_sqrtm", "Longitudinal diffusion coefficient", note=f"{rel('drift.electron_longitudinal_diffusion_relative_error')}, Magboltz"),
    I("drift.electron_transverse_diffusion_mm_sqrtbar_sqrtcm", "Transverse diffusion, card units", note="D[mm·√bar/√cm] = D[mm/√m] √P / 10"),
    I("drift.electron_longitudinal_diffusion_mm_sqrtbar_sqrtcm", "Longitudinal diffusion, card units", note="D[mm·√bar/√cm] = D[mm/√m] √P / 10"),
    D("drift.electron_transverse_spread_mm", "Electron transverse spread at full drift"),
    D("drift.electron_longitudinal_spread_mm", "Electron longitudinal spread at full drift"),
    D("drift.electron_longitudinal_spread_us", "Electron longitudinal spread in time"),
    I("drift.ion_species", "Ion"),
    I("drift.ion_drift_velocity_cm_s", "Ion drift velocity", note="reduced mobility 0.65 cm²/(V·s) in xenon, raised 17 % for helium by Blanc's law"),
    D("drift.ion_reduced_mobility_cm2_V_s", "Reduced ion mobility implied by the velocity"),
    D("drift.ion_max_drift_time_s", "Maximum ion drift time"),
    D("drift.ion_thermal_spread_mm", "Ion thermal spread at full drift"),
])}
"""))

    # 4. EL region
    sections.append(("el", "4. EL region", f"""
{table([
    I("el.depth_mm", "EL depth"),
    I("el.voltage_kV", "Voltage across the EL"),
    D("el.field_kV_cm", "EL field"),
    D("el.reduced_field_kV_cm_bar", "Reduced EL field"),
    I("el.yield_gain_photons_per_electron_kV", "Yield gain", note="Fernandes et al. (2020), Table 1, 15 % He"),
    I("el.yield_threshold_kV_cm_bar", "Yield threshold", note="Fernandes et al. (2020), Table 1, 15 % He"),
    D("el.threshold_voltage_kV", "EL voltage at zero yield"),
    D("el.photons_per_electron", "Photons per electron"),
    I("el.variance_ratio_J", "EL variance ratio J"),
])}
<p>Photons per electron = gain × (V<sub>EL</sub> − threshold × P × d), with V<sub>EL</sub> in kV, P in bar and d in cm. The reduced EL field applies inside the holes and is a different quantity from the reduced drift field of Section 3.</p>
"""))

    # 5. ELM and SiPM plane
    sections.append(("elm", "5. Electroluminescent module and SiPM plane", f"""
{table([
    I("elm.hole_pitch_mm", "Hole pitch"),
    I("elm.hole_diameter_mm", "Hole diameter"),
    D("elm.hole_open_fraction", "Open fraction of the plane", pct),
    D("elm.hole_count", "Holes, one SiPM each", lambda x: num(x, 3)),
    I("elm.transfer_gap_mm", "Transfer gap"),
    I("elm.transfer_gap_field_V_cm", "Transfer-gap field"),
    I("elm.guide_length_mm", "Light guide length", note=f"policy {inp['elm.guide_length_policy']}"),
    I("elm.guide_refractive_index", "Guide refractive index"),
    I("elm.gas_gap_radial_mm", "Radial gas gap around the guides"),
    D("elm.grading_field_kV_cm", "Average grading field along the guides"),
    D("elm.grading_reduced_field_kV_cm_bar", "Reduced grading field"),
    I("elm.sipm_side_mm", "SiPM side"),
    D("elm.sipm_coverage", "SiPM area over cell area", pct),
    I("elm.tpb_conversion_efficiency", "TPB conversion efficiency"),
    I("elm.forward_fraction", "Forward fraction after TPB"),
    I("el.sipm_plane_V", "SiPM plane potential"),
])}
<p>The SiPM device, its response model and the readout electronics are selected by the cards under <code>config/specs/response/sipm/</code> and <code>config/specs/response/sipm_electronics/</code>.</p>
"""))

    # 6. Ion readout
    sections.append(("ion", "6. Ion readout", f"""
{table([
    I("ion_readout.tpc_face_ifg_V", "TPC-facing IFG face"),
    I("ion_readout.pad_face_ifg_V", "Pad-facing IFG face"),
    I("ion_readout.pads_V", "Pads"),
    I("ion_readout.ifg_thickness_mm", "IFG thickness"),
    D("ion_readout.ifg_field_kV_cm", "Field across the IFG"),
    I("ion_readout.pad_gap_mm", "Pad gap"),
    D("ion_readout.pad_gap_field_kV_cm", "Field in the pad gap"),
])}
<p>The pad pitch, the IFG hole geometry, the charge-sensitive amplifier and the converter are selected by the run and response cards under <code>config/</code>.</p>
"""))

    # 7. Energy
    sections.append(("energy", "7. Energy", f"""
{table([
    I("energy.W_eV", "Energy per electron–ion pair, W"),
    I("energy.fano_factor", "Fano factor"),
    D("energy.qbb_keV", "Q<sub>ββ</sub> of Xe-136"),
    D("energy.electrons_at_qbb", "Electrons at Q<sub>ββ</sub>", lambda x: num(x, 6)),
    D("energy.fano_fwhm_at_qbb_keV", "Fano-limited FWHM at Q<sub>ββ</sub>"),
    D("energy.fano_fwhm_at_qbb_relative", "Fano-limited FWHM, relative", lambda x: pct(x, 3)),
])}
"""))

    # 8. References
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

    nav = "".join(f'<a href="#{sid}">{html.escape(title.split(". ", 1)[1])}</a>' for sid, title, _ in sections)
    body = "".join(f'<section id="{sid}"><h2>{title}</h2>{content}</section>' for sid, title, content in sections)
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>ITACA detector and operating conditions</title><style>body{{font:17px/1.65 system-ui;color:#203343;background:#f3f6f8;margin:0}}header{{background:#17394b;color:white;padding:30px max(24px,calc((100vw - 1100px)/2))}}header a{{color:#d4ecf6}}main{{max-width:1100px;margin:auto;padding:28px 24px}}section{{background:white;border:1px solid #d8e3e8;border-radius:10px;padding:25px;margin-bottom:24px;scroll-margin-top:15px}}h1,h2,h3{{line-height:1.25}}h2{{color:#185b83}}a{{color:#185b83}}code{{font-family:ui-monospace,monospace;font-size:.9em}}table{{border-collapse:collapse;width:100%;font-size:15px}}th,td{{text-align:left;padding:10px;border-bottom:1px solid #d8e3e8;vertical-align:top}}th{{background:#edf3f6}}td:nth-child(2){{white-space:nowrap}}nav{{display:flex;flex-wrap:wrap;gap:10px 20px;margin-bottom:25px}}footer{{padding:20px;color:#526776}} @media(max-width:650px){{main{{padding:15px}}section{{padding:18px}}table{{display:block;overflow:auto}}}}</style></head><body><header><a href="index.html">Documentation index</a><h1>ITACA detector and operating conditions</h1><p>Gas, fiducial volume, drift and EL regions, ELM, ion readout and energy constants, with the quantities derived from them.</p></header><main><nav>{nav}</nav>
<p>Inputs: <code>{html.escape(meta['card_path'])}</code> (dated {html.escape(meta['card_date'])}, SHA-256 <code>{html.escape(meta['card_sha256'][:16])}…</code>). Derived values: <code>scripts/documentation/calculate_itaca_conditions.jl</code> → <code>results/itaca_detector/itaca_conditions.h5</code>; this page is written by <code>scripts/documentation/present_itaca_detector.py</code> (<code>make itaca-detector</code>). Detector options under optimisation (pad pitch, SiPM device, electronics) are in the run and response cards, not here.</p>
{body}</main><footer>ITACA · Detector and operating conditions</footer></body></html>
"""
    TARGET.write_text(page)
    print(f"Wrote {TARGET.relative_to(ROOT)} ({len(page)} bytes, {len(sections)} sections)")


if __name__ == "__main__":
    main()
