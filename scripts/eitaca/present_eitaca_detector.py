#!/usr/bin/env python3
"""eITACA entry page: detector values, voltage chain, optics and the pilot.

  python scripts/eitaca/present_eitaca_detector.py

Reads results/eitaca/detector/conditions.h5 (scripts/eitaca/calculate_eitaca_conditions.jl)
and writes notes/html/eitaca/index.html. Every number comes from that file.
"""
from pathlib import Path
import html
import h5py

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'results/eitaca/detector/conditions.h5'
TARGET = ROOT / 'notes/html/eitaca/index.html'
STYLE = ("body{font:17px/1.65 system-ui;color:#203343;background:#f3f6f8;margin:0}header{background:#17394b;color:white;padding:30px max(24px,calc((100vw - 1100px)/2))}"
         "header a{color:#d4ecf6}main{max-width:1100px;margin:auto;padding:28px 24px}section{background:white;border:1px solid #d8e3e8;border-radius:10px;padding:25px;margin-bottom:24px}"
         "h1,h2,h3{line-height:1.25}h2{color:#185b83}a{color:#185b83}code{font-family:ui-monospace,monospace;font-size:.9em;overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;font-size:15px}"
         "th,td{text-align:left;padding:9px;border-bottom:1px solid #d8e3e8;vertical-align:top}th{background:#edf3f6}td:nth-child(2){white-space:nowrap}pre{background:#f3f7f9;padding:12px;overflow-x:auto}"
         "footer{padding:20px;color:#526776}li{margin:6px 0}@media(max-width:650px){main{padding:16px}section{padding:18px}table{display:block;overflow:auto}}")


def text(v):
    return v.decode() if isinstance(v, bytes) else v


def table(header, rows):
    return ('<table><tr>' + ''.join(f'<th>{h}</th>' for h in header) + '</tr>'
            + ''.join('<tr>' + ''.join(f'<td>{c}</td>' for c in r) + '</tr>' for r in rows) + '</table>')


def g(x, d=4):
    return f'{x:,.{d}g}' if abs(x) < 1e6 else f'{x:.{d}g}'


def render():
    with h5py.File(SOURCE) as f:
        I = {k: text(v[()]) for k, v in f['inputs'].items()}
        D = {k: (float(v['value'][()]), text(v['unit'][()]), text(v['note'][()]), int(v['order'][()])) for k, v in f['derived'].items()}
        variant, date = text(f['variant'][()]), text(f['date'][()])
    d = lambda k: D[k][0]
    inputs = [
        ('Gas', f"enriched xenon with {100*I['gas.helium_mole_fraction']:.1f} % helium (mole fraction), {I['gas.pressure_bar']:g} bar, {I['gas.temperature_K']:g} K", 'gas.toml'),
        ('Fiducial volume', f"radius {I['geometry.fiducial_diameter_cm']/200:g} m, drift length {I['geometry.fiducial_length_cm']/100:g} m (cathode to gate)", 'geometry.toml'),
        ('Drift field', f"{I['drift.field_V_cm']:g} V/cm ({d('drift.reduced_field_V_cm_bar'):g} V/cm/bar)", 'drift.toml'),
        ('Electron transport', f"v = {I['drift.electron_drift_velocity_mm_us']:g} ± {I['drift.electron_drift_velocity_error_mm_us']:g} mm/µs; "
         f"D*_T = {I['drift.electron_reduced_transverse_diffusion_um_sqrtbar_sqrtcm']:g} ± {I['drift.electron_reduced_transverse_diffusion_error_um_sqrtbar_sqrtcm']:g}, "
         f"D*_L = {I['drift.electron_reduced_longitudinal_diffusion_um_sqrtbar_sqrtcm']:g} ± {I['drift.electron_reduced_longitudinal_diffusion_error_um_sqrtbar_sqrtcm']:g} µm√bar/√cm",
         'drift.toml; papers/xehe_paper.pdf, Table 1, 15.8 % He, 25 V/cm/bar'),
        ('Electron lifetime', f"{I['drift.electron_lifetime_ms']:g} ms", 'drift.toml'),
        ('EL gap', f"{I['geometry.el_depth_mm']:g} mm", 'geometry.toml'),
        ('EL yield', f"{I['el.yield_gain_photons_per_electron_kV']:g} × (E/p − {I['el.yield_threshold_kV_cm_bar']:g} kV/cm/bar) × P × d photons per electron; J = {I['el.variance_ratio_J']:g}",
         'el.toml; papers/xehe_paper.pdf, Section 5.3'),
        ('EL plane', f"holes of {I['el.hole_diameter_mm']:g} mm on a {I['el.hole_pitch_mm']:g} mm pitch, one SiPM per hole, no light guides", 'el.toml'),
        ('SiPM', f"{I['sipm.device']}, {I['sipm.side_mm']:g} mm, against the EL anode, TPB on its window, PDE {I['sipm.pde']:g}; {1000*I['sipm.electron_bucket_us']:g} ns buckets", 'sipm.toml'),
        ('SiPM response and ASIC', f"<code>{I['sipm.response_card']}</code>, <code>{I['sipm.asic_card']}</code>, the cards of the hybrid production", 'sipm.toml'),
        ('Fast trigger', f"{I['trigger.frame_side_mm']:g} mm cube, {I['trigger.simulation_margin_mm']:g} mm margin, {I['trigger.gamma_energy_cut_keV']:g} keV gamma veto, "
         f"main track ≥ {I['trigger.el_exclusion_mm']:g} mm from the EL entry", 'trigger.toml'),
        ('Ionisation', f"W = {I['energy.W_eV']:g} eV, F = {I['energy.fano_factor']:g}", 'energy.toml'),
    ]
    ladder = [('Cathode', d('voltages.potential_cathode_kV'), 'gate − E_d L'), ('Gate (EL entry, end of the drift)', d('voltages.potential_gate_kV'), 'card'),
              ('EL anode', d('voltages.potential_el_anode_kV'), 'card'), ('SiPM plane', d('voltages.potential_sipm_plane_kV'), 'card, against the anode')]
    groups = {}
    for k, (v, u, note, order) in sorted(D.items(), key=lambda kv: kv[1][3]):
        groups.setdefault(k.split('.')[0], []).append((k.split('.', 1)[1].replace('_', ' '), f'{g(v)} {html.escape(u)}', html.escape(note)))
    derived = ''.join(f'<h3>{grp}</h3>' + table(['Quantity', 'Value', 'Definition'], rows) for grp, rows in groups.items())
    a, ax = d('el.blue_photons_per_el_photon'), d('el.axial_hole_acceptance')
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>eITACA · ITACA</title><style>{STYLE}</style></head><body>
<header><p><a href="../index.html">NAUSICAA documentation</a></p><h1>eITACA</h1>
<p>The electron-reading ITACA detector in xenon–helium (also ITACAe, IE), without ion readout. Default set: variant {variant}; user design inputs, {date}.</p></header><main>
<section><h2>Pilot</h2><p><a href="../../../results/eitaca/pilot/index.html">eITACA pilot: does the SiPM readout saturate?</a> bb0nu and single electrons of the Xe/He pilot, with the eITACA optics and, as reference, the optics and yield of the hybrid production.</p></section>
<section><h2>Design values</h2>{table(['Item', 'Value', 'Card (config/eitaca/)'], inputs)}
<p>¹³⁶Xe mass: {d('geometry.xe136_mass_kg'):,.0f} kg (xenon density {d('gas.xenon_density_kg_m3'):.1f} kg/m³, density model of <code>gas.toml</code>). D = D*/√P at {I['gas.pressure_bar']:g} bar:
D_T = {d('drift.electron_transverse_diffusion_mm_sqrtm'):.2f} mm/√m, D_L = {d('drift.electron_longitudinal_diffusion_mm_sqrtm'):.2f} mm/√m.</p>
<p>EL slope: the paper gives the slope of the mixtures equal to pure xenon within 10 % and the 15.8 % threshold 1.083 ± 0.007 kV/cm/bar (fit above 1.3 kV/cm/bar). The slope 134 follows from the measured yield ratio to pure xenon at 1.6 kV/cm/bar, 0.64 ± 0.01, with pure xenon at Y/p = 140 E/p − 116 (Monteiro et al., JINST 2 (2007) P05001, the normalisation of the paper): 0.64 × 108 / (1.6 − 1.083) = 134. The pure-xenon slope 140 would give 4.5 % more light.</p></section>
<section><h2>Voltage chain</h2>{table(['Electrode', 'Potential (kV)', 'From'], [(n, f'{v:g}', s) for n, v, s in ladder])}
<p>Drift: {d('drift.drift_voltage_kV'):g} kV over {I['geometry.fiducial_length_cm']:g} cm, {I['drift.field_V_cm']:g} V/cm. EL: {d('el.el_voltage_kV'):g} kV over {I['geometry.el_depth_mm']:g} mm,
{d('el.field_kV_cm'):g} kV/cm = {d('el.reduced_field_kV_cm_bar'):g} kV/cm/bar, {d('el.photons_per_electron'):.0f} photons per electron. The drift ends at the gate: no ion plane and no transfer gap.</p></section>
<section><h2>Optics without light guides</h2>
<p>Blue photons entering the SiPM per EL photon (<code>EItaca.direct_light_acceptance</code>):</p>
<pre>acceptance = axial_light_acceptance(r, d) × TPB conversion × forward fraction
           = {ax:.4f} × {I['el.tpb_conversion_efficiency']:g} × {I['el.forward_fraction']:g} = {a:.4f}</pre>
<ol><li>VUV emission isotropic and uniform along the hole axis; absorbing hole walls (<code>axial_light_acceptance</code>, src/event_analysis/grid_readout.jl; derivation in docs/elm.tex).</li>
<li>The SiPM sits against the hole exit and its {I['sipm.side_mm']:g} mm square covers the {I['el.hole_diameter_mm']:g} mm hole: every VUV photon leaving the hole reaches the TPB on the SiPM window.</li>
<li>TPB converts with probability {I['el.tpb_conversion_efficiency']:g} and re-emits isotropically; half of the blue light goes towards the SiPM.</li>
<li>The PDE ({I['sipm.pde']:g}) is applied by the SiPM response.</li></ol>
<p>The ITACA module with light guides (<code>elm_light_acceptance</code>) multiplies the same factors by the guide capture 1 − n_gap/n_guide = 1/3. Without guides the light per EL photon is three times larger (test <code>tests/eitaca/test_eitaca.jl</code>). With the eITACA yield ({d('el.photons_per_electron'):.0f} photons per electron) the SiPMs see {d('el.photoelectrons_per_electron'):.2f} photoelectrons per drift electron.</p></section>
<section><h2>Derived quantities</h2>{derived}</section>
<section><h2>Files</h2><ul>
<li>Cards: <code>config/eitaca/</code> (index <code>eitaca.toml</code>, variant <code>variants/IE_L1500_D400.toml</code>, run cards <code>runs/</code>).</li>
<li>Module: <code>src/eitaca/EItaca.jl</code> (resolver, checks, optics, derived quantities). Tests: <code>tests/eitaca/test_eitaca.jl</code>.</li>
<li>Driver: <code>scripts/eitaca/run_eitaca_batch.jl</code> (fast trigger, electron response, SiPM and ASIC, reconstruction; no ions).</li>
<li>This page: <code>scripts/eitaca/calculate_eitaca_conditions.jl</code> then <code>scripts/eitaca/present_eitaca_detector.py</code>.</li></ul></section>
</main><footer>Generated from results/eitaca/detector/conditions.h5.</footer></body></html>"""
    TARGET.parent.mkdir(parents=True, exist_ok=True); TARGET.write_text(page); print(TARGET)


if __name__ == '__main__':
    render()
