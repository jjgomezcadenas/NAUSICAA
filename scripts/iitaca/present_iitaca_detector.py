#!/usr/bin/env python3
"""Page of the iItaca detector and its variants.

  python scripts/iitaca/present_iitaca_detector.py

Reads results/iitaca/detector/conditions.h5 (scripts/iitaca/calculate_iitaca_conditions.jl)
and writes notes/html/iitaca_detector.html. Every number of the page comes from
that file: the card values after the overrides of each variant and the derived
quantities. Nothing is recalculated here.
"""
from pathlib import Path
import html
import h5py

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'results/iitaca/detector/conditions.h5'
PAGE = ROOT / 'notes/html/iitaca_detector.html'
STYLE = ("body{font:17px/1.65 system-ui;color:#203343;background:#f3f6f8;margin:0}header{background:#17394b;color:white;padding:30px max(24px,calc((100vw - 1100px)/2))}"
         "header a{color:#d4ecf6}main{max-width:1100px;margin:auto;padding:28px 24px}section{background:white;border:1px solid #d8e3e8;border-radius:10px;padding:25px;margin-bottom:24px;scroll-margin-top:15px}"
         "h1,h2,h3{line-height:1.25}h2{color:#185b83}a{color:#185b83}code{font-family:ui-monospace,monospace;font-size:.9em}table{border-collapse:collapse;width:100%;font-size:15px}"
         "th,td{text-align:left;padding:10px;border-bottom:1px solid #d8e3e8;vertical-align:top}th{background:#edf3f6}td:nth-child(2){white-space:nowrap}table.variants td:nth-child(3){white-space:nowrap}"
         "nav{display:flex;flex-wrap:wrap;gap:10px 20px;margin-bottom:25px}footer{padding:20px;color:#526776}li{margin:6px 0} @media(max-width:650px){main{padding:15px}section{padding:18px}table{display:block;overflow:auto}}")


def text(v):
    return v.decode() if isinstance(v, bytes) else v


def load():
    """Variants in card order: name, card path, inputs, overrides, derived."""
    with h5py.File(SOURCE) as f:
        out = []
        for name, g in sorted(f['variants'].items(), key=lambda kv: kv[1].attrs['order']):
            out.append(dict(name=name, card=text(g.attrs['card']),
                            inputs={k: text(v[()]) for k, v in g['inputs'].items()},
                            overrides={k: text(v[()]) for k, v in g['overrides'].items()},
                            derived={k: float(v[()]) for k, v in g['derived'].items()}))
        return out, text(f.attrs['detector_index']), text(f.attrs['detector_date'])


def num(x, digits=4):
    """Number with thousands separators; integers without decimals."""
    if isinstance(x, str):
        return html.escape(x)
    if abs(x) >= 1e6:
        mantissa, exponent = f'{x:.{digits - 1}e}'.split('e')
        return f'{mantissa} × 10<sup>{int(exponent)}</sup>'
    s = f'{x:,.{digits}g}' if abs(x) < 1000 else f'{round(x):,}'
    return s


def table(header, rows, variants=False):
    cls = ' class="variants"' if variants else ''
    head = ''.join(f'<th>{h}</th>' for h in header)
    body = ''.join('<tr>' + ''.join(f'<td>{c}</td>' for c in r) + '</tr>' for r in rows)
    return f'<table{cls}><tr>{head}</tr>{body}</table>'


def render():
    V, index, date = load()
    names = [v['name'] for v in V]
    I0, D0 = V[0]['inputs'], V[0]['derived']
    differs = lambda key, store: len({str(v[store][key]) for v in V}) > 1

    def one(key, digits=4, store='inputs'):
        """Value common to all variants; 'see Section 3' when the variants differ."""
        return 'by variant, Section 3' if differs(key, store) else num(V[0][store][key], digits)

    def each(key, unit='', digits=4, store='inputs', scale=1.0):
        return [f"{num(v[store][key] * scale if not isinstance(v[store][key], str) else v[store][key], digits)}{' ' + unit if unit else ''}" for v in V]

    d = lambda key, digits=4: one(key, digits, 'derived')
    de = lambda key, unit='', digits=4, scale=1.0: each(key, unit, digits, 'derived', scale)
    quantity = ['Quantity', 'Value', 'Unit', 'Definition']
    by_variant = ['Quantity'] + names + ['Definition']
    listed = ', '.join(names[:-1]) + ' and ' + names[-1] if len(names) > 1 else names[0]
    spread_mm = [f"{v['derived']['drift.ion_thermal_spread_mm']:.2f} mm" for v in V]
    spreads = [x[:-3] for x in spread_mm]
    roi = f"{D0['energy.roi_lower_keV']:.2f} to {D0['energy.roi_upper_keV']:.2f}"
    fwhm_keV = I0['energy.resolution_fwhm_at_qbb_keV']
    resolution = f"{D0['energy.resolution_fwhm_at_qbb_percent']:.2f} %"
    pitch = I0['ions.pad_pitch_mm']

    def ladder(v):
        p = lambda k: v['derived'][f'voltages.potential_{k}_kV']
        k = lambda x: f'{x:.4g}'
        return {'sipm': f"{k(p('sipm_plane'))} kV", 'guides': f"{k(p('el_exit'))} → {k(p('sipm_plane'))} kV", 'exit': f"{k(p('el_exit'))} kV",
                'el': f"{k(p('el_entry'))} → {k(p('el_exit'))} kV", 'entry': f"{k(p('el_entry'))} kV", 'gap': f"{k(p('drift_end'))} → {k(p('el_entry'))} kV",
                'end': f"{k(p('drift_end'))} kV", 'drift': f"{k(p('tpc_face_ifg'))} → {k(p('drift_end'))} kV", 'ifg_tpc': f"{k(p('tpc_face_ifg'))} kV",
                'ifg': f"{k(p('pad_face_ifg'))} → {k(p('tpc_face_ifg'))} kV", 'ifg_pad': f"{k(p('pad_face_ifg'))} kV", 'padgap': f"{k(p('pads'))} → {k(p('pad_face_ifg'))} kV",
                'pads': f"{k(p('pads'))} kV"}
    L = [ladder(v) for v in V]
    row = lambda label, key, note='': [label] + [l[key] for l in L] + [note]
    drift_note = '; '.join(f"{num(v['inputs']['geometry.fiducial_length_cm'] * 10)} mm at {num(v['inputs']['drift.field_V_cm'])} V/cm" for v in V)

    overrides = ''.join(
        f"<h3>{v['name']}</h3><p>Card <code>{v['card']}</code>. " +
        ('No override: it is the default set.' if not v['overrides'] else 'Overrides of the default set:') + '</p>' +
        (table(['Key', 'Value'], [[f'<code>{k}</code>', num(x)] for k, x in sorted(v['overrides'].items())]) if v['overrides'] else '')
        for v in V)

    page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>iItaca detector and operating conditions</title><style>{STYLE}</style></head><body><header><a href="index.html">Documentation index</a> · <a href="iitaca/index.html">iItaca</a><h1>iItaca detector and operating conditions</h1><p>The ion-reading version of ITACA in pure xenon (iItaca, also written ITACAi or II): common design values, the variants {listed}, and the quantities derived from them.</p></header><main><nav><a href="#principle">Principle</a><a href="#common">Common design values</a><a href="#variants">Variants</a><a href="#derived">Derived quantities</a><a href="#simulation">Simulation</a><a href="#cards">Cards</a><a href="#open">Open items</a><a href="#sources">Sources</a></nav><p>Design values: user design inputs, detector set of {date}. The page is generated from the cards of <code>config/iitaca/</code> (<code>make iitaca-detector</code>): the values after the overrides of each variant and the quantities derived by <code>src/iitaca/IItaca.jl</code>.</p>

<section id="principle"><h2>1. Purpose and principle</h2>
<ol>
<li>The ions give the track. The Xe<sub>2</sub><sup>+</sup> molecular ions (written Xe2+, charge +e) drift to the TopMetal pads with a thermal spread of {' to '.join(sorted(spreads))} mm at full drift and are read on a cubic {num(I0['ions.z_cell_mm'])} mm lattice.</li>
<li>The electrons give the start time, the fast trigger and the event energy. They drift to the EL plane, read by SiPMs at {num(I0['el.documentation.hole_pitch_mm'])} mm pitch.</li>
<li>The gas is pure xenon enriched in Xe-136, without helium. The electron spread at full drift is about 10 mm, so the electron plane gives a coarse image and the track topology comes from the ions.</li>
<li>The simulation covers the ion side: Geant4, ion transport, VH1 electronics, reconstruction, hits and tracks. The electron side is emulated: the fast trigger from the Geant4 deposits and the energy from their sum smeared with the energy resolution.</li>
</ol></section>

<section id="common"><h2>2. Common design values</h2>
<p>Values of the default set. A quantity that a variant changes is marked and listed in Section 3. The simulation reads the gas, the fiducial radius and length, the drift field, the ion velocity, the ion plane (pitch, frame, VH1 card), the energy and the trigger. The other values, kept under <code>[documentation]</code> tables in the cards, give the potential ladder and the light budget of this page: thicknesses and potentials of the IFG and the gaps, the whole EL module and the SiPMs, the reference ion mobility and the electron diffusion.</p>
<h3>Gas</h3>
{table(quantity, [
    ['Gas', 'pure xenon' if I0['gas.helium_mole_fraction'] == 0 else f"xenon with {num(100 * I0['gas.helium_mole_fraction'])} % helium", '', 'enriched in Xe-136'],
    ['Xenon isotopic fractions', f"Xe136 = {dict(l.split(' = ') for l in I0['gas.xenon_isotopic_fractions'].strip().splitlines())['Xe136']}", '', 'NEXUS GXeEnriched'],
    ['Pressure', one('gas.pressure_bar'), 'bar', ''],
    ['Temperature', one('gas.temperature_K'), 'K', '']])}
<h3>Geometry</h3>
{table(quantity, [
    ['Fiducial radius', d('geometry.fiducial_radius_cm'), 'cm', ''],
    ['Fiducial length (drift)', one('geometry.fiducial_length_cm'), 'cm', 'TPC-facing IFG face to the end of the drift'],
    ['Origin of z', 'TopMetal pad plane', '', 'z along the drift; an event at z has ion drift length z']])}
<h3>Ion plane</h3>
{table(quantity, [
    ['Readout chip', one('ions.readout_chip'), '', ''],
    ['Pad pitch', one('ions.pad_pitch_mm'), 'mm', ''],
    ['Collecting electrode side', one('ions.pad_side_mm'), 'mm', ''],
    ['Ion focusing grid (IFG) thickness', one('geometry.documentation.ifg_thickness_mm'), 'mm', ''],
    ['Pad gap', one('geometry.documentation.pad_gap_mm'), 'mm', 'pad-facing IFG face to the pads'],
    ['Pads', one('voltages.documentation.pads_V'), 'V', 'ground'],
    ['Pad-facing IFG face', one('voltages.documentation.pad_face_ifg_V'), 'V', f"{d('voltages.pad_gap_field_kV_cm')} kV/cm in the pad gap"],
    ['TPC-facing IFG face', one('voltages.documentation.tpc_face_ifg_V'), 'V', f"{d('voltages.ifg_field_kV_cm')} kV/cm across the IFG"],
    ['Ion', one('drift.ion_species'), '', 'Xe<sub>2</sub><sup>+</sup>'],
    ['Length drifted in one frame', one('ions.z_cell_mm'), 'mm', 'ion velocity × frame; checked when the cards are resolved'],
    ['Ion electronics', 'VH1, 20 fF capacitors', '', '<a href="event_analysis/sim_vh1_readout.html">VH1 readout</a>; the card of each variant has its frame']])}
<h3>EL plane and SiPMs (present in the detector, not simulated)</h3>
{table(quantity, [
    ['Role', 'start time, fast trigger, energy', '', ''],
    ['EL hole and SiPM pitch', one('el.documentation.hole_pitch_mm'), 'mm', 'one SiPM per hole'],
    ['Transfer gap', one('geometry.documentation.transfer_gap_mm'), 'mm', 'end of the drift to the EL entry electrode'],
    ['Transfer-gap field', one('voltages.documentation.transfer_gap_field_V_cm'), 'V/cm', ''],
    ['EL gap', one('geometry.documentation.el_depth_mm'), 'mm', ''],
    ['Voltage across the EL', one('voltages.documentation.el_voltage_kV'), 'kV', f"{d('el.field_kV_cm')} kV/cm, {d('el.reduced_field_kV_cm_bar')} kV/cm/bar"],
    ['SiPM plane', one('voltages.documentation.sipm_plane_V'), 'V', 'ground, behind the light guides'],
    ['Hole and light-guide diameter', one('el.documentation.hole_diameter_mm'), 'mm', f"open fraction {100 * D0['el.hole_open_fraction']:.0f} % of the cell"],
    ['Xenon gap around each guide', one('el.documentation.gas_gap_radial_mm'), 'mm', ''],
    ['Light-guide length', one('geometry.documentation.guide_length_mm'), 'mm', 'EL exit electrode to the SiPM plane'],
    ['Guide and gap refractive indices', f"{one('el.documentation.guide_refractive_index')}, {one('el.documentation.gap_refractive_index')}", '', ''],
    ['TPB conversion efficiency, forward fraction', f"{one('el.documentation.tpb_conversion_efficiency')}, {one('el.documentation.forward_fraction')}", '', 'TPB on the guide entrance'],
    ['SiPM', one('el.documentation.sipm_device'), '', f"{one('el.documentation.sipm_side_mm')} mm side, {num(I0['el.documentation.sipm_microcell_count'])} microcells; {100 * D0['el.sipm_coverage']:.0f} % of the cell"],
    ['SiPM PDE', one('el.documentation.sipm_pde'), '', 'peak PDE at 460 nm, used for the TPB spectrum'],
    ['EL variance ratio J', one('el.documentation.variance_ratio_J'), '', 'as in the ITACA module']])}
<h3>Energy</h3>
{table(quantity, [
    ['Energy per electron–ion pair, W', one('energy.W_eV'), 'eV', ''],
    ['Fano factor', one('energy.fano_factor'), '', ''],
    ['Q<sub>ββ</sub> of Xe-136', f"{I0['energy.qbb_keV']:.2f}", 'keV', ''],
    ['Energy resolution, FWHM at Q<sub>ββ</sub>', f'{fwhm_keV:.2f}', 'keV', f"{resolution}; Gaussian σ = {D0['energy.resolution_sigma_at_qbb_keV']:.2f} keV; scales as √E. Width obtained in the full simulation of the Xe/He electron readout, kept unchanged"],
    ['ROI', roi, 'keV', f"{num(I0['energy.roi_width_fwhm'])} FWHM around Q<sub>ββ</sub>"]])}
<h3>Event frame and fast trigger</h3>
{table(quantity, [
    ['Frame', one('trigger.frame_side_mm'), 'mm', 'side of the fixed cube around the main-track box; larger tracks are rejected'],
    ['Frame placement', 'inside the detector', '', 'cube centred on the track, moved inward at the detector boundary'],
    ['Margin of the ion response', one('trigger.simulation_margin_mm'), 'mm', 'around the track and the contained deposits, inside the cube'],
    ['Gamma veto', one('trigger.gamma_energy_cut_keV'), 'keV', 'a gamma deposit above this outside the frame rejects the event'],
    ['Distance to the EL entry', one('trigger.el_exclusion_mm'), 'mm', 'the main track stays this far from the EL entry']])}
<h3>Hits and track selection</h3>
{table(quantity, [
    ['Saved hits', one('tracks.saved_hit_cut_keV'), 'keV', 'hits above this are written to the hit file'],
    ['Hit cut', one('tracks.hit_cut_keV'), 'keV', 'hits above this enter the selection and the track'],
    ['Radius of a component', one('tracks.segment_connection_mm'), 'mm', 'touching cells of the lattice'],
    ['Track segment', one('tracks.segment_min_hits'), 'hits', 'a component with at least this many hits'],
    ['Detached deposit', one('tracks.deposit_min_keV'), 'keV', 'a smaller component with at least this energy rejects the event'],
    ['Segments', one('tracks.max_segments'), '', f"at most this many, joined at {one('tracks.connection_mm')} mm"],
    ['Cells of the track graph', one('tracks.voxel_width_mm'), 'mm', ''],
    ['Endpoint sphere', one('tracks.blob_radius_mm'), 'mm', 'radius for the endpoint energies E₁ and E₂']])}
<p>Method: <a href="event_analysis/track_selection.html">track selection on the readout lattice</a>; code in <code>src/track_selection/TrackSelection.jl</code>.</p></section>

<section id="variants"><h2>3. The variants</h2>
{table(by_variant, [
    ['Drift length'] + each('geometry.fiducial_length_cm', 'cm') + [''],
    ['Drift field'] + each('drift.field_V_cm', 'V/cm') + [''],
    ['Reduced drift field'] + de('drift.reduced_field_V_cm_bar', 'V/cm/bar') + ['plain pressure'],
    ['Fiducial volume'] + de('geometry.fiducial_volume_m3', 'm³') + ['π R² L'],
    ['Xenon mass'] + de('geometry.xenon_mass_kg', 'kg') + ['density × volume'],
    ['Xe-136 mass'] + de('geometry.xe136_mass_kg', 'kg') + [f"Xe-136 mass fraction {100 * D0['gas.xe136_mass_fraction']:.1f} %"],
    ['Drift voltage'] + de('drift.drift_voltage_kV', 'kV') + ['E<sub>d</sub> L'],
    ['Highest potential (EL exit)'] + de('voltages.potential_el_exit_kV', 'kV') + [''],
    ['Ion drift velocity'] + each('drift.ion_drift_velocity_cm_s', 'cm/s') + ['user design input'],
    ['Frame'] + each('ions.frame_ms', 'ms') + [f"one frame spans {num(I0['ions.z_cell_mm'])} mm"],
    ['Longest ion drift time'] + de('drift.ion_max_drift_time_s', 's') + ['L / v<sub>ion</sub>'],
    ['Frames in the longest drift'] + de('drift.ion_frames_in_longest_drift', '', 3) + [''],
    ['Ion thermal spread at full drift'] + spread_mm + ['√(2 k T L / (e E))']], True)}
<h3>Potentials along the axis, anode at top</h3>
<p>The pads and the SiPM plane are at ground. The potential rises from the TPC-facing face of the IFG through the drift region, the transfer gap and the EL region, and is graded back to ground along the light guides.</p>
{table(['Electrode or region'] + names + ['Length and field'], [
    row('SiPM plane', 'sipm', 'sensors at ground'),
    row('Light guides', 'guides', f"{one('geometry.documentation.guide_length_mm')} mm; grading field {' and '.join(de('voltages.grading_field_kV_cm', 'kV/cm', 3))}"),
    row('EL exit (highest potential)', 'exit'),
    row('EL region', 'el', f"{one('geometry.documentation.el_depth_mm')} mm, {one('voltages.documentation.el_voltage_kV')} kV across"),
    row('EL entry', 'entry'),
    row('Transfer gap', 'gap', f"{one('geometry.documentation.transfer_gap_mm')} mm, {one('voltages.documentation.transfer_gap_field_V_cm')} V/cm"),
    row('End of drift', 'end'),
    row('Drift region', 'drift', drift_note),
    row('TPC-facing IFG face', 'ifg_tpc'),
    row('Ion focusing grid (IFG)', 'ifg', f"{one('geometry.documentation.ifg_thickness_mm')} mm, {d('voltages.ifg_field_kV_cm')} kV/cm"),
    row('Pad-facing IFG face', 'ifg_pad'),
    row('Pad gap', 'padgap', f"{one('geometry.documentation.pad_gap_mm')} mm, {d('voltages.pad_gap_field_kV_cm')} kV/cm"),
    row('Pads (z = 0)', 'pads', 'ground')], True)}
<h3>Ion drift velocity</h3>
<p>The reduced mobility K<sub>0</sub> = {num(I0['drift.documentation.ion_reference_reduced_mobility_cm2_V_s'])} cm²/(V·s) of Xe<sub>2</sub><sup>+</sup> in xenon (Neves et al., via Cortez et al.) and N/N<sub>0</sub> = {D0['gas.number_density_relative_to_STP']:.2f} give v = K<sub>0</sub> E N<sub>0</sub>/N. The adopted velocities make one frame span exactly {num(I0['ions.z_cell_mm'])} mm.</p>
{table(by_variant, [
    ['Ion velocity from the reference mobility'] + de('drift.ion_reference_velocity_cm_s', 'cm/s', 3) + ['K<sub>0</sub> E N<sub>0</sub>/N'],
    ['Adopted ion velocity'] + each('drift.ion_drift_velocity_cm_s', 'cm/s') + ['user design input'],
    ['Reduced mobility implied by the adopted velocity'] + de('drift.ion_reduced_mobility_cm2_V_s', 'cm²/(V·s)', 3) + ['v N / (N<sub>0</sub> E)']], True)}</section>

<section id="derived"><h2>4. Derived quantities</h2>
<h3>Gas</h3>
{table(quantity, [
    ['Molar mass of the enriched xenon', d('gas.xenon_molar_mass_g_mol'), 'g/mol', 'isotopic fractions of the card'],
    ['Xe-136 mass fraction', f"{100 * D0['gas.xe136_mass_fraction']:.2f} %", '', ''],
    ['Xenon compressibility z', d('gas.xenon_compressibility'), '', 'second virial, at the card pressure and temperature'],
    ['Gas density', d('gas.density_kg_m3'), 'kg/m³', 'P M / (z R T)'],
    ['Number density', d('gas.number_density_m3'), 'm⁻³', ''],
    ['Number density relative to STP, N/N₀', d('gas.number_density_relative_to_STP'), '', 'N₀ at 273.15 K and 101325 Pa']])}
<p>Density model: z = 1 + (a − b / T<sub>r</sub><sup>1.6</sup>) P<sub>r</sub> / T<sub>r</sub>, with a = {num(I0['gas.density_model.xenon_virial_a'])}, b = {num(I0['gas.density_model.xenon_virial_b'])}, T<sub>r</sub> = T / {num(I0['gas.density_model.xenon_critical_temperature_K'])} K, P<sub>r</sub> = P / {num(I0['gas.density_model.xenon_critical_pressure_bar'])} bar; ρ = P M<sub>Xe</sub> / (z R T). NEXT DocDB 1183.</p>
<h3>Readout planes</h3>
{table(quantity, [
    ['Area of each end', d('geometry.end_area_m2'), 'm²', 'π R²'],
    ['Pads', d('ions.pad_count', 3), '', 'end area over the pad pitch squared'],
    ['Electrode area over cell area', f"{100 * D0['ions.pad_fill_fraction']:.1f} %", '', 'collecting electrode over cell area'],
    ['EL holes and SiPMs', f"{round(D0['el.channel_count']):,}" if not differs('el.channel_count', 'derived') else 'by variant', '', 'end area over the EL pitch squared']])}
<h3>Ions</h3>
<p>The ion spread is thermal and does not depend on the ion velocity: σ = √(2 k T L / (e E)), the same in x, y and z (ITACA revisited, Eq. (2.1); implementation in <code>src/event_analysis/production.jl</code>).</p>
{table(by_variant, [
    ['Ion thermal spread at full drift'] + spread_mm + ['uniform field, singly charged ions'],
    ['Spread over pad pitch'] + de('drift.ion_spread_over_pad_pitch', '', 2) + [f'σ / {num(pitch)} mm'],
    ['Longest ion drift time'] + de('drift.ion_max_drift_time_s', 's') + ['L / v<sub>ion</sub>'],
    ['Pads to SiPM plane'] + de('geometry.pads_to_sipm_plane_mm', 'mm') + ['pad gap + IFG + drift + transfer gap + EL + guides']], True)}
<h3>Electrons (design values, not simulated)</h3>
<p>Diffusion in pure xenon from Table 1 of the NEXT-DEMO++ paper (<code>papers/xehe_paper.pdf</code>), measured at about 8.7 bar. Conversion: D = D* / √P, assuming that the reduced coefficients depend only on E/p; 1,000 µm·√bar/√cm is 2.5 mm/√m at 16 bar.</p>
{table(['E/p (V/cm/bar)', 'D*<sub>L</sub> (µm·√bar/√cm)', 'D*<sub>T</sub> (µm·√bar/√cm)', 'D<sub>L</sub> at 16 bar (mm/√m)', 'D<sub>T</sub> at 16 bar (mm/√m)'], [
    ['25', '1752 ± 122', '3834 ± 119', '4.38', '9.59'], ['25', '1611 ± 82', '4124 ± 442', '4.03', '10.31'],
    ['35', '1191 ± 58', '3677 ± 249', '2.98', '9.19'], ['40', '1082 ± 21', '4121 ± 200', '2.71', '10.30']])}
<p>The cards hold one value of each coefficient per variant: at 31.25 V/cm/bar the interpolation between the measured points (the two data sets give 9.3 and 10.3 mm/√m transverse, 3.4 to 3.7 longitudinal); at 21.9 V/cm/bar, below the lowest measured point, the mean of the two measurements at 25 V/cm/bar. The longitudinal coefficient rises towards low field.</p>
{table(by_variant, [
    ['Reduced drift field'] + de('drift.reduced_field_V_cm_bar', 'V/cm/bar') + [''],
    ['Transverse diffusion coefficient'] + each('drift.documentation.electron_transverse_diffusion_mm_sqrtm', 'mm/√m') + ['card value'],
    ['Longitudinal diffusion coefficient'] + each('drift.documentation.electron_longitudinal_diffusion_mm_sqrtm', 'mm/√m') + ['card value'],
    ['Transverse spread at full drift'] + de('drift.electron_transverse_spread_mm', 'mm', 3) + ['D<sub>T</sub> √L'],
    ['Longitudinal spread at full drift'] + de('drift.electron_longitudinal_spread_mm', 'mm', 3) + ['D<sub>L</sub> √L'],
    ['Electron drift velocity'] + each('drift.documentation.electron_drift_velocity_mm_us', 'mm/µs') + ['about 0.9 mm/µs in pure xenon between 20 and 40 V/cm/bar (NEXT detectors)'],
    ['Electron lifetime'] + each('drift.documentation.electron_lifetime_ms', 'ms') + ['nominal'],
    ['Longest electron drift time'] + de('drift.electron_max_drift_time_ms', 'ms', 3) + ['L / v<sub>e</sub>'],
    ['Largest attachment loss'] + [f"{100 * v['derived']['drift.electron_max_attachment_loss']:.1f} %" for v in V] + ['1 − exp(−t<sub>max</sub> / τ)']], True)}
<p>With an electron spread of about 10 mm and {num(I0['el.documentation.hole_pitch_mm'])} mm pitch, the electron plane gives the energy, the start time and a coarse image.</p>
<h3>EL yield and energy</h3>
{table(quantity, [
    ['Reduced EL field', d('el.reduced_field_kV_cm_bar'), 'kV/cm/bar', 'V<sub>EL</sub> / (d P)'],
    ['Yield law', f"{num(I0['el.documentation.yield_gain_photons_per_electron_kV'])} × (E/p − {num(I0['el.documentation.yield_threshold_kV_cm_bar'])})", 'photons/(electron·cm·bar)', 'Y/p = 140 E/p − 116 (Monteiro et al.), the normalisation of the NEXT-DEMO++ paper, whose own fit gives a threshold of 0.869 kV/cm/bar'],
    ['Photons per electron', d('el.photons_per_electron'), '', 'gain × (V<sub>EL</sub> − threshold × P × d)'],
    ['Electrons at Q<sub>ββ</sub>', d('energy.electrons_at_qbb'), '', 'Q / W'],
    ['EL photons at Q<sub>ββ</sub>', d('el.photons_at_qbb', 2), '', 'before attachment and optical collection'],
    ['EL voltage at zero yield', d('el.threshold_voltage_kV'), 'kV', 'threshold × P × d']])}
<h3>Light budget and expected energy resolution</h3>
<p>Optical model of the ITACA module (<code>docs/elm.tex</code>): isotropic emission along the hole axis and absorbing walls give the fraction of the EL photons reaching the TPB at the guide entrance; the TPB converts them and sends half forward; the guide keeps the fraction 1 − n<sub>gap</sub>/n<sub>guide</sub>. The {num(I0['el.documentation.sipm_side_mm'])} mm SiPM covers part of the exit of the {num(I0['el.documentation.hole_diameter_mm'])} mm guide; for uniform illumination of the exit it sees the ratio of the areas.</p>
{table(quantity, [
    ['Blue photons at the guide exit per EL photon', f"{D0['el.guide_exit_photons_per_el_photon']:.4f}", '', 'hole acceptance × TPB conversion × forward fraction × guide capture'],
    ['Fraction of the guide exit on the SiPM', f"{D0['el.sensor_fraction_of_guide_exit']:.3f}", '', 'side² / (π r²)'],
    ['Photoelectrons per EL photon', f"{D0['el.photoelectrons_per_el_photon']:.5f}", '', 'guide exit × sensor fraction × PDE'],
    ['Photoelectrons per electron', f"{D0['el.photoelectrons_per_electron']:.2f}", '', ''],
    ['Photoelectrons at Q<sub>ββ</sub>', d('el.photoelectrons_at_qbb', 3), '', 'before attachment'],
    ['Fano term, FWHM at Q<sub>ββ</sub>', f"{D0['energy.fano_fwhm_at_qbb_percent']:.2f} %", '', f"2.355 √(F / N); {D0['energy.fano_fwhm_at_qbb_keV']:.2f} keV"],
    ['EL fluctuation term', f"{D0['energy.el_fluctuation_fwhm_at_qbb_percent']:.2f} %", '', '2.355 √(J / N)'],
    ['Photoelectron statistics term', f"{D0['energy.photoelectron_fwhm_at_qbb_percent']:.2f} %", '', '2.355 / √N<sub>pe</sub>'],
    ['Expected FWHM at Q<sub>ββ</sub>', f"{D0['energy.expected_fwhm_at_qbb_percent']:.2f} %", '', 'the three terms in quadrature; no attachment, SiPM excess noise or calibration term'],
    ['FWHM used by the simulation', f"{D0['energy.resolution_fwhm_at_qbb_percent']:.2f} %", '', 'energy card']])}</section>

<section id="simulation"><h2>5. What the simulation does</h2>
<h3>Simulated</h3>
<ol>
<li>Geant4 in the gas of the cards, with the radius and length of each variant.</li>
<li>Ion transport to the pad plane with the thermal spread of the variant and its ion velocity, on the {num(I0['ions.z_cell_mm'])} mm lattice with the frame of the variant.</li>
<li>VH1 readout with 20 fF capacitors and linear reconstruction of the ion signal (<a href="event_analysis/sim_vh1_readout.html">VH1 readout</a>).</li>
<li>Hits and tracks from the reconstructed ion signal.</li>
</ol>
<h3>Emulated</h3>
<ol>
<li>Fast trigger and event frame from the Geant4 deposits (<a href="event_analysis/mc_fast_trigger.html">MC fast trigger</a>).</li>
<li>Event energy: sum of the Geant4 deposits inside the box given to the ion response (<code>IItaca.frame_energy_keV</code>), smeared with a Gaussian of {resolution} FWHM at Q<sub>ββ</sub>, scaling as √E (<code>IItaca.emulated_energy_keV</code>).</li>
</ol>
<h3>Order of the cuts</h3>
<ol>
<li>Fast trigger: energy window, fiducial selection, gamma veto on deposits outside the frame.</li>
<li>ROI cut on the emulated energy: {roi} keV.</li>
<li>Track selection on the ion lattice: components at 5.5 mm, detached deposits, segments rejoined at 8 mm (<a href="event_analysis/track_selection.html">track selection</a>).</li>
<li>Endpoint-energy cut or classification network (<a href="network/classification_network.html">classification network</a>).</li>
</ol></section>

<section id="cards"><h2>6. Cards</h2>
<p>The detector set is <code>{index}</code>: nine cards (gas, geometry, drift, voltages, ions, energy, el, trigger, tracks) with design values only. Values under a <code>[documentation]</code> table are not passed to the simulation (<code>IItaca.simulation_inputs</code>); changing one changes this page and no result. A variant is a card with overrides, dotted keys into the set. A run card names the set, optionally a variant, and may add its own overrides, applied after those of the variant: radius, length, pressure, fields, potentials, pitch, frame and energy resolution are all changed that way.</p>
<p>When the cards are resolved (<code>IItaca.resolve_detector</code>) the result is checked: lengths, fields, velocity and frame positive; ion velocity × frame equal to the length drifted in one frame of the ions card; the frame of the VH1 card equal to the frame of the detector; potentials rising from the pads to the IFG.</p>
{overrides}</section>

<section id="open"><h2>7. Open items</h2>
<ul>
<li>None.</li>
</ul>
<p>Closed: VH1 with 15 ms frames gives the same noise, saturation and hits as with 20 ms (<a href="../../results/readout/bb0nu_itaca_pilot_10x100_ion20_15ms/vh1_frame_comparison/index.html">comparison on the pilot events</a>).</p></section>

<section id="sources"><h2>8. Sources</h2>
<ul>
<li>Design values: user design inputs of 2026-11-17, including the ion drift velocities.</li>
<li>Gas density: NEXT DocDB 1183 (xenon compressibility); isotopic fractions of NEXUS <code>GXeEnriched</code>.</li>
<li>Ion mobility: P. N. B. Neves et al., J. Chem. Phys. 133 (2010) 124316, as plotted in A. F. V. Cortez et al., JINST 12 (2017) P06012, Fig. 3.</li>
<li>Ion thermal spread: ITACA revisited, arXiv:2604.02443, Eq. (2.1); implementation in <code>src/event_analysis/production.jl</code>.</li>
<li>Electron diffusion in pure xenon: NEXT-DEMO++ paper, <code>papers/xehe_paper.pdf</code>, Table 1. EL yield: C. M. B. Monteiro et al., JINST 2 (2007) P05001, as used in Section 4.4 of that paper.</li>
<li>Potentials of the ion plane: <code>manuscripts/itaca3/itaca3.tex</code> (IFG at 8 kV, sensors at ground).</li>
<li>Energy resolution: full simulation of the Xe/He electron readout (<a href="../../results/readout/itaca_100k_ion15_comparison/index.html">comparison of bb0nu and single electrons</a>). Q<sub>ββ</sub>: M. Redshaw et al., PRL 98 (2007) 053003.</li>
</ul></section></main><footer>iItaca · Detector and operating conditions</footer></body></html>
'''
    PAGE.write_text(page)
    print(PAGE.relative_to(ROOT))


if __name__ == '__main__':
    render()
