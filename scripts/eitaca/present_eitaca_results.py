#!/usr/bin/env python3
"""Results page of eITACA IE_L1500_D400 (etrack) beside iItaca II_L1000_D500_P2 (itrack, 2 mm pads).

  python scripts/eitaca/present_eitaca_results.py [notes/html/eitaca/results.html]

Reads, without recalculating any physical quantity:
  results/eitaca/IE_L1500_D400_results/results.h5          cut flows, network, pooled gamma window, hit noise
  results/eitaca/IE_L1500_D400_comparison/comparison.h5    cut on E2
  results/eitaca/IE_L1500_D400_bi214_truth/summary.h5      truth classes of the selected 214Bi events
  results/iitaca/II_L1000_D500_P2_<sample>/statistics/events.h5, II_L1000_D500_P2_tl208_w (208Tl),
  results/network/II_L1000_D500_P2/trainings/itrack_sparse_e100/ (training.h5, scores_bi214.h5),
  results/iitaca/II_L1000_D500_P2_gamma_window/gamma_window.h5, II_L1000_D500_P2_comparison/comparison.h5,
  results/iitaca/II_L1000_D500_P2_bi214_truth/summary.h5   the iItaca 2 mm numbers of notes/html/iitaca/results.html

Cut flow with the ROI last, as the iItaca page: saved by Geant4, fast trigger,
track selection, network at 80 % signal, ROI. The network fraction is measured
on selected events and applied to all selected events; the ROI fraction is the
fraction of the selected events whose energy is in the ROI.
"""
import hashlib
import sys
import tomllib
from pathlib import Path
import h5py
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'iitaca'))
from present_iitaca_detector import STYLE, table, text

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'notes/html/eitaca/results.html'
SAMPLES = {'bb0nu': 'bb0nu', 'single_electrons': 'single electrons', 'bi214': '²¹⁴Bi', 'tl208': '²⁰⁸Tl'}
UNIT = {'bb0nu': 'event', 'single_electrons': 'event', 'bi214': 'gamma', 'tl208': 'gamma'}
I, E = 'iItaca 2 mm', 'eITACA'
INPUTS = []


def h5(path):
    p = ROOT / path
    INPUTS.append((path, hashlib.sha256(p.read_bytes()).hexdigest()))
    return h5py.File(p, 'r')


def n(x):
    return f'{int(round(x)):,}'


def pct(x, digits=2, e=None):
    return f'{100 * x:.{digits}f} ± {100 * e:.{digits}f} %' if e is not None else f'{100 * x:.{digits}f} %'


def sci(x):
    if x == 0:
        return '0'
    m, ex = f'{x:.3e}'.split('e')
    return f'{m} × 10<sup>{int(ex)}</sup>'


def rej(x):
    return f'{1 / x:,.4g}' if x > 0 else '—'


def binom(k, m):
    a = k / m
    return a, np.sqrt(a * (1 - a) / m)


def i80(values):
    return int(np.argmin(np.abs(np.asarray(values) - 0.8)))


def ev(x):
    return f'{x:,.1f}' if abs(x - round(x)) > 1e-6 else n(x)


# ---- iItaca 2 mm, as notes/html/iitaca/results.html -----------------------------------------------------
D = {I: {}, E: {}}
run = 'II_L1000_D500_P2'
tdir = f'results/network/{run}/trainings/itrack_sparse_e100'
tr = h5(f'{tdir}/training.h5')
sep = tr['test/separation']
j = i80(sep['signal_efficiency'][()])
net = {'bb0nu': (int(round(sep['signal_efficiency'][j] * sep['signal_events'][()])), int(sep['signal_events'][()])),
       'single_electrons': (int(round(sep['background_acceptance'][j] * sep['background_events'][()])), int(sep['background_events'][()]))}
g = h5(f'{tdir}/scores_bi214.h5')
net['bi214'] = (int(g['accepted_events'][i80(g['target_signal_efficiency'][()])]), len(g['events/score']))
gw = h5(f'results/iitaca/{run}_gamma_window/gamma_window.h5')
jw = i80(gw['target_signal_efficiency'][()])
net['tl208'] = (int(gw['tl208/accepted_events'][jw]), int(gw['tl208/window_events'][()]))
D[I]['area'] = float(sep['area'][()])
for s in SAMPLES:
    f = h5(f'results/iitaca/{run}_{s}/statistics/events.h5' if s != 'tl208' else f'results/iitaca/{run}_tl208_w/statistics/events.h5')
    flow = dict(zip([text(x) for x in f['cut_flow/step'][()]], f['cut_flow/events'][()]))
    sel, roi = f['events/selected'][()].astype(bool), f['events/in_roi'][()].astype(bool)
    k, m = net[s]
    D[I][s] = dict(generated=int(flow['generated']), read=int(flow['read']), trigger=int(flow['trigger']), selected=int(sel.sum()),
                   sel_roi=int((sel & roi).sum()), k=k, m=m)
D[I]['notes'] = {'bi214': 'selected events in the ROI', 'tl208': 'selected events with box energy in 2420–2490 keV (40 M and 100 M runs)'}

# ---- eITACA ------------------------------------------------------------------------------------------------
r = h5('results/eitaca/IE_L1500_D400_results/results.h5')
jt = i80(r['target_signal_efficiency'][()])
t = r['network']
D[E]['area'] = float(t['area'][()])
pk, pm = int(r['pooled/accepted'][jt]), int(r['pooled/events'][()])
neteI = {'bb0nu': (int(round(t['signal_efficiency'][jt] * t['signal_events'][()])), int(t['signal_events'][()])),
         'single_electrons': (int(t['background_accepted'][jt]), int(t['background_events'][()])),
         'bi214': (pk, pm), 'tl208': (pk, pm)}
for s in SAMPLES:
    gs = r[s]
    k, m = neteI[s]
    D[E][s] = dict(generated=int(gs['flow/generated'][()]), read=int(gs['flow/read'][()]), trigger=int(gs['flow/trigger'][()]),
                   selected=int(gs['selected'][()]), sel_roi=int(gs['selected_in_roi'][()]), k=k, m=m)
    assert D[E][s]['sel_roi'] == int(gs['flow/track'][()])
WIN = {s: (int(r[f'{s}/window/accepted'][jt]), int(r[f'{s}/window/events'][()])) for s in ('bi214', 'tl208')}
D[E]['notes'] = {'bi214': 'selected ²¹⁴Bi and ²⁰⁸Tl events with box energy in 2420–2470 keV, pooled',
                 'tl208': 'selected ²¹⁴Bi and ²⁰⁸Tl events with box energy in 2420–2470 keV, pooled'}
NOISE = {k: r[f'bb0nu/noise/{k}'][()] for k in ('median_keV', 'max_keV', 'dark_per_bucket_median', 'asic_rms_pe', 'rebin_factor', 'keV_per_pe', 'hit_cut_keV', 'hit_cut_sigma')}
HPM = float(r['bb0nu/hits_per_mm_median'][()])
with h5py.File(ROOT / 'results/eitaca/IE_L1500_D400_bb0nu/statistics/events.h5') as f:
    CARD = tomllib.loads(text(f['settings'][()]))
TK = CARD['tracks']


# ---- cut flow, ROI last -------------------------------------------------------------------------------------
def flow_rows(d):
    a = d['k'] / d['m']
    steps = [('Generated in Geant4', d['generated']), ('Saved by Geant4', d['read']), ('Fast trigger', d['trigger']),
             ('Track selection', d['selected']), ('Network at 80 % signal', d['selected'] * a), ('ROI', d['sel_roi'] * a)]
    rows, prev = [], None
    for name, x in steps:
        rows.append([name, n(x), '', '', '1'] if prev is None else [name, ev(x), pct(x / prev), rej(x / prev), sci(x / d['generated'])])
        prev = x
    return rows, d['sel_roi'] * a


def note(arm, s):
    k, m = D[arm][s]['k'], D[arm][s]['m']
    if s == 'bb0nu':
        return f'network: {pct(k / m)} of the {m:,} bb0nu test events (threshold set at 80 %)'
    if s == 'single_electrons':
        return f'network: {k} of the {m:,} single-electron test events, {pct(*binom(k, m)[:1], 2, binom(k, m)[1])}'
    return f'network: {k} of the {m:,} {D[arm]["notes"][s]}, {pct(binom(k, m)[0], 2, binom(k, m)[1])}'


HEAD = ['Cut', 'Events', 'Fraction of the previous cut', 'Rejection factor', 'Fraction of generated']
flow_html, FINAL = [], {}
for arm, title in ((I, 'iItaca 2 mm pads (II_L1000_D500_P2), itrack'), (E, 'eITACA (IE_L1500_D400), etrack')):
    parts = []
    for s, name in SAMPLES.items():
        rows, final = flow_rows(D[arm][s])
        FINAL[(arm, s)] = final
        parts.append(f'<h3>{name}</h3><p class="meta">{note(arm, s)}</p><div class="scroll">{table(HEAD, rows).replace("<table>", "<table class=num>")}</div>')
    flow_html.append(f'<h3 class="pitch">{title}</h3>' + ''.join(parts))

# ---- track selection × network --------------------------------------------------------------------------------
tsn = []
for s in ('bb0nu', 'single_electrons', 'bi214'):     # 208Tl: the saved windows of the two productions differ
    for arm in (I, E):
        d = D[arm][s]
        ts, a = d['selected'] / d['trigger'], d['k'] / d['m']
        tsn.append([SAMPLES[s], arm, f'{ts:.3f}', f'{a:.4f}' + (f' ({d["k"]} of {d["m"]:,})' if s != 'bb0nu' else ''), f'{ts * a:.4f}', rej(ts * a)])

# ---- final numbers ---------------------------------------------------------------------------------------------
fin = [[name] + [f'{ev(FINAL[(arm, s)])} of {D[arm][s]["generated"]:,}<br>{sci(FINAL[(arm, s)] / D[arm][s]["generated"])} per {UNIT[s]}' for arm in (I, E)]
       for s, name in SAMPLES.items()]

# ---- cut on E2 -------------------------------------------------------------------------------------------------
def e2_rows(path):
    c = h5(path)
    k = i80(c['targets'][()])
    out = {'threshold': float(c['thresholds_keV'][k])}
    for b in c['backgrounds']:
        g = c[f'backgrounds/{b}']
        out[text(g['sample'][()])] = (float(g['acceptance'][k]), float(g['acceptance_error'][k]), int(g['accepted_events'][k]), len(g['E2_keV']))
    return out


E2 = {I: e2_rows('results/iitaca/II_L1000_D500_P2_comparison/comparison.h5'), E: e2_rows('results/eitaca/IE_L1500_D400_comparison/comparison.h5')}
e2 = [['E₂ threshold at 80 % signal'] + [f'{E2[arm]["threshold"]:.1f} keV' for arm in (I, E)]]
for s in ('single_electrons', 'bi214', 'tl208'):
    e2.append([SAMPLES[s]] + [(lambda a: f'{a[2]} of {a[3]:,}: {pct(a[0], 2, a[1])}')(E2[arm][s]) for arm in (I, E)])

# ---- truth -----------------------------------------------------------------------------------------------------
T = {}
for arm, path in ((I, 'results/iitaca/II_L1000_D500_P2_bi214_truth/summary.h5'), (E, 'results/eitaca/IE_L1500_D400_bi214_truth/summary.h5')):
    f = h5(path)
    T[arm] = {g: {x: f[f'bi214/{g}/{x}'][()] for x in ('events', 'passed')} for g in ('by_class', 'by_history', 'by_second_distance', 'all')}
    T[arm]['classes'] = [text(c) for c in f['classes'][()]]
    T[arm]['se'] = (int(f['single_electrons/all/passed'][0]), int(f['single_electrons/all/events'][0]))
CLASS = {'electron_only': 'no secondary deposit of 10 keV or more', 'fluorescence': 'fluorescence photon of the gamma',
         'second_electron': 'second electron (Compton class)', 'photon': 'other photon (main-track bremsstrahlung, …)'}


def acc(G, i):
    k, m = int(G['passed'][i]), int(G['events'][i])
    return f'{k} of {m:,}: {pct(k / m, 2, np.sqrt(k / m * (1 - k / m) / m))}' if m else '—'


cls = [[CLASS[c]] + [acc(T[arm]['by_class'], i) for arm in (I, E)] for i, c in enumerate(T[I]['classes'])]
cls.append(['all selected ²¹⁴Bi in the ROI'] + [acc(T[arm]['all'], 0) for arm in (I, E)])
cls.append(['single electrons (test sample)'] + [f'{T[arm]["se"][0]} of {T[arm]["se"][1]:,}: {pct(binom(*T[arm]["se"])[0], 2, binom(*T[arm]["se"])[1])}' for arm in (I, E)])
hist = [[h] + [acc(T[arm]['by_history'], i) for arm in (I, E)] for i, h in ((1, 'photoelectric absorption only'), (2, 'Compton scattering(s), then absorption'))]
dist = [[lab] + [acc(T[arm]['by_second_distance'], i) for arm in (I, E)] for i, lab in enumerate(('&lt; 4 mm', '4–8 mm', '&gt; 8 mm'))]

# ---- page ------------------------------------------------------------------------------------------------------
up = '../../../'
settings = table(['', 'iItaca 2 mm (II_L1000_D500_P2)', 'eITACA (IE_L1500_D400)'], [
    ['Readout', 'ions on 2 mm pads, 10 ms frames: cubic 2 mm cells', 'electrons on SiPMs of 4 mm pitch, 20 × 250 ns buckets at 0.49 mm/µs: 4 × 4 × 2.45 mm cells'],
    ['Transverse spread of the carriers', 'ion, thermal, up to 0.9 mm', 'electron, up to σ<sub>T</sub> 3.4 mm, σ<sub>L</sub> 4.6 mm (1.5 m drift)'],
    ['Hit cut', '3 keV', f'{TK["hit_cut_keV"]:g} keV, the saved cut: {NOISE["hit_cut_sigma"]:.0f} σ of the noise of one hit'],
    ['Components (touching cells)', '3.7 mm', f'{TK["segment_connection_mm"]:g} mm'],
    ['Segment: hits of 18 mm of track', '9 hits', f'{TK["segment_min_hits"]} hits ({HPM:.2f} hits per mm)'],
    ['Segments joined (across one empty cell)', '5.3 mm', f'{TK["connection_mm"]:g} mm'],
    ['Detached deposit; segments', '10 keV; at most 3', f'{TK["deposit_min_keV"]:g} keV; at most {TK["max_segments"]}'],
    ['Track graph cells; end-point spheres', '2 mm; 7 mm', f'{TK["voxel_width_mm"]:g} mm; {TK["blob_radius_mm"]:g} mm'],
    ['Network', f'itrack, sparse, 100 epochs; area {D[I]["area"]:.4f}', f'etrack, sparse, 100 epochs; area {D[E]["area"]:.4f}'],
    ['Gamma samples', '²¹⁴Bi 40 M; ²⁰⁸Tl 100 M saved in 2420–2490 keV', '²¹⁴Bi 40 M and ²⁰⁸Tl 40 M, saved in 2300–2470 keV'],
])
page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>eITACA results</title><style>{STYLE}
table.num td:not(:first-child){{text-align:right;white-space:nowrap}}table.num th:not(:first-child){{text-align:right}}.meta{{color:#526776;font-size:14px;margin:4px 0 8px}}
h3.pitch{{margin-top:30px;padding-top:14px;border-top:2px solid #185b83}}.scroll{{overflow-x:auto}}td{{font-size:14px}}</style></head><body>
<header><a href="{up}notes/html/index.html">Documentation index</a> · <a href="index.html">eITACA</a> · <a href="../iitaca/results.html">iItaca results</a><h1>eITACA IE_L1500_D400 against iItaca 2 mm: results</h1>
<p>The same analysis on the electron track of eITACA (Xe with 15.8 % He, 16 bar, 1.5 m drift at 400 V/cm, SiPMs of 4 mm pitch) and on the ion track of iItaca II_L1000_D500_P2 (pure Xe, 16 bar, 1.0 m drift, 2 mm pads). The event energy is the energy in the simulated box smeared with 0.472 % FWHM at Q<sub>ββ</sub> in both, so that the comparison isolates the topology.</p></header><main>
<nav><a href="#settings">Settings</a><a href="#flow">Cut flow</a><a href="#product">Track selection × network</a><a href="#final">Final numbers</a><a href="#e2">Cut on E₂</a><a href="#truth">²¹⁴Bi truth</a><a href="#links">Pages</a></nav>

<section id="settings"><h2>Settings</h2>
{settings}
<p>The eITACA track settings are those of <code>config/eitaca/tracks.toml</code>, set by the rules of the iItaca cards in cells of the etrack lattice. The components radius is just above the cell diagonal, 6.17 mm; on these flat cells it also joins cells two buckets apart in z (4.9 mm). The minimum segment counts the hits of 18 mm of track, measured on the selected bb0nu tracks in the ROI: diffusion spreads the track over several SiPMs, so a segment needs more hits than on the pad lattice. The join radius is 8 mm × 4 mm / 3 mm. Noise of one hit (one SiPM over {int(NOISE["rebin_factor"])} buckets: ASIC {NOISE["asic_rms_pe"]:g} PE per bucket, median {NOISE["dark_per_bucket_median"]:.3f} dark pulses per bucket, {1000 * NOISE["keV_per_pe"]:.3f} eV per PE): median {1000 * NOISE["median_keV"]:.1f} eV, at most {1000 * NOISE["max_keV"]:.1f} eV; the noise sets no cut above the saved one. The ROI is one FWHM around Q<sub>ββ</sub> (2452.0–2463.6 keV). The network threshold gives 80 % of the bb0nu test events of each training.</p>
<p>eITACA gamma backgrounds: the network acceptance is measured on the selected ²¹⁴Bi and ²⁰⁸Tl events with box energy in 2420–2470 keV, pooled (²¹⁴Bi alone {WIN["bi214"][0]} of {WIN["bi214"][1]:,}, {pct(binom(*WIN["bi214"])[0], 2, binom(*WIN["bi214"])[1])}; ²⁰⁸Tl alone {WIN["tl208"][0]} of {WIN["tl208"][1]:,}, {pct(binom(*WIN["tl208"])[0], 2, binom(*WIN["tl208"])[1])}), and applied to the selected events of each isotope. The ²⁰⁸Tl files hold no event above 2470 keV of deposited energy; those can reach the ROI only through the smearing, more than 1.3 σ, and are missing.</p></section>

<section id="flow"><h2>Cut flow, ROI last</h2>
<p>The network fraction is measured on the selected events and applied to all selected events: over this energy range the topology does not depend on the energy. For bb0nu and single electrons it is that of the test sample of the training. The ROI fraction is the fraction of the selected events whose energy is in the ROI.</p>
{''.join(flow_html)}</section>

<section id="product"><h2>Track selection × network</h2>
<p>Track selection as a fraction of the events accepted by the fast trigger; network at 80 % signal.</p>
<div class="scroll">{table(['Sample', 'Readout', 'Track selection', 'Network', 'Product', 'Rejection'], tsn).replace('<table>', '<table class=num>')}</div></section>

<section id="final"><h2>Final numbers</h2>
<p>Events after all cuts (saved, fast trigger, track selection, network at 80 % signal, ROI) and per generated event or gamma.</p>
<div class="scroll">{table(['Sample', I, E], fin).replace('<table>', '<table class=num>')}</div></section>

<section id="e2"><h2>Cut on E₂</h2>
<p>Selected events in the ROI; an event is accepted when the energy E₂ of the less energetic end-point sphere exceeds the threshold that keeps 80 % of the bb0nu events.</p>
<div class="scroll">{table(['', I, E], e2).replace('<table>', '<table class=num>')}</div></section>

<section id="truth"><h2>²¹⁴Bi truth</h2>
<p>Selected ²¹⁴Bi events in the ROI, classified from the Monte Carlo truth (<code>summarize_iitaca_truth.jl</code>); network at 80 % signal. Second electron: the most energetic Compton electron, other photoelectron or positron with at least 10 keV in the box; distance: its smallest distance to the main track.</p>
<h3>Network acceptance by truth class</h3>
<div class="scroll">{table(['Class of the selected ²¹⁴Bi event', I, E], cls).replace('<table>', '<table class=num>')}</div>
<h3>By interaction history of the gamma</h3>
<div class="scroll">{table(['History', I, E], hist).replace('<table>', '<table class=num>')}</div>
<h3>Compton class by distance of the second electron to the main track</h3>
<div class="scroll">{table(['Distance', I, E], dist).replace('<table>', '<table class=num>')}</div></section>

<section id="links"><h2>Pages</h2><ul>
<li>Cut on E₂: <a href="{up}results/eitaca/IE_L1500_D400_comparison/index.html">eITACA</a> · <a href="{up}results/iitaca/II_L1000_D500_P2_comparison/index.html">iItaca 2 mm</a></li>
<li>Trainings: <a href="{up}results/network/IE_L1500_D400/trainings/index.html">eITACA</a> · <a href="{up}results/network/II_L1000_D500_P2/trainings/index.html">iItaca 2 mm</a></li>
<li><a href="../iitaca/results.html">iItaca results</a> · <a href="../iitaca/analysis.html">iItaca analysis</a> (the same steps)</li></ul></section>

<section id="inputs"><h2>Inputs</h2><p class="meta">Numbers: <code>julia --project=. scripts/eitaca/prepare_eitaca_results.jl</code>; page: <code>python scripts/eitaca/present_eitaca_results.py</code>. SHA-256 of the files read:</p>
<div class="scroll">{table(['File', 'SHA-256'], [[f'<code>{p}</code>', f'<code>{h[:16]}…</code>'] for p, h in INPUTS])}</div></section>
</main><footer>eITACA · results · IE_L1500_D400</footer></body></html>
'''
OUT.write_text(page)
print('wrote', OUT)
