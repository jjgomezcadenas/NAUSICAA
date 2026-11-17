#!/usr/bin/env python3
"""Results page of iItaca II_L1000_D500 on 3 mm and on 2 mm pads.

  python scripts/iitaca/present_iitaca_results.py [notes/html/iitaca/results.html]

Reads, without recalculating any physical quantity:
  results/iitaca/<run>_<sample>/statistics/events.h5     cut flow, columns selected and in_roi
  results/network/<run>/trainings/itrack_sparse_e100/     training.h5 (bb0nu, single electrons),
                                                          scores_bi214.h5, scores_tl208.h5
  results/iitaca/<run>_bi214_truth/summary.h5             truth classes of the selected 214Bi events
  results/iitaca/II_L1000_D500_bi214_floor/floor.h5       drift thirds and voxelised Monte Carlo (3 mm)
  results/iitaca/II_L1000_D500_pitch_comparison/comparison.h5   cross-check of the counts
for <run> = II_L1000_D500 (3 mm) and II_L1000_D500_P2 (2 mm, radii scaled), and writes the page.

Cut flow with the ROI last: saved by Geant4, fast trigger, track selection,
network at 80 % signal, ROI. The network fraction is measured on the selected
events in the ROI and applied to all selected events.
"""
import hashlib
import sys
from pathlib import Path
import h5py
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from present_iitaca_detector import STYLE, table, text

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'notes/html/iitaca/results.html'
RUNS = {'II_L1000_D500': '3 mm', 'II_L1000_D500_P2': '2 mm'}
SAMPLES = {'bb0nu': 'bb0nu', 'single_electrons': 'single electrons', 'bi214': '²¹⁴Bi', 'tl208': '²⁰⁸Tl'}
UNIT = {'bb0nu': 'event', 'single_electrons': 'event', 'bi214': 'gamma', 'tl208': 'gamma'}
POISSON_UL90_ZERO = -np.log(0.10)          # 2.303 events: 90 % CL upper limit for 0 observed, no background
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
    m, e = f'{x:.3e}'.split('e')
    return f'{m} × 10<sup>{int(e)}</sup>'


def rej(x):
    return f'{1 / x:,.4g}' if x > 0 else '—'


def binom(k, m):
    a = k / m
    return a, np.sqrt(a * (1 - a) / m)


# ---- read ---------------------------------------------------------------------------------------
pitch = h5('results/iitaca/II_L1000_D500_pitch_comparison/comparison.h5')
D = {}
for run in RUNS:
    tdir = f'results/network/{run}/trainings/itrack_sparse_e100'
    tr = h5(f'{tdir}/training.h5')
    sep = tr['test/separation']
    i80 = int(np.argmin(np.abs(sep['signal_efficiency'][()] - 0.8)))
    lab, score = tr['test/label'][()], tr['test/score'][()]
    net = {'bb0nu': (int(round(sep['signal_efficiency'][i80] * sep['signal_events'][()])), int(sep['signal_events'][()])),
           'single_electrons': (int(round(sep['background_acceptance'][i80] * sep['background_events'][()])), int(sep['background_events'][()]))}
    for s in ('bi214', 'tl208'):
        g = h5(f'{tdir}/scores_{s}.h5')
        assert abs(g['target_signal_efficiency'][i80] - 0.8) < 1e-9
        net[s] = (int(g['accepted_events'][i80]), len(g['events/score']))
    D[run] = {'net': net, 'area': float(sep['area'][()])}
    for s in SAMPLES:
        f = h5(f'results/iitaca/{run}_{s}/statistics/events.h5')
        flow = dict(zip([text(x) for x in f['cut_flow/step'][()]], f['cut_flow/events'][()]))
        sel, roi = f['events/selected'][()].astype(bool), f['events/in_roi'][()].astype(bool)
        assert len(sel) == flow['trigger'] and roi.sum() == flow['roi'] and (sel & roi).sum() == flow['track']
        cmp = dict(zip([text(x) for x in pitch[f'{run}/flow/{s}/step'][()]], pitch[f'{run}/flow/{s}/events'][()]))
        assert cmp == flow, (run, s)
        k, m = net[s]
        assert s in ('bb0nu', 'single_electrons') or m == flow['track'], (run, s)
        D[run][s] = dict(generated=flow['generated'], read=flow['read'], trigger=flow['trigger'], selected=int(sel.sum()),
                         sel_roi=int((sel & roi).sum()), roi=flow['roi'], k=k, m=m)
    t = h5(f'results/iitaca/{run}_bi214_truth/summary.h5')
    D[run]['truth'] = {g: {x: t[f'bi214/{g}/{x}'][()] for x in ('events', 'passed', 'fraction', 'error')}
                       for g in ('by_class', 'by_history', 'by_second_distance', 'all')}
    D[run]['classes'] = [text(c) for c in t['classes'][()]]
    D[run]['se_truth'] = (int(t['single_electrons/all/passed'][0]), int(t['single_electrons/all/events'][0]))
    assert D[run]['truth']['all']['passed'][0] == D[run]['bi214']['k']
# 208Tl at 2 mm: the 100 M production saved in 2420-2490 keV (II_L1000_D500_P2_tl208_w) for the cut flow, and the
# network acceptance measured on the selected events with box energy in 2420-2490 keV, pooled over the 40 M and 100 M runs.
gw = h5('results/iitaca/II_L1000_D500_P2_gamma_window/gamma_window.h5')
i80w = int(np.argmin(np.abs(gw['target_signal_efficiency'][()] - 0.8)))
f = h5('results/iitaca/II_L1000_D500_P2_tl208_w/statistics/events.h5')
flow = dict(zip([text(x) for x in f['cut_flow/step'][()]], f['cut_flow/events'][()]))
sel, roi = f['events/selected'][()].astype(bool), f['events/in_roi'][()].astype(bool)
assert int((sel & roi).sum()) == int(gw['tl208/roi_selected_events'][()])
D['II_L1000_D500_P2']['tl208'] = dict(generated=flow['generated'], read=flow['read'], trigger=flow['trigger'], selected=int(sel.sum()),
    sel_roi=int((sel & roi).sum()), roi=flow['roi'], k=int(gw['tl208/accepted_events'][i80w]), m=int(gw['tl208/window_events'][()]), window=True)
fl = h5('results/iitaca/II_L1000_D500_bi214_floor/floor.h5')
TH = {x: fl[f'thirds/II_L1000_D500/{x}'][()] for x in fl['thirds/II_L1000_D500'].keys()}
BT = {x: fl[f'bi214_thirds/{x}'][()] for x in fl['bi214_thirds'].keys()}
VX = {x: fl[f'voxelised/{x}'][()] for x in fl['voxelised'].keys()}
tr_ev, tr_cls, tr_ps = fl['truth/event_number'][()], fl['truth/class'][()], fl['truth/passed'][()]


# ---- cut flow, ROI last ----------------------------------------------------------------------------
def flow_rows(r, s):
    """Steps, events, fraction of the previous step, cumulative fraction of the generated events."""
    d = D[r][s]
    k, m = d['k'], d['m']
    upper = (k == 0)
    a = (POISSON_UL90_ZERO / m) if upper else k / m
    net_events = d['selected'] * a
    final = d['sel_roi'] * a                        # = selected × network × (selected in ROI / selected)
    steps = [('Generated in Geant4', d['generated']), ('Saved by Geant4', d['read']), ('Fast trigger', d['trigger']),
             ('Track selection', d['selected']), ('Network at 80 % signal', net_events), ('ROI', final)]
    rows, prev = [], None
    for i, (name, ev) in enumerate(steps):
        lt = '&lt; ' if upper and i >= 4 else ''
        if prev is None:
            rows.append([name, n(ev), '', '', '1'])
        else:
            frac = ev / prev
            evs = lt + (f'{ev:,.1f}' if abs(ev - round(ev)) > 1e-6 else n(ev))
            rows.append([name, evs, lt + pct(frac), ('&gt; ' if lt else '') + rej(frac), lt + sci(ev / d['generated'])])
        prev = ev
    return rows, final, upper


def network_note(r, s):
    k, m = D[r][s]['k'], D[r][s]['m']
    if s == 'bb0nu':
        return f'network: {pct(k / m)} of the {m:,} bb0nu test events (threshold set at 80 %)'
    if s == 'single_electrons':
        a, e = binom(k, m)
        return f'network: {k} of the {m:,} single-electron test events, {pct(a, 2, e)}'
    if D[r][s].get('window'):
        a, e = binom(k, m)
        return (f'100 M gammas saved in 2420–2490 keV; network: {k} of the {m:,} selected events with energy in 2420–2490 keV '
                f'(40 M and 100 M runs), {pct(a, 2, e)}')
    if k == 0:
        return f'network: 0 of the {m:,} selected events in the ROI; 90 % CL upper limit {POISSON_UL90_ZERO:.2f} events (Poisson, no background)'
    a, e = binom(k, m)
    return f'network: {k} of the {m:,} selected events in the ROI, {pct(a, 2, e)}'


FLOW_HEAD = ['Cut', 'Events', 'Fraction of the previous cut', 'Rejection factor', 'Fraction of generated']
flow_html, FINAL = [], {}
for r, label in RUNS.items():
    parts = []
    for s, sname in SAMPLES.items():
        rows, final, upper = flow_rows(r, s)
        FINAL[(r, s)] = (final, upper)
        extra = ''
        if s == 'tl208' and upper:
            extra = '<p class="meta">At 2 mm the network acceptance is measured on 2,150 selected events in 2420–2490 keV.</p>'
        parts.append(f'<h3>{sname}</h3><p class="meta">{network_note(r, s)}</p>'
                     f'<div class="scroll">{table(FLOW_HEAD, rows).replace("<table>", "<table class=num>")}</div>{extra}')
    flow_html.append(f'<h3 class="pitch">{label} pads ({r})</h3>' + ''.join(parts))


# ---- track selection × network ---------------------------------------------------------------------
tsn_rows = []
for s in ('bb0nu', 'single_electrons', 'bi214'):
    for r, label in RUNS.items():
        d = D[r][s]
        ts = d['selected'] / d['trigger']
        a = d['k'] / d['m']
        tsn_rows.append([SAMPLES[s], label, f'{ts:.3f}', f'{a:.4f}' + (f' ({d["k"]} of {d["m"]:,})' if s != 'bb0nu' else ''),
                         f'{ts * a:.4f}', rej(ts * a)])

# ---- final numbers -------------------------------------------------------------------------------------
fin_rows = []
for s, sname in SAMPLES.items():
    cells = [sname]
    for r in RUNS:
        final, upper = FINAL[(r, s)]
        g = D[r][s]['generated']
        if upper:
            cells.append(f'&lt; {final:.2f} of {g:,} (90 % CL)<br>&lt; {sci(final / g)} per {UNIT[s]}')
        else:
            ev = f'{final:,.1f}' if abs(final - round(final)) > 1e-6 else n(final)
            cells.append(f'{ev} of {g:,}<br>{sci(final / g)} per {UNIT[s]}')
    fin_rows.append(cells)

# ---- truth ---------------------------------------------------------------------------------------------------
CLASS_NAME = {'electron_only': 'no secondary deposit of 10 keV or more', 'fluorescence': 'fluorescence photon of the gamma',
              'second_electron': 'second electron (Compton class)', 'photon': 'other photon (main-track bremsstrahlung, …)'}


def acc(T, i):
    k, m = int(T['passed'][i]), int(T['events'][i])
    return f'{k} of {m:,}: {pct(k / m, 2, np.sqrt(k / m * (1 - k / m) / m))}' if m else '—'


r3, r2 = list(RUNS)
cls_rows = [[CLASS_NAME[c], acc(D[r3]['truth']['by_class'], i), acc(D[r2]['truth']['by_class'], i)] for i, c in enumerate(D[r3]['classes'])]
cls_rows.append(['all selected ²¹⁴Bi', acc(D[r3]['truth']['all'], 0), acc(D[r2]['truth']['all'], 0)])
for r in RUNS:
    k, m = D[r]['net']['single_electrons']
    assert D[r]['se_truth'] == (k, m)
cls_rows.append(['single electrons (test sample)'] + [f"{D[r]['net']['single_electrons'][0]} of {D[r]['net']['single_electrons'][1]:,}: "
                                                      f"{pct(*(lambda a: (a[0], 2, a[1]))(binom(*D[r]['net']['single_electrons'])))}" for r in RUNS])
hist_rows = [[h, acc(D[r3]['truth']['by_history'], i), acc(D[r2]['truth']['by_history'], i)]
             for i, h in ((1, 'photoelectric absorption only'), (2, 'Compton scattering(s), then absorption'))]
dist_rows = [[lab, acc(D[r3]['truth']['by_second_distance'], i), acc(D[r2]['truth']['by_second_distance'], i)]
             for i, lab in enumerate(('&lt; 4 mm', '4–8 mm', '&gt; 8 mm'))]
TD = {r: D[r]['truth']['by_second_distance'] for r in RUNS}
TH2 = {r: D[r]['truth']['by_history'] for r in RUNS}

third_rows = []
for i in range(3):
    third_rows.append([('first', 'second', 'third')[i], f"{TH['edges_cm'][i]:.0f}–{TH['edges_cm'][i + 1]:.0f}",
                       f"{TH['sigma_median_mm'][i]:.2f}",
                       f"{TH['passed'][i]} of {TH['events'][i]:,}: {pct(TH['fraction'][i], 2, TH['error'][i])}",
                       f"{BT['edges_cm'][i]:.0f}–{BT['edges_cm'][i + 1]:.0f}",
                       f"{BT['passed'][i]} of {BT['selected'][i]:,}: {pct(*(lambda a: (a[0], 2, a[1]))(binom(BT['passed'][i], BT['selected'][i])))}",
                       f"{1000 * BT['passed'][i] / BT['roi'][i]:.2f}"])
vs, rs = VX['selected'] == 1, VX['readout_selected'] == 1
idx = dict(zip(VX['event_number'], range(len(vs))))
cv = vs[np.array([idx[e] for e in tr_ev])]
comp = tr_cls == 3
vox_rows = [['²¹⁴Bi events in the ROI', f'{len(vs):,}', f'{len(vs):,}'],
            ['selected', f'{rs.sum():,}', f'{vs.sum():,}'],
            ['Compton class among the full-readout selection, kept', f'{comp.sum():,}', f'{cv[comp].sum():,} ({100 * cv[comp].mean():.0f} %)'],
            ['… of which passing the network, kept', f'{(comp & (tr_ps == 1)).sum():,}', f'{cv[comp & (tr_ps == 1)].sum():,}']]

# ---- page --------------------------------------------------------------------------------------------------------
up = '../../../'
pc = f'{up}results/iitaca/II_L1000_D500_pitch_comparison/index.html'
stat = lambda r, s: f'<a href="{up}results/iitaca/{r}_{s}/statistics/index.html">{SAMPLES[s]}</a>'
a3 = TD[r3]['fraction'][0]; a2 = TD[r2]['fraction'][0]
page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>iItaca results</title><style>{STYLE}
table.num td:not(:first-child){{text-align:right;white-space:nowrap}}table.num th:not(:first-child){{text-align:right}}.meta{{color:#526776;font-size:14px;margin:4px 0 8px}}
h3.pitch{{margin-top:30px;padding-top:14px;border-top:2px solid #185b83}}.pending{{color:#a3261b;font-weight:650}}.scroll{{overflow-x:auto}}td{{font-size:14px}}</style></head><body>
<header><a href="{up}notes/html/index.html">Documentation index</a> · <a href="index.html">iItaca</a><h1>iItaca II_L1000_D500: results</h1>
<p>Pure xenon at 16 bar, 1.0 m drift at 500 V/cm. Ions read on 3 mm pads with 15 ms frames, and the same Geant4 events, seeds and ion drift on 2 mm pads with 10 ms frames.</p></header><main>
<nav><a href="#settings">Settings</a><a href="#flow">Cut flow</a><a href="#product">Track selection × network</a><a href="#final">Final numbers</a><a href="#floor">The ²¹⁴Bi floor</a><a href="#links">Pages</a></nav>

<section id="settings"><h2>Settings</h2>
{table(['Pitch', 'Run', 'Hit cut', 'Components', 'Segments', 'Segments joined', 'Network'], [
    ['3 mm', 'II_L1000_D500', '5 keV', '5.5 mm', '6 hits', '8 mm', f'itrack, sparse, 100 epochs; area {D[r3]["area"]:.4f}'],
    ['2 mm', 'II_L1000_D500_P2', '3 keV', '3.7 mm', '9 hits', '5.3 mm', f'itrack, sparse, 100 epochs; area {D[r2]["area"]:.4f}']])}
<p>Both: detached deposits of 10 keV rejected, at most three segments, end-point spheres of 7 mm. The ROI is one FWHM around Q<sub>ββ</sub> (2452.0–2463.6 keV). ²¹⁴Bi and ²⁰⁸Tl: one gamma per generated event (2447.86 and 2614.5 keV). The network threshold gives 80 % of the bb0nu test events of each training.</p></section>

<section id="flow"><h2>Cut flow, ROI last</h2>
<p>The network fraction is measured on the selected events in the ROI and applied to all selected events: over this energy range the topology does not depend on the energy. For bb0nu and single electrons it is that of the test sample of the training. The ROI fraction is the fraction of the selected events whose energy is in the ROI.</p>
{''.join(flow_html)}</section>

<section id="product"><h2>Track selection × network</h2>
<p>Track selection as a fraction of the events accepted by the fast trigger; network at 80 % signal.</p>
<div class="scroll">{table(['Sample', 'Pitch', 'Track selection', 'Network', 'Product', 'Rejection'], tsn_rows).replace('<table>', '<table class=num>')}</div></section>

<section id="final"><h2>Final numbers</h2>
<p>Events after all cuts (saved, fast trigger, track selection, network at 80 % signal, ROI) and per generated event or gamma.</p>
<div class="scroll">{table(['Sample', '3 mm', '2 mm'], fin_rows).replace('<table>', '<table class=num>')}</div></section>

<section id="floor"><h2>The ²¹⁴Bi floor</h2>
<p>What passes is mostly Compton events whose scattered photon was absorbed within about 1 cm of the main track: two electrons start from almost one point, the topology of a double-beta decay. Compton events are {pct(TH2[r3]['events'][2] / TH2[r3]['events'].sum(), 1)} of the selected ²¹⁴Bi at 3 mm and give {TH2[r3]['passed'][2]} of the {D[r3]['bi214']['k']} passing events; at 2 mm, {pct(TH2[r2]['events'][2] / TH2[r2]['events'].sum(), 1)} and {TH2[r2]['passed'][2]} of {D[r2]['bi214']['k']}. Photoelectric events behave like single electrons: at 3 mm and 2 mm they pass at {pct(TH2[r3]['fraction'][1])} and {pct(TH2[r2]['fraction'][1])}, against {pct(binom(*D[r3]['net']['single_electrons'])[0])} and {pct(binom(*D[r2]['net']['single_electrons'])[0])} for single electrons; those with no secondary deposit pass at {pct(D[r3]['truth']['by_class']['fraction'][0])} and {pct(D[r2]['truth']['by_class']['fraction'][0])}. Details, mechanism and figures: <a href="../../physics/bi214_compton_floor.pdf">note on the ²¹⁴Bi Compton floor</a> (PDF).</p>
<h3>Network acceptance by truth class (80 % signal)</h3>
<div class="scroll">{table(['Class of the selected ²¹⁴Bi event', '3 mm', '2 mm'], cls_rows).replace('<table>', '<table class=num>')}</div>
<h3>By interaction history of the gamma</h3>
<div class="scroll">{table(['History', '3 mm', '2 mm'], hist_rows).replace('<table>', '<table class=num>')}</div>
<h3>Compton class by distance of the second electron to the main track</h3>
<div class="scroll">{table(['Distance', '3 mm', '2 mm'], dist_rows).replace('<table>', '<table class=num>')}</div>
<p>The 2 mm selection joins segments at 5.3 mm and removes most second electrons beyond 8 mm ({TD[r3]['events'][2]} → {TD[r2]['events'][2]} events) and many at 4–8 mm ({TD[r3]['events'][1]} → {TD[r2]['events'][1]}). Below 4 mm the acceptance does not improve, {pct(a3, 1)} at 3 mm and {pct(a2, 1)} at 2 mm: the second electron touches the main track and the lattice records one object. This class is the floor.</p>
<h3>Diffusion is not the driver (3 mm)</h3>
<p>Network acceptance in thirds of the drift length (equal numbers of events; threshold of all test events). σ: thermal ion spread at the median drift of the third.</p>
<div class="scroll">{table(['Third', 'Single electrons: drift (cm)', 'σ (mm)', 'accepted', '²¹⁴Bi: drift (cm)', 'selected events accepted', 'passing per 1000 ROI events'], third_rows).replace('<table>', '<table class=num>')}</div>
<h3>Voxelised Monte Carlo (3 mm)</h3>
<p>Geant4 hits summed in 3 mm cells (fixed lattice, cells above 3 keV), then the same track selection. It misses the deposits merged into the track by the ion readout.</p>
<div class="scroll">{table(['', 'Full readout', 'Voxelised Monte Carlo'], vox_rows).replace('<table>', '<table class=num>')}</div></section>

<section id="links"><h2>Pages</h2><ul>
<li>Statistics, 3 mm: {' · '.join(stat(r3, s) for s in SAMPLES)}</li>
<li>Statistics, 2 mm: {' · '.join(stat(r2, s) for s in SAMPLES)}</li>
<li>Comparisons, selection and cut on E₂: <a href="{up}results/iitaca/II_L1000_D500_comparison/index.html">3 mm</a> · <a href="{up}results/iitaca/II_L1000_D500_P2_comparison/index.html">2 mm</a></li>
<li><a href="{pc}">3 mm against 2 mm pads</a>: cut flows, E₂ and network, truth, rejection factors</li>
<li><a href="{up}results/iitaca/II_L1000_D500_bi214_truth/index.html">²¹⁴Bi passing the network: Monte Carlo truth (3 mm)</a></li>
<li>Trainings: <a href="{up}results/network/II_L1000_D500/trainings/index.html">3 mm</a> · <a href="{up}results/network/II_L1000_D500_P2/trainings/index.html">2 mm</a></li>
<li><a href="../../physics/bi214_compton_floor.pdf">Note: the ²¹⁴Bi background that passes the topological selection</a></li>
<li><a href="../eitaca/results.html">eITACA against iItaca 2 mm</a>: the same analysis on the electron track</li></ul></section>

<section id="inputs"><h2>Inputs</h2><p class="meta">Page: <code>python scripts/iitaca/present_iitaca_results.py</code>. SHA-256 of the files read:</p>
<div class="scroll">{table(['File', 'SHA-256'], [[f'<code>{p}</code>', f'<code>{h[:16]}…</code>'] for p, h in INPUTS])}</div></section>
</main><footer>iItaca · results · II_L1000_D500</footer></body></html>
'''
OUT.write_text(page)
print('wrote', OUT)
