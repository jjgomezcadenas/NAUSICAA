#!/usr/bin/env python3
"""Page of the same Monte Carlo events read on 3 mm and on 2 mm pads.

  python scripts/iitaca/present_iitaca_pitch_comparison.py OUT_DIR

Reads OUT_DIR/comparison.h5 (scripts/iitaca/prepare_iitaca_pitch_comparison.jl)
and writes OUT_DIR/index.html with its plots. Nothing is recalculated.
"""
import os
import sys
import tomllib
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
from present_iitaca_detector import STYLE, table, text
from present_iitaca_statistics import SAMPLE, FLOW

ROOT = Path(__file__).resolve().parents[2]
COLORS = ('#2a78d6', '#eb6834', '#1baf7a', '#eda100')       # categorical slots 1-4, fixed order of the settings
NAME = {'II_L1000_D500': '3 mm', 'II_L1000_D500_P2': '2 mm, radii scaled (A)', 'II_L1000_D500_P2_R3': '2 mm, 3 mm radii (B)'}
HISTORY = {'photoelectric': 'photoelectric absorption only', 'compton': 'Compton scattering(s), then absorption or escape', 'other': 'other'}
DISTANCE = {'lt4mm': 'second electron below 4 mm', '4to8mm': 'second electron 4 to 8 mm', 'gt8mm': 'second electron beyond 8 mm',
            'no_second_electron': 'no second electron of 10 keV or more'}
BACKGROUNDS = ('single_electrons', 'bi214', 'tl208')
STEPS = ('generated', 'read', 'trigger', 'roi', 'track')


def pct(x, e=None, digits=2):
    if x is None or not np.isfinite(x):
        return '—'
    return f'{100 * x:.{digits}f} ± {100 * e:.{digits}f} %' if e is not None else f'{100 * x:.{digits}f} %'


def factor(before, after):
    return f'{before / after:,.3g}' if after > 0 else '—'


def load(path):
    with h5py.File(path) as f:
        S = [text(s) for s in f['settings'][:]]; targets = f['targets'][:]
        hist = [text(s) for s in f['history_groups'][:]]; dist = [text(s) for s in f['distance_groups'][:]]
        common = {k: int(f['bi214_common'][k][()]) for k in f['bi214_common']} if 'bi214_common' in f else {}
        R = {}
        for s in S:
            g = f[s]; r = dict(tracks=tomllib.loads(text(g['track_settings'][()]))['tracks'])
            r['flow'] = {k: dict(zip([text(x) for x in g['flow'][k]['step'][:]], [int(n) for n in g['flow'][k]['events'][:]])) for k in g['flow']}
            r['e2_threshold'] = g['e2/threshold_keV'][:]; r['e2_efficiency'] = g['e2/signal_efficiency'][:]
            r['e2'] = {b: dict(acceptance=g['e2'][b]['acceptance'][:], error=g['e2'][b]['error'][:], accepted=g['e2'][b]['accepted'][:]) for b in BACKGROUNDS}
            n = g['network']; r['net_threshold'] = n['threshold'][:]; r['net_efficiency'] = n['signal_efficiency'][:]
            r['area'] = float(n['area'][()]); r['signal_test'] = int(n['signal_test_events'][()])
            r['net'] = {b: dict(acceptance=n[b]['acceptance'][:], error=n[b]['error'][:], accepted=n[b]['accepted'][:], events=int(n[b]['events'][()])) for b in BACKGROUNDS}
            for name in ('history', 'distance'):
                h = g[f'bi214_{name}']
                r[name] = {k: h[k][()] for k in h}        # 2-D arrays: [target, group]
            R[s] = r
    return S, targets, hist, dist, common, R


def plots(out, S, targets, dist, R):
    (out / 'plots').mkdir(exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2), sharey=False)
    for k, b in enumerate(BACKGROUNDS):
        for i, s in enumerate(S):
            a = R[s]['net'][b]
            ax[k].errorbar(100 * targets, 100 * a['acceptance'], yerr=100 * a['error'], marker='o', ms=8, lw=2, capsize=3, color=COLORS[i], label=NAME.get(s, s))
            j = list(np.isclose(targets, 0.8)).index(True)
            ax[k].annotate(f"{100 * a['acceptance'][j]:.2f} %", (80, 100 * a['acceptance'][j]), textcoords='offset points', xytext=(8, 4 - 12 * i), fontsize=8, color='#203343')
        ax[k].set_title(f"{SAMPLE[b]}: network acceptance\n(ROI and track selection; {'test sample' if b == 'single_electrons' else 'all selected events'})", fontsize=10)
        ax[k].set_xlabel('bb0nu signal efficiency (%)'); ax[k].set_ylabel('accepted (%)'); ax[k].grid(alpha=0.25); ax[k].invert_xaxis()
    ax[0].legend(fontsize=9)
    fig.tight_layout(); fig.savefig(out / 'plots/network.png', dpi=130); plt.close(fig)

    j = list(np.isclose(targets, 0.8)).index(True); x = np.arange(len(dist)); w = 0.8 / len(S)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    for k, (key, label) in enumerate((('net', 'network'), ('e2', 'E₂ cut'))):
        for i, s in enumerate(S):
            d = R[s]['distance']; frac = d[f'{key}_fraction'][j]; err = d[f'{key}_error'][j]
            bars = ax[k].bar(x + (i - (len(S) - 1) / 2) * w, 100 * frac, w * 0.92, yerr=100 * err, capsize=2, color=COLORS[i], label=NAME.get(s, s))
            for xb, f_, e_, n in zip(x + (i - (len(S) - 1) / 2) * w, frac, err, d['events']):
                ax[k].text(xb, 100 * (f_ + e_) + 1, f'{100 * f_:.1f}\nn={int(n)}', ha='center', va='bottom', fontsize=7, color='#203343')
        ax[k].set_xticks(x, ['no second electron' if g == 'no_second_electron' else DISTANCE[g].replace('second electron ', '') for g in dist])
        ax[k].set_title(f'²¹⁴Bi selected events: {label} acceptance at 80 % signal,\nby distance of the second electron to the main track (% and events)', fontsize=10)
        ax[k].set_ylim(0, ax[k].get_ylim()[1] * 1.12)
        ax[k].set_ylabel('accepted (%)'); ax[k].grid(axis='y', alpha=0.25)
    ax[0].legend(fontsize=9)
    fig.tight_layout(); fig.savefig(out / 'plots/bi214_distance.png', dpi=130); plt.close(fig)


def render(out):
    out = Path(out).resolve(); S, targets, hist, dist, common, R = load(out / 'comparison.h5')
    names = [NAME.get(s, s) for s in S]; plots(out, S, targets, dist, R)
    j80 = list(np.isclose(targets, 0.8)).index(True)

    settings = table(['Setting', 'Run prefix', 'Hit cut', 'Components', 'Segment', 'Segments joined', 'Graph cells', 'Detached deposit', 'Endpoint sphere'],
        [[NAME.get(s, s), f'<code>{s}</code>', f"{R[s]['tracks']['hit_cut_keV']:g} keV", f"{R[s]['tracks']['segment_connection_mm']:g} mm", f"{R[s]['tracks']['segment_min_hits']} hits",
          f"{R[s]['tracks']['connection_mm']:g} mm", f"{R[s]['tracks']['voxel_width_mm']:g} mm", f"{R[s]['tracks']['deposit_min_keV']:g} keV", f"{R[s]['tracks']['blob_radius_mm']:g} mm"] for s in S])

    rows = []
    for k in ('bb0nu',) + BACKGROUNDS:
        for st in STEPS:
            rows.append([SAMPLE[k] if st == 'generated' else '', FLOW[st]] + [f"{R[s]['flow'][k][st]:,}" for s in S])
    flow = table(['Sample', 'Step'] + names, rows)

    rows = []
    for i, t in enumerate(targets):
        rows.append([f'<b>{100 * t:.0f} %</b>', 'threshold'] + [f"E₂ > {R[s]['e2_threshold'][i]:.0f} keV · score > {R[s]['net_threshold'][i]:.2f}" for s in S])
        for b in BACKGROUNDS:
            rows.append(['', SAMPLE[b]] + [f"E₂ {pct(R[s]['e2'][b]['acceptance'][i], R[s]['e2'][b]['error'][i])} · network {pct(R[s]['net'][b]['acceptance'][i], R[s]['net'][b]['error'][i])}" for s in S])
    accept = table(['Signal efficiency', 'Sample'] + names, rows)

    def grouped(name, labels, keys):
        rows = []
        for g, key in enumerate(keys):
            if not any(R[s][name]['events'][g] for s in S):
                continue                                   # a group without events in any setting
            for i, t in enumerate(targets):
                rows.append([labels[key] if i == 0 else '', f'{100 * t:.0f} %'] +
                            [f"{int(R[s][name]['events'][g])} events · E₂ {pct(R[s][name]['e2_fraction'][i][g], R[s][name]['e2_error'][i][g], 1)} · network {pct(R[s][name]['net_fraction'][i][g], R[s][name]['net_error'][i][g], 1)}" for s in S])
        return table(['²¹⁴Bi group', 'Signal efficiency'] + names, rows)
    by_history = grouped('history', HISTORY, hist); by_distance = grouped('distance', DISTANCE, dist)

    def rejection(s):
        r = R[s]; fl = r['flow']; rows = []
        prev = {k: fl[k]['generated'] for k in fl}
        for st in STEPS[1:]:
            rows.append([FLOW[st], f"{100 * fl['bb0nu'][st] / prev['bb0nu']:.1f} % · {100 * fl['bb0nu'][st] / fl['bb0nu']['generated']:.1f} %"] +
                        [f"{factor(prev[b], fl[b][st])} · {factor(fl[b]['generated'], fl[b][st])}" for b in BACKGROUNDS])
            prev = {k: fl[k][st] for k in fl}
        for label, key in (('E₂ cut at 80 % signal', 'e2'), ('network at 80 % signal', 'net')):
            eff = r['e2_efficiency'][j80] if key == 'e2' else r['net_efficiency'][j80]
            cells = []
            for b in BACKGROUNDS:
                a = r[key][b]['acceptance'][j80]
                cells.append(f"{factor(1, a)} · {factor(fl[b]['generated'], fl[b]['track'] * a)}" if a > 0 else '— (no event accepted)')
            rows.append([label, f"{100 * eff:.1f} % · {100 * eff * fl['bb0nu']['track'] / fl['bb0nu']['generated']:.1f} %"] + cells)
        return f'<h3>{NAME.get(s, s)}</h3>' + table(['Step', 'bb0nu kept: step · cumulative'] + [f'{SAMPLE[b]} rejection: step · cumulative' for b in BACKGROUNDS], rows)
    rejections = ''.join(rejection(s) for s in S)

    head = ['3 mm', '2 mm A', '2 mm B']
    summary = table(['At 80 % signal'] + names,
        [['bb0nu after ROI and track selection'] + [f"{R[s]['flow']['bb0nu']['track']:,} ({100 * R[s]['flow']['bb0nu']['track'] / R[s]['flow']['bb0nu']['roi']:.1f} % of the ROI)" for s in S]] +
        [[f'{SAMPLE[b]}, network'] + [pct(R[s]['net'][b]['acceptance'][j80], R[s]['net'][b]['error'][j80]) + f" ({int(R[s]['net'][b]['accepted'][j80])} of {R[s]['net'][b]['events']:,})" for s in S] for b in BACKGROUNDS] +
        [[f'{SAMPLE[b]}, E₂ cut'] + [pct(R[s]['e2'][b]['acceptance'][j80], R[s]['e2'][b]['error'][j80]) for s in S] for b in BACKGROUNDS] +
        [['²¹⁴Bi over single electrons, network'] + [f"{R[s]['net']['bi214']['acceptance'][j80] / R[s]['net']['single_electrons']['acceptance'][j80]:.2f}" for s in S]] +
        [['Network: area under the curve (test sample)'] + [f"{R[s]['area']:.4f}" for s in S]])
    common_text = '; '.join(f"{NAME.get(a, a)} and {NAME.get(b, b)}: {n:,}" for k, n in common.items() for a, b in [k.split('__')])

    depth = len(os.path.relpath(out, ROOT).split(os.sep)); up = '../' * depth
    page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>iItaca pad pitch</title><style>{STYLE}img{{max-width:100%}}td{{font-size:14px}}</style></head><body><header><a href="{up}notes/html/index.html">Documentation index</a> · <a href="{up}notes/html/iitaca/index.html">iItaca</a><h1>II_L1000_D500 on 3 mm and on 2 mm pads</h1><p>The same Geant4 events, seeds and ion drift, read on 3 mm pads with 15 ms frames and on 2 mm pads with 10 ms frames; for 2 mm, two track selections from the same hit files.</p></header><main>
<section><h2>1. Settings</h2><p>Detector of II_L1000_D500 (1.0 m drift, 500 V/cm, 20.0 cm/s). The 2 mm readout uses the VH1 card with 10 ms frames; hits are saved above 3 keV in both pitches. Setting A scales the radii of the selection to the 2 mm lattice; setting B keeps the radii of the 3 mm selection. Both 2 mm settings use a hit cut of 3 keV: at 5 keV the smaller cells break the track.</p>{settings}</section>
<section><h2>2. Summary at 80 % signal</h2><p>Acceptances are fractions of the events in the ROI that pass the track selection of each setting; single electrons for the network are those of its test sample; errors are binomial.</p>{summary}</section>
<section><h2>3. Cut flow</h2>{flow}</section>
<section><h2>4. Cut on E₂ and network</h2><p>The E₂ threshold and the network threshold of each setting give the stated bb0nu efficiency on its own selected events (the network on its test sample). Each network is a sparse residual network on the itrack, trained for 100 epochs on bb0nu against single electrons of that setting (validation batches 8–13, test batches 0–7).</p>{accept}<img src="plots/network.png" alt="Network acceptance of single electrons, 214Bi and 208Tl against signal efficiency"></section>
<section><h2>5. ²¹⁴Bi by Monte Carlo truth</h2><p>Truth of the selected ²¹⁴Bi events of each setting (<code>prepare_iitaca_truth.jl</code>, <code>summarize_iitaca_truth.jl</code>). The second electron is the most energetic Compton electron, photoelectron other than the leading electron, or positron with at least 10 keV in the simulated box; its distance is the smallest distance between its hits and the main-track hits. Selected events common to two settings: {common_text}.</p>
<h3>Interaction history of the gamma</h3>{by_history}
<h3>Distance of the second electron</h3>{by_distance}<img src="plots/bi214_distance.png" alt="214Bi acceptance by distance of the second electron"></section>
<section><h2>6. Rejection factors</h2><p>For each step: bb0nu kept at the step and since generation; for the backgrounds, events before the step over events after it, and generated events over events after it. The network and E₂ rows follow the track selection; for single electrons the network acceptance is that of the test sample.</p>{rejections}</section>
<p>Numerical results: <a href="comparison.h5">comparison.h5</a> (<code>scripts/iitaca/prepare_iitaca_pitch_comparison.jl</code>); page: <code>scripts/iitaca/present_iitaca_pitch_comparison.py</code>.</p>
</main><footer>iItaca · pad pitch · II_L1000_D500</footer></body></html>
'''
    (out / 'index.html').write_text(page); print(out / 'index.html')


if __name__ == '__main__':
    render(sys.argv[1])
