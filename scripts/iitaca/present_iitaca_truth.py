#!/usr/bin/env python3
"""Page of the Monte Carlo truth of the selected 214Bi events against the network score.

  python scripts/iitaca/present_iitaca_truth.py TRUTH_DIR

Reads TRUTH_DIR/truth_bi214.h5, truth_single_electrons.h5
(scripts/iitaca/prepare_iitaca_truth.jl) and summary.h5
(scripts/iitaca/summarize_iitaca_truth.jl) and writes TRUTH_DIR/index.html
with its plots. Nothing is recalculated: the page histograms the saved
columns and prints the saved counts.
"""
import os
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
from present_iitaca_detector import STYLE, table, text

ROOT = Path(__file__).resolve().parents[2]
GROUPS = (('bi214', 1, '²¹⁴Bi passing', '#b23a48'), ('bi214', 0, '²¹⁴Bi rejected', '#e0a458'),
          ('single_electrons', 1, 'single electrons passing', '#185b83'), ('single_electrons', 0, 'single electrons rejected', '#8fb8cf'))
CLASS_NAMES = {'electron_only': 'no secondary deposit', 'fluorescence': 'fluorescence photon of the gamma',
               'second_electron': 'second electron (Compton electron, photoelectron or positron)', 'photon': 'photon of the main track or other photon'}
XRAY_NAMES = {'no_photon': 'no fluorescence photon', 'below_10keV_in_box': 'photon with < 10 keV in the box',
              'lt4mm': 'deposit < 4 mm from the main track', '4to8mm': 'deposit 4–8 mm', 'gt8mm': 'deposit > 8 mm'}
HISTORY_NAMES = {'photoelectric': 'photoelectric absorption only', 'compton_photoelectric': 'Compton scattering(s), then photoelectric absorption',
                 'compton_only': 'Compton scattering only', 'pair': 'pair production', 'other': 'other', 'generated_electron': 'generated electron'}
DEPOSIT_NAMES = ['Compton electron', 'photoelectron (not leading)', 'relaxation electron', 'fluorescence photon', 'other photon of the gamma',
                 'positron', 'photon of the main track', 'other', 'gamma itself (local)']


def load_truth(path):
    with h5py.File(path) as f:
        ev = {k: f['events'][k][:] for k in f['events']}; dep = {k: f['deposits'][k][:] for k in f['deposits']}
        meta = dict(threshold=float(f['threshold_80'][()]), label=text(f['run_label'][()]), dE=float(f['max_box_energy_difference_keV'][()]))
    return ev, dep, meta


def stats(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return f'{len(x):,} ev · mean {x.mean():.4g} · median {np.median(x):.4g}' if len(x) else '0 ev'


def pct(k, n, err=None):
    if n == 0:
        return '—'
    e = err if err is not None else np.sqrt(max(k, 1) * (1 - k / n)) / n
    return f'{100 * k / n:.2f} ± {100 * e:.2f} % ({k:,} of {n:,})'


def hist_groups(ax, data, key, bins, xlabel, title, groups=GROUPS, logy=False, transform=None):
    for sample, p, name, color in groups:
        ev = data[sample][0]; x = ev[key][ev['passed'] == p].astype(float)
        if transform is not None:
            x = transform(ev, p)
        x = x[np.isfinite(x)]
        if len(x):
            ax.hist(np.clip(x, bins[0], bins[-1]), bins=bins, density=True, histtype='step', lw=1.7, color=color, label=f'{name}: {stats(x)}')
    ax.set_xlabel(xlabel); ax.set_ylabel('Events, normalised'); ax.set_title(title, fontsize=9.5); ax.legend(fontsize=7.2)
    if logy:
        ax.set_yscale('log')


def render(directory):
    out = Path(directory).resolve(); plots = out / 'plots'; plots.mkdir(exist_ok=True)
    data = {s: load_truth(out / f'truth_{s}.h5') for s in ('bi214', 'single_electrons')}
    with h5py.File(out / 'summary.h5') as f:
        classes = [text(c) for c in f['classes'][:]]; xclasses = [text(c) for c in f['xray_classes'][:]]; histories = [text(c) for c in f['histories'][:]]
        edges = f['start_edges_keV'][:]; e2_thr = float(f['e2_threshold_keV'][()])
        expected = float(f['expected_bi214_fraction'][()]); expected_err = float(f['expected_bi214_fraction_error'][()]); expected_e2 = float(f['expected_bi214_e2_fraction'][()])
        groups_ = ('by_class', 'by_xray', 'by_history', 'by_second_distance', 'by_second_energy', 'by_start', 'all')
        S = {s: {g: {k: f[s][g][k][:] for k in f[s][g]} for g in groups_} for s in ('bi214', 'single_electrons')}
        dedges = f['second_distance_edges_mm'][:]; eedges = f['second_energy_edges_keV'][:]
        second = {k: f['bi214'][k][:] for k in ('second_distance_mm', 'second_keV', 'passed')}
        cls = {s: f[s]['class'][:] for s in S}
    bi, be = data['bi214'][0], data['single_electrons'][0]; thr = data['bi214'][2]['threshold']
    sel = 'Events in the ROI passing the track selection'

    # 1. acceptance per class
    fig, ax = plt.subplots(figsize=(11, 4.3)); x = np.arange(len(classes)); wbar = .38
    for j, (s, name, color) in enumerate((('bi214', '²¹⁴Bi', '#b23a48'), ('single_electrons', 'single electrons', '#185b83'))):
        a = S[s]['by_class']; ok = a['events'] > 0
        ax.bar(x[ok] + (j - .5) * wbar, 100 * a['fraction'][ok], wbar, yerr=100 * a['error'][ok], color=color, alpha=.8, label=name, capsize=3)
        for i in np.flatnonzero(ok):
            ax.text(x[i] + (j - .5) * wbar, 100 * a['fraction'][i] + 100 * a['error'][i] + .15, f"{a['passed'][i]}/{a['events'][i]:,}", ha='center', fontsize=7.5)
    for s, ls in (('bi214', '--'), ('single_electrons', ':')):
        ax.axhline(100 * S[s]['all']['fraction'][0], color='#203343', ls=ls, lw=1, label=f"all {'²¹⁴Bi' if s == 'bi214' else 'single electrons'}: {100 * S[s]['all']['fraction'][0]:.2f} %")
    ax.set_xticks(x, [CLASS_NAMES[c].replace(' (', '\n(') for c in classes], fontsize=8); ax.set_ylabel('Network acceptance at 80 % signal (%)'); ax.legend(fontsize=8)
    ax.set_title(f'{sel}: network acceptance per class of secondary deposit (≥ 10 keV in the box)', fontsize=9.5)
    fig.tight_layout(); fig.savefig(plots / 'acceptance_class.png', dpi=130); plt.close(fig)

    # 1b. second electron
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4)); m = np.isfinite(second['second_distance_mm'])
    bins = np.linspace(0, 16, 33)
    for p, name, color in ((1, '²¹⁴Bi passing', '#b23a48'), (0, '²¹⁴Bi rejected', '#e0a458')):
        x = second['second_distance_mm'][m & (second['passed'] == p)]
        ax[0].hist(np.clip(x, 0, 16), bins=bins, histtype='step', lw=1.7, color=color, label=f'{name}: {stats(x)}')
        y = second['second_keV'][m & (second['passed'] == p)]
        ax[1].hist(np.clip(y, 0, 2000), bins=np.linspace(0, 2000, 41), histtype='step', lw=1.7, color=color, label=f'{name}: {stats(y)}')
    for v in (4, 8):
        ax[0].axvline(v, color='#6b7480', ls='--', lw=1)
    ax[0].set_xlabel('Smallest distance of the second electron to the main track (mm)'); ax[1].set_xlabel('Energy of the second electron in the box (keV)')
    for a_ in ax:
        a_.set_ylabel('Events'); a_.legend(fontsize=7.5); a_.set_title('²¹⁴Bi selected events with a second electron ≥ 10 keV (the most energetic one)', fontsize=9.5)
    fig.tight_layout(); fig.savefig(plots / 'second_electron.png', dpi=130); plt.close(fig)

    # 2. fluorescence photon
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4)); g2 = GROUPS[:2]
    bins = np.linspace(0, 40, 41)
    hist_groups(ax[0], data, 'xray_distance_mm', bins, 'Smallest distance of the fluorescence deposit to the main track (mm)',
                '²¹⁴Bi with a fluorescence photon depositing ≥ 10 keV in the box', groups=g2,
                transform=lambda ev, p: ev['xray_distance_mm'][(ev['passed'] == p) & (ev['xray_deposit_keV'] >= 10)])
    for v in (4, 8):
        ax[0].axvline(v, color='#6b7480', ls='--', lw=1)
    hist_groups(ax[1], data, 'xray_start_mm', bins, 'Distance of the fluorescence deposit centre to the start of the leading electron (mm)',
                '²¹⁴Bi with a fluorescence photon depositing ≥ 10 keV in the box', groups=g2,
                transform=lambda ev, p: ev['xray_start_mm'][(ev['passed'] == p) & (ev['xray_deposit_keV'] >= 10)])
    fig.tight_layout(); fig.savefig(plots / 'fluorescence.png', dpi=130); plt.close(fig)

    # 3. energy at the start of the leading electron
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    hist_groups(ax[0], data, 'start_all_keV', np.linspace(0, 800, 41), 'Energy within 7 mm of the start of the leading electron, every hit in the box (keV)', sel)
    mid = np.arange(len(edges) - 1)
    for s, name, color in (('bi214', '²¹⁴Bi', '#b23a48'), ('single_electrons', 'single electrons', '#185b83')):
        a = S[s]['by_start']; ok = a['events'] > 0
        ax[1].errorbar(mid[ok] + (.08 if s == 'bi214' else -.08), 100 * a['fraction'][ok], yerr=100 * a['error'][ok], fmt='o', color=color, capsize=3, label=name)
    ax[1].set_xticks(mid, [f'{edges[i]:.0f}–{edges[i + 1]:.0f}' for i in mid], rotation=35, fontsize=8)
    ax[1].set_xlabel('Energy within 7 mm of the start (keV)'); ax[1].set_ylabel('Network acceptance at 80 % signal (%)'); ax[1].set_yscale('log'); ax[1].legend(fontsize=8)
    ax[1].set_title('Acceptance per bin of the start energy (same selection)', fontsize=9.5); ax[1].grid(alpha=.3, which='both')
    fig.tight_layout(); fig.savefig(plots / 'start_energy.png', dpi=130); plt.close(fig)

    # 4. secondary energy near the start, from the truth
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    hist_groups(ax[0], data, 'start_all_keV', np.linspace(0, 120, 49), 'Energy of hits outside the main track within 7 mm of the start (keV)', sel,
                transform=lambda ev, p: (ev['start_all_keV'] - ev['start_main_keV'])[ev['passed'] == p], logy=True)
    hist_groups(ax[1], data, 'start_main_keV', np.linspace(0, 800, 41), 'Energy of main-track hits within 7 mm of the start (keV)', sel)
    fig.tight_layout(); fig.savefig(plots / 'start_parts.png', dpi=130); plt.close(fig)

    # 5. main track
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.4))
    hist_groups(ax[0, 0], data, 'leading_keV', np.linspace(1000, 2460, 74), 'Initial kinetic energy of the leading electron (keV)', sel, logy=True)
    hist_groups(ax[0, 1], data, 'E2_keV', np.linspace(0, 700, 36), 'E₂, the smaller endpoint energy of the itrack (keV)', sel)
    ax[0, 1].axvline(e2_thr, color='#203343', ls='--', lw=1)
    hist_groups(ax[1, 0], data, 'length_mm', np.linspace(0, 400, 41), 'Length of the itrack (mm)', sel)
    hist_groups(ax[1, 1], data, 'ion_drift_cm', np.linspace(0, 100, 21), 'Ion drift length (cm)', sel)
    fig.tight_layout(); fig.savefig(plots / 'main_track.png', dpi=130); plt.close(fig)

    # 6. secondary energy by distance
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for a_, key, name in zip(ax, ('secondary_lt4_keV', 'secondary_4to8_keV', 'secondary_gt8_keV'), ('below 4 mm', '4 to 8 mm', 'beyond 8 mm')):
        hist_groups(a_, data, key, np.linspace(0, 400, 41), f'Energy of secondary hits {name} from the main track (keV)', 'Selected events, in the box', logy=True)
    fig.tight_layout(); fig.savefig(plots / 'secondary_distance.png', dpi=130); plt.close(fig)

    # tables
    def group_rows(rows):
        return [[name] + [fn(data[s][0], data[s][0]['passed'] == p) for fn in rows] for s, p, name, _ in GROUPS]
    n_of = lambda ev, m: f'{m.sum():,}'
    frac = lambda cond: (lambda ev, m: f'{100 * np.mean(cond(ev)[m]):.1f} %')
    mean = lambda key: (lambda ev, m: f'{np.nanmean(ev[key][m]):.1f}')
    ab_passed = int(S['bi214']['all']['passed'][0])
    t_hist = table(['Interaction history of the gamma', '²¹⁴Bi passing', '²¹⁴Bi rejected', 'network acceptance'],
                   [[HISTORY_NAMES[h], f"{S['bi214']['by_history']['passed'][i]:,}", f"{S['bi214']['by_history']['events'][i] - S['bi214']['by_history']['passed'][i]:,}",
                     pct(S['bi214']['by_history']['passed'][i], S['bi214']['by_history']['events'][i], S['bi214']['by_history']['error'][i])]
                    for i, h in enumerate(histories) if S['bi214']['by_history']['events'][i] > 0])
    lead_compton = lambda ev: ev['leading_process'] == 2
    t_main = table(['Group', 'events', 'leading electron from Compton', 'leading electron (keV), mean', 'main-track energy (keV), mean', 'secondary energy in the box (keV), mean',
                    'E₂ (keV), mean', 'itrack length (mm), mean', 'start energy, 7 mm (keV), mean'],
                   group_rows([n_of, frac(lead_compton), mean('leading_keV'), mean('main_keV_truth'), mean('secondary_keV'), mean('E2_keV'), mean('length_mm'), mean('start_all_keV')]))
    t_sec = table(['Group', '≥ 1 secondary deposit ≥ 10 keV < 4 mm', '4–8 mm', '> 8 mm', 'secondary energy < 4 mm (keV), mean', '4–8 mm', '> 8 mm'],
                  group_rows([frac(lambda ev: ev['deposits_ge10_lt4'] > 0), frac(lambda ev: ev['deposits_ge10_4to8'] > 0), frac(lambda ev: ev['deposits_ge10_gt8'] > 0),
                              mean('secondary_lt4_keV'), mean('secondary_4to8_keV'), mean('secondary_gt8_keV')]))
    t_class = table(['Class of the event', '²¹⁴Bi: share of events', '²¹⁴Bi: network', '²¹⁴Bi: E₂ cut', 'single electrons: share', 'single electrons: network', 'single electrons: E₂ cut'],
                    [[CLASS_NAMES[c]] + sum(([f"{100 * S[s]['by_class']['events'][i] / S[s]['by_class']['events'].sum():.1f} %",
                                               pct(S[s]['by_class']['passed'][i], S[s]['by_class']['events'][i], S[s]['by_class']['error'][i]),
                                               pct(S[s]['by_class']['e2_passed'][i], S[s]['by_class']['events'][i], S[s]['by_class']['e2_error'][i])]
                                              for s in ('bi214', 'single_electrons')), []) for i, c in enumerate(classes)])
    xb = S['bi214']['by_xray']
    t_xray = table(['Fluorescence photon of the gamma (photoelectric absorption only)', '²¹⁴Bi: share of events', 'share of passing events', 'network acceptance', 'E₂ cut'],
                   [[XRAY_NAMES[c], f"{100 * xb['events'][i] / xb['events'].sum():.1f} %", f"{100 * xb['passed'][i] / xb['passed'].sum():.1f} %",
                     pct(xb['passed'][i], xb['events'][i], xb['error'][i]), pct(xb['e2_passed'][i], xb['events'][i], xb['e2_error'][i])] for i, c in enumerate(xclasses)])
    sd_, se_ = S['bi214']['by_second_distance'], S['bi214']['by_second_energy']
    t_second = table(['Most energetic second electron', 'events', 'share of passing ²¹⁴Bi', 'network acceptance', 'E₂ cut'],
                     [[f'distance {dedges[i]:.0f}–{dedges[i + 1]:.0f} mm' if dedges[i + 1] < 1e6 else f'distance > {dedges[i]:.0f} mm', f"{sd_['events'][i]:,}",
                       f"{100 * sd_['passed'][i] / ab_passed:.1f} %", pct(sd_['passed'][i], sd_['events'][i], sd_['error'][i]), pct(sd_['e2_passed'][i], sd_['events'][i], sd_['e2_error'][i])]
                      for i in range(len(dedges) - 1)] +
                     [[f'energy {eedges[i]:.0f}–{eedges[i + 1]:.0f} keV', f"{se_['events'][i]:,}", f"{100 * se_['passed'][i] / ab_passed:.1f} %",
                       pct(se_['passed'][i], se_['events'][i], se_['error'][i]), pct(se_['e2_passed'][i], se_['events'][i], se_['e2_error'][i])] for i in range(len(eedges) - 1)])
    sb, se = S['bi214']['by_start'], S['single_electrons']['by_start']
    t_start = table(['Start energy (keV)', '²¹⁴Bi: share', '²¹⁴Bi: network', 'single electrons: share', 'single electrons: network'],
                    [[f'{edges[i]:.0f}–{edges[i + 1]:.0f}', f"{100 * sb['events'][i] / sb['events'].sum():.1f} %", pct(sb['passed'][i], sb['events'][i], sb['error'][i]),
                      f"{100 * se['events'][i] / se['events'].sum():.1f} %", pct(se['passed'][i], se['events'][i], se['error'][i])] for i in range(len(edges) - 1)])
    ab, ae = S['bi214']['all'], S['single_electrons']['all']
    t_all = table(['Sample', 'network at 80 % signal', f'E₂ > {e2_thr:.0f} keV (80 % signal)'],
                  [['²¹⁴Bi', pct(ab['passed'][0], ab['events'][0]), pct(ab['e2_passed'][0], ab['events'][0])],
                   ['single electrons (test events of the training)', pct(ae['passed'][0], ae['events'][0]), pct(ae['e2_passed'][0], ae['events'][0])],
                   ['²¹⁴Bi expected from its start energy, with the single-electron acceptance per bin', f'{100 * expected:.2f} ± {100 * expected_err:.2f} %', f'{100 * expected_e2:.2f} %']])
    # deposit types of the secondary deposits >= 10 keV, Bi214
    db = data['bi214'][1]; m10 = db['deposit_in_box_keV'] >= 10
    passed_ev = set(bi['event_number'][bi['passed'] == 1].tolist()); pd_ = np.array([e in passed_ev for e in db['event_number']])
    t_dep = table(['Type of deposit ≥ 10 keV in the box (²¹⁴Bi)', 'in passing events', 'per passing event', 'in rejected events', 'per rejected event', 'median smallest distance, all (mm)'],
                  [[DEPOSIT_NAMES[t - 1], f'{np.sum(m10 & pd_ & (db["type"] == t)):,}', f'{np.sum(m10 & pd_ & (db["type"] == t)) / max(len(passed_ev), 1):.3f}',
                    f'{np.sum(m10 & ~pd_ & (db["type"] == t)):,}', f'{np.sum(m10 & ~pd_ & (db["type"] == t)) / max(len(bi["passed"]) - len(passed_ev), 1):.3f}',
                    f'{np.median(db["min_distance_mm"][m10 & (db["type"] == t)]):.1f}'] for t in range(1, 10) if np.any(m10 & (db['type'] == t))])

    depth = len(os.path.relpath(out, ROOT).split(os.sep)); up = '../' * depth
    page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>iItaca 214Bi truth</title><style>{STYLE}img{{max-width:100%}}</style></head><body><header><a href="{up}notes/html/index.html">Documentation index</a> · <a href="{up}notes/html/iitaca/index.html">iItaca</a><h1>²¹⁴Bi events passing the network · Monte Carlo truth · II_L1000_D500</h1><p>Why selected ²¹⁴Bi events pass the bb0nu/single-electron network more often than single electrons.</p></header><main>
@@CONCLUSION@@
<section><h2>Samples</h2><p>²¹⁴Bi: the 7,640 events in the ROI passing the track selection, all scored by the network <code>itrack_sparse_e100</code>. Single electrons: the 12,523 test events of its training (the same selection). Passing: score above {thr:.3f}, the threshold for 80 % of the bb0nu test events. E₂ cut: E₂ above {e2_thr:.1f} keV, the threshold for 80 % of the signal.</p>{t_all}</section>
<section><h2>Method</h2><ol>
<li>Read the Geant4 hits (reduced hit table, mm, keV) and particles (NEXUS <code>MC/particles</code>) of every event.</li>
<li>Main track as in the production: the leading electron (for ²¹⁴Bi the most energetic electron produced directly by the gamma; for single electrons the generated one) and its descendants not passing through a photon.</li>
<li>Secondary deposit: a particle produced directly by the gamma, or a photon emitted by the main track, with its descendants. Its distance to the main track is the smallest distance of its hits to a main-track hit (Geant4 steps of at most 1 mm).</li>
<li>Energies are summed inside the simulated box of the event (main track and gamma deposits inside the 160 mm frame, with a 20 mm margin). The box energy recomputed from the hits equals the production value for every event (largest difference {data['bi214'][2]['dE']:.2g} keV).</li>
<li>Start energy: energy of the hits within 7 mm (the blob radius) of the creation point of the leading electron.</li>
<li>Classes: a deposit counts when it has ≥ 10 keV in the box (the deposit threshold of the track selection). The relaxation electrons and the local energy of the gamma at the absorption point are always at the start of the photoelectron and are not counted as deposits.</li></ol>
<p>Code: <code>scripts/iitaca/prepare_iitaca_truth.jl</code> (truth, Hyperion), <code>scripts/iitaca/summarize_iitaca_truth.jl</code> (classes and acceptances), <code>scripts/iitaca/present_iitaca_truth.py</code> (this page).</p></section>
<section><h2>Interaction history</h2>{t_hist}<p>Main track of the four groups:</p>{t_main}</section>
<section><h2>Secondary deposits</h2><p>Fraction of events with at least one secondary deposit of ≥ 10 keV in the box by smallest distance to the main track, and mean energy of the secondary hits by their distance.</p>{t_sec}{t_dep}<img src="plots/secondary_distance.png" alt="Energy of secondary hits by distance to the main track"></section>
<section><h2>Acceptance per class</h2>{t_class}<img src="plots/acceptance_class.png" alt="Network acceptance per class of secondary deposit"></section>
<section><h2>Second electron</h2><p>²¹⁴Bi events of the class “second electron”: all of them have a Compton scattering followed by the photoelectric absorption.</p>{t_second}<img src="plots/second_electron.png" alt="Distance and energy of the second electron"></section>
<section><h2>Fluorescence photon</h2>{t_xray}<img src="plots/fluorescence.png" alt="Distance of the fluorescence deposit to the main track and to its start"></section>
<section><h2>Energy at the start of the leading electron</h2><p>The network and E₂ respond to the energy near the two ends of the track. The start energy of the truth separates it into main-track hits and the rest.</p>{t_start}<img src="plots/start_energy.png" alt="Start energy and acceptance per bin"><img src="plots/start_parts.png" alt="Start energy of secondary and main-track hits"></section>
<section><h2>Main track</h2><img src="plots/main_track.png" alt="Leading electron energy, E2, itrack length and ion drift"></section>
<p>Numerical results: <a href="truth_bi214.h5">truth_bi214.h5</a> · <a href="truth_single_electrons.h5">truth_single_electrons.h5</a> · <a href="summary.h5">summary.h5</a>.</p>
</main><footer>iItaca · ²¹⁴Bi truth · II_L1000_D500</footer></body></html>
'''
    # conclusion, from the saved counts
    hb, cb, xb_ = S['bi214']['by_history'], S['bi214']['by_class'], S['bi214']['by_xray']
    pe = ae['passed'][0] / ae['events'][0]; nb = int(ab['events'][0]); kb = int(ab['passed'][0])
    i_ph, i_cp = histories.index('photoelectric'), histories.index('compton_photoelectric'); i_se = classes.index('second_electron')
    excess = kb - nb * pe; excess_c = cb['passed'][i_se] - cb['events'][i_se] * pe
    i_no, i_4 = xclasses.index('no_photon'), xclasses.index('lt4mm')
    lc = (bi['passed'] == 1) & (bi['leading_process'] == 2); lc_pass = int(lc.sum()); lc_keV = float(np.median(bi['leading_keV'][lc]))
    conclusion = f'''<section><h2>Conclusion</h2><ol>
<li>Of the {kb} passing ²¹⁴Bi events, {cb['passed'][i_se]} ({100 * cb['passed'][i_se] / kb:.0f} %) have a second electron of ≥ 10 keV in the box: the gamma scattered (Compton) and was then absorbed, and the two electrons lie within about 15 mm of each other, so the track selection keeps them as one itrack. These events are {100 * cb['events'][i_se] / nb:.1f} % of the selected ²¹⁴Bi and pass at {pct(cb['passed'][i_se], cb['events'][i_se], cb['error'][i_se])}; the E₂ cut accepts {100 * cb['e2_fraction'][i_se]:.0f} % of them, close to the 80 % of the signal: two electrons from nearby points look like bb0nu. In {lc_pass} of the {cb['passed'][i_se]} passing ones the leading electron is the Compton electron (median {lc_keV:.0f} keV): the gamma scattered at a large angle and the scattered photon, of a few hundred keV, was absorbed within about 1 cm.</li>
<li>They make {100 * excess_c / excess:.0f} % of the excess over the single-electron acceptance ({excess_c:.0f} of {excess:.0f} events). The acceptance falls with the distance of the second electron: {pct(sd_['passed'][0], sd_['events'][0], sd_['error'][0])} below 4 mm, {pct(sd_['passed'][1], sd_['events'][1], sd_['error'][1])} at 4–8 mm, {pct(sd_['passed'][2], sd_['events'][2], sd_['error'][2])} beyond.</li>
<li>Events absorbed without Compton scattering pass at {pct(hb['passed'][i_ph], hb['events'][i_ph], hb['error'][i_ph])}, against {100 * pe:.2f} % for single electrons. Without a fluorescence photon they pass as single electrons ({pct(xb_['passed'][i_no], xb_['events'][i_no], xb_['error'][i_no])}); with the photon absorbed within 4 mm of the main track, sharing its cells, {pct(xb_['passed'][i_4], xb_['events'][i_4], xb_['error'][i_4])}. The merged fluorescence photon is a small effect.</li>
<li>Weighting the single-electron acceptance per bin of the energy within 7 mm of the start of the leading electron with the ²¹⁴Bi distribution of that energy gives {100 * expected:.2f} ± {100 * expected_err:.2f} %, against the {100 * kb / nb:.2f} % measured: the network accepts events with energy near the start of the track, from a second electron, a fluorescence photon or the track folding back, and responds to it in about the same way in both samples. Passing single electrons and passing photoelectric events have a compact main track whose start region holds a large energy.</li>
<li>The second electrons below 4 mm ({sd_['passed'][0]} of the passing events) touch the main track within one or two cells and cannot be separated on the 3 mm lattice. Those at 4–8 mm and beyond ({sd_['passed'][1] + sd_['passed'][2]} events) start at a point separate from the main track; training against ²¹⁴Bi, or a feature giving the number of track ends or the separation of the two starts, could address them.</li></ol></section>
'''
    (out / 'index.html').write_text(page.replace('@@CONCLUSION@@', conclusion)); print(out / 'index.html')


if __name__ == '__main__':
    render(sys.argv[1])
