#!/usr/bin/env python3
"""Page of the diffusion-limited readout bound (prepare_diffusion_limit.jl).

  python scripts/iitaca/present_diffusion_limit.py RESULT_DIR

Reads RESULT_DIR/diffusion_limit.h5 and writes RESULT_DIR/index.html with its
plots. Nothing is recalculated: the page prints the saved fractions.
"""
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
SAMPLES = (('bb0nu', 'bb0nu (signal loss)', '#2a9d8f'), ('single_electrons', 'single electrons', '#185b83'), ('bi214', '²¹⁴Bi', '#b23a48'))
CLASS_NAMES = {'lt4mm': 'largest deposit < 4 mm from the main track', '4to8mm': '4–8 mm', 'gt8mm': '> 8 mm', 'none': 'no deposit of 3 keV or more'}
DRIFTS = (('production', 'II_L1000_D500: drift of the event, 1 m at 500 V/cm'), ('long', 'II_L1500_D350: drift uniform over 1.5 m at 350 V/cm'))


def pct(p, e=None):
    return f'{100 * p:.1f} %' if e is None else f'{100 * p:.1f} ± {100 * e:.1f} %'


def main(folder):
    out = Path(folder); plots = out / 'plots'; plots.mkdir(parents=True, exist_ok=True)
    with h5py.File(out / 'diffusion_limit.h5') as f:
        K = f['k_grid'][:]; E = f['E_grid_keV'][:]; S = f['sigma_grid_mm'][:]; classes = [text(c) for c in f['classes'][:]]
        # Julia arrays (k, E) and (sigma, k, E) come out of h5py with the axes reversed
        kE = lambda a: a[:].T; skE = lambda a: a[:].transpose(2, 1, 0)
        d = {s: {k: f[s][k][()] for k in ('selected', 'passed')} for s, _, _ in SAMPLES}
        for s, _, _ in SAMPLES:
            for name, _ in DRIFTS:
                d[s][name] = dict(flagged=kE(f[s][name]['flagged']), error=kE(f[s][name]['flagged_error']), sigma=f[s][name]['sigma_mm'][:])
            d[s]['fixed'] = dict(flagged=skE(f[s]['fixed']['flagged']), error=skE(f[s]['fixed']['flagged_error']))
        cls = {name: {c: dict(events=int(f['bi214'][name]['class'][c]['events'][()]), flagged=kE(f['bi214'][name]['class'][c]['flagged']),
                              n=kE(f['bi214'][name]['class'][c]['flagged_events'])) for c in classes} for name, _ in DRIFTS}
        drift = {name: dict(variant=text(f['drift'][name]['variant'][()]), smax=float(f['drift'][name]['sigma_max_mm'][()]),
                            E=float(f['drift'][name]['E_V_cm'][()]), L=float(f['drift'][name]['L_cm'][()])) for name, _ in DRIFTS}
        folder_truth = text(f['truth_folder'][()])
    k3 = int(np.argmin(np.abs(K - 3.0))); e10 = int(np.argmin(np.abs(E - 10.0)))

    # flagged fraction against k, 10 keV, both drifts
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, (name, title) in zip(axes, DRIFTS):
        for s, label, color in SAMPLES:
            ax.errorbar(K, 100 * d[s][name]['flagged'][:, e10], yerr=100 * d[s][name]['error'][:, e10], marker='o', color=color, label=label)
        ax.set_title(title, fontsize=10); ax.set_xlabel('k, distance threshold in units of σ'); ax.grid(alpha=0.3)
    axes[0].set_ylabel('events with ≥ 10 keV beyond k σ (%)'); axes[0].legend(fontsize=9)
    fig.tight_layout(); fig.savefig(plots / 'flagged_k.png', dpi=130); plt.close(fig)
    # against sigma, k = 3, 10 keV
    fig, ax = plt.subplots(figsize=(6, 4))
    for s, label, color in SAMPLES:
        ax.errorbar(S, 100 * d[s]['fixed']['flagged'][:, k3, e10], yerr=100 * d[s]['fixed']['error'][:, k3, e10], marker='o', color=color, label=label)
    ax.set_xlabel('σ (mm), the same for every event'); ax.set_ylabel('events with ≥ 10 keV beyond 3 σ (%)'); ax.grid(alpha=0.3); ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(plots / 'flagged_sigma.png', dpi=130); plt.close(fig)
    # sigma distributions
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for name, title in DRIFTS:
        ax.hist(d['bi214'][name]['sigma'], bins=30, histtype='step', label=title)
    ax.set_xlabel('σ of the passing ²¹⁴Bi events (mm)'); ax.set_ylabel('events'); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(plots / 'sigma.png', dpi=130); plt.close(fig)

    rows = []
    for name, title in DRIFTS:
        for b, e in enumerate(E):
            for a, k in enumerate(K):
                r = [title if (a == 0 and b == 0) else '', f'{e:.0f}', f'{k:.1f}']
                for s, _, _ in SAMPLES:
                    r.append(pct(d[s][name]['flagged'][a, b], d[s][name]['error'][a, b]))
                fb = d['bi214'][name]['flagged'][a, b]; fe = d['single_electrons'][name]['flagged'][a, b]; fs = d['bb0nu'][name]['flagged'][a, b]
                r += [f'{1 / (1 - fe):.2f}' if fe < 1 else '∞', f'{1 / (1 - fb):.2f}' if fb < 1 else '∞', f'{100 * 0.8 * (1 - fs):.1f} %']
                rows.append(r)
    main_table = table(['Drift', 'E_thr (keV)', 'k'] + [l for _, l, _ in SAMPLES] + ['single-electron rejection ×', '²¹⁴Bi rejection ×', 'signal efficiency after'], rows)
    class_rows = []
    for name, title in DRIFTS:
        for c in classes:
            x = cls[name][c]
            class_rows.append([title if c == classes[0] else '', CLASS_NAMES[c], x['events'], f"{x['n'][k3, e10]} ({pct(x['flagged'][k3, e10])})",
                               f"{x['n'][k3, 0]} ({pct(x['flagged'][k3, 0])})", f"{x['n'][1, e10]} ({pct(x['flagged'][1, e10])})"])
    class_table = table(['Drift', '²¹⁴Bi class (passing events)', 'events', 'flagged, k = 3, 10 keV', 'k = 3, 5 keV', 'k = 2, 10 keV'], class_rows)
    counts = table(['Sample', 'selected events with a score', 'passing the network at 80 % signal'],
                   [[l, d[s]['selected'], d[s]['passed']] for s, l, _ in SAMPLES])
    # orientative: the copper-shield rates of iItaca at the bound (k = 3, 10 keV, the 1.5 m drift)
    rates = ROOT / 'results/background/rates.h5'; orientative = ''
    if rates.exists():
        with h5py.File(rates) as f:
            g = f['detectors/iitaca']; bi = float(g['bi214/events_per_year'][()]); tl = float(g['tl208/events_per_year'][()])
            mxe = float(g['xenon_mass_kg'][()]) / 1000
        fb = d['bi214']['long']['flagged'][k3, e10]; fe = d['single_electrons']['long']['flagged'][k3, e10]; fs = d['bb0nu']['long']['flagged'][k3, e10]
        rows_o = [['²¹⁴Bi', f'{bi:.3f}', f'{bi * (1 - fb):.3f}', f'{bi / mxe:.4f}', f'{bi * (1 - fb) / mxe:.4f}'],
                  ['²⁰⁸Tl (not measured: between the single-electron and the ²¹⁴Bi fractions)', f'{tl:.3f}', f'{tl * (1 - fb):.3f} – {tl * (1 - fe):.3f}', f'{tl / mxe:.4f}', f'{tl * (1 - fb) / mxe:.4f} – {tl * (1 - fe) / mxe:.4f}'],
                  ['Total', f'{bi + tl:.3f}', f'{(bi + tl) * (1 - fb):.3f} – {bi * (1 - fb) + tl * (1 - fe):.3f}', f'{(bi + tl) / mxe:.4f}', f'{(bi + tl) * (1 - fb) / mxe:.4f} – {(bi * (1 - fb) + tl * (1 - fe)) / mxe:.4f}'],
                  ['Signal efficiency', '80 %', f'{100 * 0.8 * (1 - fs):.0f} %', '', '']]
        orientative = f'''<section class="panel"><h2>Orientative: the copper-shield background of iItaca at the bound</h2><p>The <a href="../../../results/background/index.html">rates of the 20 cm copper shield</a> for iItaca (1.5 m, 350 V/cm; 2 mm pads, network at 80 % signal) multiplied by the fraction of survivors a diffusion-limited readout would keep at k = 3 and 10 keV: ²¹⁴Bi by 1 − {100 * fb:.0f} %, the signal by 1 − {100 * fs:.0f} %. The ²⁰⁸Tl survivors were not examined; the range takes the single-electron fraction ({100 * fe:.0f} %) and the ²¹⁴Bi one as limits. A bound, not a prediction: no readout noise, no threshold per pixel and no network.</p>{table(['', 'events/yr, 2 mm pads', 'events/yr at the bound', 'events / (t yr), 2 mm pads', 'events / (t yr) at the bound'], rows_o)}</section>'''
    up = '../../../'
    page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>iItaca diffusion limit</title><style>{STYLE}img{{max-width:100%}}</style></head><body><header><a href="{up}notes/html/index.html">Documentation index</a> · <a href="{up}notes/html/iitaca/index.html">iItaca</a><h1>Diffusion-limited ion readout · bound from the Monte Carlo truth · II_L1000_D500_P2</h1><p>What a readout limited only by ion diffusion could reject among the events that pass the 2 mm pads, the track selection and the network.</p></header><main>
<section class="panel"><h2>Method</h2><p>The events are those of the 2 mm study (II_L1000_D500_P2) that pass the network at 80 % signal: the test bb0nu events of the training, its test single electrons, and every selected ²¹⁴Bi event. For each one the Monte Carlo truth gives the energy of the secondary deposits (every hit not on the main track, in the simulated box) at a distance of at least d from the nearest main-track hit, on a grid of 0.25 mm (<code>scripts/iitaca/prepare_iitaca_truth.jl</code>, <code>{folder_truth}</code>).</p>
<p>An event is <b>flagged</b> when it holds at least E<sub>thr</sub> of secondary energy beyond k σ from the main track, σ being the thermal spread of the ions after the drift of the event, σ = √(2 k<sub>B</sub> T L / e E): a readout blurred only by diffusion would show a detached deposit there. Two drifts: the one of the event in the production ({drift['production']['L']:.0f} cm at {drift['production']['E']:.0f} V/cm, σ up to {drift['production']['smax']:.2f} mm) and a drift drawn uniformly over {drift['long']['L']:.0f} cm at {drift['long']['E']:.0f} V/cm ({drift['long']['variant']}, σ up to {drift['long']['smax']:.2f} mm). The energy beyond k σ is read at the first grid point not below k σ, so it is never overestimated. The flagged fraction of the backgrounds is the most a diffusion-limited readout could add to the rejection of the 2 mm pads and the network; the flagged fraction of bb0nu is the signal it would cost. No readout noise, no threshold per pixel and no network are involved: this is a bound, not a prediction.</p>{counts}</section>
<section class="panel"><h2>Flagged fractions</h2><img src="plots/flagged_k.png" alt="flagged fraction against k"><p>Rejection factors are 1 / (1 − flagged fraction) on top of the current selection; the signal efficiency after is 80 % × (1 − flagged fraction of bb0nu).</p>{main_table}</section>
<section class="panel"><h2>²¹⁴Bi by class</h2>{class_table}</section>
{orientative}
<section class="panel"><h2>Fixed σ</h2><img src="plots/flagged_sigma.png" alt="flagged fraction against sigma"><img src="plots/sigma.png" alt="sigma of the passing events"></section>
</main></body></html>'''
    (out / 'index.html').write_text(page); print(out / 'index.html')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / 'results/iitaca/II_L1000_D500_P2_diffusion_limit')
