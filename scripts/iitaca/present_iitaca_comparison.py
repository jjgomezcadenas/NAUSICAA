#!/usr/bin/env python3
"""Page of the comparison of the signal with the backgrounds of an iItaca variant.

  python scripts/iitaca/present_iitaca_comparison.py COMPARISON_DIR

Reads COMPARISON_DIR/comparison.h5 (scripts/iitaca/prepare_iitaca_comparison.jl)
and writes COMPARISON_DIR/index.html with its plots. Nothing is recalculated.
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
COLORS = ('#b16d31', '#7a4fa3', '#b23a48', '#6b7480')


def load(g):
    r = dict(label=text(g['run_label'][()]), sample=text(g['sample'][()]), E1=g['E1_keV'][:], E2=g['E2_keV'][:],
             flow=dict(zip([text(s) for s in g['flow_step'][:]], [int(n) for n in g['flow_events'][:]])))
    if 'acceptance' in g:
        r.update(acceptance=g['acceptance'][:], error=g['acceptance_error'][:], accepted=g['accepted_events'][:])
    if 'curve' in g:
        r.update(curve=(g['curve/signal_efficiency'][:], g['curve/background_acceptance'][:]), area=float(g['curve/area'][()]))
    return r


def render(directory):
    out = Path(directory).resolve(); plots = out / 'plots'; plots.mkdir(exist_ok=True)
    with h5py.File(out / 'comparison.h5') as f:
        variant = text(f['variant'][()]); settings = tomllib.loads(text(f['settings'][()])); roi = (float(f['roi_lower_keV'][()]), float(f['roi_upper_keV'][()]))
        targets = f['targets'][:]; thresholds = f['thresholds_keV'][:]; efficiency = f['signal_efficiency'][:]
        signal = load(f['signal']); backgrounds = [load(f['backgrounds'][k]) for k in sorted(f['backgrounds'])]
    name = lambda r: SAMPLE.get(r['sample'], r['sample']); tracks = settings['tracks']

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    top = max([signal['E2'].max(initial=1)] + [b['E2'].max(initial=1) for b in backgrounds]); bins = np.linspace(0, top, 51)
    ax[0].hist(signal['E2'], bins=bins, density=True, histtype='stepfilled', color='#287ba0', alpha=.45, label=f"{name(signal)}: {len(signal['E2']):,} events, mean {signal['E2'].mean():.0f} keV")
    for b, c in zip(backgrounds, COLORS):
        if len(b['E2']):
            ax[0].hist(b['E2'], bins=bins, density=True, histtype='step', color=c, lw=1.6, label=f"{name(b)}: {len(b['E2']):,} events, mean {b['E2'].mean():.0f} keV")
    ax[0].axvline(thresholds[list(targets).index(.8)], color='#203343', ls='--', lw=1, label=f'threshold for 80 % of the signal')
    ax[0].set_xlabel('E₂, the smaller endpoint energy (keV)'); ax[0].set_ylabel('Events, normalised'); ax[0].legend(fontsize=8)
    ax[0].set_title(f'{variant or "default set"}: events in the ROI passing the track selection', fontsize=10)
    for b, c in zip(backgrounds, COLORS):
        if 'curve' in b:
            ax[1].plot(b['curve'][0], np.clip(b['curve'][1], 1e-5, 1), color=c, lw=1.6, label=f"{name(b)}: area {b['area']:.4f}")
    ax[1].set_yscale('log'); ax[1].set_xlim(.4, 1); ax[1].set_ylim(1e-3, 1); ax[1].set_xlabel('Signal efficiency of the E₂ cut'); ax[1].set_ylabel('Background acceptance'); ax[1].grid(alpha=.3, which='both')
    ax[1].legend(fontsize=8); ax[1].set_title('Cut on E₂: accepted when E₂ exceeds the threshold', fontsize=10)
    fig.tight_layout(); fig.savefig(plots / 'e2.png', dpi=130); plt.close(fig)

    runs = [signal] + backgrounds; steps = [s for s in ('generated', 'read', 'trigger', 'roi', 'track') if any(r['flow'].get(s, 0) for r in runs)]
    def fraction(r, s):
        base = r['flow'].get('generated') or r['flow']['read']
        return f"{r['flow'][s]:,} ({100 * r['flow'][s] / base:.3g} %)" if r['flow'].get(s) is not None and base else ''
    flow = table(['Step'] + [name(r) for r in runs], [[FLOW[s]] + [fraction(r, s) for r in runs] for s in steps])
    cut = table(['Signal efficiency of the cut', 'Threshold on E₂ (keV)'] + [f'{name(b)} accepted' for b in backgrounds],
                [[f'{100 * e:.1f} %', f'{t:.0f}'] + [f"{100 * b['acceptance'][i]:.2f} ± {100 * b['error'][i]:.2f} % ({int(b['accepted'][i])} of {len(b['E2']):,})" if len(b['E2']) else 'no events' for b in backgrounds]
                 for i, (e, t) in enumerate(zip(efficiency, thresholds))])
    depth = len(os.path.relpath(out, ROOT).split(os.sep)); up = '../' * depth
    links = ' · '.join(f"<code>{r['label']}</code>" for r in runs)
    page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>iItaca · signal and backgrounds</title><style>{STYLE}img{{max-width:100%}}</style></head><body><header><a href="{up}notes/html/index.html">Documentation index</a> · <a href="{up}notes/html/iitaca/index.html">iItaca</a><h1>{name(signal)} against {', '.join(name(b) for b in backgrounds)} · {variant or 'default set'}</h1><p>Selection step by step, and the cut on the endpoint energy for the selected events.</p></header><main>
<section><h2>Selection</h2><p>Events after each step; in brackets, the fraction of the first row. ROI: {roi[0]:.2f} to {roi[1]:.2f} keV.</p>{flow}</section>
<section><h2>Cut on the endpoint energy</h2><p>E₁ ≥ E₂ are the energies in spheres of {tracks['blob_radius_mm']:g} mm around the two ends of the main track. An event is accepted when E₂ exceeds the threshold. Efficiencies and acceptances are fractions of the events in the ROI that pass the track selection; errors are binomial.</p>{cut}<img src="plots/e2.png" alt="Distribution of E2 and background acceptance against signal efficiency"></section>
<p>Runs: {links}. Numerical results: <a href="comparison.h5">comparison.h5</a>.</p>
</main><footer>iItaca · signal and backgrounds · {variant or 'default set'}</footer></body></html>
'''
    (out / 'index.html').write_text(page); print(out / 'index.html')


if __name__ == '__main__':
    render(sys.argv[1])
