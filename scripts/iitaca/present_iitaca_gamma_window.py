#!/usr/bin/env python3
"""Page of the network acceptance of the gamma backgrounds of II_L1000_D500_P2, measured on a box-energy window.

  python scripts/iitaca/present_iitaca_gamma_window.py

Reads results/iitaca/II_L1000_D500_P2_gamma_window/gamma_window.h5
(scripts/iitaca/prepare_iitaca_gamma_window.jl) and writes index.html beside it.
Nothing is recalculated.
"""
import sys
from pathlib import Path
import h5py
sys.path.insert(0, str(Path(__file__).parent))
from present_iitaca_detector import STYLE, table, text

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/iitaca/II_L1000_D500_P2_gamma_window'
NAME = {'tl208': '²⁰⁸Tl', 'bi214': '²¹⁴Bi'}


def sci(x, e):
    """Value and error with the exponent of the value, for instance (3.1 ± 0.4) × 10⁻⁹."""
    if x == 0:
        return f'0 ± {e:.1e}'
    p = int(f'{x:e}'.split('e')[1])
    return f'({x / 10**p:.2f} ± {e / 10**p:.2f}) × 10<sup>{p}</sup>'


def big(x):
    return '∞' if x == float('inf') else f'{x:.2e}'.replace('e+0', ' × 10<sup>').replace('e+', ' × 10<sup>') + '</sup>'


def render():
    with h5py.File(OUT / 'gamma_window.h5') as f:
        lo, hi = float(f['window_lower_keV'][()]), float(f['window_upper_keV'][()]); training = text(f['training'][()])
        targets = f['target_signal_efficiency'][:]; data = {}
        for iso in ('tl208', 'bi214'):
            g = f[iso]
            data[iso] = {k: (g[k][()] if g[k].shape == () else g[k][:]) for k in g if isinstance(g[k], h5py.Dataset)}
            data[iso]['flow'] = dict(zip([text(s) for s in g['cut_flow/step'][:]], g['cut_flow/events'][:]))
    sections = []
    rows = []
    for iso, d in data.items():
        runs = [text(r) for r in d['window_runs']]
        rows.append([NAME[iso], ' + '.join(f'<code>{r}</code> ({n:,})' for r, n in zip(runs, d['window_events_per_run'])), f"{int(d['window_events']):,}"])
    sections.append(f'<section><h2>Window samples</h2><p>Events passing the track selection (<code>selected</code>) with box energy <code>energy_true_keV</code> in [{lo:.0f}, {hi:.0f}] keV. Within this window the event topology does not depend on energy; the network acceptance measured there applies to the events in the ROI. Network: <code>{training}</code>; thresholds from the signal test events of the training.</p>'
                    + table(['Isotope', 'Runs (events in the window)', 'Events'], rows) + '</section>')
    rows = []
    for j, t in enumerate(targets):
        r = [f"{100 * t:.0f} % ({100 * data['tl208']['signal_efficiency'][j]:.1f} %)", f"{data['tl208']['threshold'][j]:.4f}"]
        for iso, d in data.items():
            r += [f"{int(d['accepted_events'][j])} / {int(d['window_events']):,}", f"{100 * d['acceptance'][j]:.3f} ± {100 * d['acceptance_error'][j]:.3f}"]
        rows.append(r)
    sections.append('<section><h2>Network acceptance in the window</h2><p>Events above the threshold over the events of the window; binomial errors.</p>'
                    + table(['Signal efficiency (achieved on the test events)', 'Threshold'] + sum(([f'{NAME[i]} accepted', f'{NAME[i]} acceptance (%)'] for i in data), []), rows) + '</section>')
    rows = []
    for iso, d in data.items():
        fl = d['flow']
        rows.append([NAME[iso], f"<code>{text(d['roi_run'])}</code>", f"{int(d['generated_gammas']):,}", f"{int(d['saved_events']):,}", f"{int(fl['trigger']):,}", f"{int(fl['roi']):,}",
                     f"{int(d['roi_selected_events']):,}", sci(d['roi_selected_per_gamma'], d['roi_selected_per_gamma_error'])])
    roi = (float(data['tl208']['roi_lower_keV']), float(data['tl208']['roi_upper_keV']))
    sections.append(f'<section><h2>ROI events passing the track selection</h2><p>ROI {roi[0]:.1f} to {roi[1]:.1f} keV on the event energy. ²⁰⁸Tl: the 100 M run saved in 2420–2490 keV, which holds every event that can reach the ROI after the energy smearing; ²¹⁴Bi: the 40 M run saved in 2300–2470 keV. Binomial error.</p>'
                    + table(['Isotope', 'Run', 'Generated gammas', 'Saved by Geant4', 'Trigger', 'ROI', 'ROI and track selection', 'Per generated gamma'], rows) + '</section>')
    rows = []
    for j, t in enumerate(targets):
        r = [f'{100 * t:.0f} %']
        for iso, d in data.items():
            r += [sci(d['background_per_gamma'][j], d['background_per_gamma_error'][j]), big(d['rejection_saved'][j]), big(d['rejection_generated'][j])]
        rows.append(r)
    sections.append('<section><h2>Background after selection and network</h2><p>Background per generated gamma: ROI events passing the track selection per generated gamma × network acceptance in the window; relative errors added in quadrature. Rejection factors: events saved by Geant4 (in the saved window of the ROI run) and generated gammas, over the background after selection and network.</p>'
                    + table(['Signal efficiency'] + sum(([f'{NAME[i]} per gamma', f'{NAME[i]} rejection (saved)', f'{NAME[i]} rejection (generated)'] for i in data), []), rows) + '</section>')
    up = '../../../'
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>iItaca gamma window</title><style>{STYLE}</style></head><body>'
            f'<header><a href="{up}notes/html/index.html">Documentation index</a> · <a href="{up}notes/html/iitaca/index.html">iItaca</a><h1>Gamma backgrounds on 2 mm pads · II_L1000_D500_P2</h1>'
            f'<p>Network acceptance of ²⁰⁸Tl and ²¹⁴Bi measured on the box-energy window {lo:.0f}–{hi:.0f} keV, and the background per generated gamma in the ROI.</p></header><main>'
            + ''.join(sections) + '<section><h2>Reproduce</h2><p><code>julia --project=. scripts/iitaca/prepare_iitaca_gamma_window.jl</code>, then <code>python scripts/iitaca/present_iitaca_gamma_window.py</code>. Inputs: the score files and cut flows listed in <code>gamma_window.h5</code>, with their SHA-256.</p></section></main></body></html>')
    (OUT / 'index.html').write_text(page); print(OUT / 'index.html')


if __name__ == '__main__':
    render()
