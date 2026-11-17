#!/usr/bin/env python3
"""Page of the statistics of an iItaca run.

  python scripts/iitaca/present_iitaca_statistics.py RUN_DIR

Reads RUN_DIR/statistics/events.h5 (scripts/iitaca/prepare_iitaca_statistics.jl)
and writes RUN_DIR/statistics/index.html with its plots. Every number comes
from that file; the plots are histograms of its columns.
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

ROOT = Path(__file__).resolve().parents[2]
BLUE, GREEN, GREY = '#287ba0', '#3d8f5a', '#8a97a0'
SAMPLE = {'bb0nu': 'bb0nu', 'single_electrons': 'single electrons', 'bi214': '²¹⁴Bi', 'tl208': '²⁰⁸Tl'}
DECISIONS = (('outside_detector_events', 'track or energetic gamma deposit outside the detector'), ('main_track_near_el_events', 'main track near the EL entry'),
             ('track_larger_than_frame_events', 'track larger than the frame'), ('frame_outside_detector_events', 'frame outside the detector'),
             ('outside_frame_gamma_events', 'gamma deposit outside the frame'))
FLOW = {'generated': 'Generated in Geant4', 'read': 'In the input files', 'trigger': 'Accepted by the fast trigger', 'roi': 'Event energy in the ROI', 'track': 'Track selection'}


def numerics(x):
    return f'entries {len(x):,}\nmean {np.mean(x):.4g}\nstd {np.std(x):.3g}\nmedian {np.median(x):.4g}' if len(x) else 'entries 0'


def histogram(ax, x, bins, label, unit, title, lines=(), color=BLUE):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    ax.hist(x, bins=bins, color=color, alpha=.85)
    for v in lines:
        ax.axvline(v, color='#b23a48', lw=1.2, ls='--')
    ax.set_xlabel(f'{label} ({unit})' if unit else label); ax.set_ylabel('Events'); ax.set_title(title, fontsize=10)
    ax.text(.98, .97, numerics(x), transform=ax.transAxes, ha='right', va='top', fontsize=8, family='monospace', bbox=dict(fc='white', ec='#d8e3e8'))


def render(run_dir):
    run_dir = Path(run_dir).resolve(); out = run_dir / 'statistics'; plots = out / 'plots'; plots.mkdir(exist_ok=True)
    with h5py.File(out / 'events.h5') as f:
        e = {k: v[:] for k, v in f['events'].items()}
        counts = {k: int(v[()]) for k, v in f['selection/counts'].items()}
        flow = list(zip([text(s) for s in f['cut_flow/step'][:]], [int(n) for n in f['cut_flow/events'][:]]))
        label, sample, variant = text(f['run_label'][()]), text(f['sample'][()]), text(f['variant'][()])
        settings = tomllib.loads(text(f['settings'][()])); roi = (float(f['roi_lower_keV'][()]), float(f['roi_upper_keV'][()]))
        batches = len(f['batches'])
    name = SAMPLE.get(sample, sample); tracks = settings['tracks']; qbb = settings['energy']['qbb_keV']
    in_roi = e['in_roi'] == 1; selected = in_roi & (e['selected'] == 1); all_events = np.ones(len(in_roi), bool)
    where = f'{name}, {variant or "default set"}'

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    fwhm = settings['energy']['resolution_fwhm_at_qbb_keV']
    histogram(ax[0], e['energy_true_keV'], np.linspace(qbb - 160, qbb + 15, 71), 'Energy deposited in the box', 'keV', f'{where}: events accepted by the trigger', roi)
    histogram(ax[1], e['energy_keV'], np.linspace(qbb - 5 * fwhm, qbb + 5 * fwhm, 81), 'Event energy (emulated)', 'keV', 'Events accepted by the trigger; lines: ROI', roi)
    histogram(ax[2], e['ion_energy_keV'], 60, 'Ion energy (reconstructed ions × W)', 'keV', 'Events accepted by the trigger')
    fig.tight_layout(); fig.savefig(plots / 'energy.png', dpi=130); plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    top = int(max(4, e['segments'][in_roi].max(initial=0), e['deposits'][in_roi].max(initial=0))) + 1
    histogram(ax[0], e['segments'][in_roi], np.arange(-.5, top + .5), 'Track segments', '', f'{where}: events in the ROI', ())
    histogram(ax[1], e['deposits'][in_roi], np.arange(-.5, top + .5), 'Detached deposits', '', 'Events in the ROI', ())
    d = e['deposit_keV'][in_roi & (e['deposits'] > 0)]
    histogram(ax[2], d, np.linspace(0, max(100, d.max(initial=0)), 51), 'Largest detached deposit', 'keV', f"Events in the ROI with a detached deposit (at least {tracks['deposit_min_keV']:g} keV)", ())
    fig.tight_layout(); fig.savefig(plots / 'selection.png', dpi=130); plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    s = selected
    ax[0].scatter(e['E1_keV'][s], e['E2_keV'][s], s=4, color=BLUE, alpha=.5)
    ax[0].set_xlabel('E₁ (keV)'); ax[0].set_ylabel('E₂ (keV)'); ax[0].set_title(f"{where}: selected events; spheres of {tracks['blob_radius_mm']:g} mm", fontsize=10)
    ax[0].text(.98, .97, f'entries {s.sum():,}', transform=ax[0].transAxes, ha='right', va='top', fontsize=8, family='monospace', bbox=dict(fc='white', ec='#d8e3e8'))
    histogram(ax[1], e['E2_keV'][s], 50, 'E₂, the smaller endpoint energy', 'keV', 'Selected events', ())
    histogram(ax[2], e['length_mm'][s], 50, 'Length of the main track', 'mm', 'Selected events', ())
    fig.tight_layout(); fig.savefig(plots / 'track.png', dpi=130); plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    histogram(ax[0], e['ion_drift_cm'], 40, 'Ion drift length', 'cm', f'{where}: events accepted by the trigger', ())
    histogram(ax[1], e['hit_count'], 50, f"Hits above {tracks['saved_hit_cut_keV']:g} keV", '', 'Events accepted by the trigger', ())
    histogram(ax[2], e['ion_noise_ions'], 40, 'Noise per pad and frame after reconstruction', 'ions', 'Events accepted by the trigger', ())
    fig.tight_layout(); fig.savefig(plots / 'readout.png', dpi=130); plt.close(fig)

    rows = []; previous = None
    for step, n in flow:
        if step == 'generated' and n == 0:
            continue
        first = rows[0][3] if rows else n
        rows.append([FLOW[step], f'{n:,}', '' if previous is None else f'{100 * n / previous:.2f} %' if previous else '', n])
        previous = n
    reference = rows[0][3]
    flow_table = table(['Step', 'Events', 'Fraction of the previous step', f'Fraction of “{rows[0][0].lower()}”'],
                       [[r[0], r[1], r[2], f'{100 * r[3] / reference:.3g} %' if reference else ''] for r in rows])
    rejected = table(['Rejected by the fast trigger', 'Events', 'Fraction of the events read'],
                     [[what, f'{counts[k]:,}', f"{100 * counts[k] / counts['input_events']:.2f} %"] for k, what in DECISIONS])
    invalid = int((e['invalid_ion_buckets'] > 0).sum()); saturated = int((e['saturated_pads'] > 0).sum())
    depth = len(os.path.relpath(out, ROOT).split(os.sep)); up = '../' * depth
    page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>iItaca · {name} · statistics</title><style>{STYLE}img{{max-width:100%}}</style></head><body><header><a href="{up}notes/html/index.html">Documentation index</a> · <a href="{up}notes/html/iitaca/index.html">iItaca</a><h1>{name} · {variant or 'default set'}</h1><p>Run <code>{label}</code>, {batches} batches. Selection and statistics of the events accepted by the fast trigger.</p></header><main>
<section><h2>Selection</h2>
{flow_table}
<p>ROI: event energy between {roi[0]:.2f} and {roi[1]:.2f} keV. Track selection: components at {tracks['segment_connection_mm']:g} mm, no detached deposit of {tracks['deposit_min_keV']:g} keV or more, one segment or at most {tracks['max_segments']} joined at {tracks['connection_mm']:g} mm (<a href="{up}notes/html/event_analysis/track_selection.html">track selection</a>).</p>
{rejected}
<p>Accepted events with a pad-frame without valid reconstruction: {invalid:,}; with a pad outside the linear window of VH1: {saturated:,}.</p></section>
<section><h2>Energy</h2><p>The event energy is the energy deposited in the simulated box smeared with {100 * fwhm / qbb:.2f} % FWHM at Q<sub>ββ</sub>; the ion energy is the sum of the reconstructed ions.</p><img src="plots/energy.png" alt="Energy deposited, emulated event energy and ion energy"></section>
<section><h2>Track selection</h2><img src="plots/selection.png" alt="Segments, detached deposits and largest deposit of the events in the ROI"></section>
<section><h2>Selected tracks</h2><img src="plots/track.png" alt="Endpoint energies and length of the selected tracks"></section>
<section><h2>Readout</h2><img src="plots/readout.png" alt="Drift length, hits and noise"></section>
<p>Numerical results: <a href="events.h5">events.h5</a>. Chain and columns: <a href="{up}notes/html/iitaca/simulation.html">simulation chain</a>.</p>
</main><footer>iItaca · {label}</footer></body></html>
'''
    (out / 'index.html').write_text(page); print(out / 'index.html')


if __name__ == '__main__':
    render(sys.argv[1])
