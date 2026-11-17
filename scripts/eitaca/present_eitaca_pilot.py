#!/usr/bin/env python3
"""Page of the eITACA pilot: does the SiPM readout saturate?

  python scripts/eitaca/present_eitaca_pilot.py [results/eitaca/pilot]

Reads <pilot>/pilot_summary.h5 (scripts/eitaca/prepare_eitaca_pilot.jl) and
writes <pilot>/index.html with its plots under <pilot>/plots/. Every number
comes from that file.
"""
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
BLUE, ORANGE, RED, GREY = '#287ba0', '#c46a1b', '#b23a48', '#8a97a0'
SAMPLES = (('bb0nu', 'bb0nu'), ('single_electrons', 'single electrons'))
STYLE = ("body{font:17px/1.65 system-ui;color:#203343;background:#f3f6f8;margin:0}header{background:#17394b;color:white;padding:30px max(24px,calc((100vw - 1100px)/2))}"
         "header a{color:#d4ecf6}main{max-width:1100px;margin:auto;padding:28px 24px}section{background:white;border:1px solid #d8e3e8;border-radius:10px;padding:25px;margin-bottom:24px}"
         "h1,h2,h3{line-height:1.25}h2{color:#185b83}a{color:#185b83}code{font-family:ui-monospace,monospace;font-size:.9em}table{border-collapse:collapse;width:100%;font-size:15px}"
         "th,td{text-align:left;padding:9px;border-bottom:1px solid #d8e3e8;vertical-align:top}th{background:#edf3f6}img{width:100%;height:auto}"
         ".answer{border-left:5px solid #185b83;background:#eef5f9;padding:14px 18px}footer{padding:20px;color:#526776}li{margin:6px 0}"
         "@media(max-width:650px){main{padding:16px}section{padding:18px}table{display:block;overflow:auto}}")


def text(v):
    return v.decode() if isinstance(v, bytes) else v


def load(path):
    out = {}
    with h5py.File(path) as f:
        out['peak_edges'] = f['peak_edges_pe'][()]
        for s, _ in SAMPLES:
            for arm in ('', '_guides', '_no_exclusion', '_guides_no_exclusion'):
                g = f[s + arm]
                out[s + arm] = dict(summary={k: v[()] for k, v in g['summary'].items()}, events={k: v[()] for k, v in g['events'].items()},
                                    cell_edges=g['cell_histogram/edges_pe'][()], cell_counts=g['cell_histogram/counts'][()],
                                    peak_counts=g['peak_histogram/counts'][()], optics=text(g['optics'][()]), exclusion=float(g['el_exclusion_mm'][()]),
                                    paired={k: v[()] for k, v in g['paired'].items()} if 'paired' in g else None)
            h = f['hybrid_production/' + s]
            out['hybrid_' + s] = {k: h[k][()] for k in ('events', 'events_with_invalid_electron_buckets', 'events_with_nan_energy')}
    return out


def stats(x):
    x = np.asarray(x, float)
    return f'entries {len(x):,}\nmean {np.mean(x):.4g}\nmedian {np.median(x):.4g}\nmax {np.max(x):.4g}' if len(x) else 'entries 0'


def box(ax, s, x=.02, y=.97, ha='left'):
    ax.text(x, y, s, transform=ax.transAxes, ha=ha, va='top', fontsize=7.5, family='monospace', bbox=dict(fc='white', ec='#d8e3e8', alpha=.9))


def plots(D, out):
    out.mkdir(exist_ok=True)
    limit = float(D['bb0nu']['summary']['asic_limit_pe'])
    names = {}
    # 1. Event peak
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    for ax, (s, title) in zip(axes, SAMPLES):
        e = D['peak_edges']
        for arm, color, label in (('', BLUE, 'eITACA, no guides'), ('_guides', ORANGE, 'hybrid optics and yield')):
            ax.stairs(D[s + arm]['peak_counts'], e, color=color, lw=1.6, label=label)
        ax.axvline(limit, color=RED, ls='--', lw=1.2, label=f'ASIC limit {limit:.0f} PE')
        ax.set_xscale('log'); ax.set_xlim(100, 2e4); ax.set_xlabel('Event peak: largest expected photoelectrons in one SiPM and one 250 ns bucket (PE)'); ax.set_ylabel('Events')
        ax.set_title(f'{title}: events accepted by the fast trigger, 50 mm EL exclusion', fontsize=10); ax.legend(fontsize=8, loc='upper right')
        box(ax, 'eITACA\n' + stats(D[s]['events']['peak_expected_pe']) + '\n\nhybrid optics\n' + stats(D[s + '_guides']['events']['peak_expected_pe']))
    fig.tight_layout(); names['peak'] = 'event_peak.png'; fig.savefig(out / names['peak'], dpi=130); plt.close(fig)
    # 2. All cells with photons
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    for ax, (s, title) in zip(axes, SAMPLES):
        for arm, color, label in (('', BLUE, 'eITACA, no guides'), ('_guides', ORANGE, 'hybrid optics and yield')):
            r = D[s + arm]; ax.stairs(r['cell_counts'], r['cell_edges'], color=color, lw=1.4, label=label)
        ax.axvline(limit, color=RED, ls='--', lw=1.2, label='ASIC limit')
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e-1, 1e4)
        ax.set_xlabel('Expected photoelectrons in one SiPM and one 250 ns bucket (PE)'); ax.set_ylabel('SiPM-bucket cells with photons')
        ax.set_title(f'{title}: all cells with at least one blue photon, 50 mm EL exclusion', fontsize=10); ax.legend(fontsize=8, loc='lower left')
        a, b = D[s]['summary'], D[s + '_guides']['summary']
        box(ax, f"eITACA: cells {int(a['cells_with_photons']):,}\n  clipped {int(a['cells_clipped_with_photons']):,}\n"
                f"hybrid optics: cells {int(b['cells_with_photons']):,}\n  clipped {int(b['cells_clipped_with_photons']):,}", x=.98, ha='right')
    fig.tight_layout(); names['cells'] = 'cell_pe.png'; fig.savefig(out / names['cells'], dpi=130); plt.close(fig)
    # 3. Peak against the distance of the track to the EL, without exclusion
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, (s, title) in zip(axes, SAMPLES):
        for arm, color, label in (('_no_exclusion', BLUE, 'eITACA, no guides'), ('_guides_no_exclusion', ORANGE, 'hybrid optics and yield')):
            e = D[s + arm]['events']; c = e['clipped_buckets'] > 0
            ax.scatter(e['track_distance_to_el_mm'][~c], e['peak_expected_pe'][~c], s=6, color=color, alpha=.6, label=label)
            ax.scatter(e['track_distance_to_el_mm'][c], e['peak_expected_pe'][c], s=16, facecolor='none', edgecolor=color, lw=1, label=label + ', clipped')
        ax.axhline(limit, color=RED, ls='--', lw=1.2, label='ASIC limit'); ax.axvline(50, color=GREY, ls=':', lw=1.2, label='50 mm EL exclusion')
        ax.set_xscale('symlog', linthresh=10); ax.set_xlim(0, 1500); ax.set_yscale('log'); ax.set_ylim(30, 4e4)
        ax.set_xlabel('Distance from the main track to the EL entry (mm)'); ax.set_ylabel('Event peak (PE per SiPM and 250 ns)')
        a, b = D[s + '_no_exclusion']['summary'], D[s + '_guides_no_exclusion']['summary']
        ax.set_title(f'{title}: events accepted without the EL exclusion', fontsize=10)
        ax.legend(fontsize=7, loc='lower left', title_fontsize=7.5,
                  title=f"entries {int(a['accepted_events'])}; clipped events: eITACA {int(a['events_with_clipping'])}, hybrid optics {int(b['events_with_clipping'])}")
    fig.tight_layout(); names['distance'] = 'peak_distance.png'; fig.savefig(out / names['distance'], dpi=130); plt.close(fig)
    # 4. Energy
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    ax = axes[0]
    for (s, title), color in zip(SAMPLES, (BLUE, ORANGE)):
        e = D[s]['events']; ok = e['clipped_buckets'] == 0; r = e['electron_energy_valid_keV'][ok] / e['energy_true_keV'][ok]
        ax.hist(r, bins=np.linspace(0.985, 1.015, 61), histtype='step', lw=1.5, color=color, label=title)
        box(ax, '\n'.join(f"{t}: entries {int(D[s2]['summary']['energy_ratio_unclipped_entries'])}\n  mean {D[s2]['summary']['energy_ratio_unclipped_mean']:.5f}"
                          f"  std {D[s2]['summary']['energy_ratio_unclipped_std']:.5f}" for s2, t in SAMPLES))
    ax.set_xlabel('Reconstructed electron energy / energy deposited in the box'); ax.set_ylabel('Events')
    ax.set_title('eITACA, 50 mm EL exclusion, events without clipped buckets', fontsize=10); ax.legend(fontsize=8, loc='upper right')
    ax = axes[1]
    for (s, title), marker in zip(SAMPLES, ('o', 's')):
        for arm, color in (('_no_exclusion', BLUE), ('_guides_no_exclusion', ORANGE)):
            e = D[s + arm]['events']; c = e['clipped_buckets'] > 0
            if c.any():
                loss = 100 * (1 - e['electron_energy_clipped_kept_keV'][c] / e['electron_energy_unclipped_keV'][c])
                ax.scatter(e['track_distance_to_el_mm'][c], loss, s=18, marker=marker, color=color,
                           label=f"{title}, {'eITACA' if arm == '_no_exclusion' else 'hybrid optics'}: {int(c.sum())} clipped")
    ax.axvline(50, color=GREY, ls=':', lw=1.2)
    ax.set_xlabel('Distance from the main track to the EL entry (mm)'); ax.set_ylabel('Energy loss by clipping (%)')
    ax.set_title('Clipped events without the EL exclusion: 1 - E(clipped buckets at the limit) / E(no limit)', fontsize=9)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=7)
    fig.tight_layout(); names['energy'] = 'energy.png'; fig.savefig(out / names['energy'], dpi=130); plt.close(fig)
    return names


def fmt(x, d=3):
    return f'{x:.{d}g}'


def render(pilot):
    pilot = Path(pilot).resolve(); D = load(pilot / 'pilot_summary.h5'); names = plots(D, pilot / 'plots')
    a0 = D['bb0nu']['summary']; limit = float(a0['asic_limit_pe'])
    rows = []
    for s, title in SAMPLES:
        for arm, label in (('', 'eITACA'), ('_guides', 'hybrid optics and yield'), ('_no_exclusion', 'eITACA, no EL exclusion'), ('_guides_no_exclusion', 'hybrid optics, no EL exclusion')):
            a = D[s + arm]['summary']
            rows.append([title, label, f"{int(a['input_events'])}", f"{int(a['accepted_events'])}", f"{a['seconds_per_accepted_event']:.2f}",
                         f"{a['detected_pe_mean']/1e6:.3f}", f"{a['peak_pe_median']:.0f}", f"{a['peak_pe_max']:.0f}", f"{a['peak_pe_over_limit_max']:.2f}",
                         f"{100*a['cells_clipped_fraction']:.2g}", f"{int(a['events_with_clipping'])} ({100*a['events_with_clipping_fraction']:.1f} %)",
                         f"{100*(1-a['sample_charge_retained_fraction']):.2g}", f"{100*a['clipped_event_energy_loss_max']:.2g}",
                         f"{100*a['recovery_loss_max']:.2g}", f"{100*a['peak_occupancy_max']:.2g}"])
    head = ['Sample', 'Optics', 'Input events', 'Accepted', 's per accepted event', 'Detected PE per event (10⁶)', 'Peak median (PE)', 'Peak max (PE)',
            'Peak max / ASIC limit', 'Cells with photons clipped (%)', 'Events with clipping', 'Sample charge lost (%)', 'Largest energy loss of a clipped event (%)',
            'Largest recovery loss (%)', 'Largest peak occupancy (%)']
    table = '<div style="overflow-x:auto"><table><tr>' + ''.join(f'<th>{h}</th>' for h in head) + '</tr>' + ''.join('<tr>' + ''.join(f'<td>{c}</td>' for c in r) + '</tr>' for r in rows) + '</table></div>'
    ratios = []
    for s, title in SAMPLES:
        for arm, label in (('', '50 mm EL exclusion'), ('_no_exclusion', 'no EL exclusion')):
            p = D[s + arm]['paired']
            ratios.append([title, label, f"{len(p['event_number'])}", f"{p['peak_ratio_mean']:.3f}", f"{p['peak_ratio_median']:.3f}", f"{p['detected_ratio_mean']:.3f}"])
    ratio_table = ('<table><tr><th>Sample</th><th>Selection</th><th>Paired events</th><th>Event peak ratio, mean</th><th>Event peak ratio, median</th><th>Detected PE ratio, mean</th></tr>'
                   + ''.join('<tr>' + ''.join(f'<td>{c}</td>' for c in r) + '</tr>' for r in ratios) + '</table>')
    hb = [D['hybrid_' + s] for s, _ in SAMPLES]
    se = D['single_electrons']['summary']
    acc_e, acc_g = float(a0['light_acceptance']), float(D['bb0nu_guides']['summary']['light_acceptance'])
    y_e, y_g = float(a0['photons_per_electron']), float(D['bb0nu_guides']['summary']['photons_per_electron'])
    clipped_default = int(a0['events_with_clipping']) + int(se['events_with_clipping'])
    nb, ns = D['bb0nu_no_exclusion']['summary'], D['single_electrons_no_exclusion']['summary']
    answer = (f"<p class='answer'><b>With the 50 mm EL exclusion, the eITACA SiPM readout {'does not saturate' if clipped_default == 0 else 'saturates in a few events'}:</b> "
              f"{clipped_default} of {int(a0['accepted_events']) + int(se['accepted_events'])} accepted events have a clipped bucket. The largest event peak is "
              f"{a0['peak_pe_max']:.0f} PE (bb0nu) and {se['peak_pe_max']:.0f} PE (single electrons) per SiPM and 250 ns bucket, "
              f"{a0['peak_pe_over_limit_max']:.2f} and {se['peak_pe_over_limit_max']:.2f} of the ASIC limit of {limit:.0f} PE. "
              f"Microcell recovery removes at most {100*max(a0['recovery_loss_max'], se['recovery_loss_max']):.2g} % of the charge of an event, and the reconstruction corrects it. "
              f"Without the exclusion, {int(nb['events_with_clipping'])} of {int(nb['accepted_events'])} bb0nu and {int(ns['events_with_clipping'])} of {int(ns['accepted_events'])} "
              f"single-electron events clip; the largest energy loss of a clipped event is {100*max(nb['clipped_event_energy_loss_max'], ns['clipped_event_energy_loss_max']):.2g} %. "
              f"With the exclusion, the clipped events lose at most {100*max(a0['clipped_event_energy_loss_max'], se['clipped_event_energy_loss_max']):.2g} % of their energy, and their main tracks lie "
              f"at most {max(a0['clipped_event_distance_to_el_max_mm'], se['clipped_event_distance_to_el_max_mm']):.0f} mm from the EL entry; "
              f"{int(a0['events_within_clipped_distance'])} bb0nu and {int(se['events_within_clipped_distance'])} single-electron accepted events lie that close. "
              f"With the hybrid optics and yield, {int(D['bb0nu_guides']['summary']['events_with_clipping']) + int(D['single_electrons_guides']['summary']['events_with_clipping'])} events clip with the exclusion "
              f"and {int(D['bb0nu_guides_no_exclusion']['summary']['events_with_clipping']) + int(D['single_electrons_guides_no_exclusion']['summary']['events_with_clipping'])} without it.</p>")
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>eITACA pilot · SiPM saturation</title><style>{STYLE}</style></head><body>
<header><p><a href="../../../notes/html/eitaca/index.html">eITACA</a></p><h1>eITACA pilot: does the SiPM readout saturate?</h1>
<p>Electron readout of the eITACA detector without light guides, on the Xe/He pilot events of the Mac; the same events with the optics and EL yield of the hybrid production as reference.</p></header><main>
<section><h2>Answer</h2>{answer}
<p>Without guides, each EL photon gives {acc_e:.4f} blue photons at the SiPM, {acc_e/acc_g:.3f} times the {acc_g:.4f} of the guides. The eITACA yield is {y_e:.0f} photons per electron against
{y_g:.0f} of the hybrid card, so the light per electron is {acc_e*y_e/(acc_g*y_g):.3f} times the hybrid one.</p></section>
<section><h2>Inputs and chain</h2><ol>
<li>Events: the Xe/He pilot <code>$NAUSICAA_DATA/pilot_10x100/{{bb0nu,electron}}/</code>, 10 batches of 100 generated events of each sample, generated in Xe–15 % He at 16 bar in a chamber of 1.8 m radius and 1.3 m length. The run cards override the eITACA radius (2.0 m) and length (1.5 m) to 1.8 m and 1.3 m so that the fast trigger and the drift lengths match the files. Light per SiPM does not depend on them. Transport and yield are those of 15.8 % He; the tracks were simulated in 15 % He.</li>
<li>Fast trigger: fixed 160 mm cube, 20 mm simulation margin, 3 keV gamma veto, main track at least 50 mm from the EL entry (<code>config/eitaca/trigger.toml</code>). The runs marked “no EL exclusion” set that distance to 0.</li>
<li>Electron response (<code>sim_electron_detector_response</code>): Fano pairs, attachment with 100 ms lifetime, diffusion D_T = 2.735 and D_L = 3.72 mm/√m (D*/√P at 16 bar), drift velocity 0.49 mm/µs, 4 mm SiPM lattice with a random offset per event, 250 ns buckets, EL photons and blue photons at the SiPM.</li>
<li>SiPM and ASIC (<code>sim_sipm_readout</code>), with the hybrid cards unchanged: <code>config/specs/response/sipm/elm_pitch4mm_3015ps_tau12p5ns_recovered.toml</code> (PDE 0.32, 21,808 illuminated microcells, 12.5 ns recovery, fully recovered start) and <code>config/specs/response/sipm_electronics/readout_v1.toml</code> (0.5 PE noise, 1 mA upper limit = {limit:.0f} PE per 250 ns bucket).</li>
<li>Reconstruction (<code>reco_electron_signal</code>, <code>rebin_electron_signal</code>, <code>reconstruct_hits</code>): mean dark subtracted, recovery inverse m/(1 − αm), clipped buckets invalid, attachment corrected, 20 buckets per hit bucket.</li>
<li>Random streams per event (<code>EItaca.event_seed</code>), with the seeds of the hybrid pilot: the primary pairs and diffusion draws are those of the hybrid production, and the eITACA and reference runs differ only in the light.</li></ol>
<p>Saturation criteria, as in the hybrid clipping pilots (<code>results/sipm_statistical_production_pilot/clipping_comparison.html</code>): ASIC clipping (cells with photons whose input exceeds the 1 mA limit, events with a clipped bucket, charge retained) and the microcell recovery loss, 1 − SiPM charge / (detected PE + dark avalanches). The event peak is the largest blue-photon count in one SiPM and one 250 ns bucket times the PDE. The peak occupancy is α × peak, α = τ/(M Δt): the fraction of the bucket's recovered-cell equivalents the peak uses. For an event with clipped buckets the readout is repeated with the same random streams and no ASIC limit, which gives the energy without clipping.</p></section>
<section><h2>Saturation</h2>{table}
<p>Sample charge lost: ASIC charge removed by clipping over the SiPM charge of all accepted events. Energy loss of a clipped event: 1 − E(clipped buckets kept at the limit) / E(no limit). In the reconstruction, as in the hybrid production, a clipped bucket is invalid and the event energy undefined.</p>
<img src="plots/{names['peak']}" alt="Event peak distributions"><img src="plots/{names['cells']}" alt="Photoelectrons per SiPM and bucket">
<img src="plots/{names['distance']}" alt="Event peak against the distance to the EL"></section>
<section><h2>eITACA against the hybrid optics, event by event</h2>{ratio_table}
<p>Ratios eITACA / reference for the same event, same pairs and diffusion draws. The light per electron gives the expected ratio {acc_e*y_e/(acc_g*y_g):.3f}.</p>
<p>Hybrid production of the same events (<code>results/readout/{{bb0nu,single_electrons}}_itaca_pilot_10x100</code>, 350 V/cm, Magboltz transport, guides): events with invalid electron buckets
{int(hb[0]['events_with_invalid_electron_buckets'])} of {int(hb[0]['events'])} (bb0nu) and {int(hb[1]['events_with_invalid_electron_buckets'])} of {int(hb[1]['events'])} (single electrons).</p></section>
<section><h2>Energy</h2><img src="plots/{names['energy']}" alt="Reconstructed energy"></section>
<section><h2>Files</h2><ul><li>Batches: <code>results/eitaca/pilot/IE_pilot_*/batches/b*/electrons.h5</code> (<code>scripts/eitaca/run_eitaca_batch.jl</code>, run cards <code>config/eitaca/runs/IE_pilot_*.toml</code>).</li>
<li>Summary: <code>results/eitaca/pilot/pilot_summary.h5</code> (<code>scripts/eitaca/prepare_eitaca_pilot.jl</code>); this page: <code>scripts/eitaca/present_eitaca_pilot.py</code>.</li></ul></section>
</main><footer>Generated from results/eitaca/pilot/pilot_summary.h5.</footer></body></html>"""
    (pilot / 'index.html').write_text(page)
    print(pilot / 'index.html')


if __name__ == '__main__':
    render(sys.argv[1] if len(sys.argv) > 1 else ROOT / 'results/eitaca/pilot')
