#!/usr/bin/env python3
"""Present the saved 128 versus 64 substep timing comparison."""
from pathlib import Path
import h5py
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/readout/bb0nu_full_ion3mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42/electronics_sipm_substeps64'
with h5py.File(OUT/'controls.h5') as f:
    s=f['summary'];e=f['events']
    rows=''.join(f'<tr><td>{int(event)}</td><td>{1000*e["before_s"][i]:.2f}</td><td>{1000*e["after_s"][i]:.2f}</td></tr>' for i,event in enumerate(e['event_number']))
    accuracy=''
    if (OUT/'accuracy.h5').exists():
        with h5py.File(OUT/'accuracy.h5') as a:
            q=a['summary'];cases=a['cases']
            details=''.join(f'<tr><td>{cases["pattern"][i].decode()}</td><td>{int(cases["bucket"][i])}</td><td>{cases["mean_difference_pe"][i]:.4g}</td><td>{cases["variance_difference_pe2"][i]:.4g}</td><td>{100*cases["clipping_difference"][i]:.3g}</td></tr>' for i in range(len(cases['bucket'])))
            accuracy=f'''<h2>64-substep accuracy</h2><ul>
<li>Reference: individual-microcell simulation, starting fully recovered.</li>
<li>Six patterns: dark, weak, steady, clipping edge, isolated pulse, rising and falling signal. {a['repetitions'][()]} repetitions per pattern.</li>
<li>Both calculations use the same sampled optical counts and ASIC noise. Cell assignments and dark avalanches are independent.</li>
<li>Failures: charge mean {q['mean_failures'][()]}; variance {q['variance_failures'][()]}; clipping probability {q['clipping_failures'][()]}.</li>
<li>Largest relative mean difference for buckets with at least 3,000 incident photons: {q['maximum_bright_mean_difference_percent'][()]:.4f}%. Largest relative variance difference: {q['maximum_bright_variance_difference_percent'][()]:.3f}%.</li>
<li>Largest clipping-probability difference: {q['maximum_clipping_difference_pp'][()]:.2f} percentage points, within sampling uncertainty.</li></ul>
<p>Acceptance limits: mean difference ≤ max(0.05 PE, 0.5% of reference mean, 4 standard errors); variance difference ≤ max(0.01 PE², 5% of reference variance, 4 standard errors); clipping-probability difference ≤ max(0.5 percentage points, 4 standard errors).</p>
<details><summary>Per-bucket differences: 64 substeps minus reference</summary><table><tr><th>Pattern</th><th>Bucket</th><th>Mean (PE)</th><th>Variance (PE²)</th><th>Clipping (percentage points)</th></tr>{details}</table></details>'''
    page=f''' <!doctype html><html lang="en"><head><meta charset="utf-8"><title>ITACA · SiPM substep timing</title><style>body{{font:17px/1.6 system-ui;max-width:900px;margin:30px auto;color:#203343}}a,h1{{color:#185b83}}td,th{{padding:8px 18px;text-align:left}}</style></head><body>
<a href="../../../../notes/html/index.html">ITACA index</a><h1>SiPM timing: 128 versus 64 substeps</h1>
<ul><li>Same ten photon events; both start fully recovered.</li><li>Three repetitions per event after warm-up, alternating execution order. Each row shows median times.</li><li>Times include the full SiPM and ASIC calculation for all channels.</li><li>The timing comparison uses identical events. The separate accuracy test below checks 64 substeps; the active production card now selects 64.</li></ul>
<p>Mean: <strong>{1000*s['before_s'][()]:.2f} → {1000*s['after_s'][()]:.2f} ms/event</strong>. Speedup: {s['speedup'][()]:.3f}×.</p>
<table><tr><th>Event</th><th>128 substeps (ms)</th><th>64 substeps (ms)</th></tr>{rows}</table>{accuracy}</body></html>'''
    (OUT/'index.html').write_text(page)
print(OUT/'index.html')
