#!/usr/bin/env python3
"""Compare the three saved VH1 pilots in one HTML table."""
from pathlib import Path
import h5py,numpy as np

ROOT=Path(__file__).resolve().parents[2]
base=ROOT/'results/readout'
ion4=base/'bb0nu_full_ion4mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42'
ion3=base/'bb0nu_full_ion3mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42'
runs=[('4 mm · 3 fF',ion4/'vh1_pilot_100'),('4 mm · 30 fF',ion4/'vh1_pilot_100_cf30f'),('3 mm · 20 fF',ion3/'vh1_pilot_100_cf20f')]
out=base/'vh1_pilots_100';out.mkdir(exist_ok=True)
results=[]
for label,path in runs:
 with h5py.File(path/'controls.h5') as f:
  g=f['events'];v={k:g[k][:] for k in g}
  row={'label':label,'path':path,'ids':v['event_number'],'pads':sum(v['ion_pads']),'saturated':sum(v['vh1_saturated_pads']),
       'invalid':sum(v['vh1_invalid_pads']),'loss':100*sum(v['vh1_ions_on_invalid_pads'])/sum(v['ions_input']),
       'noise':np.mean(v['vh1_predicted_bucket_noise']),'rms':np.mean(v['vh1_bucket_residual_rms']),
       'time':np.mean(v['vh1_simulation_s']+v['vh1_plan_s']+v['vh1_reconstruction_s'])}
  row['snr']=sum(v['signal_snr_sum'])/sum(v['signal_pad_buckets']) if 'signal_snr_sum' in v else None
  row['input']=sum(v['signal_input_ions'])/sum(v['signal_pad_buckets']) if 'signal_input_ions' in v else None
  results.append(row)
assert all(np.array_equal(results[0]['ids'],r['ids']) for r in results[1:])
rows=''
for r in results:
 href=(r['path']/'index.html').relative_to(out) if False else __import__('os').path.relpath(r['path']/'index.html',out)
 snr='—' if r['snr'] is None else f"{r['snr']:.2f}"
 input_ions='—' if r['input'] is None else f"{r['input']:.1f}"
 rows+=f'<tr><th><a href="{href}">{r["label"]}</a></th><td>{r["saturated"]}</td><td>{r["invalid"]}</td><td>{r["loss"]:.2f}%</td><td>{r["noise"]:.2f}</td><td>{r["rms"]:.2f}</td><td>{input_ions}</td><td>{snr}</td><td>{r["time"]:.3f}</td></tr>'
page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>VH1 pitch and capacitance pilots</title><style>body{{font:17px/1.5 system-ui;max-width:1200px;margin:2em auto;color:#203343}}table{{border-collapse:collapse}}th,td{{padding:9px;border-bottom:1px solid #ddd;text-align:left}}</style><a href="../../../notes/html/index.html">ITACA index</a><h1>VH1 · same 100 bb0nu events</h1><p>Ion pitch and feedback capacitance change together in the last row. All runs use the same selected event numbers and positions. Each report links to its controls and plots.</p><table><tr><th>Configuration</th><th>Saturated pads</th><th>Invalid pads</th><th>Input ions on invalid pads</th><th>Predicted bucket noise (ions)</th><th>Bucket residual RMS (ions)</th><th>Mean input ions / signal cell</th><th>Signal S/N</th><th>VH1 time (s/event)</th></tr>{rows}</table><p>Signal S/N is the mean of reconstructed ions divided by predicted noise per pad-bucket, selecting valid cells with input ions. The earlier 3 fF run did not save this diagnostic. The 3 mm and 4 mm rows use identical selection criteria. Differences in signal S/N reflect both pad sampling and the VH1 capacitance.</p><h2>Frame length</h2><p><a href="../bb0nu_itaca_pilot_10x100_ion20_15ms/vh1_frame_comparison/index.html">VH1 with 15 ms frames against 20 ms</a>: the same 130 bb0nu events of the ITACA pilot (batches 0 and 1), 3 mm pads and 20 fF, read at 20.0 cm/s with 15 ms frames (two 7.5 ms conversions) and at 15.0 cm/s with 20 ms frames; noise after reconstruction, S/N, saturation, residuals, energy and hits.</p><p><a href="../bb0nu_itaca_pilot_10x100_ion30_10ms/vh1_frame_comparison/index.html">VH1 with 10 ms frames against 20 ms</a>: the same events, read at 30.0 cm/s with 10 ms frames (two 5 ms conversions).</p></html>'''
(out/'index.html').write_text(page)
print(out/'index.html')
for r in results:print(r['label'],'saturated',r['saturated'],'noise',round(r['noise'],2),'SNR',r['snr'])
