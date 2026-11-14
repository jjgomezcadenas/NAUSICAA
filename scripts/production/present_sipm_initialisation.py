#!/usr/bin/env python3
"""Present saved initial-state timing comparisons."""
from pathlib import Path
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/readout/bb0nu_full_ion3mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42/electronics_sipm_initialisation'
def main():
 with h5py.File(OUT/'controls.h5') as f:
  e=f['events'];s=f['summary']
  fig,ax=plt.subplots(figsize=(9,4.5),layout='constrained')
  ax.plot(e['event_number'][:],1000*e['before_s'][:],'o-',label='Sample each initial microcell')
  ax.plot(e['event_number'][:],1000*e['after_s'][:],'o-',label='Start fully recovered')
  ax.set_xlabel('Event number');ax.set_ylabel('SiPM + ASIC time (ms/event)');ax.grid(alpha=.2)
  ax.text(.03,.97,f"Events = {s['events'][()]}\nMean before = {1000*s['before_s'][()]:.1f} ms\nMean after = {1000*s['after_s'][()]:.1f} ms\nSpeedup = {s['speedup'][()]:.2f}×",transform=ax.transAxes,va='top',fontsize=10)
  ax.legend(loc='upper right');fig.savefig(OUT/'timing.png',dpi=130);plt.close(fig)
  page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><title>ITACA · SiPM initial recovery</title><style>body{{font:17px/1.6 system-ui;max-width:1000px;margin:30px auto;padding:0 20px;color:#203343}}h1,a{{color:#185b83}}img{{max-width:100%}}td,th{{padding:10px;text-align:left}}</style></head><body>
<a href="../../../../notes/html/index.html">ITACA index</a><h1>SiPM initial recovery</h1>
<ul><li>Start fully recovered: R = S = M.</li><li>Dark activity, PDE sampling, occupancy, recovery, 128 substeps and ASIC response remain enabled.</li><li>Ten identical input photon events for both settings; three timed repetitions per event, after warm-up. Each plotted point is the median of those repetitions.</li><li>Times cover all SiPMs and buckets in the event, including channels with zero photon input.</li></ul>
<img src="timing.png" alt="SiPM and ASIC time before and after changing initial recovery">
<table><tr><th></th><th>Sampled initial state</th><th>Fully recovered</th></tr>
<tr><td>Mean SiPM + ASIC time (s/event)</td><td>{s['before_s'][()]:.4f}</td><td>{s['after_s'][()]:.4f}</td></tr>
<tr><td>Mean allocations (MiB/event)</td><td>{s['before_MiB'][()]:.2f}</td><td>{s['after_MiB'][()]:.2f}</td></tr></table>
<p>Removed initial microcell draws: {s['mean_initial_draws_removed'][()]/1e6:.2f} million per event on this sample.</p>
<p>Historical cards retain stationary-dark initialization. New framed runs use the fully recovered card.</p>
</body></html>'''
  (OUT/'index.html').write_text(page)
 print(OUT/'index.html')
if __name__=='__main__':main()
