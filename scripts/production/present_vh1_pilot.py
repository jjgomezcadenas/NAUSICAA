#!/usr/bin/env python3
"""Present saved VH1/present-CSA comparison without rerunning detector physics."""
from pathlib import Path
import sys,html
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def render(path):
 out=path.parent;plots=out/'plots';plots.mkdir(exist_ok=True)
 with h5py.File(path) as f:
  g=f['events'];v={k:g[k][:] for k in g}
  n=len(v['event_number']);input_total=int(np.sum(v['ions_input']))
  invalid=int(np.sum(v['vh1_ions_on_invalid_pads']))
  qn=int(np.sum(v['vh1_quiet_pads']));qm=np.sum(v['vh1_quiet_sum'])/qn
  qvar=np.sum(v['vh1_quiet_sumsq'])/qn-qm*qm
  predicted=np.sum(v['vh1_quiet_pads']*v['vh1_predicted_quiet_variance'])/qn
  old_valid=v['ions_input']-v['old_ions_on_invalid_pads'];vh_valid=v['ions_input']-v['vh1_ions_on_invalid_pads']
  fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
  for a,key,title in [(ax[0,0],'vh1_bucket_residual_rms','VH1 bucket residual'),(ax[0,1],'old_bucket_residual_rms','Present CSA bucket residual')]:
   a.hist(v[key],bins=25,color='#287ba0');a.set(xlabel='RMS (ions)',ylabel='Events',title=title)
  ax[1,0].hist(v['vh1_saturated_pads'],bins=np.arange(-.5,max(v['vh1_saturated_pads'])+1.5),color='#c37b34');ax[1,0].set(xlabel='Saturated pads/event',ylabel='Events')
  ax[1,1].hist(v['vh1_ions_on_invalid_pads']/v['ions_input'],bins=25,color='#c37b34');ax[1,1].set(xlabel='Fraction of input ions on invalid pads',ylabel='Events')
  for a in ax.flat:a.grid(alpha=.15)
  fig.savefig(plots/'comparison.png',dpi=140);plt.close(fig)
  def row(name,val):return '<tr><th>'+html.escape(name)+'</th><td>'+html.escape(str(val))+'</td></tr>'
  rows=''.join([
   row('Events',n),row('Mean pads/event',f"{np.mean(v['ion_pads']):.1f}"),
   row('Mean bucket RMS: VH1 / present',f"{np.mean(v['vh1_bucket_residual_rms']):.2f} / {np.mean(v['old_bucket_residual_rms']):.2f} ions"),
   row('VH1 predicted bucket noise',f"{np.mean(v['vh1_predicted_bucket_noise']):.2f} ions"),
   row('Mean signal pad-bucket S/N',f"{np.sum(v['signal_snr_sum'])/np.sum(v['signal_pad_buckets']):.2f} (reconstructed ions / predicted noise; {np.sum(v['signal_pad_buckets']):,} cells with input ions)"),
   row('Mean input ions per signal pad-bucket',f"{np.sum(v['signal_input_ions'])/np.sum(v['signal_pad_buckets']):.1f} ions"),
   row('Quiet-pad sum variance: observed / predicted',f'{qvar:.1f} / {predicted:.1f} ions² (ratio {qvar/predicted:.3f})'),
   row('Saturated VH1 pads',f"{np.sum(v['vh1_saturated_pads'])} of {np.sum(v['ion_pads'])} ({100*np.sum(v['vh1_saturated_pads'])/np.sum(v['ion_pads']):.2f}%)"),
   row('VH1 ADC-limit pads',int(np.sum(v['vh1_adc_limit_pads']))),
   row('VH1 invalid pads / present invalid pads',f"{np.sum(v['vh1_invalid_pads'])} / {np.sum(v['old_invalid_pads'])}"),
   row('Input ions on invalid VH1 pads',f'{invalid:,} of {input_total:,} ({100*invalid/input_total:.2f}%)'),
   row('Mean reconstructed minus valid input: VH1 / present',f"{np.mean(v['vh1_ions_reconstructed']-vh_valid):.1f} / {np.mean(v['old_ions_reconstructed']-old_valid):.1f} ions/event"),
   row('Mean runtime: detector / VH1 / present',f"{np.mean(v['detector_s']):.3f} / {np.mean(v['vh1_simulation_s']+v['vh1_plan_s']+v['vh1_reconstruction_s']):.3f} / {np.mean(v['old_chain_s']):.3f} s/event"),
   row('Total pilot runtime',f"{f['total_s'][()]:.2f} s")])
  page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>VH1 100-event pilot</title><style>body{{font:17px/1.5 system-ui;max-width:1050px;margin:2em auto;color:#203343}}table{{border-collapse:collapse}}th,td{{text-align:left;padding:8px 14px;border-bottom:1px solid #ddd}}img{{max-width:100%}}</style><a href="../../../../../notes/html/index.html">ITACA index</a><h1>VH1 and present CSA: {n} bb0nu events</h1><p>Both circuits receive identical simulated ion counts and frames. Each uses its own noise and ADC stream. Residuals compare reconstructed bucket ions with input bucket ions on valid pads.</p><table>{rows}</table><img src="plots/comparison.png" alt="Ion residuals and VH1 saturation"><p><a href="controls.h5">Numerical results and example pad vectors</a></p></html>'''
  (out/'index.html').write_text(page)
 print(out/'index.html')
if __name__=='__main__':render(Path(sys.argv[1]))
