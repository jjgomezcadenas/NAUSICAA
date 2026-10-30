#!/usr/bin/env python3
"""Present the saved Julia longer-retention comparison without repeating calculations."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
import argparse,csv,hashlib,os,tomllib
import h5py,numpy as np
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/retention_matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def scalar(x):return x.decode() if isinstance(x,bytes) else x

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=ROOT/'config/retention_scan.toml');a=p.parse_args()
 c=tomllib.loads(a.config.read_text());data=ROOT/c['output'];out=ROOT/c['report'];out.mkdir(parents=True,exist_ok=True)
 rows=[];budgets=[];compression=[];validation=[];baseline_checks=[];cal_csv=[];comp_csv=[];noise_csv=[]
 with h5py.File(data/'presentation.h5') as h:
  names=[scalar(x) for x in h['variant_names'][:]]
  fig,axes=plt.subplots(1,3,figsize=(13.5,3.7))
  for frame,color in ((20,'#2563a6'),(40,'#b65228')):
   groups=[h[n][f'{frame}ms'] for n in names]
   tau=[float(g.attrs['tau_s']) for g in groups]
   enc=[float(g['ENC_e_rms'][0]) for g in groups]
   oldband=[float(g['ENC_by_low_cut_e'][2]) for g in groups]
   axes[0].plot(tau,enc,'o-',color=color,label=f'{frame} ms, 0.0001 Hz lower bound')
   axes[0].plot(tau,oldband,':',color=color,label=f'{frame} ms, 0.1 Hz lower bound')
   for sample,style in (('bb0nu','-'),('electrons','--')):
    rate=[]
    for n in names:
     g=h[n][f'histories/p{2 if frame==20 else 4}_{sample}'];j=int(np.flatnonzero(g['loss_fraction'][:]==.5)[0])
     rate.append(100*float(g['event_fraction'][j]))
    axes[1 if frame==20 else 2].plot(tau,rate,'o'+style,color=color,label='Double beta' if sample=='bb0nu' else 'Electron')
  axes[0].set_ylabel('Frame ENC (electrons RMS)')
  for ax,frame in zip(axes[1:],(20,40)):
   ax.set_ylabel('Events with a confirmed ≥50% loss (%)');ax.set_title(f'{frame} ms frames, {2 if frame==20 else 4} mm pitch')
  for ax in axes:
   ax.set_xlabel('Fitted discharge time (s)');ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=7)
  fig.tight_layout()
  for suffix in ('png','pdf'):fig.savefig(out/f'retention_comparison.{suffix}',dpi=180)
  plt.close(fig)
  for name in names:
   v=h[name];x=v['20ms'];y=v['40ms']
   length=1e6*float(v.attrs['length_m']);gate=float(v.attrs['m4_gate_V'])
   rows.append(f"{length:g} & {gate:.4f} & {float(v.attrs['baseline_V']):.3f} & {x.attrs['tau_s']:.2f} & {x['ENC_e_rms'][0]:.2f} & {y['ENC_e_rms'][0]:.2f} & {1e6*x['gain_V_per_ion'][0]:.2f} & {1e6*y['gain_V_per_ion'][0]:.2f} "+r'\\')
   for frame in ('20ms','40ms'):
    g=v[frame];devices=[scalar(s) for s in g['device_names'][:]]
    cal_csv.append([length,gate,float(v.attrs['baseline_V']),float(g.attrs['tau_s']),int(frame[:-2]),float(g['ENC_e_rms'][0]),1e6*float(g['gain_V_per_ion'][0]),float(g.attrs['maximum_frame_error_e'])])
    for device,enc,frac in zip(devices,g['device_enc_e'][:],g['device_variance_fraction'][:]):
     budgets.append(f"{length:g} & {frame[:-2]} & {device[2:-1]} & {enc:.2f} & {100*frac:.2f} "+r'\\')
     noise_csv.append([length,int(frame[:-2]),device[2:-1],float(enc),float(frac)])
    validation.append(f"{length:g} & {frame[:-2]} & {g.attrs['charge_error']:.2g} & {g.attrs['transient_error']:.2g} & {g.attrs['tau_refinement_error']:.2g} & {g.attrs['spectral_error']:.2g} "+r'\\')
    baseline_checks.append(f"{length:g} & {frame[:-2]} & {1e9*g.attrs['maximum_baseline_error_V']:.2g} & {g.attrs['maximum_frame_error_e']:.2g} "+r'\\')
   for sn in ('p2_bb0nu','p2_electrons','p4_bb0nu','p4_electrons'):
    g=v['histories/'+sn];rates=g['event_fraction'][:];loss=g['loss_fraction'][:]
    def percent(value):return 100*float(rates[np.flatnonzero(loss==value)[0]])
    label=sn.replace('_',r'\_')
    compression.append(f"{length:g} & {label} & {g.attrs['completed']} & {g.attrs['events']} & {percent(.25):.3f} & {percent(.5):.3f} & {percent(.75):.3f} "+r'\\')
    comp_csv.append([length,sn,int(g.attrs['completed']),int(g.attrs['events']),percent(.25),percent(.5),percent(.75),int(g.attrs['refined_histories']),float(g.attrs['maximum_refinement_error'])])
  reference=h[names[0]]
  baseline=f"The 400 $\\mu$m feedback device gives {reference['20ms']['ENC_e_rms'][0]:.2f} and {reference['40ms']['ENC_e_rms'][0]:.2f} electrons RMS for 20 and 40 ms frames in the extended noise band."
  last=h[names[-1]]
  def loss_rate(v,sample):
   g=v['histories/'+sample];j=int(np.flatnonzero(g['loss_fraction'][:]==.5)[0]);return 100*float(g['event_fraction'][j])
  m2=[]
  for frame in ('20ms','40ms'):
   g=last[frame];dev=[scalar(x) for x in g['device_names'][:]]
   j=dev.index('V(m2)');m2.append(100*float(g['device_variance_fraction'][j]))
  interpretation=(f"Extending the discharge time from {reference['20ms'].attrs['tau_s']:.2f} to {last['20ms'].attrs['tau_s']:.2f} s "
    f"lowers the 20 ms frame ENC from {reference['20ms']['ENC_e_rms'][0]:.2f} to {last['20ms']['ENC_e_rms'][0]:.2f} electrons RMS, "
    f"and the 40 ms frame ENC from {reference['40ms']['ENC_e_rms'][0]:.2f} to {last['40ms']['ENC_e_rms'][0]:.2f}. "
    f"M2 supplies {m2[0]:.1f}\\% and {m2[1]:.1f}\\% of the noise variance at the longest retention. "
    f"For double-beta events, the confirmed 50\\% local-loss rate changes from {loss_rate(reference,'p2_bb0nu'):.3f}\\% to {loss_rate(last,'p2_bb0nu'):.3f}\\% "
    f"with 2 mm pads, and from {loss_rate(reference,'p4_bb0nu'):.2f}\\% to {loss_rate(last,'p4_bb0nu'):.2f}\\% with 4 mm pads.")
 tex=r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm]{geometry}
\usepackage{amsmath,booktabs,longtable,graphicx,hyperref}
\title{Longer CSA feedback retention for NAUSICA}
\date{\today}
\begin{document}\maketitle
\section{Circuit and physical inputs}
The scan extends the existing 180 nm circuit by varying the feedback transistor
length. Its width is 0.36 $\mu$m and its gate is grounded. The external input
capacitance is 1.25 pF, the feedback capacitance is 3 fF and the supply is 1.8 V
at 20 $^\circ$C. The amplifier dimensions and other biases are held fixed.
The M4 bias is selected to restore the reference quiescent output and rounded
to four decimal places. These are transistor-model calculations.

The physical input is the established on-axis benchmark: a 0.5 mm square pad in
an infinite grounded plane, with pad weighting potential one and all other
boundaries, including infinity, zero. The dielectric is homogeneous. Ions
approach at constant speed 3 m/s through a finite 5 mm gap and are collected
uniformly within the specified frame. The input integral per ion remains
$e[1-\phi_w(g)]$; no extra collection impulse or charge renormalization is added.
The stored event histories use the corrected direct carrier collection.

\section{Discharge calibration and measurement}
Independent zero-current transients establish the baseline. Ten- and hundred-ion
packets provide six-second tails. Exponential fits are checked against a finer
solver grid, three tail intervals and the two charge amplitudes.
For the fitted $\tau$, $\alpha=\exp(-\Delta t/\tau)$ and the measured signal is
$D_n=-(\overline V_n-\alpha\overline V_{n-1})$ relative to the baseline.
Each voltage uses a centred 5 ms average. The small-signal gain is calibrated
with 100 ions and checked through 2000 ions.

SPICE supplies the one-sided output noise spectrum. Julia integrates the spectrum
with the averaging and frame-difference transfer function from 0.0001 Hz to
1 MHz. Lower bounds of 0.001 and 0.1 Hz are also saved. Two spectral densities,
two quadrature orders and two transient resolutions check convergence.
'''+baseline+r'''
\begin{center}\small
\begin{tabular}{rrrrrrrr}\toprule
$L_f$ & M4 bias & $V_0$ & $\tau$ & ENC (20 ms) & ENC (40 ms) & $g_{20}$ & $g_{40}$\\
($\mu$m) & (V) & (V) & (s) & ($e^-$ RMS) & ($e^-$ RMS) & ($\mu$V/ion) & ($\mu$V/ion)\\\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule\end{tabular}\end{center}
\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{retention_comparison.pdf}
\caption{Frame noise and confirmed local compression versus fitted retention.
The lower-frequency comparison uses the same spectrum and frame estimator.
Compression flags require a missing signal of at least one nominal frame ENC.}
\end{figure}

\section{Full-history compression}
All stored pad histories, including empty buckets, are screened with the
calibrated linear response. Histories with predicted output below 0.35 V are
simulated nonlinearly in LTspice, together with 12 fixed-seed controls in each
higher-voltage stratum. Each sample contains 8000 events. The 2 mm histories
use 20 ms frames and the 4 mm histories use 40 ms frames.
An event is counted once when a tested pad has a positive-count bucket whose
nonlinear response loses the selected fraction and at least one nominal ENC.
These fractions count confirmed local distortions within the screening study.
Electronic noise is used to set the deficit criterion; random noise is not
added to the mean transients. Event usefulness is assessed with track reconstruction.
\begin{center}\scriptsize
\begin{longtable}{rlrrrrr}\toprule
$L_f$ ($\mu$m) & Sample & Pads simulated & Events & $\geq25\%$ & $\geq50\%$ & $\geq75\%$\\\midrule
'''+ '\n'.join(compression)+r'''
\bottomrule\end{longtable}\end{center}

\section{Noise and compression}
'''+interpretation+r'''

\section{Transistor noise budget}
The entries are contributions to the recovered charge; independent transistor
noise variances add. Percentages refer to the total variance.
\begin{center}\scriptsize\begin{longtable}{rrlrr}\toprule
$L_f$ ($\mu$m) & Frame (ms) & Device & ENC contribution & Variance (\%)\\\midrule
'''+ '\n'.join(budgets)+r'''
\bottomrule\end{longtable}\end{center}

\section{Numerical validation}
Zero-input controls are also expressed as recovered charge. A dedicated study
worker allows a recorded 100 nV baseline bound when required, limited to less
than 0.005 ion by the calibrated gain. Independent control measurements check
that the frame error stays below 0.001 ion and the drift below 0.01\% of the
noise RMS. The default viewer keeps its original baseline policy.

\begin{center}\small\begin{longtable}{rrrr}\toprule
$L_f$ ($\mu$m) & Frame (ms) & Baseline drift (nV) & Frame error (ions)\\\midrule
'''+ '\n'.join(baseline_checks)+r'''
\bottomrule\end{longtable}\end{center}

The errors below compare independent charge integration, transient resolution,
tail-fit resolution and spectral sampling. The HDF5 results also retain the
fit residual, local tail constants, amplitude dependence, gain linearity and
quadrature convergence. Nonlinear event statistics include their refinement
and event-accounting checks.
\begin{center}\scriptsize\begin{longtable}{rrrrrr}\toprule
$L_f$ ($\mu$m) & Frame (ms) & Charge & Transient & Tail time & Spectrum\\\midrule
'''+ '\n'.join(validation)+r'''
\bottomrule\end{longtable}\end{center}
\section{Reproduction}
From the repository root:
\begin{verbatim}
python scripts/electronics/feedback/run_retention_scan.py \
  --simulator /path/to/LTspice
python scripts/electronics/feedback/present_retention_scan.py
\end{verbatim}
The pipeline writes input currents, acquisitions, analysis, operating points,
history cards and nonlinear statistics under \texttt{results/ltspice/retention\_scan}.
The presentation reads \texttt{presentation.h5}. Source snapshots and hashes
preserve acquisition and calculation provenance. Numerical HDF5 files, raw acquisitions, full screens
and case caches remain local. Report compilation is a separate explicit action.
\end{document}
'''
 if 'study_cards' in c:
  start=tex.index('From the repository root:')
  stop=tex.index('The pipeline writes',start)
  commands="From the repository root:\n\\begin{verbatim}\npython scripts/electronics/feedback/run_retention_comparison.py \\\n  --simulator /path/to/LTspice\npython scripts/electronics/feedback/present_retention_scan.py \\\n  --config config/retention_comparison.toml\n\\end{verbatim}\n"
  tex=tex[:start]+commands+tex[stop:]
  tex=tex.replace('results/ltspice/retention\\_scan','results/ltspice/retention\\_comparison')
 for filename,header,values in (
  ('calibration.csv',['feedback_length_um','M4_gate_V','baseline_V','tau_s','frame_ms','ENC_e_rms','gain_uV_per_ion','maximum_zero_input_frame_error_ions'],cal_csv),
  ('compression.csv',['feedback_length_um','sample','pads_simulated','events','event_percent_loss25','event_percent_loss50','event_percent_loss75','refined_histories','maximum_refinement_error'],comp_csv),
  ('noise_budget.csv',['feedback_length_um','frame_ms','device','ENC_contribution_e_rms','variance_fraction'],noise_csv)):
  with (data/filename).open('w',newline='') as stream:
   writer=csv.writer(stream,lineterminator="\n");writer.writerow(header);writer.writerows(values)
 (out/'report.tex').write_text(tex)
 (out/'sources.txt').write_text('\n'.join(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+str(path.relative_to(ROOT)) for path in (data/'analysis.h5',data/'presentation.h5',a.config.resolve(),Path(__file__).resolve()))+'\n')
 print(out/'report.tex')

if __name__=='__main__':main()
