#!/usr/bin/env python3
"""Present Julia's historical/direct history comparison without recalculating it."""
from pathlib import Path
import tomllib,hashlib
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
LABELS={'p2_bb0nu':r'2 mm / 20 ms, $0\nu\beta\beta$','p2_electrons':'2 mm / 20 ms, electrons','p4_bb0nu':r'4 mm / 40 ms, $0\nu\beta\beta$','p4_electrons':'4 mm / 40 ms, electrons'}
def main():
 c=tomllib.loads((ROOT/'config/pad_history_direct_comparison.toml').read_text());source=ROOT/c['output']/'comparison.h5';out=ROOT/c['report'];out.mkdir(parents=True,exist_ok=True)
 with h5py.File(source) as h:
  rates=[];pulses=[];coverage=[];eventrows=[]
  for name,label in LABELS.items():
   g=h[name]
   for version in ('historical','direct'):
    group=g[version];a=group.attrs;loss=group['criterion_loss_fraction'][:];noise=group['criterion_noise_multiple'][:]
    selected=[i for i,(r,s) in enumerate(zip(loss,noise)) if s==1 and r in (.1,.25,.5,.75)]
    rates.append(label+' & '+version+' & '+' & '.join(f'{100*group["compression_event_fraction"][i]:.3f}' for i in selected)+r' \\')
    pulses.append(label+' & '+version+f" & {a['mean_ions_per_occupied_cell']:.1f} & {a['median_ions_per_occupied_cell']:.0f} & {a['mean_event_mean_ions']:.1f} & {a['maximum_bucket_ions']:,} "+r'\\')
    coverage.append(label+' & '+version+f" & {a['events']:,} & {a['occupied_cells']:,} & {a['simulated_pads']:,} & {a['control_pads']} "+r'\\')
   if name=='p2_electrons':
    fig,axes=plt.subplots(2,1,figsize=(8,5.4),sharex=True)
    for ax,version in zip(axes,('historical','direct')):
     e=g[version]['event_3790'];ax.stairs(e['ions_per_frame'][:],e['time_edges_ms'][:]/1000)
     ax.set_ylabel('Ions per 20 ms bucket');ax.set_title(version.capitalize());ax.grid(alpha=.2)
     eventrows.append(f"{version} & {e.attrs['original_buckets']} & {e.attrs['occupied_buckets']} "+r'\\')
    axes[-1].set_xlabel("Time relative to each dataset's first retained slice (s)")
    fig.tight_layout();fig.savefig(out/'event_3790.pdf');fig.savefig(out/'event_3790.png',dpi=160);plt.close(fig)
  fig,ax=plt.subplots(figsize=(7.4,4.4))
  for name,label in LABELS.items():
   for version,style in (('historical','--'),('direct','-')):
    g=h[name][version];cuts=list(zip(g['criterion_loss_fraction'][:],g['criterion_noise_multiple'][:]));indices=[i for i,(_,s) in enumerate(cuts) if s==1]
    ax.plot([100*cuts[i][0] for i in indices],[100*g['compression_event_fraction'][i] for i in indices],style,label=label+', '+version)
  ax.set(xlabel='Minimum relative response loss (%)',ylabel='Events with confirmed loss (%)',ylim=(0,None));ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=7,ncols=2);fig.tight_layout();fig.savefig(out/'compression_comparison.pdf');fig.savefig(out/'compression_comparison.png',dpi=160);plt.close(fig)
  body=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=22mm]{geometry}
\usepackage{amsmath,booktabs,graphicx,hyperref}
\title{CSA performance with direct carrier collection}
\date{\today}
\begin{document}\maketitle
\section{Inputs and unchanged circuit}
The corrected ion histories assign sampled Gaussian carrier positions directly
to physical pixels and time buckets. The historical histories assigned occupied
diffusion-cell centres to those grids. Both versions use 8000 events per sample.
Checks confirm identical source event IDs, event seeds, sampled drift distances,
diffusion widths and produced carrier counts. Each file passes its independent
collected/retained/outside ledger checks. The configurations are 2 mm pitch with
20 ms frames and 4 mm pitch with 40 ms frames.

The electronics are unchanged: 180 nm optimized CSA, external capacitance
1.25 pF, $C_f=3$ fF, grounded feedback gate with $M_f$ width 0.36 $\mu$m and
length 400 $\mu$m, and $M_4$ gate at 0.6164 V. Five-ms centred averages are
compared using the same fixed low-charge gain and fitted discharge correction.
Nominal quiescent ENC remains ENC20 electrons at 20 ms and ENC40 at 40 ms.

Input currents use uniform collection times within each saved bucket. The
single-ion current remains the constant-speed final approach at 3 m/s over a
5 mm gap to a 0.5 mm square pad in a grounded plane. Its integral is
$e[1-\phi_w(g)]$ per ion. Event simulations start at the quiescent operating
point. No random noise is added to the mean transient responses.

\section{Coverage}
\begin{center}\scriptsize\begin{tabular}{llrrrr}\toprule
Sample & Input & Events & Occupied cells & Simulated pads & Controls\\\midrule
COVERAGE
\bottomrule\end{tabular}\end{center}
The complete linear screen retains every pad and empty frame. The nonlinear
batch simulates every history with a predicted minimum below 0.35 V, together
with twelve fixed-seed controls in each higher-voltage stratum per sample.
The separate twelve-case study retains matched baselines and two resolutions.

\section{Input pulses}
\begin{center}\small\begin{tabular}{llrrrr}\toprule
Sample & Input & Mean ions/cell & Median & Mean event mean & Maximum\\\midrule
PULSES
\bottomrule\end{tabular}\end{center}
These means and medians refer to positive pad/frame cells across the sample.
They differ from an average of per-event means, which is also saved in the
comparison HDF5. No charge threshold is applied.

\section{Confirmed event compression}
Relative response loss compares the nonlinear and linear frame signals for the
same pad history. A bucket is flagged when it contains collected ions and
exceeds both the selected fractional loss and one nominal ENC of missing
signal. An event is counted once if any tested pad contains a flagged bucket.

\begin{center}\scriptsize\begin{tabular}{llrrrr}\toprule
Sample & Input & $\geq10\%$ & $\geq25\%$ & $\geq50\%$ & $\geq75\%$\\\midrule
RATES
\bottomrule\end{tabular}\end{center}
Entries are percentages of each full 8000-event sample. They count confirmed
cases in the conservative screened study and its controls. A flag identifies
local distortion; event usefulness will be evaluated with reconstructed tracks.
Zero observed flags in the tested 2 mm cases is an empirical result for this
sample and screening coverage.

\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{compression_comparison.pdf}
\caption{Historical and direct-collected inputs, evaluated with the same circuit
and one-ENC missing-signal criterion.}
\end{figure}

\section{Event 3790}
\begin{center}\begin{tabular}{lrr}\toprule
Electron input & Event buckets & Occupied buckets across all pads\\\midrule
EVENTROWS
\bottomrule\end{tabular}\end{center}
The extended arrival interval follows the source event geometry. Direct
collection fills the numerical gaps introduced by coarsening the diffusion
cloud before readout. Each plotted axis begins at that dataset's first retained
slice; the reference shifts slightly when the binned approximation is removed.
\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{event_3790.pdf}
\caption{Full-event collected counts before and after the numerical correction.}
\end{figure}

\section{Reproduction}
From the repository root:
\begin{verbatim}
python scripts/run_pad_history_pipeline.py \
  --config config/pad_history_readout_direct.toml \
  --simulator /path/to/LTspice
python scripts/run_pad_compression.py \
  --config config/pad_history_viewer_direct.toml --workers 4
julia --project=. scripts/electronics/histories/compare_inputs.jl
python scripts/electronics/histories/present_input_comparison.py
python scripts/view_pad_histories.py
\end{verbatim}
Julia writes \texttt{results/ltspice/pad\_history\_direct\_comparison/comparison.h5};
Python presents its saved values. The source files, configurations, screen
outputs and nonlinear statistics have recorded hashes. Historical results are
preserved, and the corrected viewer uses the direct-input configuration by
default. The report source is generated without compiling a PDF.
\end{document}
'''
  for key,value in {'COVERAGE':'\n'.join(coverage),'PULSES':'\n'.join(pulses),'RATES':'\n'.join(rates),'EVENTROWS':'\n'.join(eventrows),'ENC20':f"{h['p2_bb0nu'].attrs['nominal_ENC_e']:.1f}",'ENC40':f"{h['p4_bb0nu'].attrs['nominal_ENC_e']:.1f}"}.items():body=body.replace(key,value)
  (out/'report.tex').write_text(body)
  (out/'sources.txt').write_text(hashlib.sha256(source.read_bytes()).hexdigest()+'  '+str(source.relative_to(ROOT))+'\n')
 print(out/'report.tex')
if __name__=='__main__':main()
