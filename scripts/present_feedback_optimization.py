#!/usr/bin/env python3
"""Present saved Julia feedback/averaging results and their validation."""
from pathlib import Path
import hashlib,json,os,tomllib
import h5py,numpy as np
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 c=tomllib.loads((ROOT/'config/feedback_optimization.toml').read_text());source=ROOT/c['output']/'analysis.h5';out=ROOT/c['report'];out.mkdir(parents=True,exist_ok=True)
 selected='length_400um_gate0_rounded';rows=[];params=[];scans=[];rejections=[]
 with h5py.File(source) as h:
  assert selected in h
  sg=h[selected];sa=dict(sg.attrs)
  summary=[]
  for frame in ('20ms','40ms'):
   ref=h['reference'][frame];best=sg[frame];j=int(best.attrs['best_index'])-1
   for label,g,index in [('Reference',ref,0),('Reference, longer average',ref,int(ref.attrs['best_index'])-1),('Feedback candidate',best,0),('Feedback candidate, longer average',best,j)]:
    summary.append(f"{frame[:-2]} & {label} & {g['width_s'][index]*1e3:g} & {g['ENC_e_rms'][index]:.2f} & {g['noise_rms_V'][index]*1e3:.3f} & {g['mean_signal_V'][index]*1e3:.3f} & {g['mean_snr'][index]:.3f}"+r' \\')
  for name in sorted(h):
   g=h[name];a=dict(g.attrs);safe=name.replace('_',r'\_')
   if 'rejected_operating_point' in a:rejections.append(safe+f" & DC output {a['baseline_V']:.3f} V"+r' \\');continue
   if 'incomplete' in a:rejections.append(safe+r' & Incomplete transient acquisition'+r' \\');continue
   if '40ms' not in g:continue
   p=g['40ms'];pa=dict(p.attrs)
   params.append(f"{safe} & {a['width_m']*1e6:g} & {a['length_m']*1e6:g} & {a['gate_V']*1e3:g} & {a.get('m4_gate_V',.556803417206):.4f} & {a['baseline_V']:.3f} & {pa.get('tau_s',float('nan')):.3f}"+r' \\')
   if 'best_index' not in pa:rejections.append(safe+r' & Fails readout or linearity criteria'+r' \\');continue
   j=int(pa['best_index'])-1
   rows.append(f"{safe} & {p['ENC_e_rms'][0]:.2f} & {p['width_s'][j]*1e3:g} & {p['ENC_e_rms'][j]:.2f} & {p['mean_snr'][j]:.3f} & {100*p['max_gain_change'][j]:.3f}"+r' \\')
   scans.append((name,p['ENC_e_rms'][j]))
  p=sg['40ms'];j=int(p.attrs['best_index'])-1
  v40={k:p[k][j] for k in ('gain_V_per_ion','noise_rms_V','ENC_e_rms','mean_snr','mean_signal_V','previous_frame_fraction','following_frame_fraction','max_gain_change','transient_error','spectral_error','quadrature_error')}
  tau=float(p.attrs['tau_s']);v20=sg['20ms'];j20=int(v20.attrs['best_index'])-1
  v20_values={k:v20[k][j20] for k in ('ENC_e_rms','mean_snr')}
  reference_power=float(h['reference'].attrs['power_W']);selected_power=float(sa['power_W'])
  validation='\n'.join(label+f" & {value:.3g}"+r' \\' for label,value in [('40 ms transient refinement',v40['transient_error']),('40 ms spectral refinement',v40['spectral_error']),('40 ms noise quadrature',v40['quadrature_error']),('40 ms charge conservation',float(p.attrs['charge_error'])),('40 ms tail-fit residual',float(p.attrs['fit_error'])),('40 ms decay-time refinement',float(p.attrs['tau_refinement_error']))])
  budget='\n'.join(f"{name} & {value:.1f}"+r' \\' for name,value in zip(p['device_names'].asstr()[:],100*p['device_variance_fraction'][:]))
  with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
   fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
   names=['reference','gate_0mV_retuned','length_200um_retuned','length_400um_retuned',selected]
   for name in names:
    g=h[name]['40ms'];axes[0].plot(g['width_s'][:]*1e3,g['ENC_e_rms'][:],'o-',label=name.replace('_',' '))
   axes[0].set(xlabel='Averaging window (ms)',ylabel='ENC (elementary charges RMS)',xscale='log');axes[0].legend(frameon=False,fontsize=7);axes[0].grid(alpha=.2)
   names0=h['reference/40ms/device_names'].asstr()[:];names1=p['device_names'].asstr()[:]
   assert list(names0)==list(names1)
   x=np.arange(len(names0));axes[1].bar(x-.18,h['reference/40ms/device_enc_e'][:],width=.36,label='Reference, 5 ms')
   axes[1].bar(x+.18,p['device_enc_e'][:],width=.36,label='Candidate, 5 ms')
   axes[1].set(xticks=x,xticklabels=[n[2:-1] for n in names0],ylabel='Contribution to ENC (RMS)',xlabel='Transistor');axes[1].legend(frameon=False,fontsize=8)
   fig.savefig(out/'noise_comparison.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(out/'noise_comparison.png',dpi=160);plt.close(fig)
 text=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=20mm]{geometry}
\usepackage{amsmath,booktabs,graphicx,longtable,hyperref}
\title{NAUSICA: feedback transistor and averaging optimization}
\author{Circuit study}\date{\today}
\begin{document}\maketitle
\section{Result}
The selected tested candidate has feedback transistor $M_f$ with
$W=0.36\,\mu$m, $L=400\,\mu$m and its gate connected to ground. The
M4 gate bias is 0.6164 V. The other transistor dimensions and biases retain
the paper's values. The external pad capacitance is 1.25 pF, feedback
capacitance 3 fF, supply 1.8 V and temperature 20 degrees Celsius.

With 5 ms averaging and 40 ms buckets, the candidate gives ENC=@ENC@ elementary
charges RMS and S/N=@SNR@ for the assumed mean of 400 ions. The corresponding
20 ms result is ENC=@ENC20@ and S/N=@SNR20@ for 120 ions. The 400-ion voltage
measurement is @SIGNAL@ mV and its noise is @NOISE@ mV RMS. The output decay
is @TAU@ s. CSA quiescent power is @POWER@ $\mu$W, compared with
@REFPOWER@ $\mu$W for the reference. Five milliseconds is the largest averaging window in this scan;
the selected value is the best tested setting within the stated range and criteria.

\begin{center}\small\begin{tabular}{rlrrrrr}\toprule
Frame & Setting & Average & ENC & Noise & Signal & S/N\\
(ms) & & (ms) & ($e$ RMS) & (mV RMS) & (mV) &\\\midrule
@SUMMARY@
\bottomrule\end{tabular}\end{center}
\includegraphics[width=\linewidth]{noise_comparison.pdf}

\section{Physical input and readout}
The inputs retain the central 0.5 mm square pad in an infinite grounded plane,
in a homogeneous gas half-space. The pad weighting potential is one; the
surrounding plane and infinity are at zero. Ions follow a central straight
trajectory through the final 5 mm at constant 3 m/s. Collection times have
uniform mean density over each 20 or 40 ms bucket. The finite induced charge
$e[1-\phi_w(g)]$ is preserved. The assumed mean populations, 120 and 400 ions,
are design inputs for 2 and 4 mm pitch. The local benchmark and external
capacitance stay fixed; the focusing field and event redistribution are not
recalculated by this circuit study.

Each voltage reading is a centred average over $T_a$. The measured frame
signal is $D_n=-(\overline V_n-\alpha\overline V_{n-1})$, with
$\alpha=\exp(-\Delta t/\tau_0)$. A matched zero-input transient defines the
baseline. An exponential fit to a 100-ion current-free tail, from 1 ms after
collection ends to 120 ms after first collection, supplies $\tau_0$.
The signal gain is calibrated with that same readout and a 100-ion packet.
The assumed mean packet is independently simulated in SPICE.

Stationary thermal and flicker noise from the transistor models is integrated
from 0.1 Hz to 1 MHz after the averaging and corrected-difference filters:
\[
 \sigma_D^2=\int S_V(f)\operatorname{sinc}^2(fT_a)
 |1-\alpha e^{-2\pi i f\Delta t}|^2\,df,\qquad
 \mathrm{ENC}=\sigma_D/g.
\]
The noiseless signal is zero before the finite-gap current starts. This exact
pre-pulse baseline extends the archived prehistory when a 5 ms centred window
starts before the simulation's original 2 ms prehistory.

\section{Scan and operating point}
The initial scan changes only the feedback gate or geometry. Several settings
move the DC output close to a rail. A second scan adjusts the M4 gate to restore
the reference output near 0.879 V, then repeats the signal and noise calculation.
The combined scan tests longer feedback channels with the feedback gate at ground.
Finally the selected M4 bias is rounded to four decimal places and its complete
signal and noise runs are repeated. The table below shows the actual simulated
DC baseline. Its bias entries are displayed with four decimal places.

Candidates are ranked by calibrated ENC. The accepted settings have DC output
between 0.2 and 1.6 V, less than 1\% gain change between the 100-ion calibration
and the tested mean/2000-ion packets, and less than 5\% response in either
neighbouring frame. They also satisfy the numerical checks below. These are
scan criteria, not a final reconstruction threshold or a routed-layout result.

\small
\begin{longtable}{lrrrrrr}
\toprule
Variant & $W_f$ & $L_f$ & $V_f$ & $V_{M4}$ & Baseline & $\tau_0$\\
 & ($\mu$m) & ($\mu$m) & (mV) & (V) & (V) & (s)\\\midrule\endhead
@PARAMETERS@
\bottomrule\end{longtable}
\begin{longtable}{lrrrrr}
\toprule
Variant & ENC, 0.3 ms & Best $T_a$ & Best ENC & S/N, 400 ions & Gain change\\
 & ($e$ RMS) & (ms) & ($e$ RMS) & & (\%)\\\midrule\endhead
@SCAN@
\bottomrule\end{longtable}
\normalsize

\section{Time response and retention}
For the selected 40 ms setting, a packet contributes @PREVIOUS@\% of its
own-bucket signal to the following empty measurement. An equal packet collected
in the following bucket contributes @FOLLOWING@\% to the selected measurement,
through induction before collection and the centred averaging window. These
coefficients describe the response to the uniform packets used here.
The 5 ms averaging window corresponds to 0.5 mm drift at 0.1 m/s, compared with
a 4 mm nominal longitudinal bucket. The signal calibration includes the averaging.

Weaker feedback increases the retention time to @TAU@ s. Charge from a long
track or earlier events therefore remains on the output for longer. The
2000-ion transients test local amplitude response; complete pad histories
are needed to establish cumulative headroom and the eventual reset schedule.
The feedback channel area is $WL=144\,\mu\mathrm{m}^2$; its packing and parasitic
capacitances must be included when the layout is updated.

\section{Noise sources and numerical validation}
The selected 40 ms variance budget is:
\begin{center}\begin{tabular}{lr}\toprule
Device & Fraction of variance (\%)\\\midrule
@BUDGET@
\bottomrule\end{tabular}\end{center}
All completed accepted variants compare 5 and 2.5 $\mu$s maximum transient
steps, 200 and 400 spectral points per decade, and noise quadrature orders
8 and 16. Charge integrals agree with the finite-path input. The 10- and
100-ion tail fits check amplitude dependence. The reference runs reproduce
the saved 20 and 40 ms results. For the selected 40 ms case:
\begin{center}\begin{tabular}{lr}\toprule
Check & Relative error\\\midrule
@VALIDATION@
\bottomrule\end{tabular}\end{center}
Rejected fixed-bias points remain recorded:
\begin{center}\small\begin{tabular}{ll}\toprule
Variant & Reason\\\midrule
@REJECTIONS@
\bottomrule\end{tabular}\end{center}

\section{Reproduction}
From the repository root:
\begin{verbatim}
PYTHONPATH=src python scripts/run_feedback_pipeline.py \
  --julia /path/to/julia --simulator /path/to/LTspice
\end{verbatim}
The configuration files are \texttt{config/feedback\_optimization.toml} and
\texttt{config/feedback\_combined.toml}. Julia performs the operating-point search,
reduction and checks. Python acquires SPICE output and presents the saved HDF5
results. Presentation alone uses \texttt{scripts/present\_feedback\_optimization.py}.
The selected exact netlists, model card and input currents are under
\texttt{results/ltspice/feedback\_optimization/decks/length\_400um\_gate0\_rounded/}.
Relative PWL references point to the shared sibling \texttt{inputs/} directory.
\end{document}
'''
 values={'ENC':f"{v40['ENC_e_rms']:.2f}",'SNR':f"{v40['mean_snr']:.3f}",'ENC20':f"{v20_values['ENC_e_rms']:.2f}",'SNR20':f"{v20_values['mean_snr']:.3f}",'SIGNAL':f"{v40['mean_signal_V']*1e3:.3f}",'NOISE':f"{v40['noise_rms_V']*1e3:.3f}",'TAU':f'{tau:.3f}','POWER':f'{selected_power*1e6:.1f}','REFPOWER':f'{reference_power*1e6:.1f}','SUMMARY':'\n'.join(summary),'PARAMETERS':'\n'.join(params),'SCAN':'\n'.join(rows),'BUDGET':budget,'VALIDATION':validation,'REJECTIONS':'\n'.join(rejections),'PREVIOUS':f"{100*v40['previous_frame_fraction']:.3f}",'FOLLOWING':f"{100*v40['following_frame_fraction']:.3f}"}
 for key,value in values.items():text=text.replace('@'+key+'@',value)
 (out/'report.tex').write_text(text)
 (out/'provenance.json').write_text(json.dumps({'source':str(source.relative_to(ROOT)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'generator':str(Path(__file__).relative_to(ROOT)),'selected':selected},indent=2)+'\n')
 print('Report:',out/'report.tex')
if __name__=='__main__':main()
