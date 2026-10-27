#!/usr/bin/env python3
"""Present the saved bias sweep without repeating the physical calculation."""
from pathlib import Path
import os,hashlib,json
import h5py
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def main():
 source=ROOT/'results/ltspice/bias_sensitivity/analysis.h5'
 out=ROOT/'latex_generated/bias_sensitivity';out.mkdir(parents=True,exist_ok=True)
 with h5py.File(source) as h:
  data={name:dict(h[name].attrs) for name in h if name!='ranges'}
  def label(a):
   family=a['family'].decode() if isinstance(a['family'],bytes) else a['family']
   return 'Paper' if family=='paper' else 'v2 retuned'
  chosen=['paper_m4_6','paper_rounded_3','paper_rounded_4','paper_rounded_5','v2_m4_6','v2_rounded_3','v2_rounded_4']
  rows=[]
  for name in chosen:
   a=data[name]
   precision=int(name[-1]) if 'rounded' in name else 9
   rows.append(f"{label(a)} & {a['m3_gate_V']:.{precision}f} & {a['m4_gate_V']:.{precision}f} & {a['baseline_V']:.3f} & {a['gain_V_per_ion']*1e6:.2f} & {a['ENC_e_rms']:.2f} & {a['gain1000_relative_change']*100:.2f}"+r' \\')
  allrows=[]
  for name,a in data.items():
   allrows.append(f"{name.replace('_',r'\_')} & {a['m3_gate_V']:.7f} & {a['m4_gate_V']:.7f} & {a['baseline_V']:.3f} & {a['gain_V_per_ion']*1e6:.2f} & {a['ENC_e_rms']:.2f} & {a['gain1000_relative_change']*100:.2f}"+r' \\')
  with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
   fig,axes=plt.subplots(1,3,figsize=(11,3.7),layout='constrained')
   for prefix,legend in [('paper_m3_','Paper: M3 gate'),('paper_m4_','Paper: M4 gate'),('v2_m4_','v2: M4 gate')]:
    points=sorted([a for name,a in data.items() if name.startswith(prefix)],key=lambda a:a['offset_V'])
    x=[a['offset_V']*1e3 for a in points]
    for ax,key,scale in zip(axes,['baseline_V','gain_V_per_ion','ENC_e_rms'],[1,1e6,1]):
     ax.plot(x,[a[key]*scale for a in points],'o-',ms=3,label=legend)
   for ax,ylabel in zip(axes,['Output baseline (V)','Gain (µV/ion)','ENC (electrons RMS)']):
    ax.set(xlabel='Gate-voltage offset (mV)',ylabel=ylabel);ax.grid(alpha=.2)
   axes[0].legend(fontsize=8)
   fig.savefig(out/'bias_sweep.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
  text=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}
\usepackage[margin=18mm]{geometry}
\usepackage{amsmath,booktabs,longtable,graphicx,hyperref}
\title{CSA sensitivity to bias-voltage precision}
\author{NAUSICA calculation record}\date{\today}
\begin{document}\maketitle
The long decimal strings specify the numerical operating point used in the
simulation. Rounding the paper circuit's M4 and M3 gate voltages to
0.5568 and 0.5385 V preserves its gain, noise and near-linear response.
The same is true for the retuned v2 circuit with its M4 gate rounded to
0.6514 V. The exact output baseline of 0.879 V is a convenient reference
for comparing circuits; the sensitivity results establish a useful range
around that operating point.

\section*{Calculation}
The paper circuit and the retuned v2 circuit use their previously recorded
transistor dimensions. Temperature is $20\,^{\circ}$C, external input
capacitance 1.25 pF, feedback capacitance 3 fF and supply 1.8 V. One gate
voltage is varied at a time. Rounded settings change both tuned voltages
together in the paper circuit; v2 retains its M3 gate at 0.6 V.

The input is the saved mean uniform-arrival current for a 20 ms frame.
It comes from the central square 0.5 mm pad benchmark in an infinite
grounded plane, with homogeneous gas and weighting potential zero at
infinity. Ions approach through the final 5 mm at 3 m/s. Finite-path charge
$e\Delta\phi_w$ is retained without renormalization or a collection impulse.
The readout averages over centred 300 $\mu$s windows. Each bias point has
its own fitted tail constant and calibrated gain. Noise is the quiescent
output spectrum over 0.1 Hz--1 MHz after the same averaging and frame
difference. These choices match the previous circuit comparison.

Gain is obtained from the 100-ion signal. A 1000-ion transient tests its
amplitude dependence. The final table column is the relative change of
gain per ion between these two signals. This definition makes rounding
errors in calibration visible alongside the input-equivalent noise.

\section*{Rounded voltages}
\begin{center}\small\begin{tabular}{lrrrrrr}\toprule
Circuit & M3 gate (V) & M4 gate (V) & Baseline (V) & $g$ ($\mu$V/ion) & ENC ($e$) & Gain change (\%)\\\midrule
@ROWS@
\bottomrule\end{tabular}\end{center}

The four-decimal paper settings change ENC by @ROUNDENC@\% and retain
a gain change of @ROUNDLIN@\% between 100 and 1000 ions. Three-decimal
rounding shifts the baseline to @ROUND3BASE@ V and introduces a
@ROUND3LIN@\% gain change over the same charge range. The intrinsic ENC
remains similar, but the signal calibration becomes more amplitude-dependent.

\section*{Bias ranges}
Figure~\ref{fig:scan} shows the full scan. The following tested ranges
retain ENC within 5\% of nominal, gain within 10\% of nominal, and less
than 1\% gain change between 100 and 1000 ions:
\begin{itemize}
\item Paper M4 gate: offsets from @P4LOW@ to @P4HIGH@ mV.
\item Paper M3 gate: offsets from @P3LOW@ to @P3HIGH@ mV.
\item Retuned v2 M4 gate: offsets from @V4LOW@ to @V4HIGH@ mV.
\end{itemize}
The paper M3 gate is more sensitive: offsets of $\pm0.5$ mV produce
several-percent changes in gain between 100 and 1000 ions. Across the wider
scan, the output noise can fall together with the signal gain, leaving ENC
nearly unchanged. Therefore ENC alone does not define a useful bias range.
Sufficient signal voltage and a stable charge calibration also matter for
the digitisation chain.

\begin{figure}[htbp]\centering
\includegraphics[width=\linewidth]{bias_sweep.pdf}
\caption{One gate varied at a time, with the other biases fixed.}
\label{fig:scan}\end{figure}

The four-decimal values are adequate for reporting the present operating
point and reproducing its performance. A practical circuit can trim or
calibrate its bias to a suitable operating region. This scan uses the supplied
model at one temperature; process and temperature scans would determine
the adjustment range needed in hardware.

\section*{Checks and reproduction}
The rounded settings and representative offset points use paired maximum
steps of 5 and 2.5 $\mu$s with tightened solver tolerances, noise spectra
at 200 and 400 points per decade, and quadrature orders 8 and 16.
Input-charge conservation is checked for every point. The stored analyses
include the output extrema for 1000 ions and the input-node voltage shift.

Configuration is \texttt{config/bias\_sensitivity.toml}; Julia reduction is
\texttt{src/ion\_signal/bias\_sensitivity.jl}. Reproduction commands and the
list of refined points are in \texttt{results/ltspice/bias\_sensitivity/README.md}.
Python generates this report from \texttt{analysis.h5}.

\section*{Complete scan}
Indices 1--11 correspond to gate offsets
$-5,-2,-1,-0.5,-0.1,0,0.1,0.5,1,2,5$ mV.
\small
\begin{longtable}{lrrrrrr}\toprule
Case & M3 (V) & M4 (V) & Baseline (V) & $g$ ($\mu$V/ion) & ENC ($e$) & Gain change (\%)\\\midrule\endhead
@ALLROWS@
\bottomrule\end{longtable}
\end{document}
'''
  a=data['paper_rounded_4'];b=data['paper_rounded_3']
  values={'ROWS':'\n'.join(rows),'ALLROWS':'\n'.join(allrows),
          'ROUNDENC':f"{(a['ENC_relative_to_nominal']-1)*100:.2f}",'ROUNDLIN':f"{a['gain1000_relative_change']*100:.2f}",
          'ROUND3BASE':f"{b['baseline_V']:.3f}",'ROUND3LIN':f"{b['gain1000_relative_change']*100:.2f}"}
  for prefix,name in [('P4','paper_m4'),('P3','paper_m3'),('V4','v2_m4')]:
   ra=h['ranges/'+name].attrs
   values[prefix+'LOW']=f"{ra['low_offset_V']*1e3:g}"
   values[prefix+'HIGH']=f"{ra['high_offset_V']*1e3:g}"
  for k,v in values.items():text=text.replace('@'+k+'@',v)
  (out/'report.tex').write_text(text)
 (out/'provenance.json').write_text(json.dumps({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,Path(__file__))},indent=2)+'\n')
if __name__=='__main__':main()
