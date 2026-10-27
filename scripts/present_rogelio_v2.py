#!/usr/bin/env python3
"""Create a standalone report from the saved Julia v2 circuit assessment."""
from pathlib import Path
import os,json,hashlib
import h5py
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 source=ROOT/'results/ltspice/rogelio_v2/analysis.h5';diagnostic=source.parent/'diagnostic_analysis.h5'
 out=ROOT/'latex_generated/rogelio_v2';out.mkdir(parents=True,exist_ok=True)
 labels={'paper_reference':'Paper circuit','v2_native_corrected':'v2, native conditions',
         'v2_matched':'v2, matched conditions','v2_m3_large':'v2 with enlarged M3'}
 tuned=ROOT/'results/ltspice/rogelio_v2_retuned/analysis.h5'
 with h5py.File(source) as h,h5py.File(diagnostic) as diag,h5py.File(tuned) as ht:
  rows=[]
  for name,label in labels.items():
   a=dict(h[name].attrs)
   rows.append(f"{label} & {a['gain_V_per_ion']*1e6:.3f} & {a['noise_rms_V']*1e3:.3f} & {a['ENC_e_rms']:.1f} & {a.get('verified_N_snr5','--')}"+r' \\')
  ta=dict(ht['v2_retuned'].attrs)
  rows.append(f"v2, retuned gate & {ta['gain_V_per_ion']*1e6:.3f} & {ta['noise_rms_V']*1e3:.3f} & {ta['ENC_e_rms']:.1f} & {ta.get('verified_N_snr5','--')}"+r' \\')
  biasrows=[]
  for name,label in [('paper','Paper circuit'),('v2_native','v2, native conditions'),('v2_matched','v2, matched conditions'),('v2_retuned','v2, retuned gate')]:
   a=dict(diag[name].attrs)
   biasrows.append(f"{label} & {a['baseline_V']:.3f} & {a['gain_magnitude_1Hz']:.1f} & {a['M4_ro_Ohm']*1e-6:.1f}"+r' \\')
  checks=[]
  for name in list(labels)[1:]:
   a=dict(h[name].attrs)
   checks.append(f"{labels[name]} & {a['tau_s']:.4f} & {a['transient_relative_error']:.2g} & {a['spectral_relative_error']:.2g} & {a['tail_fit_max_relative_error']:.2g}"+r' \\')
  with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
   fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
   for name in ['paper_reference','v2_matched','v2_m3_large']:
    g=h[name];axes[0].plot(g['time_s'][:]*1e3,-g['voltage_V_per_ion'][:]*1e6,label=labels[name])
    axes[1].loglog(g['frequency_Hz'][:],g['output_asd_V_sqrtHz'][:]*1e6,label=labels[name])
   g=ht['v2_retuned'];axes[0].plot(g['time_s'][:]*1e3,-g['voltage_V_per_ion'][:]*1e6,label='v2, retuned gate')
   axes[1].loglog(g['frequency_Hz'][:],g['output_asd_V_sqrtHz'][:]*1e6,label='v2, retuned gate')
   axes[0].set(xlim=(-2,45),xlabel='Time relative to frame start (ms)',ylabel='Output drop per ion (µV)')
   axes[1].set(xlabel='Frequency (Hz)',ylabel='Output noise ASD (µV/√Hz)')
   for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.2)
   fig.savefig(out/'response_noise.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
  repro=h['reproduction'];native=dict(h['v2_native_corrected'].attrs);matched=dict(h['v2_matched'].attrs)
  template=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[margin=22mm]{geometry}
\usepackage{amsmath,booktabs,graphicx,hyperref}
\title{Assessment of Rogelio Palomo's v2 CSA}
\author{NAUSICA calculation record}\date{\today}
\begin{document}\maketitle
The v2 circuit reduces output noise but also reduces the voltage produced by
an arriving ion. Under the paper's input and readout conditions, its
input-equivalent noise is @ENC@ electrons RMS, compared with @REFERENCE@ for
the circuit currently used in the paper. The calculation retains the supplied
v2 transistor dimensions and gate voltages.

\section*{Circuits and conditions}
The v2 M4 has $W=0.54\,\mu$m, $L=2.4\,\mu$m and gate voltage 0.64 V.
Its M3 has $W=50\,\mu$m, $L=0.5\,\mu$m and gate voltage 0.6 V.
The paper circuit uses M4 with $W=2.88\,\mu$m, $L=1.44\,\mu$m and
gate voltage 0.556803 V; M3 has $W=200\,\mu$m, $L=2\,\mu$m and gate
voltage 0.538477 V. All cases use a 3 fF feedback capacitor and 1.8 V supply.

The native v2 conditions are 0.25 pF external input capacitance and
$27\,^{\circ}$C. The matched comparison uses 1.25 pF and $20\,^{\circ}$C
for both circuits. A third calculation combines the v2 M4 and its gate bias
with the paper's enlarged M3 and its adjusted gate bias.

\section*{Reproduction and model syntax}
LTspice 17.2.4 for macOS reproduces the supplied Windows 26.0.1 spectrum
over 1 Hz--1 GHz with maximum relative difference @REPRO@, using the
supplied model card unchanged. The card contains \texttt{NOIBE} where
\texttt{NOIB} is intended. The supplied log reports an unknown parameter.
The corrected calculation replaces this spelling and closes the model
parentheses, using the same corrected card as the paper. At 1 Hz the output
noise amplitude density changes from @OLDASD@ to @NEWASD@
$\mu$V/$\sqrt{\mathrm{Hz}}$.

\section*{Operating point and voltage gain}
\begin{center}\begin{tabular}{lrrr}\toprule
Circuit & Baseline (V) & $|V_{\rm out}/V_{\rm in}|$ at 1 Hz & M4 $1/g_{ds}$ (M$\Omega$)\\\midrule
@BIASROWS@
\bottomrule\end{tabular}\end{center}
The voltage ratio is measured by an AC current excitation with the feedback
network connected, preserving the operating point. The resistance is M4's
intrinsic small-signal drain resistance, read from the operating-point log
at its printed precision. It is distinct from the output resistance of the
complete amplifier. The matched v2 output lies near the positive supply;
M5 has only @M5VDS@ $\mu$V between drain and source. The increased M4
resistance therefore accompanies a much lower amplifier voltage gain at
the supplied bias.

Changing only the v2 M4 gate voltage to @RETUNEGATE@ V restores an output
baseline of @RETUNEBASE@ V. The transistor dimensions, M3 and other biases
are unchanged. The per-ion gain then rises to @RETUNEGAIN@ $\mu$V/ion,
demonstrating that the small signal at the supplied bias comes from the
operating point. The ENC is @RETUNEENC@ electrons RMS: output noise rises
along with signal gain. This variant is a separate bias-adjustment study.

\section*{Ion current, readout and sensitivity}
The current is that of a central square pad of 0.5 mm side in an infinite
grounded plane, with homogeneous gas in the upper half-space. The weighting
potential is one on the pad, zero on the surrounding plane, and zero at
infinity. Ions approach centrally at 3 m/s through the final 5 mm, with mean
collection times uniformly distributed over 20 ms. Each ion contributes
$e\,\Delta\phi_w$ on this finite path. The input waveform is reused from the
saved Julia calculation without renormalization or an extra collection pulse.

The voltage is averaged over centred 300 $\mu$s windows at consecutive frame
boundaries. The measured signal is
$D=-(\overline V_1-\alpha\overline V_0)$, where
$\alpha=\exp(-20\,\mathrm{ms}/\tau)$.
Each circuit's $\tau$ is obtained from its 100-ion, current-free tail over
21--120 ms. Noise uses the quiescent SPICE one-sided output spectrum over
0.1 Hz--1 MHz, filtered by the same averaging and frame difference.
The one-ion gain $g$ defines $\mathrm{ENC}=\sigma_D/g$. The final column
is checked using finite-count transients on both sides of $S/N=5$.

\begin{center}\begin{tabular}{lrrrr}\toprule
Circuit & $g$ ($\mu$V/ion) & $\sigma_D$ (mV) & ENC ($e$ RMS) & Minimum ions\\\midrule
@ROWS@
\bottomrule\end{tabular}\end{center}

\begin{figure}[htbp]\centering
\includegraphics[width=\linewidth]{response_noise.pdf}
\caption{Matched input capacitance and temperature. Left: voltage response
to the uniform-arrival current, normalized per ion. Right: quiescent output
noise spectrum.}
\end{figure}

\section*{Numerical checks}
Both 5 and 2.5 $\mu$s maximum time steps are evaluated, with transient
tolerances tightened from $(10^{-5},10^{-17}\,\mathrm A,10^{-8}\,\mathrm V)$
to $(5\times10^{-6},5\times10^{-18}\,\mathrm A,5\times10^{-9}\,\mathrm V)$
for relative, current and voltage tolerances. The more stringent inherited
voltage tolerance of $10^{-10}$ V caused the v2 zero-input transient to stall
in this simulator. The table gives the successful paired-refinement errors.
Frequency densities of 200 and 400 points per decade and independent
Gauss quadrature orders of 8 and 16 check the noise integral. Charge
conservation is checked against the finite-path input charge.

\begin{center}\small\begin{tabular}{lrrrr}\toprule
Circuit & $\tau$ (s) & Signal error & Noise error & Fit residual\\\midrule
@CHECKS@
\bottomrule\end{tabular}\end{center}
The gain is also checked for 10, 100, 564, 1000 and 2000 ions, and for the
counts that bracket the threshold. In the matched v2 case its largest change
relative to the one-ion gain is @NONLINEAR@\% over the tested range.
This weak nonlinearity is included in the threshold calculation.

\section*{Reproduction}
Configuration: \texttt{config/rogelio\_v2.toml}. Acquisition:
\texttt{scripts/run\_rogelio\_v2.py}. Julia analysis:
\texttt{src/ion\_signal/rogelio\_v2.jl}. Presentation:
\texttt{scripts/present\_rogelio\_v2.py}.
The pipeline command and full data locations are in
\texttt{results/ltspice/rogelio\_v2/README.md}. Supplied files and the
paper's historical results are preserved. This report is generated from HDF5.
\end{document}
'''
  values={'ENC':f"{matched['ENC_e_rms']:.1f}",'REFERENCE':f"{h['paper_reference'].attrs['ENC_e_rms']:.1f}",
          'REPRO':f"${repro['original'].attrs['max_relative_difference_supplied']:.2g}$",
          'OLDASD':f"{repro['original'].attrs['asd_1Hz_V_sqrtHz']*1e6:.2f}",
          'NEWASD':f"{repro['corrected'].attrs['asd_1Hz_V_sqrtHz']*1e6:.2f}",
          'M5VDS':f"{matched['M5_Vds_V']*1e6:.1f}",'NONLINEAR':f"{matched['max_gain_relative_change_to_2000']*100:.2f}",
          'RETUNEGATE':f"{ta['m4_gate_V']:.9f}",'RETUNEBASE':f"{ta['baseline_V']:.6f}",
          'RETUNEGAIN':f"{ta['gain_V_per_ion']*1e6:.3f}",'RETUNEENC':f"{ta['ENC_e_rms']:.1f}",
          'ROWS':'\n'.join(rows),'BIASROWS':'\n'.join(biasrows),'CHECKS':'\n'.join(checks)}
  for key,value in values.items():template=template.replace('@'+key+'@',value)
  (out/'report.tex').write_text(template)
  summary={name:{k:(v.item() if hasattr(v,'item') else v) for k,v in h[name].attrs.items()
                 if isinstance(v,(int,float)) or getattr(v,'dtype',None) is not None and v.dtype.kind in 'fiu'} for name in labels}
  summary['v2_retuned']={k:v.item() if hasattr(v,'item') else v for k,v in ta.items() if isinstance(v,(int,float)) or getattr(v,'dtype',None) is not None and v.dtype.kind in 'fiu'}
  (out/'values.json').write_text(json.dumps(summary,indent=2)+'\n')
 (out/'provenance.json').write_text(json.dumps({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (source,diagnostic,tuned,Path(__file__))},indent=2)+'\n')
if __name__=='__main__':main()
