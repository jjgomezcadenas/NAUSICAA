#!/usr/bin/env python3
"""Render saved SPICE spectra and Julia reductions; no physics recomputation."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[4]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import os
import h5py
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
DATA=ROOT/'results/ltspice/spice2_reproduction'
OUT=ROOT/'latex_generated/spice2_reproduction';OUT.mkdir(parents=True,exist_ok=True)
PREVIEW=ROOT/'notes/physics/work/spice2-docs'
with h5py.File(DATA/'raw.h5') as h:
    new={k:h['updated_0p1Hz'][k][:] for k in h['updated_0p1Hz']}
    old={k:h['previous_0p1Hz'][k][:] for k in h['previous_0p1Hz']}
with h5py.File(DATA/'analysis.h5') as h:
    numbers={k:h['updated_0p1Hz'][k][:] for k in h['updated_0p1Hz']}
    checks=dict(h.attrs)
    op={k:h['operating_point'][k][0] for k in h['operating_point']}
    closure=h['updated_0p1Hz'].attrs['device_variance_closure']
f=new['frequency']
with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
    fig,axes=plt.subplots(6,1,figsize=(8,10),sharex=True,layout='constrained')
    panels=[(['V(onoise)'],1e3,'Output noise (mV/√Hz)'),(['V(m1)'],1e6,'M1 total (µV/√Hz)'),
            (['V(m1.1overf)'],1e6,'M1 flicker (µV/√Hz)'),(['V(m1.id)','V(m1.rd)'],1e9,'M1 thermal (nV/√Hz)'),
            (['gain'],1e-12,'Gain (TΩ)'),(['inoise'],1e12,'Input noise (pA/√Hz)')]
    for ax,(names,scale,ylabel) in zip(axes,panels):
        for name in names:ax.semilogx(f,new[name]*scale,label=name)
        ax.set(ylabel=ylabel.replace(' (', '\n('),xlim=(.1,1e9),ylim=(0,None))
        ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=8,loc='upper right')
    axes[-1].set_xlabel('Frequency (Hz)')
    fig.savefig(OUT/'noise_reproduction.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(PREVIEW/'noise_reproduction.png',dpi=140);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
    ax.loglog(f,new['V(onoise)'],label='Updated model (NOIMOD=2)')
    ax.loglog(old['frequency'],old['V(onoise)'],label='Previous model')
    ax.set(xlabel='Frequency (Hz)',ylabel='Output noise (V/√Hz)',xlim=(.1,1e9))
    ax.grid(alpha=.2);ax.legend(frameon=False)
    fig.savefig(OUT/'noise_comparison.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(PREVIEW/'noise_comparison.png',dpi=140);plt.close(fig)
rows='\n'.join(f'{f:g} & {a*1e3:.6g} & {b*1e6:.6g} & {100*c:.4g} '+r'\\' for f,a,b,c in zip(numbers['frequency_Hz'],numbers['output_ASD_V_sqrtHz'],numbers['M1_flicker_ASD_V_sqrtHz'],numbers['flicker_variance_fraction']))
text=r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[margin=22mm]{geometry}
\usepackage{graphicx,booktabs,amsmath}
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt}
\begin{document}
\section*{Reproduction of Rogelio's CSA noise simulation}
The received Spice2 schematic is simulated with its original device dimensions,
biases, 250 fF external input capacitance and 3 fF feedback capacitance.
LTspice 17.2.4 uses the default temperature of 27 degrees Celsius.
A derived model copy replaces NMOS NOIBE with NOIB and closes both model
parameter lists. All numerical coefficients are preserved.
The noise sweep is 1 Hz to 1 GHz, with 1000 points per decade. A second sweep
starts at 0.1 Hz to match the supplied screenshot. All results below use the
extended sweep. The supplied BSIM3 noise parameters are design inputs from
Rogelio; this reproduction does not establish their process calibration.

\begin{center}
\begin{tabular}{rrrr}\toprule
$f$ (Hz) & Total (mV/$\sqrt{\mathrm{Hz}}$) & M1 flicker ($\mu$V/$\sqrt{\mathrm{Hz}}$) & Flicker variance (\%)\\\midrule
"""+rows+r"""
\bottomrule\end{tabular}
\end{center}
Flicker contributions are nonzero for all seven transistors. The last column
includes all transistor flicker sources. M1 alone does not represent the
whole circuit. The input-noise trace is equivalent current noise; output
traces are amplitude spectral densities.

At 0.1 Hz, the total output noise is @TOTAL@ mV/$\sqrt{\mathrm{Hz}}$,
and M1 flicker is @M1@ $\mu$V/$\sqrt{\mathrm{Hz}}$.
These scales agree with the supplied screenshot (about 20 mV and 500 $\mu$V,
respectively). The received plot file contains presentation settings rather
than numerical data, so this comparison is limited by the screenshot.

The operating-point change is @OPERR@, and the maximum relative gain change is
@GAINERR@. The device variances sum to the output variance with maximum
relative error @CLOSURE@. The independently started noise sweeps agree at
common frequencies to @SWEEPERR@ relative error. Each device's individual
noise mechanisms also sum to its total variance within $10^{-10}$.
The solver reaches the operating point using Gmin stepping. It retains the
original warnings about unspecified diffusion perimeters and the isolated
charge-monitor node; no unknown model parameters are reported.

\begin{center}\includegraphics[width=.9\linewidth]{noise_comparison.pdf}\end{center}
\clearpage
\begin{center}\includegraphics[height=.91\textheight]{noise_reproduction.pdf}\end{center}
\end{document}
"""
values={'TOTAL':f"{numbers['output_ASD_V_sqrtHz'][0]*1e3:.4g}",'M1':f"{numbers['M1_flicker_ASD_V_sqrtHz'][0]*1e6:.4g}",
'OPERR':f"{checks['maximum_op_absolute_change']:.3g}",'GAINERR':f"{checks['maximum_gain_relative_change']:.3g}",
'CLOSURE':f'{closure:.3g}','SWEEPERR':f"{checks['overlapping_sweep_relative_difference']:.3g}"}
for key,val in values.items():text=text.replace('@'+key+'@',val)
(OUT/'report.tex').write_text(text)
print(values)
print({k:v for k,v in op.items() if k in ('V(out)','V(in)','V(s)','V(med)','V(vd6)')})
