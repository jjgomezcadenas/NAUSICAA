#!/usr/bin/env python3
"""Generate the 4 mm / 40 ms exercise report using saved Julia results."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import hashlib,json,os,tomllib
import h5py
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    c=tomllib.loads((ROOT/'config/pitch_4mm_40ms.toml').read_text());source=ROOT/c['output']/'analysis.h5';out=ROOT/c['report'];out.mkdir(parents=True,exist_ok=True)
    with h5py.File(source) as h:
        a=dict(h.attrs);old=dict(h['reference_20ms'].attrs)
        t=h['time_s'][:];v=h['voltage_V'][:];i=h['current_A'][:]
        targets=h['snr_targets'][:];ions=h['linear_ions_at_target'][:]
        counts=h['counts'][:];snr=h['snr'][:];ratios=h['gain_relative_to_small_signal'][:]
        names=h['device_names'].asstr()[:];fractions=h['device_variance_fraction'][:]
    with h5py.File(ROOT/c['reference']) as h:
        ot=h['time_s'][:];ov=h['voltage_V_per_ion'][:]*c['reference_mean_ions'];oi=h['current_A_per_ion'][:]*c['reference_mean_ions']
    with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
        fig,axes=plt.subplots(1,2,figsize=(9,3.6),layout='constrained')
        for time,current,voltage,label in ((ot,oi,ov,'2 mm, 20 ms, 120 ions'),(t,i,v,'4 mm, 40 ms, 400 ions')):
            axes[0].plot(time*1e3,current*1e15,label=label)
            axes[1].plot(time*1e3,-voltage*1e3,label=label)
        axes[0].set(xlabel='Time since first collection (ms)',ylabel='Input current (fA)',xlim=(-2,50))
        axes[1].set(xlabel='Time since first collection (ms)',ylabel='CSA voltage magnitude relative to baseline (mV)',xlim=(-2,100))
        for ax in axes:ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.2)
        fig.savefig(out/'response_comparison.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(out/'response_comparison.png',dpi=160);plt.close(fig)
    fields=[('Mean ions per bucket','mean_ions',1,'.0f'),('Gain ($\\mu$V/ion)','gain_V_per_ion',1e6,'.3f'),('Readout noise (mV RMS)','noise_rms_V',1e3,'.3f'),('ENC ($e$ RMS)','ENC_e_rms',1,'.2f'),('Mean signal (mV)','mean_signal_V',1e3,'.3f'),('Mean signal-to-noise ratio','mean_snr',1,'.3f'),('Fitted discharge time (s)','tau_s',1,'.5f')]
    rows='\n'.join(label+' & '+format(old[key]*factor,fmt)+' & '+format(a[key]*factor,fmt)+r' \\' for label,key,factor,fmt in fields)
    threshold='\n'.join(f'{target:g} & {value:.1f}'+r' \\' for target,value in zip(targets,ions))
    checks='\n'.join(f'{key.replace("_"," ")} & {a[key]:.3g}'+r' \\' for key in ('transient_error','spectral_error','quadrature_error','tail_fit_relative_error','tau_timestep_error','noise_variance_closure_error'))
    devices='\n'.join(f'{name} & {100*fraction:.1f}'+r' \\' for name,fraction in zip(names,fractions))
    text=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=22mm]{geometry}
\usepackage{amsmath,booktabs,graphicx,siunitx,hyperref}
\title{NAUSICA: 4 mm pitch and 40 ms buckets}
\author{Circuit exercise}\date{\today}
\begin{document}\maketitle
The 4 mm pitch case assumes a mean of 400 ions per pad per 40 ms bucket,
provided by J. J. G\'omez-Cadenas as a design input. The comparison uses the
updated mean of 120 ions for 2 mm pitch and 20 ms buckets.

\section{Circuit and input}
The exercise retains the paper's 180 nm transistor circuit, 1.8 V supply,
20 degree Celsius temperature, external capacitance 1.25 pF and feedback
capacitance 3 fF. The collecting pad remains a 0.5 mm square. M3 has
$W=200\,\mu$m and $L=2\,\mu$m; M4 has $W=2.88\,\mu$m and $L=1.44\,\mu$m.
The archived gate settings and complete model card are recorded in the netlists.
Thermal and flicker noise are evaluated at the DC operating point.

The input uses uniform mean collection times over 40 ms. Each ion follows a
central straight trajectory at 3 m/s through the final 5 mm gap above the pad.
The weighting problem has unit potential on the pad, zero on the surrounding
infinite grounded plane and at infinity, in a homogeneous gas half-space.
The integrated signal is $e[1-\phi_w(g)]$, retained without renormalization.
The 4 mm pitch affects the specified charge population and spatial sampling;
the calculation retains the same local pad benchmark and circuit capacitance.
An electrostatic focusing calculation and event redistribution are separate
inputs to this circuit exercise.

For the 400-ion packet the integrated input charge is @CHARGE@ fC and the
peak current is @CURRENT@ fA. Collection begins at simulation time 2 ms and
ends at 42 ms; induction starts before collection.

\section{Readout and results}
At frame boundaries the voltage is averaged over centred 300 $\mu$s windows.
After subtraction of a matched zero-input baseline, the measured difference is
\[
D_n=-[\overline V_n-\alpha\overline V_{n-1}],\qquad
\alpha=e^{-\Delta t/\tau_0}.
\]
The decay time is fitted from the 100-ion packet between 41 and 120 ms after
first collection. The output-noise spectrum is the same as the paper's circuit;
its filtered integral is recalculated for the 40 ms separation:
\[
\sigma_D^2=\int_{0.1\,\mathrm{Hz}}^{1\,\mathrm{MHz}}S_V(f)
\left|\frac{\sin(\pi fT_a)}{\pi fT_a}\right|^2
\left|1-\alpha e^{-2\pi i f\Delta t}\right|^2\,df.
\]
The small-signal ENC is $\sigma_D/g$. The 400-ion signal and S/N are obtained
from a finite-charge SPICE transient. The 20 ms mean signal uses the saved
small-signal calibration multiplied by 120.

\begin{center}\begin{tabular}{lrr}\toprule
Quantity & 2 mm / 20 ms & 4 mm / 40 ms\\\midrule
@ROWS@
\bottomrule\end{tabular}\end{center}
\includegraphics[width=\linewidth]{response_comparison.pdf}

The following populations correspond to illustrative S/N values in the linear
model. Reconstruction can combine information from neighbouring measurements;
these values specify individual-bucket significance.
\begin{center}\begin{tabular}{rr}\toprule
S/N & Ions per 40 ms bucket\\\midrule
@THRESHOLDS@
\bottomrule\end{tabular}\end{center}
The nominal spatial sampling changes from 2 to 4 mm transversely and from 2 to
4 mm along drift, using 0.1 m/s bulk ion speed. The assumed mean population
must be evaluated together with track and blob preservation at that sampling.

\section{Noise contributions and validation}
\begin{center}\begin{tabular}{lr}\toprule
Device & Fraction of filtered noise variance (\%)\\\midrule
@DEVICES@
\bottomrule\end{tabular}\end{center}
The calculation compares 5 and 2.5 $\mu$s transient steps, 200 and 400 frequency
points per decade, and quadrature orders 8 and 16. It checks charge conservation,
independent analytical-current quadrature, 10- and 100-ion discharge fits,
local fit windows, and linearity through 400 ions. Twenty saved-result checks pass.
\begin{center}\begin{tabular}{lr}\toprule
Check & Relative error\\\midrule
@CHECKS@
\bottomrule\end{tabular}\end{center}

\section{Reproduction}
From the repository root, with the pinned Python and Julia dependencies:
\begin{verbatim}
PYTHONPATH=src python scripts/electronics/readout/run_pitch_exercise.py \
  --julia /path/to/julia --simulator /path/to/LTspice
\end{verbatim}
Julia performs the input calculation, reduction and checks. HDF5 stores inputs,
raw transients, the reused noise spectrum, numerical results and provenance.
Python acquires LTspice output and generates this report and figure. Presentation
alone is regenerated with \texttt{python scripts/electronics/readout/present\_pitch\_exercise.py}.
\end{document}
'''
    for key,value in {'ROWS':rows,'THRESHOLDS':threshold,'CHECKS':checks,'DEVICES':devices,'CHARGE':f"{a['mean_input_charge_C']*1e15:.5f}",'CURRENT':f"{a['mean_peak_current_A']*1e15:.4f}"}.items():text=text.replace('@'+key+'@',value)
    (out/'report.tex').write_text(text)
    (out/'provenance.json').write_text(json.dumps({'source':str(source.relative_to(ROOT)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'generator':str(Path(__file__).relative_to(ROOT))},indent=2)+'\n')
    (ROOT/c['output']/'README.md').write_text('''# 4 mm pitch / 40 ms circuit exercise

User design inputs: 4 mm pitch and a mean of 400 ions per bucket; comparison
uses 2 mm / 20 ms with the updated mean of 120 ions. Pad geometry, external
capacitance, circuit and operating point stay those of the paper. The local
square-pad finite-gap benchmark and uniform mean arrivals are retained; no
charge renormalization. Noise is stationary at the DC operating point.

Reproduce from the repository root:

```sh
PYTHONPATH=src python scripts/electronics/readout/run_pitch_exercise.py --julia /path/to/julia --simulator /path/to/LTspice
```

Julia input generation and reduction: src/ion_signal/pitch_exercise.jl.
HDF5 inputs, raw acquisitions and analysis retain hashes and configuration.
The existing noise spectrum is reused only after its circuit/model hashes match.
Both transient resolutions, current PWLs, netlists and readable logs are saved.
Native LTspice raw/log scratch is ignored; raw HDF5 retains the numerical data.
Report source and figure: latex_generated/pitch_4mm_40ms/.
No manuscript values or previous 20 ms results are changed.
''')
    print('Report:',out/'report.tex')
if __name__=='__main__':main()
