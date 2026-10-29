#!/usr/bin/env python3
"""Present the Julia/HDF5 ion-frame results for the revised flicker model."""

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
DATA=ROOT/'results/ltspice/spice2_enc'
OUT=ROOT/'latex_generated/spice2_enc';OUT.mkdir(parents=True,exist_ok=True)
PREVIEW=ROOT/'notes/physics/work/spice2-docs'
with h5py.File(DATA/'analysis.h5') as h:
    a=dict(h.attrs);w=h['averaging_s'][:]*1e6;enc=h['scan_enc_electrons_rms'][:]
    gain=h['scan_gain_V_per_ion'][:];noise=h['scan_noise_rms_V'][:]
    names=[x.decode().replace('V(','').replace(')','').upper() for x in h['device_names'][:]]
    fractions=h['device_variance_fraction'][:]*100
    counts=h['checked_counts'][:];snrs=h['checked_snr'][:]
    refinement=h['count_refinement_error'][:]
assert a['verified_threshold_ions']>0
with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
    fig,axes=plt.subplots(1,2,figsize=(10,3.9),layout='constrained')
    axes[0].plot(w,enc,'o-',color='#245e8b')
    j=int(np.argmin(abs(w-300)));axes[0].plot(w[j],enc[j],'o',color='#b46618',ms=8)
    axes[0].set(xlabel='Averaging time (µs)',ylabel='ENC (electrons RMS)',xticks=w)
    axes[0].grid(alpha=.2)
    axes[1].bar(names,fractions,color='#245e8b')
    axes[1].set(xlabel='Transistor',ylabel='Share of filtered noise variance (%)',ylim=(0,100))
    axes[1].grid(axis='y',alpha=.2)
    fig.savefig(OUT/'enc_and_sources.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(PREVIEW/'enc_and_sources.png',dpi=160);plt.close(fig)
fmt=lambda v,p=3:f'{v:.{p}f}'.rstrip('0').rstrip('.')
rows='\n'.join(f'{x:g} & {fmt(g*1e6)} & {fmt(s*1e3)} & {fmt(e,2)} '+r'\\' for x,g,s,e in zip(w,gain,noise,enc))
threshold=int(a['verified_threshold_ions']);bracket={int(n):float(s) for n,s in zip(counts,snrs)}
text=r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[margin=22mm]{geometry}
\usepackage{graphicx,booktabs,amsmath}
\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}
\begin{document}
\section*{NAUSICA ion-frame noise with the revised transistor model}
The calculation uses the supplied Spice2 transistor card with NOIMOD=2.
The derived copy corrects NMOS NOIBE to NOIB and closes both model parameter
lists. Device dimensions, biases and numerical model parameters are unchanged.
The coefficients are treated as inputs provided by Rogelio, without independent
process calibration.

The external pad capacitance is @PAD@ pF, with intrinsic MOS capacitances
added by SPICE. The feedback capacitance is @CF@ fF and the temperature is
@TEMP@ K. Uniform mean ion arrivals extend over @FRAME@ ms. The input is
the saved current for a 0.5 mm square pad in a grounded plane, integrated along
a finite 5 mm approach without charge renormalization.

Centered voltage averages of @AVERAGE@ $\mu$s are taken at the frame boundaries.
The corrected voltage difference is $D=-(\overline V_n-\alpha\overline V_{n-1})$,
with $\alpha=\exp(-\Delta t/\tau_0)$ and $\tau_0=@TAU@$ s.
The circuit's actual discharge is determined by the feedback transistor.
The one-sided stationary output noise is integrated from @LOW@ Hz to @HIGH@ MHz
with the same averaging and subtraction filter. ADC noise is excluded.

The signal calibration is @GAIN@ $\mu$V per ion, unchanged from the previous
model. The filtered noise is @SIGMA@ mV RMS, giving
\[
 \mathrm{ENC}=@ENC@\ e\ \mathrm{RMS}.
\]
The previous model gave @PREVIOUS@ electrons RMS. Flicker contributes
@FLICKER@ electrons RMS and thermal sources contribute @THERMAL@ electrons RMS,
combined in quadrature. Flicker accounts for @FRACTION@\% of the filtered
noise variance. M4 contributes @MFOUR@\% of the total variance.

Direct LTspice transients give $S/N=@BELOWSNR@$ for @BELOW@ ions and
$S/N=@ABOVESNR@$ for @ABOVE@ ions. Thus the verified integer threshold is
@ABOVE@ ions for $S/N>5$. The linear calibration predicts @LINEAR@ ions;
the direct test includes the small finite-amplitude change of response.

\begin{center}
\begin{tabular}{rrrr}\toprule
Average ($\mu$s) & Gain ($\mu$V/ion) & Noise (mV RMS) & ENC ($e$ RMS)\\\midrule
"""+rows+r"""
\bottomrule\end{tabular}
\end{center}
\begin{center}\includegraphics[width=\linewidth]{enc_and_sources.pdf}\end{center}

\subsection*{Validation}
Independent finite-path charge, halved transient steps, doubled spectral density
and quadrature orders 8 and 16 are checked. Relative errors are respectively
@CHARGE@, @TRANSIENT@, @SPECTRAL@ and @QUADRATURE@. The two adjacent integer
counts remain on opposite sides of $S/N=5$ at both transient resolutions;
the combined numerical refinement estimate is less than one tenth of their
distance from the boundary. The total noise variance agrees with the sum of
the transistor variances. Separately integrated flicker and thermal variances
close to the integrated total within @CLOSURE@ relative error.

The reproducible pipeline is \texttt{make spice2-enc}; HDF5 results, input hashes,
source waveforms, derived decks and solver logs are saved in the repository.
\end{document}
"""
values={'PAD':fmt(a['external_pad_F']*1e12),'CF':fmt(a['feedback_F']*1e15),'TEMP':fmt(a['temperature_K']),
'FRAME':fmt(a['frame_s']*1e3),'AVERAGE':fmt(a['averaging_s']*1e6),'TAU':fmt(a['correction_tau_s']),
'LOW':fmt(a['frequency_low_Hz']),'HIGH':fmt(a['frequency_high_Hz']*1e-6),
'GAIN':fmt(a['gain_V_per_ion']*1e6),'SIGMA':fmt(a['noise_rms_V']*1e3),'ENC':fmt(a['enc_electrons_rms'],2),
'PREVIOUS':fmt(a['previous_enc_electrons_rms'],2),'FLICKER':fmt(a['flicker_enc_electrons_rms'],2),
'THERMAL':fmt(a['thermal_enc_electrons_rms'],2),'FRACTION':fmt(a['flicker_variance_fraction']*100,2),
'MFOUR':fmt(fractions[names.index('M4')],2),'BELOW':str(threshold-1),'ABOVE':str(threshold),
'BELOWSNR':fmt(bracket[threshold-1],6),'ABOVESNR':fmt(bracket[threshold],6),'LINEAR':str(a['linear_threshold_ions']),
'CHARGE':f"{a['charge_error']:.3g}",'TRANSIENT':f"{max(a['transient_error'],max(refinement)):.3g}",
'SPECTRAL':f"{a['spectral_error']:.3g}",'QUADRATURE':f"{a['quadrature_error']:.3g}",'CLOSURE':f"{a['filtered_variance_closure']:.3g}"}
for key,value in values.items():text=text.replace('@'+key+'@',value)
(OUT/'report.tex').write_text(text)
print(values)
