"""Present the Julia/HDF5 scan without repeating the physical calculation."""
from pathlib import Path
import json
import os
import shutil
import h5py
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'notes/physics/generated'
REPORT=ROOT/'latex_generated/spice_averaging'
REPORT.mkdir(parents=True,exist_ok=True)
with h5py.File(ROOT/'results/ltspice/averaging_scan.h5') as h:
    widths=h['averaging_s'][:]*1e6; enc=h['enc_electrons_rms'][:]
    gains=h['gain_V_per_ion'][:]; noises=h['noise_rms_V'][:]
    errors={name:float(max(h[name][:])) for name in ('quadrature_relative_error','spectral_relative_error','transient_relative_error')}
    assumptions=h.attrs['assumptions']
selected=int(np.argmin(abs(widths-300)))
reference=json.loads((ROOT/'results/ltspice/capacitance/pad_1.25pF/study.json').read_text())['ions']
assert np.isclose(gains[selected],reference['signal_V_per_ion'],rtol=1e-9)
assert np.isclose(noises[selected],reference['output_noise_rms_V'],rtol=1e-9)
with plt.rc_context({'font.size':11,'axes.spines.top':False,'axes.spines.right':False}):
    fig,ax=plt.subplots(figsize=(6.7,3.8),layout='constrained')
    ax.plot(widths,enc,'o-',color='#245e8b',lw=1.5)
    ax.plot(widths[selected],enc[selected],'o',color='#b46618',ms=8,zorder=3)
    ax.annotate(f'300 µs: {enc[selected]:.1f} electrons RMS',xy=(widths[selected],enc[selected]),
                xytext=(235,enc.max()-.5),arrowprops={'arrowstyle':'->','color':'#b46618'},color='#925317')
    ax.set(xlabel='Averaging time (µs)',ylabel='ENC (electrons RMS)',xticks=widths)
    ax.grid(alpha=.18)
    fig.savefig(OUT/'spice_averaging.pdf',metadata={'CreationDate':None,'ModDate':None})
    preview=ROOT/'notes/physics/work/spice-paper/spice_averaging.png'
    preview.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(preview,dpi=170)
    plt.close(fig)
shutil.copyfile(OUT/'spice_averaging.pdf',REPORT/'spice_averaging.pdf')
values={'SpiceScanFirstEnc':f'{enc[0]:.1f}','SpiceScanLastEnc':f'{enc[-1]:.1f}'}
(OUT/'spice_averaging_values.tex').write_text('% Generated from results/ltspice/averaging_scan.h5\n'+''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k,v in values.items()))
display=lambda x, places: f'{x:.{places}f}'.rstrip('0').rstrip('.')
rows='\n'.join(f'{w:g} & {display(g*1e6,2)} & {display(s*1e3,3)} & {display(e,2)} '+r'\\' for w,g,s,e in zip(widths,gains,noises,enc))
report=r"""\documentclass[11pt]{article}
\usepackage[margin=24mm]{geometry}
\usepackage{graphicx,booktabs,amsmath}
\begin{document}
\section*{CSA averaging scan}
External pad capacitance: 1.25 pF; feedback: 3 fF; temperature: 293.15 K.
Uniform ion arrivals last 10 ms. Voltage averages are centered on the frame
boundaries; the discharge correction is $\alpha=\exp(-0.01/0.3)$.
Signal gain is recalculated for each averaging duration from the saved LTspice
transient. The stationary output noise is integrated from 0.1 Hz to 1 MHz
with the same averaging and frame subtraction. Intrinsic MOS capacitances
are included. The supplied model gives zero flicker contribution.
The input uses the central square pad in a grounded plane and a finite
5 mm approach, without charge renormalization.
\begin{center}
\includegraphics[width=.9\linewidth]{spice_averaging.pdf}
\begin{tabular}{rrrr}\toprule
Average ($\mu$s) & Gain ($\mu$V/ion) & Noise (mV RMS) & ENC ($e$ RMS)\\\midrule
""" + rows + r"""
\bottomrule\end{tabular}
\end{center}
The 300 $\mu$s point reproduces the independent saved Python calculation
within $10^{-9}$ relative tolerance. Checks compare quadrature orders 8 and 16,
noise spectra with 200 and 400 points per decade, and two transient resolutions.
The white-noise integral agrees with the non-overlapping box-average identity
to 0.2\% over the finite integration band.
""" + '\n\nMaximum relative differences: '+', '.join(f'{k.replace("_"," ")}: {v:.3g}' for k,v in errors.items()) + r""".
\end{document}
"""
(REPORT/'report.tex').write_text(report)
print('Generated averaging figure, value macros and standalone report source.')
