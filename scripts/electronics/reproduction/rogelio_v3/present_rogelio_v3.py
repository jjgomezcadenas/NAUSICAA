#!/usr/bin/env python3
"""Present the v3 input-charge and pulse-amplitude comparison."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[4]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import hashlib,json,os
import h5py
ROOT=Path(__file__).resolve().parents[4]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def main():
 source=ROOT/'results/ltspice/rogelio_v3/analysis.h5';out=ROOT/'latex_generated/rogelio_v3';out.mkdir(parents=True,exist_ok=True)
 with h5py.File(source) as h:
  a=dict(h.attrs);ex=dict(h['exp_r2'].attrs);un=dict(h['uniform_1000_r2'].attrs)
  with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
   fig,axes=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
   for name,label,shift in [('exp_r2','Supplied exponential pulse',.2),('uniform_1000_r2','1000 ions, uniform arrivals',.002)]:
    g=h[name];t=(g['time_s'][:]-shift)*1e3
    axes[0].plot(t,g['current_A'][:]*1e15,label=label)
    axes[1].plot(t,g['voltage_drop_V'][:]*1e3,label=label)
   for ax in axes:ax.set(xlim=(-2,30),xlabel='Time relative to pulse/frame start (ms)');ax.grid(alpha=.2)
   axes[0].set_ylabel('Input current (fA)');axes[1].set_ylabel('Output voltage drop (mV)');axes[1].legend(fontsize=8)
   fig.savefig(out/'pulse_comparison.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
  text=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}
\usepackage[margin=22mm]{geometry}\usepackage{amsmath,booktabs,graphicx,hyperref}
\title{Rogelio v3: input charge and the 14 mV signal}
\author{NAUSICA calculation record}\date{\today}
\begin{document}\maketitle
The supplied v3 pulse produces an output drop of @EXPPEAK@ mV. Its
integrated input charge is @EXPQ@ fC, equivalent to @EXPE@ elementary
charges. The paper's 1000-ion current applied to exactly the same circuit
produces @UNPEAK@ mV. The approximately 14 mV result is therefore consistent
with the paper: the supplied test injects about half as much charge.

\section*{Reproduction}
The circuit is read from
\texttt{LTSpice/Spice2/v3/CSAfolded180nm-2\_Astronio.asc}, with its supplied
model card. It uses the paper's transistor dimensions, 3 fF feedback,
1.25 pF external input capacitance and $20\,^{\circ}$C. Its gate voltages
are 0.53847 V for M3 and 0.556803 V for M4. The output baseline is
@BASE@ V. The source is
\begin{center}\texttt{EXP(0 48f 200m 0.05n 0 1.7m)}.\end{center}
The simulated current begins at 0.2 s, reaches 48 fA and decays with a
fitted time constant of @TAU@ ms. Its integral is consistent with
$48\,\mathrm{fA}\times1.7\,\mathrm{ms}=0.0816\,\mathrm{fC}$.
The schematic annotation ``0.5 ke = 0.8 fC'' has a factor-ten error;
500 elementary charges correspond to approximately 0.08 fC.

\begin{center}\begin{tabular}{lrrr}\toprule
Input & Charge (fC) & Charge ($e$) & Output drop (mV)\\\midrule
Supplied exponential & @EXPQ@ & @EXPE@ & @EXPPEAK@\\
1000 ions, uniform over 20 ms & @UNQ@ & @UNE@ & @UNPEAK@\\\bottomrule
\end{tabular}\end{center}
The uniform-ion current retains the original finite-path normalization.
It uses the central 0.5 mm square pad in a grounded plane, homogeneous gas,
and a central 5 mm approach at 3 m/s. The weighting potential is one on
the pad, zero on the remaining plane and zero at infinity. Each ion
contributes $e\Delta\phi_w$ rather than being renormalized to $e$.

\begin{figure}[htbp]\centering
\includegraphics[width=\linewidth]{pulse_comparison.pdf}
\caption{Both input waveforms simulated with the supplied v3 circuit.}
\end{figure}

\section*{Charge integration and convergence}
The calculation uses maximum timesteps of 5 and 2.5 $\mu$s, tight solver
tolerances and disabled waveform compression. Integrating the saved input
current agrees with the independent behavioral \texttt{idt(I(Is1))}
accumulator. Step refinement changes the charge by @QERR@\% and the
output peak by @VERR@\%.
The original 40 s run reproduces the approximately 14 mV pulse, but its
compressed current trace is unsuitable for charge integration. The supplied
\texttt{.meas} result must be checked with resolved, uncompressed current data.

\section*{Resistance label}
The schematic's resistance label evaluates \texttt{V(out)/Id(M4)}. At this
operating point it gives @STATIC@ M$\Omega$. This is a DC ratio involving
the output node voltage. M4's incremental drain resistance is
$r_o=1/g_{ds}$, which the operating-point log gives as approximately
@RO@ M$\Omega$. These are different quantities; the displayed ratio does
not establish the amplifier's small-signal voltage gain.

\section*{Reproduction commands}
Acquisition is \texttt{scripts/electronics/reproduction/rogelio\_v3/run\_rogelio\_v3.py}; Julia analysis is
\texttt{scripts/electronics/reproduction/rogelio\_v3/analyze\_rogelio\_v3.jl}; presentation is
\texttt{scripts/electronics/reproduction/rogelio\_v3/present\_rogelio\_v3.py}. Commands and saved data are in
\texttt{results/ltspice/rogelio\_v3/README.md}. The manuscript is unchanged.
\end{document}
'''
  values={'EXPPEAK':f"{ex['peak_voltage_V']*1e3:.4f}",'EXPQ':f"{ex['charge_C']*1e15:.7f}",'EXPE':f"{ex['equivalent_electrons']:.2f}",
          'UNPEAK':f"{un['peak_voltage_V']*1e3:.4f}",'UNQ':f"{un['charge_C']*1e15:.7f}",'UNE':f"{un['equivalent_electrons']:.2f}",
          'BASE':f"{a['baseline_V']:.6f}",'TAU':f"{ex['fitted_current_tau_s']*1e3:g}",
          'QERR':f"{ex['charge_refinement_error']*100:.5f}",'VERR':f"{ex['voltage_refinement_error']*100:.4f}",
          'STATIC':f"{a['M4_display_Vout_over_Id_Ohm']*1e-6:.3f}",'RO':f"{a['M4_incremental_ro_Ohm']*1e-6:.1f}"}
  for k,v in values.items():text=text.replace('@'+k+'@',v)
  (out/'report.tex').write_text(text)
 (out/'provenance.json').write_text(json.dumps({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,Path(__file__))},indent=2)+'\n')
if __name__=='__main__':main()
