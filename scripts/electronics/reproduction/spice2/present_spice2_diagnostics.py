#!/usr/bin/env python3
"""Generate the diagnostic report from Julia HDF5 results."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[4]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import h5py
ROOT=Path(__file__).resolve().parents[4]
def main():
    result=ROOT/'results/ltspice/spice2_diagnostics'
    report=ROOT/'latex_generated/spice2_diagnostics'; report.mkdir(parents=True,exist_ok=True)
    with h5py.File(result/'analysis.h5') as h:
        rows=[]; markdown=[]
        for tag in ('cap_0.75pF','cap_1.25pF','cap_3pF'):
            a=h[tag].attrs
            cap=a['external_capacitance_F']*1e12;enc=a['frame_enc_e']
            rows.append(f'{cap:g} & {enc:.1f}'+r' \\')
            markdown.append(f'| {cap:g} | {enc:.1f} |')
        g=h['cap_1.25pF'];a=g.attrs
        inject=g['injection_300us'].attrs
        names=g['device_names'].asstr()[:]; fractions=g['device_variance_fraction'][:]
        sources='\n'.join(f"{name[2:-1].upper()} & {100*fraction:.2f}"+r' \\' for name,fraction in zip(names,fractions))
        scan='\n'.join(f"{w} & {g[f'injection_{w}us'].attrs['enc_e']:.1f}"+r' \\' for w in (180,200,250,300,350))
        checks=max(float(h[tag].attrs[key]) for tag in ('cap_0.75pF','cap_1.25pF','cap_3pF') for key in ('transient_error','spectral_error','quadrature_error','charge_error'))
        text=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=25mm]{geometry}
\usepackage{amsmath,booktabs,siunitx,hyperref}
\title{CSA noise: capacitance and injection controls}
\author{NAUSICAA calculation record}
\date{\today}
\begin{document}\maketitle
The supplied 180 nm circuit is evaluated with the corrected Spice2 card,
including its NOIMOD=2 noise parameters. The external input capacitance is
varied; intrinsic MOS capacitances, transistor sizes and biases are retained.
The three external loads bracket the design assumptions requested by the author.
Topmetal-S reports approximately 3 pF for its total input capacitance.

\section*{Conditions}
The feedback capacitor is 3 fF and the temperature is 20 degrees Celsius.
Noise is stationary about the zero-signal operating point and integrated over
0.1 Hz--1 MHz. The ion current is the existing uniform mean arrival distribution
over 10 ms, for a 0.5 mm square pad in an infinite grounded plane and a finite
5 mm approach at 3 m/s. Its charge is $e\Delta\phi_w$, without renormalization.
The frame estimator subtracts centered 300 microsecond voltage averages at
the frame boundaries, with $\alpha=\exp(-10\,\mathrm{ms}/0.3\,\mathrm{s})$.
ENC is the filtered output RMS divided by the simulated one-ion signal gain.

\begin{center}\begin{tabular}{rr}\toprule
External capacitance (pF) & Frame ENC (electrons RMS)\\\midrule
@ROWS@
\bottomrule\end{tabular}\end{center}

\section*{Injection at 1.25 pF}
A pulse with a 1 microsecond plateau and 10 ns edges injects one elementary
charge at simulation time 5 ms. Frame windows are centered at 2 and 12 ms,
so collection occurs 3 ms into this frame.
The same 10 ms estimator gives an ENC of @SAME@ electrons RMS, compared
with @ION@ for the uniform ion signal.

The trapezoidal filter is a difference of two moving averages of duration $T$,
separated by $2T$. Its rise and flat-top durations are both $T$.
The signal is evaluated at the noiseless-template peak time. Noise is
integrated with exactly the same filter; no pole-zero correction is applied.
At $T=300$ microseconds the ENC is @INJ@ electrons RMS.

\begin{center}\begin{tabular}{rr}\toprule
Rise and flat top (microseconds) & Injection ENC (electrons RMS)\\\midrule
@SCAN@
\bottomrule\end{tabular}\end{center}

\section*{Noise source diagnosis}
At 1.25 pF, flicker contributes @FLICKER@\% of the frame-noise variance.
The contributions of the individual transistors to total variance are:
\begin{center}\begin{tabular}{rr}\toprule
Transistor & Variance fraction (percent)\\\midrule
@SOURCES@
\bottomrule\end{tabular}\end{center}
M4 has $W=0.72$ and $L=0.36$ micrometres in the supplied netlist.
This identifies a specific device for geometry and bias optimization.
The present calculation does not compare calibrated 130 and 180 nm process
cards, so it cannot attribute the difference from Topmetal-S to the process node.
The Topmetal-S measurement is approximately 115 electrons RMS after digital
trapezoidal filtering; its filter timing is not reported in that paper.

\section*{Validation and reproduction}
The analysis checks finite-path and injected charge, halves transient timesteps,
doubles the spectral density from 200 to 400 points per decade, increases
quadrature order from 8 to 16, and refines the injection peak-time grid.
The largest relative error among the frame checks is @CHECK@.
Each injection convergence check is below $10^{-3}$ and its quadrature
check is below $10^{-5}$. Integrated transistor variances close to the total
within $10^{-3}$. The 1.25 pF frame result reproduces the preceding calculation
to a relative tolerance of $10^{-9}$.

Run \texttt{scripts/electronics/reproduction/spice2/run\_spice2\_diagnostics.py} with \texttt{--simulator},
then \texttt{scripts/electronics/reproduction/spice2/analyze\_spice2\_diagnostics.jl}, and
\texttt{scripts/electronics/reproduction/spice2/present\_spice2\_diagnostics.py}, from the repository root.
Acquisition saves native decks and HDF5; Julia saves the numerical analysis;
Python generates this report from that analysis.

\begin{thebibliography}{9}
\bibitem{yang} Y. Yang et al., Design and preliminary test results of the
charge sensitive amplifier for gain-less charge readout in high-pressure TPC,
JINST 19 (2024) C03031, \url{https://doi.org/10.1088/1748-0221/19/03/C03031}.
\bibitem{bsim} W. Liu et al., BSIM3v3.2.2 MOSFET Model Users' Manual,
UC Berkeley (1999), Chapter 8.
\end{thebibliography}
\end{document}
'''
        for key,value in {'ROWS':'\n'.join(rows),'SAME':f"{a['injection_frame_enc_e']:.1f}",'ION':f"{a['frame_enc_e']:.1f}",'INJ':f"{inject['enc_e']:.1f}",'FLICKER':f"{100*a['flicker_variance_fraction']:.2f}",'SCAN':scan,'SOURCES':sources,'CHECK':f'{checks:.3g}'}.items():text=text.replace('@'+key+'@',value)
        (report/'report.tex').write_text(text)
        (result/'README.md').write_text('''# Revised-card capacitance and injection diagnostics

From the repository root, with the pinned Python environment and Julia project:

```sh
PYTHONPATH=src python scripts/electronics/reproduction/spice2/run_spice2_diagnostics.py --simulator /path/to/LTspice
julia --project=. scripts/electronics/reproduction/spice2/analyze_spice2_diagnostics.jl
python scripts/electronics/reproduction/spice2/present_spice2_diagnostics.py
```

Raw data and provenance: raw.h5. Numerical reduction: analysis.h5.
Generated standalone report: latex_generated/spice2_diagnostics/report.tex.
The report source is generated without compiling or changing the manuscript.

External capacitance plus intrinsic MOS capacitances; unchanged circuit and
revised NOIMOD=2 card; 20 C; 0.1 Hz–1 MHz stationary noise. Frame: uniform
finite-gap arrivals over 10 ms, 300 us centered averages, 0.3 s nominal
discharge correction. Injection: one elementary charge over 1 us plus 10 ns
edges; trapezoid rise and flat top both 300 us in the main comparison.

| External capacitance (pF) | Frame ENC (electrons RMS) |
|---:|---:|
'''+ '\n'.join(markdown)+f"\n\nAt 1.25 pF, injection with the frame estimator: {a['injection_frame_enc_e']:.1f} e RMS; with the 300 us trapezoid: {inject['enc_e']:.1f} e RMS.\n")
if __name__=='__main__':main()
