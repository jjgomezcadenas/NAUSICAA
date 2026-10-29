#!/usr/bin/env python3
"""Present saved Julia pad-history calculations; no circuit or event calculation."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import hashlib
import argparse
import tomllib
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LABELS = {'p2_bb0nu': r'2 mm, 20 ms, $0\nu\beta\beta$',
          'p2_electrons': '2 mm, 20 ms, electrons',
          'p4_bb0nu': r'4 mm, 40 ms, $0\nu\beta\beta$',
          'p4_electrons': '4 mm, 40 ms, electrons'}

def decode(value):
    return value.decode() if isinstance(value, bytes) else str(value)

def tex(value):
    return decode(value).replace('_', r'\_').replace('%', r'\%')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path,default=ROOT/'config/pad_history_readout_direct.toml');args=parser.parse_args()
    config = tomllib.loads(args.config.read_text())
    data = ROOT/config['output']
    out = ROOT/config['report']
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'savefig.bbox': 'tight'})
    with h5py.File(data/'summary.h5') as summary, h5py.File(data/'spice_analysis.h5') as results:
        names = [s['name'] for s in config['samples']]
        rows = []
        for name in names:
            a = summary[name].attrs
            rows.append(f"{LABELS[name]} & {a['events']:,} & {a['pads']:,} & "
                        f"{a['pads_below_guard']:,} & {100*a['pad_guard_fraction']:.4f} & "
                        f"{a['events_with_guard_risk']:,} & {100*a['event_guard_fraction']:.2f} "+r'\\')
        ensemble_table = '\n'.join(rows)
        rows = []
        for name, group in summary['spice'].items():
            a = group.attrs
            label = {'largest_bucket':'largest bucket', 'near_guard':'near guard',
                     'worst_voltage':'lowest linear voltage'}[decode(a['label'])]
            short = decode(a['sample']).replace('electrons','e').replace('bb0nu',r'$\beta\beta$').replace('_',' ')
            rows.append(f"{short}, {label} & {a['event_id']} & {a['minimum_output_V']:.3f} & "
                        f"{a['recovery_after_collection_s']:.2f} & "
                        f"{a['post_collection_positive_peak_e']:.1f} & "
                        f"{a['post_collection_negative_peak_e']:.1f} "+r'\\')
        spice_table = '\n'.join(rows)
        rows = []
        for name, group in summary['adc'].items():
            for i, bits in enumerate(group['bits']):
                if bits not in (10,12,13):
                    continue
                rows.append(f"{name.replace('ms','')} & {group['span_V'][i]:.0f} & {bits} & "
                            f"{1e6*group['lsb_V'][i]:.0f} & {group['quantization_ENC_e'][i]:.2f} & "
                            f"{100*group['total_noise_relative_increase'][i]:.3f} "+r'\\')
        adc_table = '\n'.join(rows)
        hashes = '\n'.join(r'\noindent '+tex(summary[name].attrs['source'])+r'\\'+
                           r'{\ttfamily\scriptsize '+decode(summary[name].attrs['source_sha256'])+r'}\par'
                           for name in names)
        for name in ('p4_electrons_worst_voltage','p4_bb0nu_largest_bucket'):
            g = results[name]
            counts = g['counts'][:]
            edges = g['edge_time_s'][:]
            positive = np.flatnonzero(counts)
            origin = edges[positive[0]]
            stop = edges[positive[-1]+1]-origin+2.5
            time = g['time_s'][:]-origin
            gain = summary[decode(g.attrs['sample'])].attrs['gain_V_per_ion']
            fig, axes = plt.subplots(2, 1, figsize=(8.2,6), sharex=True,
                                     gridspec_kw={'height_ratios':[1,1.1]})
            axes[0].plot(time,g['output_V'][:],label='Transistor circuit',lw=1.8)
            axes[0].plot(edges-origin,g['linear_output_at_edges_V'][:],
                         '--',label='Linear packet sum',lw=1.3)
            axes[0].axhline(config['range_guard_low_V'],color='0.4',ls=':',
                           label='0.2 V diagnostic guard')
            axes[0].set_ylabel('CSA output (V)')
            axes[0].legend(frameon=False,fontsize=9,loc='lower right')
            width = g.attrs['frame_s']
            axes[1].bar(edges[:-1]-origin,counts,width=width,align='edge',
                        color='0.85',label='Collected ions')
            axes[1].stairs(g['linear_frame_signal_V'][:]/gain,edges-origin,
                           color='C1',ls='--',lw=1.1,label='Linear frame estimate')
            axes[1].stairs(g['frame_signal_V'][:]/gain,edges-origin,
                           color='C0',lw=1.3,label='SPICE frame estimate')
            axes[1].set_ylabel('Ions per frame / estimated charge ($e$)')
            axes[1].set_xlabel('Time since first occupied bucket (s)')
            axes[1].legend(frameon=False,fontsize=9,loc='upper right')
            axes[1].set_xlim(-.02,stop)
            axes[0].set_title(LABELS[decode(g.attrs['sample'])]+f", event {g.attrs['event_id']}, "
                              f"pad ({g.attrs['pixel_x']}, {g.attrs['pixel_y']})")
            for ax in axes:
                ax.grid(alpha=.18)
            fig.tight_layout()
            fig.savefig(out/(name+'.pdf'))
            fig.savefig(out/(name+'.png'),dpi=160)
            plt.close(fig)
        fig, ax = plt.subplots(figsize=(7.2,4))
        for name in names:
            g=summary[name]
            ax.plot(g['voltage_threshold_V'][:],g['pad_fraction_below_threshold'][:],label=LABELS[name])
        ax.axvline(config['range_guard_low_V'],color='0.5',ls=':')
        ax.set(xlim=(-.2,.7),ylim=(1e-6,1),yscale='log',
               xlabel='Threshold on minimum linear output (V)',
               ylabel='Fraction of pad histories below threshold')
        ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.2)
        fig.tight_layout();fig.savefig(out/'linear_screen.pdf');fig.savefig(out/'linear_screen.png',dpi=160)
        plt.close(fig)
        a=summary[names[0]].attrs
        refinement=max(g.attrs['transient_refinement_error'] for g in summary['spice'].values())
        charge_error=max(g.attrs['charge_error'] for g in summary['spice'].values())
        text = r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=22mm]{geometry}
\usepackage{amsmath,graphicx,booktabs,siunitx,hyperref}
\usepackage[T1]{fontenc}
\title{Complete pad histories with the optimized NAUSICAA CSA}
\date{\today}
\begin{document}
\maketitle
\section{Inputs and calculation}
The stored uncut ion counts provide a complete time sequence for each pad.
Missing sparse entries are filled with zero. Four files contain 8000 events each:
the same double-beta and electron samples at 2 mm pitch with 20 ms buckets and
4 mm pitch with 40 ms buckets. Each event starts with the CSA at its quiescent
operating point. Empty buckets and all large deposits are retained.

The circuit is the selected optimized 180 nm CSA: $C_f=3$ fF, external input
capacitance 1.25 pF, $M_f$ width 0.36 $\mu$m and length 400 $\mu$m, grounded
feedback gate, and $M_4$ gate at 0.6164 V. The output baseline is BASELINE V
and the fitted small-signal discharge time is TAU s.
Voltages are averaged over 5 ms centred at each frame boundary. The reported
frame signal is
\begin{equation}
 D_n=-[\Delta\overline V_n-\alpha\Delta\overline V_{n-1}],\qquad
 \alpha=\exp(-\Delta t/\tau),\qquad \Delta V=V-V_{\rm baseline}.
\end{equation}
Its low-charge gain is calibrated with the uniform-arrival SPICE packet.

Arrivals are uniform within each stored bucket. Each ion follows the on-axis,
constant-speed final approach used by the paper: a 0.5 mm square pad in a grounded
plane, 5 mm gap and speed 3 m/s. The pad weighting potential has unit pad and zero
surrounding-plane boundary conditions. This is the planar benchmark for the
finite-gap induced current. The input integral remains
$e[1-\phi_w(g)]$ per ion. No collection impulse or charge renormalization is added.

Julia first propagates every history using the calibrated small-signal packet
response and an exponential late tail. A minimum output of 0.2 V is used as a
diagnostic guard to select demanding cases. It is not a calibrated boundary of
linear operation. Three histories per input file are then simulated with the
full transistor circuit: the lowest predicted voltage, largest single bucket,
and the history closest to the guard.

\section{All-history linear screen}
\begin{table}[ht]
\centering\scriptsize
\begin{tabular}{lrrrrrr}\toprule
Sample & Events & Pads & Guard pads & Pads (\%) & Guard events & Events (\%)\\\midrule
ENSEMBLE
\bottomrule\end{tabular}
\caption{Histories for which the linear output falls below the diagnostic 0.2 V
guard. These counts identify histories requiring nonlinear simulation.}
\end{table}
\begin{figure}[ht]\centering
\includegraphics[width=.9\linewidth]{linear_screen.pdf}
\caption{Cumulative fraction of pad histories below a chosen minimum linear
output. Predictions below the supply range indicate failure of the linear approximation.}
\end{figure}

\section{Transistor-level histories}
Matched zero-input references and two time resolutions are used for all twelve
histories. Current transitions use 5 and 2.5 $\mu$s knot spacing; the solver maximum
steps are 0.5 and 0.25 ms. Each history includes 20 s of recovery after the stored
event interval. The maximum change of the frame signal under refinement,
normalized to its peak, is REFINEMENT. The maximum relative input-charge error is
CHARGEERROR. The Julia accounting and zero-input checks pass (32 assertions),
as do the SPICE charge, refinement and voltage-bound checks (48 assertions).

\begin{table}[ht]\centering\scriptsize
\begin{tabular}{lrrrrr}\toprule
History & Event & $V_{\min}$ (V) & Recovery (s) & Positive ($e$) & Negative ($e$)\\\midrule
SPICETABLE
\bottomrule\end{tabular}
\caption{Actual output minima and recovery above the diagnostic guard, timed
from the last collection. The last two columns give extrema of the difference
between the nonlinear and linear frame estimates in empty buckets, starting one
full frame after the last collection. They use the low-charge gain to express
voltage in units of elementary charge.}
\end{table}

The extreme histories compress the output. A voltage minimum above the ground
rail therefore does not establish a linear response. In the selected cases the
output can remain below the guard after collection. Later positive residuals in
empty buckets remain smaller than one nominal quiescent ENC. The largest
negative residuals occur in the first evaluated empty bucket. Small residuals
in empty buckets do not establish that the gain for a new charge packet has
recovered while the output remains displaced.

Saturation flags should use a calibrated limit on CSA compression. ADC overflow
alone will not flag compressed signals that remain within its input range.
The event reconstruction can retain flagged deposits as high-charge locations;
its topology performance remains to be evaluated with these flags.

\begin{figure}[p]\centering
\includegraphics[width=\linewidth]{p4_electrons_worst_voltage.pdf}
\caption{A demanding electron pad history. The upper panel compares the full
transistor simulation with linear accumulation. The lower panel compares collected
counts and the charge inferred using the fixed low-charge gain.}
\end{figure}
\begin{figure}[p]\centering
\includegraphics[width=\linewidth]{p4_bb0nu_largest_bucket.pdf}
\caption{The largest double-beta bucket at 4 mm pitch. Empty buckets retain the
circuit state and expose the recovery after the large deposit.}
\end{figure}

\section{ADC resolution}
\begin{table}[ht]\centering\small
\begin{tabular}{rrrrrr}\toprule
Frame (ms) & Span (V) & Bits & LSB ($\mu$V) & Quantization ENC ($e$) & RMS increase (\%)\\\midrule
ADCTABLE
\bottomrule\end{tabular}
\caption{Quantization contribution for independent uniform errors in the two
averaged voltage readings: $\sigma_{D,q}=\mathrm{LSB}\sqrt{(1+\alpha^2)/12}$.
The last column combines this contribution with the saved quiescent analogue noise.}
\end{table}
The selected SPICE histories fit within a 0--1 V input range. Twelve bits add
little noise even for a 2 V span. This resolution estimate assumes independent
uniform quantization errors on the averaged readings. It includes no reduction
from ADC oversampling and no additional ADC or buffer noise. The nominal ENC is
calculated at the quiescent operating point; noise at strongly displaced
operating points has not been evaluated in this study.

\section{Reproduction and provenance}
Physical propagation, PWL generation and numerical reduction are implemented in
\texttt{src/ion\_signal/pad\_history\_readout.jl} and
\texttt{src/ion\_signal/pad\_history\_spice.jl}. Python acquires LTspice output and
renders the saved HDF5 results. From the repository root, run
\begin{verbatim}
python scripts/electronics/histories/run_pad_history_pipeline.py --simulator /path/to/LTspice
\end{verbatim}
Existing matching inputs and raw results are reused. Presentation alone is
\begin{verbatim}
python scripts/electronics/histories/present_pad_histories.py
\end{verbatim}
The report source is generated without compiling the PDF. The numerical files
are in \texttt{results/ltspice/pad\_history\_readout/}. Configuration, code, source,
model, simulator and netlist hashes are stored with the results.

Input-file SHA256 identifiers:
\small
HASHES
\end{document}
'''
        replacements = {'BASELINE':f"{a['baseline_V']:.3f}", 'TAU':f"{a['tau_s']:.2f}",
                        'ENSEMBLE':ensemble_table,'SPICETABLE':spice_table,'ADCTABLE':adc_table,
                        'REFINEMENT':rf"\num{{{refinement:.3e}}}",
                        'CHARGEERROR':rf"\num{{{charge_error:.3e}}}", 'HASHES':hashes}
        for key,value in replacements.items():
            text=text.replace(key,value)
        (out/'report.tex').write_text(text)
        manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in (data/'summary.h5',data/'spice_analysis.h5',Path(__file__))}
        (out/'sources.txt').write_text('\n'.join(f'{v}  {k}' for k,v in manifest.items())+'\n')
    print(out/'report.tex')

if __name__ == '__main__':
    main()
