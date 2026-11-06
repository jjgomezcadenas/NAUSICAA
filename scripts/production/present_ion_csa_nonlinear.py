#!/usr/bin/env python3
"""Present saved Julia nonlinear-CSA validation as figures, HTML and compiled LaTeX."""
from pathlib import Path
import argparse
import html
import json
import os
import sys
import tomllib
import h5py
import numpy as np
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[2]/'work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/documentation'))
from build_report import build_report

def text(dataset):
    value=dataset[()]
    return value.decode() if isinstance(value,bytes) else str(value)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=ROOT/'results/ion_csa_nonlinear/v1')
    parser.add_argument('--output',type=Path,default=ROOT/'latex_generated/ion_readout/ion_csa_nonlinear')
    args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    saved=tomllib.loads((args.results/'validation.toml').read_text())
    with h5py.File(args.results/'calibration.h5') as cal,h5py.File(args.results/'report.h5') as report:
        assert cal['passed'][()] and report['passed'][()]
        s={k:report['summary'][k][()] for k in report['summary']}
        frames={k:report['frames'][k][()] for k in report['frames']}
        pilot={k:report['pilot'][k][()] for k in report['pilot']}
        manifest=Path(text(report['pilot_manifest']))
        study=next((manifest.parent/'execution/pitch4/filtered').glob('*/acquisition_summary.json')).parent
        timing=json.loads((study/'acquisition_summary.json').read_text())
        nl=tomllib.loads((study/'nonlinear_summary.toml').read_text())
        with h5py.File(study/'nonlinear_control_validation.h5') as control:
            control_error=control['maximum_error_ENC'][()]
        valid=cal['validation/valid'][()].astype(bool)
        errors=cal['validation/error_ENC'][()][valid]
        minimum=cal['validation/minv'][()][valid]
        widths=cal['validation/frame_ms'][()][valid]
        knots=cal['knots_V'][()];gain=cal['gain_V_per_e'][()];recovery=cal['recovery_V_s'][()]
        bias=cal['events/charge_difference_ions'][()]
        fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
        ax[0,0].plot(knots,gain*1e6,color='#176589');ax[0,0].set(xlabel='Voltage excursion (V)',ylabel='Gain (µV / ion)')
        twin=ax[0,0].twinx();twin.plot(knots,recovery*1e3,color='#b65b32');twin.set_ylabel('Recovery (mV / s)',color='#b65b32')
        for ms in sorted(set(widths)):
            mask=widths==ms
            ax[0,1].scatter(minimum[mask],errors[mask],s=4,alpha=.5,label=f'{ms} ms')
        ax[0,1].set(xlabel='Linear predicted minimum output (V)',ylabel='Maximum corrected bucket error (ENC)',yscale='log')
        ax[0,1].axhline(.5,color='#b65b32',ls='--',label='0.5 ENC limit');ax[0,1].legend(fontsize=8)
        ax[1,0].hist(errors,bins=50,color='#176589');ax[1,0].set(xlabel='Maximum corrected bucket error (ENC)',ylabel='Held-out pad histories')
        ax[1,1].hist(bias,bins=50,color='#176589');ax[1,1].set(xlabel='Event charge difference (ions)',ylabel='Held-out event groups')
        for a in ax.flat:a.grid(alpha=.18)
        fig.savefig(out/'validation.png',dpi=170);plt.close(fig)
        fig,ax=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
        ax[0].hist(pilot['nonlinear_charge_difference_ions'],bins=25,color='#176589')
        ax[0].set(xlabel='Event difference from nonlinear replacement (ions)',ylabel='Pilot events')
        ax[1].hist(pilot['maximum_pulse_difference_ions'],bins=25,color='#b65b32')
        ax[1].set(xlabel='Event maximum-pulse difference (ions)',ylabel='Pilot events')
        fig.savefig(out/'pilot.png',dpi=170);plt.close(fig)
        rows=[]
        for j,ms in enumerate(frames['frame_ms']):
            rows.append([str(ms),str(frames['histories'][j]),str(frames['fallbacks'][j]),
                f"{frames['maximum_error_ENC'][j]:.5f}",f"{frames['median_error_ENC'][j]:.5f}",
                f"{1000*frames['median_seconds'][j]:.2f}"])
        macros={
            'TrainingCount':s['training_histories'],'ValidationCount':saved['validation_histories'],
            'AcceptedCount':saved['valid_histories'],'FallbackCount':saved['fallback_histories'],
            'WorstError':f"{saved['maximum_error_ENC']:.5f}",'WorstCharge':f"{s['maximum_event_charge_difference_ions']:.3f}",
            'RefinementError':f"{s['maximum_refinement_ENC']:.7f}",'ModelMilliseconds':f"{1000*saved['median_model_seconds']:.2f}",
            'PilotCount':s['pilot_events'],'PilotFastPads':s['pilot_nonlinear_pads'],'PilotSpicePads':s['pilot_spice_pads'],
            'PilotLinearPads':s['pilot_linear_pads'],'PilotFastError':f"{s['pilot_nonlinear_maximum_error_ENC']:.5f}",
            'PilotFastCharge':f"{s['pilot_nonlinear_maximum_charge_difference_ions']:.3f}",
            'PilotAllError':f"{s['pilot_maximum_error_ENC']:.4f}",'PilotAllCharge':f"{s['pilot_maximum_charge_difference_ions']:.3f}",
            'PilotMeanCharge':f"{s['pilot_mean_charge_difference_ions']:.3f}",
            'PilotMaxPulse':f"{s['pilot_maximum_pulse_difference_ions']:.3f}",
            'NonlinearSeconds':f"{nl['elapsed_s']:.2f}",'WorkflowSeconds':f"{timing['elapsed_s']:.1f}",
            'DomainExcursion':f"{s['maximum_excursion_V']:.6f}",'DomainCurrent':f"{s['maximum_current_A']*1e15:.3f}",
        }
        (out/'values.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+str(v)+'}' for k,v in macros.items())+'\n')
        table='\n'.join(' & '.join(row)+r' \\' for row in rows)
        (out/'frames.tex').write_text(r'\begin{tabular}{rrrrrr}\toprule'+'\n'+
            r'Frame (ms)&Histories&Fallbacks&Max (ENC)&Median (ENC)&Runtime (ms)\\\midrule'+'\n'+
            table+'\n'+r'\bottomrule\end{tabular}'+'\n')
        prose=(f"Validated on {saved['validation_histories']:,} held-out pad histories across 20, 30 and 40 ms. "
               f"{saved['valid_histories']:,} are inside the calibrated domain; {saved['fallback_histories']} requires SPICE. "
               f"Worst corrected-bucket error: {macros['WorstError']} ENC (limit 0.5 ENC). "
               f"Maximum event-charge difference: {macros['WorstCharge']} ions.")
        pilotprose=(f"The {s['pilot_events']}-event 4 mm pilot produced {s['pilot_linear_pads']:,} linear, "
                    f"{s['pilot_nonlinear_pads']:,} fast nonlinear and {s['pilot_spice_pads']:,} SPICE-control pad histories. "
                    f"Nonlinear replacement alone differs from saved SPICE by at most {macros['PilotFastError']} ENC per bucket "
                    f"and {macros['PilotFastCharge']} ions per event. "
                    f"The 607 nonlinear evaluations, including control predictions, took {macros['NonlinearSeconds']} s. "
                    f"The filtered workflow took {macros['WorkflowSeconds']} s after preparation, including Julia startup, controls, refinement and reduction.")
        distinction=(f"The 100-event pilot selects a different set of above-guard SPICE controls from the original 8,000-event run. "
                     f"Including that control-selection change gives a maximum bucket difference of {macros['PilotAllError']} ENC "
                     f"and a maximum event-total difference of {macros['PilotAllCharge']} ions "
                     f"(mean {macros['PilotMeanCharge']} ions). The event maximum-pulse difference is {macros['PilotMaxPulse']} ions. "
                     "The model-only comparison isolates the nonlinear replacement; the held-out calibration comparison uses identical histories.")
        headings=['Bucket (ms)','Histories','Fallbacks','Worst error (ENC)','Median error (ENC)','Median model (ms)']
        htmltable='<table><tr>'+''.join('<th>'+x+'</th>' for x in headings)+'</tr>'+''.join('<tr>'+''.join('<td>'+x+'</td>' for x in row)+'</tr>' for row in rows)+'</table>'
        page='''<!doctype html><html lang="en"><meta charset="utf-8"><title>ITACA nonlinear CSA validation</title>
<style>body{font:17px system-ui;line-height:1.55;max-width:1100px;margin:40px auto;padding:0 20px;color:#193c4a}img{width:100%}table{border-collapse:collapse}th,td{padding:9px;border:1px solid #ccd8df}pre{padding:16px;background:#edf4f6;overflow:auto}.status{padding:18px;background:#e6f3ed}</style>
<h1>ITACA: fast nonlinear CSA response</h1><p class="status">Validated. New production defaults to <code>hybrid_nonlinear</code>; <code>hybrid_spice</code> remains selectable.</p>'''
        page+=f'<p>{prose}</p>{htmltable}<img src="validation.png" alt="Model tables and independent validation"><h2>Production pilot</h2><p>{pilotprose}</p><p>{distinction}</p><img src="pilot.png" alt="End-to-end pilot differences">'
        page+='''<h2>Model and boundary conditions</h2><pre>x = baseline_V - Vout
dx/dt = gain(x) * I(t)/e - recovery(x)
x(0) = 0</pre><p>Continuous piecewise-linear gain and recovery tables, fitted to 240 SPICE histories with disjoint validation events. Propagation starts at equilibrium and follows the complete finite-path input current and quiet tail. RK4 stops at input knots and acquisition-window edges; maximum step is 100 µs. The existing 5 ms averaging window, three-frame inverse and stationary temporal noise covariance are retained.</p>'''
        page+=f'<p>The validated excursion extends to {macros["DomainExcursion"]} V and input current to {macros["DomainCurrent"]} fA. The same 0.63 V selection guard applies. Controls use SPICE; histories outside the calibrated domain use SPICE. The largest timestep-refinement difference was {macros["RefinementError"]} ENC. The geometry/current prescription and circuit hashes are checked before execution.</p>'
        page+='''<p>Production HDF5 retains signed estimates and full noise covariance. Model codes: 0 linear, 1 SPICE, 2 calibrated nonlinear. Run identity includes calibration contents. Identical noise is used when comparing model effects; independent production identities obtain their own deterministic noise realizations.</p>
<h2>Run and reproduce</h2><pre>make ion-csa-nonlinear-pilot
make ion-csa-nonlinear-report</pre><p>Set <code>response_model = "hybrid_spice"</code> in the production run card for the reference mode. Numerical evidence lives under <code>results/ion_csa_nonlinear/v1/</code>. <a href="ion_csa_nonlinear.pdf">PDF report</a> · <a href="../../../docs/ion_csa_nonlinear_model.html">Interfaces and development model</a></p></html>'''
        (out/'index.html').write_text(page)
        tex=r'''\documentclass[11pt]{article}
\usepackage[a4paper,margin=22mm]{geometry}
\usepackage{amsmath,graphicx,booktabs,hyperref}
\input{values.tex}
\title{ITACA: fast nonlinear CSA response}
\author{Generated from Julia validation results}\date{2026-11-05}
\begin{document}\maketitle
\section{Model and domain}
The state $x=V_0-V_{\mathrm{out}}$ follows
\[\dot{x}=g(x)I(t)/e-r(x),\qquad x(0)=0.\]
The gain $g$ (V per elementary charge) and recovery $r$ (V/s) are continuous
piecewise-linear tables fitted to \TrainingCount{} saved SPICE histories.
Recovery vanishes at equilibrium. Each prediction follows the complete current
history from baseline, including the quiet tail. RK4 steps end at current knots
and acquisition-window boundaries, with a maximum step of $100\,\mu$s.
The calibrated excursion extends to \DomainExcursion{} V and current to
\DomainCurrent{} fA, for the existing $0.5$ mm square pad, $5$ mm final gap,
$3$ m/s gap velocity and uniform mean arrivals within each frame.
The existing 5 ms voltage average, three-frame inverse and stationary
Gaussian temporal noise covariance are retained. The guard remains 0.63 V.
Controls and out-of-domain histories use SPICE.
\section{Independent validation}
Calibration and validation use disjoint event IDs. Of \ValidationCount{}
held-out histories, \AcceptedCount{} lie inside the calibration domain and
\FallbackCount{} requires SPICE. The worst corrected-bucket error is
\WorstError{} ENC against the 0.5 ENC criterion. The largest event-total
difference is \WorstCharge{} ions. Halving the integration step changes the
corrected estimates by at most \RefinementError{} ENC. Median model runtime
is \ModelMilliseconds{} ms per history on this machine.
\begin{center}\small\input{frames.tex}\end{center}
\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{validation.png}
\caption{Fitted response and held-out validation. Event totals replace only
below-guard histories inside the domain; unchanged linear and fallback SPICE
histories cancel.}\end{figure}
\clearpage\section{End-to-end production pilot}
The \PilotCount{}-event 4 mm pilot has \PilotLinearPads{} linear,
\PilotFastPads{} fast nonlinear and \PilotSpicePads{} SPICE-control histories.
The nonlinear replacement alone differs from saved SPICE by at most
\PilotFastError{} ENC per bucket and \PilotFastCharge{} ions per event.
The nonlinear evaluation stage, including control predictions, took
\NonlinearSeconds{} s. The filtered workflow took \WorkflowSeconds{} s after
preparation, including startup, SPICE controls, solver refinement and reduction.
These timings do not include input preparation or the later correction and
publication stages; they are not a full-production speedup measurement.

The smaller pilot selects different above-guard SPICE controls from the
original 8,000-event run. Including that selection change gives a maximum
bucket difference of \PilotAllError{} ENC, maximum event-charge difference
\PilotAllCharge{} ions and mean difference \PilotMeanCharge{} ions.
The maximum-pulse difference is \PilotMaxPulse{} ions. The model-only metric
isolates the nonlinear replacement. A common noise draw isolates response
differences; new production identities otherwise have independent deterministic
noise draws. Signed estimates, quiet buckets, covariance and ion ledgers pass
the existing production checks.
\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{pilot.png}
\caption{Model-only charge differences and complete-pipeline maximum-pulse differences.}\end{figure}
\section{Production and provenance}
New production defaults to \texttt{hybrid\_nonlinear}. An explicit
\texttt{response\_model = "hybrid\_spice"} selects the reference backend.
HDF5 model codes are 0 linear, 1 SPICE and 2 calibrated nonlinear.
Calibration content, circuit and model code are checked and hashed into
the resolved run identity. Numerical results and source hashes are in
\path{results/ion_csa_nonlinear/v1/}. Reproduce the pilot with
\texttt{make ion-csa-nonlinear-pilot}, and regenerate this report with
\texttt{make ion-csa-nonlinear-report}. The numerical calibration is a local
HDF5 artifact regenerated by the calibration script and its TOML card.
\end{document}
'''
        (out/'ion_csa_nonlinear.tex').write_text(tex)
    build_report(out/'ion_csa_nonlinear.tex')
    print(out/'index.html')
    print(out/'ion_csa_nonlinear.pdf')

if __name__=='__main__':main()
