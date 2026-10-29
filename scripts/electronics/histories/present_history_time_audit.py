#!/usr/bin/env python3
"""Present the saved Julia time-axis audit without repeating its calculations."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('event',nargs='?',type=int,default=3790);args=p.parse_args()
 out=ROOT/'latex_generated/history_time_audit'/f'event_{args.event}';out.mkdir(parents=True,exist_ok=True)
 with h5py.File(ROOT/'results/ltspice/history_time_audit'/f'event_{args.event}.h5') as h:
  g=h['electrons'];a=g.attrs
  fig,axes=plt.subplots(2,1,figsize=(8,5.8))
  axes[0].scatter(g['source_hit_z_mm'][:],g['source_hit_energy_MeV'][:]*1000,s=8)
  axes[0].set(xlabel='Source-hit z (mm)',ylabel='Deposited energy per hit (keV)')
  axes[1].stairs(g['collected_ions_per_frame'][:],g['time_edges_ms'][:]/1000)
  axes[1].set(xlabel='Relative event arrival time (s)',ylabel='Collected ions per 20 ms frame')
  for ax in axes:ax.grid(alpha=.2)
  fig.tight_layout();fig.savefig(out/'time_axis.pdf');fig.savefig(out/'time_axis.png',dpi=160);plt.close(fig)
  body=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=25mm]{geometry}
\usepackage{graphicx,amsmath}
\title{Event EVENT: origin and resolution of the arrival-time axis}
\date{\today}
\begin{document}\maketitle
For the single-electron sample, the saved event axis contains BUCKETS buckets
and spans 0--DURATION s. It is a relative arrival-time axis,
\[
 t=(z-z_{\mathrm{first}})/v,\qquad v=SPEED\ \mathrm{m/s}.
\]
The first longitudinal slice of the complete retained event defines zero.
The common absolute drift delay is not included. Individual pads use this same
axis, so charge reaching a late-arriving pad can appear after a long initial
interval of zero counts.

The source MC deposits extend from ZMIN to ZMAX mm along z; the generation
vertex is at VERTEX mm. Deposits separated by approximately one metre therefore
have approximately ten seconds of separation at this speed. FINAL\% of the
collected ions fall in the final second of the saved event interval. Earlier
charge belongs to other positions and pads. Only OCCUPIED of the BUCKETS
event buckets contain collected charge anywhere in the event.

The saved diffusion table has longitudinal cells WIDTH mm wide. The 20 ms
readout corresponds to 2 mm. Assigning the coarser diffusion-cell centres to
readout buckets can introduce artificial gaps and concentrate charge into
individual buckets. A finer diffusion-grid convergence check is required before
treating the derived event compression fractions as final. The circuit
calculation remains a response to the stored input counts.

\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{time_axis.pdf}
\caption{Source-hit longitudinal positions and the full collected-ion time
sequence for the single-electron event. The timing plot sums all pads.}
\end{figure}

The same numerical event ID in the double-beta sample is a different event;
its saved duration is BBDURATION s. Input arrays, grid widths, fractions and
source hashes are saved in \texttt{results/ltspice/history\_time\_audit/}.
Reproduce from the repository root:
\begin{verbatim}
julia --project=. scripts/electronics/histories/audit_history_time.jl EVENT
python scripts/electronics/histories/present_history_time_audit.py EVENT
\end{verbatim}
The report source is generated without compiling a PDF.
\end{document}
'''
  values={'EVENT':str(args.event),'BUCKETS':str(a['original_buckets']),'DURATION':f"{a['duration_s']:.2f}",'SPEED':f"{a['bulk_speed_m_s']:.1f}",'ZMIN':f"{a['source_min_z_mm']:.1f}",'ZMAX':f"{a['source_max_z_mm']:.1f}",'VERTEX':f"{a['generation_vertex_z_mm']:.1f}",'FINAL':f"{100*a['fraction_in_final_second']:.1f}",'OCCUPIED':str(a['occupied_event_buckets']),'WIDTH':f"{a['diffusion_z_cell_width_mm']:.2f}",'BBDURATION':f"{h['bb0nu'].attrs['duration_s']:.2f}"}
  # Replace the longer token first, since it contains DURATION.
  for key in sorted(values,key=len,reverse=True):body=body.replace(key,values[key])
  (out/'report.tex').write_text(body)
 print(out/'report.tex')
if __name__=='__main__':main()
