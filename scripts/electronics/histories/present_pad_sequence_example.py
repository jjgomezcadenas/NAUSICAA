#!/usr/bin/env python3
"""Render the saved ten-bucket SPICE calculation as a table and figure."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import tomllib
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
def main():
 c=tomllib.loads((ROOT/'config/pad_sequence_example.toml').read_text());out=ROOT/c['report'];out.mkdir(parents=True,exist_ok=True)
 with h5py.File(ROOT/c['output']/'analysis.h5') as h:
  g=h['requested'];a=h.attrs;edges=g['edges_s'][:];counts=g['counts'][:]
  rows=[]
  for j,n in enumerate(counts):
   rows.append(f"{j+1} & {n} & {g['output_at_edges_V'][j+1]:.3f} & {g['linear_estimated_ions'][j]:.1f} & {g['estimated_ions'][j]:.1f} "+r'\\')
  plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
  fig,axes=plt.subplots(3,1,figsize=(8,7),sharex=True)
  axes[0].stairs(counts,edges,color='0.3',fill=True,alpha=.5)
  axes[0].set_ylabel('Collected ions')
  axes[1].plot(g['time_s'][:],g['output_V'][:],label='Transistor circuit')
  axes[1].plot(edges,g['linear_output_at_edges_V'][:],'--',label='Linear prediction')
  axes[1].set_ylabel('CSA output (V)');axes[1].legend(frameon=False)
  axes[2].stairs(g['estimated_ions'][:],edges,label='SPICE estimate')
  axes[2].stairs(g['linear_estimated_ions'][:],edges,ls='--',label='Linear estimate')
  axes[2].set_ylabel('Estimated charge ($e$)');axes[2].set_xlabel('Time (s)')
  axes[2].legend(frameon=False)
  for ax in axes:ax.set_xlim(0,edges[-1]);ax.grid(alpha=.2)
  fig.tight_layout();fig.savefig(out/'sequence.pdf',bbox_inches='tight');fig.savefig(out/'sequence.png',dpi=160,bbox_inches='tight');plt.close(fig)
  body=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=23mm]{geometry}
\usepackage{amsmath,booktabs,graphicx,siunitx,hyperref}
\title{A ten-bucket example through the optimized CSA}
\date{\today}
\begin{document}\maketitle
The user-provided sequence is 200, 200, 100, 400, 10, 10, 0, 5000, 200, 0 ions.
Buckets last FRAME ms. Arrivals are uniform within each bucket, with the paper's
finite-gap current for a 0.5 mm square pad in a grounded plane, 5 mm gap and
constant speed 3 m/s. Its integrated charge is $e[1-\phi_w(g)]$ per ion.
The optimized transistor circuit starts at its quiescent baseline BASELINE V;
its low-charge fitted discharge time is TAU s. Five-ms centred averages at
frame boundaries give
\[
 D_n=-[\Delta\overline V_n-\alpha\Delta\overline V_{n-1}],\qquad
 \alpha=\exp(-\Delta t/\tau),\qquad \widehat N_n=D_n/g.
\]
The isolated uniform-packet calibration gives $g=GAIN$ $\mu$V per ion.
All following signals are mean noiseless responses. The saved quiescent
noise reference is ENC elementary charges RMS per corrected frame.

\begin{table}[ht]\centering
\begin{tabular}{rrrrr}\toprule
Bucket & Collected ions & End voltage (V) & Linear estimate ($e$) & SPICE estimate ($e$)\\\midrule
ROWS
\bottomrule\end{tabular}
\caption{Raw voltage at each boundary and charge inferred from the centred
averages with the fixed low-charge gain. Both estimates retain the response
spreading between adjacent buckets.}\end{table}

For this sequence the circuit remains close to the linear prediction.
The 5000-ion deposit does not drive it into the strongly compressed regime seen
in the extreme event histories. Earlier charge remains on the CSA as each
bucket arrives; the output therefore decreases across the occupied buckets and
recovers gradually after collection.

Bucket 7 has zero collected ions but a positive estimated signal from ions
collected in bucket 8. Induction begins before collection, and the centred
averaging window also extends across the boundary. Bucket 9 includes signal
from the large deposit as well as its own 200 ions. These deterministic adjacent
responses must be included when converting frame estimates to collected counts.

Paired controls remove the bucket-9 deposit while keeping the rest of the
history fixed. The marginal estimated contribution of its 200 ions is PROBEWITH
with the large deposit present. Removing the 5000-ion deposit and repeating the
same paired comparison gives PROBEWITHOUT. The relative reduction of the
marginal response is REDUCTION\%.

\begin{figure}[ht]\centering\includegraphics[width=.95\linewidth]{sequence.pdf}
\caption{The complete sequence retains the accumulated CSA state. The linear
and transistor-level predictions nearly coincide at this charge scale.}
\end{figure}

\paragraph{Validation and reproduction.}
Four histories (requested sequence, removal of bucket 9, removal of bucket 8,
and removal of both) are simulated at two resolutions with matched zero-input
references. Eight PWL charge checks and eight transient charge/refinement
checks pass. Configuration, code, model, simulator and netlist hashes are saved
with HDF5 outputs in \texttt{results/ltspice/pad\_sequence\_example/}.
Julia performs current generation and reduction; Python acquires SPICE output
and presents the HDF5 results. From the repository root:
\begin{verbatim}
PYTHONPATH=src python scripts/electronics/histories/run_pad_sequence_example.py \
  --simulator /path/to/LTspice
python scripts/electronics/histories/present_pad_sequence_example.py
\end{verbatim}
The report source is generated without compiling the PDF.
\end{document}
'''
  replacements={'FRAME':f"{a['frame_s']*1000:.0f}",'BASELINE':f"{a['baseline_V']:.3f}",'TAU':f"{a['tau_s']:.2f}",'GAIN':f"{a['gain_V_per_ion']*1e6:.2f}",'ENC':f"{a['ENC_e']:.1f}",'ROWS':'\n'.join(rows),'PROBEWITHOUT':f"{a['probe_estimated_without_large_e']:.1f}",'PROBEWITH':f"{a['probe_estimated_with_large_e']:.1f}",'REDUCTION':f"{-100*a['probe_relative_gain_change']:.2f}"}
  for key,value in replacements.items():body=body.replace(key,value)
  (out/'report.tex').write_text(body)
 print(out/'report.tex')
if __name__=='__main__':main()
