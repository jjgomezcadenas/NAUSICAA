#!/usr/bin/env python3
"""Present a saved six-stream propagation summary without repeating diffusion."""
from pathlib import Path
import argparse
import subprocess
import h5py

ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    with h5py.File(args.summary) as f:
        workers = int(f["workers"][()])
        events = int(f["events_per_stream"][()])
        wall = float(f["wall_seconds"][()])
        samples = [x.decode() if isinstance(x, bytes) else x for x in f["sample"][()]]
        carriers = [x.decode() if isinstance(x, bytes) else x for x in f["carrier"][()]]
        times = f["stream_seconds"][()].tolist()
        sizes = f["output_bytes"][()].tolist()
        memory = f["peak_rss_bytes"][()].tolist()
    target = (args.output or ROOT / "latex_generated" / f"three_clouds_{events}").resolve()
    stem = f"three_clouds_{events}"
    rows = []
    for sample, carrier, seconds, size, rss in zip(samples, carriers, times, sizes, memory):
        display_carrier = {"ion": "ion", "pure_xe": "pure Xe electron",
                           "xe_he": "Xe--15\\% He electron"}[carrier]
        rows.append(f"{sample} & {display_carrier} & {seconds:.2f} & "
                    f"{size/1e6:.1f} & {rss/1e9:.2f}" + r" \\")
    tex = r"""\documentclass{article}
\usepackage[margin=23mm]{geometry}
\usepackage{booktabs}
\title{Six-stream cloud propagation: """ + str(events) + r"""-event sample}
\author{NAUSICAA analysis}
\date{29 October 2026}
\begin{document}
\maketitle
\section*{Method}
For each of two MC source samples, """ + str(events) + r""" events were propagated into three
HDF5 clouds: positive ions, pure-Xe electrons and Xe--15\%~He electrons.
One event location was sampled uniformly over the 0--150~cm drift length.
Ions and electrons used complementary distances. The Julia launcher used
six independent processes and separate output files. Timing includes Julia
startup, source reads, per-carrier diffusion, histogramming and HDF5 writes.

\section*{Run results}
""" + f"Workers: {workers}. Events per stream: {events}. " + \
          f"Total wall time: {wall:.2f}~s. Total cloud storage: {sum(sizes)/1e6:.1f}~MB.\n" + r"""
\begin{center}
\begin{tabular}{llrrr}
\toprule
Sample & Cloud & Time (s) & HDF5 (MB) & Peak RSS (GB) \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{center}
""" + f"The sum of the six individual peak resident-memory measurements is " \
          f"{sum(memory)/1e9:.2f}~GB; it bounds the simultaneous peak from above.\n" + r"""
\section*{Output and checks}
Each cloud file contains occupied Float32 histogram-bin positions in mm,
Int64 carrier counts, explicit bin edges, common event location, both drift
distances, evaluated diffusion widths, source identity and random seeds.
The validation checks each event's three matching locations and produced
carrier totals, the complementary drift distances, and the charge ledger
between sampled, retained and out-of-range carriers.
\end{document}
"""
    target.mkdir(parents=True, exist_ok=True)
    (target / f"{stem}.tex").write_text(tex)
    with (target / "latexmk.stdout").open("w") as log:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
                        f"{stem}.tex"], cwd=target, check=True,
                       stdout=log, stderr=subprocess.STDOUT)
    build_log = (target / f"{stem}.log").read_text()
    assert "! LaTeX Error" not in build_log


if __name__ == "__main__":
    main()
