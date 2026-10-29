#!/usr/bin/env python3
"""Render the saved Julia cloud benchmark to a self-contained LaTeX report."""
from pathlib import Path
import subprocess
import h5py

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results/three_cloud_benchmark/benchmark.h5"
DEST = ROOT / "latex_generated/three_cloud_benchmark"


def main():
    with h5py.File(DATA) as f:
        workers = f["workers"][()].tolist()
        seconds = f["wall_seconds"][()].tolist()
        outputs = f["output_bytes"][()].tolist()
        rss = f["sum_worker_peak_rss_bytes"][()].tolist()
        events = f["events_per_stream"][()].tolist()
    baseline = seconds[workers.index(1)]
    lines = []
    for w, n, t, b in zip(workers, events, seconds, outputs):
        lines.append(f"{w} & {n} & {t:.2f} & {baseline/t:.2f} & {b/1e6:.1f}" + r" \\")
    peak = rss[workers.index(6)] / 1e9 if 6 in workers else None
    memory = (f"The sum of individual peak resident-memory measurements for six "
              f"processes is {peak:.2f} GB. This sum is a conservative ceiling "
              "for simultaneous resident memory, not a measured concurrent peak."
              if peak is not None and peak >= 0 else "Peak memory was not measured.")
    tex = r"""\documentclass{article}
\usepackage[margin=25mm]{geometry}
\usepackage{booktabs}
\usepackage{amsmath}
\title{Three-cloud propagation prototype: concurrency benchmark}
\author{NAUSICAA analysis}
\date{29 October 2026}
\begin{document}
\maketitle
\section*{Method}
For each of 100 events from each 1~mm source sample (double beta and single
electrons), the Julia prototype produced an ion cloud and electron clouds using
the pure-Xe and Xe--15\%~He Petit diffusion coefficients. The event location was
drawn once per source event; ion and electron distances sum to 150~cm. The
benchmark includes Julia process startup, HDF5 source reads, diffusion and
compressed HDF5 output. Each stream writes a separate file.

\section*{Results}
\begin{center}
\begin{tabular}{rrrrr}
\toprule
Workers & Events/stream & Wall time (s) & Speedup & Output (MB) \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{center}
""" + memory + r"""

\section*{Validation and limits}
All 100 events in each source sample have the same sampled ion drift distance
across their three clouds. Ion and electron drift distances sum to 150~cm.
Every cloud conserves its sampled carriers across retained and outside-bin
counts. Per-event arrays and counts match exactly for one, two, four and six
workers. These measurements cover 100 events per stream and include output
to the local filesystem; production-scale throughput and concurrent peak
resident memory require a separate full-run measurement.
\end{document}
"""
    DEST.mkdir(parents=True, exist_ok=True)
    (DEST / "three_cloud_benchmark.tex").write_text(tex)
    subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
                    "three_cloud_benchmark.tex"], cwd=DEST, check=True,
                   stdout=(DEST / "latexmk.stdout").open("w"), stderr=subprocess.STDOUT)
    log = (DEST / "three_cloud_benchmark.log").read_text()
    assert "! LaTeX Error" not in log


if __name__ == "__main__":
    main()
