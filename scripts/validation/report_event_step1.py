"""Render the step-one validation HDF5 as a standalone LaTeX checkpoint.

Julia supplies every measured value. This script formats tables and preserves
the scope/settings with the report. Compile with latexmk in the output folder.
"""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import h5py

ROOT = Path(__file__).resolve().parents[2]


def render():
    output = ROOT / "latex_generated/event_analysis_step1"
    output.mkdir(parents=True, exist_ok=True)
    with h5py.File(ROOT / "results/event_analysis/step1_validation.h5") as f:
        tests = "\n".join(
            f"{name} & {count} \\\\" for name, count in
            zip(f["tests/suite"].asstr()[:], f["tests/passed"][:])
        )
        rows = []
        for name, g in f["samples"].items():
            label = "Double beta" if name == "0nubb" else "Electron background"
            rows.append(
                f"{label} & {g['events'][()]} & {g['hits'][()]} & "
                f"{g['bytes'][()]/1e6:.2f} & "
                f"{g['max_energy_relative_error'][()]:.3g} & "
                f"{g['mean_indexed_read_ms'][()]:.2f} \\\\"
            )
        text = r"""\documentclass[11pt]{article}
\usepackage[margin=20mm]{geometry}
\usepackage{booktabs,hyperref}
\begin{document}
\title{Event analysis: step-one checkpoint}
\author{NAUSICAA}
\date{}
\maketitle
The implemented checkpoint covers the versioned hit reader and Petit's
deterministic binning, graph, backbone and truth-path functions. Physical arrays
and arithmetic use Float32 in mm, MeV and ns. The reduced v3 files retain MC time
and actual primary IDs; their values are verified against the original sources.

\section*{Validation}
\begin{center}\begin{tabular}{lr}\toprule
Julia suite & Passed assertions \\\midrule
TESTS
\bottomrule\end{tabular}\end{center}
The pinned Petit fixtures cover line, curved, branched, dense, lattice and
separated geometries. Tie-free fixtures require matching graph edges and paths
(up to reversal), with numerical tolerances for Float32 lengths and scores.
Tied cases check valid neighbour choices, component sizes and energies, and
path validity. The maximum accepted edge distance is explicit in every case.

\begin{center}\small\begin{tabular}{lrrrrr}\toprule
Sample & Events & Hits & MB & Max relative error & Read ms \\\midrule
SAMPLES
\bottomrule\end{tabular}\end{center}
The energy error compares Float32 source-order event sums against independent
Float64 sums over all events. The detector's reference relative energy resolution
is 0.005. Read time is the mean of five indexed event reads after warming the
reader, with the file already open; it depends on the machine and file cache.

Five events spread across each sample also exercise the deterministic stages
on raw hits: 2 mm bins, kNN with k=10, 4 mm maximum distance and zero energy
threshold. These are software checks. Scientific selection efficiencies require
the subsequent ion-production and diffusion implementation. The approved
production reconstruction prescription is b=3 sigma and d=6 sigma at every L.

\section*{Reproduction and next checkpoint}
\begin{verbatim}
make validate-event-step1 PYTHON=/path/to/python3.11-or-newer
# Presentation only, from the existing HDF5:
make report-event-step1 PYTHON=/path/to/python3.11-or-newer
\end{verbatim}
Numerical values and source hashes are saved in
\texttt{results/event\_analysis/step1\_validation.h5}.
The HTML review guide contains the implemented interfaces, docstrings and
remaining design contracts. Implementation stops at this checkpoint.
\end{document}
""".replace("TESTS", tests).replace("SAMPLES", "\n".join(rows))
    (output / "report.tex").write_text(text)
    print(output / "report.tex")


if __name__ == "__main__":
    render()
