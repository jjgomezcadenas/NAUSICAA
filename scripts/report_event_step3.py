"""Render saved selection/pad histories as final figures, tables and LaTeX.

Julia computes every count, projection and sample statistic. Python reads the
checkpoint HDF5, selects stored examples and formats presentation artifacts.
"""
from pathlib import Path
import json
import tomllib
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "latex_generated/event_analysis_step3"


def render():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    with h5py.File(ROOT / "results/event_analysis/step3_validation.h5") as f:
        config_text = f["configuration/card"].asstr()[()]
        config = tomllib.loads(config_text)
        summary = {"passed": int(f["tests/passed"][()]), "samples": [],
                   "drift_direction": config["pads"]["drift_direction"]}
        rows, figures = [], []
        for name, sg in f["samples"].items():
            label = "Double beta" if name == "0nubb" else "Electron background"
            stats = {key: value[()].item() for key, value in sg["summary"].items()}
            summary["samples"].append({"sample": label, **stats})
            rows.append(f"{label} & {stats['accepted_events']}/{stats['tested_events']} & "
                        f"{stats['accepted_retained_ions']} & {stats['pad_ions']} & "
                        f"{stats['accepted_below_graph_threshold_ions']} & {stats['mean_seconds']*1000:.1f} \\\\")
            if "example" not in sg:
                continue
            example = sg["example"]
            event = int(example["event_id"][()])
            history = sg[f"events/{event}/pad_history"]
            xedges, yedges = history["x_edges_mm"][:], history["y_edges_mm"][:]
            time = history["time_edges_ms"][:]
            # Julia stores x-first matrices; h5py reads the physical HDF5 y,x order.
            xy = example["xy_ions"][:]
            assert xy.shape == (len(yedges)-1, len(xedges)-1)
            fig = plt.figure(figsize=(11, 7), layout="constrained")
            grid = fig.add_gridspec(2, 2)
            axmap = fig.add_subplot(grid[:, 0])
            mapdata = np.ma.masked_where(xy == 0, xy)
            m = axmap.pcolormesh(xedges, yedges, mapdata, shading="flat", cmap="viridis",
                                norm=LogNorm(vmin=1, vmax=max(2, int(xy.max()))), rasterized=True)
            axmap.set(xlabel="Source x [mm]", ylabel="Source y [mm]", title="Full retained charge per pad")
            axmap.set_aspect("equal", adjustable="box")
            fig.colorbar(m, ax=axmap, label="Ions / pad", shrink=0.7)
            axz = fig.add_subplot(grid[0, 1])
            axz.stairs(example["slice_ions"][:], time, color="#2166ac", fill=True, alpha=0.6)
            axz.set(xlabel="Relative time [ms]", ylabel="Ions / slice", title="All pads, time projection")
            axpad = fig.add_subplot(grid[1, 1])
            axpad.stairs(example["strongest_pad_ions"][:], time, color="#b35806", fill=True, alpha=0.7)
            ix, iy = example["strongest_pad_indices"][:]
            axpad.set(xlabel="Relative time [ms]", ylabel="Ions / slice", title=f"Pad ({ix}, {iy}): full sequence")
            fig.suptitle(f"{label} — accepted event {event}", fontsize=14)
            filename = f"{name}_pad_history.png"
            fig.savefig(OUT / filename, dpi=170)
            plt.close(fig)
            figures.append(r"\subsection*{" + label + f"; event {event}" + "}\n" +
                           r"\includegraphics[width=\linewidth]{" + filename + "}")
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        text = r"""\documentclass[11pt]{article}
\usepackage[margin=18mm]{geometry}
\usepackage{booktabs,graphicx,hyperref}
\begin{document}
\title{Event analysis: selection and pad histories}
\author{NAUSICAA}\date{}\maketitle
The integrated checkpoint reads source-frame hits, produces rounded ion counts,
samples a common drift distance, diffuses each ion, and constructs Petit's kNN
graph. Events pass when exactly one component survives the configured ion
threshold. Accepted events supply their full retained diffusion table to pad
counting. The cumulative Julia suite passes TESTS assertions.

\section*{Selection and exact charge accounting}
\begin{center}\small\setlength{\tabcolsep}{4pt}
\begin{tabular}{lrrrrr}\toprule
Sample & Accepted/tested & Retained ions & Pad ions & Below cut, saved & Mean ms \\\midrule
ROWS
\bottomrule\end{tabular}\end{center}
These events are a reproducible validation subset spread across each input file;
the acceptance fractions describe this subset. Retained and below-cut totals in
the table refer to accepted events. Every accepted event satisfies
$N_{\rm pads}=N_{\rm diffusion,retained}$ exactly. Rejected events retain their
decision and diffusion ledger. All events satisfy
$N_{\rm sampled}=N_{\rm retained}+N_{\rm outside}$ and
$N_{\rm retained}=N_{\rm graph}+N_{\rm below\ threshold}$.
Timings cover the computational pipeline after a warm-up, including allocations;
file reads and writes are outside the timed region.

\section*{Model and conventions}
The model retains the common-event thermal diffusion and per-hit rounding from
step two. Reconstruction cells have width $3\sigma_T$, maximum graph edges
$6\sigma_T$, and the configured nearest-neighbour count. The graph threshold is
applied to summed Int64 ion counts, inclusively; energies retain Petit's Float32
weights for the backbone. This gives exact equality at the configured threshold.

The selected track's energy-weighted x/y centre and minimum z define the grid
references, following Petit. Transverse grid edges receive one independently
sampled offset per axis in [0,pitch), shared by all slices. Whole-cell extensions
cover every retained diffusion-bin centre. All ions at a diffusion-bin centre
are assigned to its containing cell: this is the geometric collection
approximation. The independent physical pitch is read from the detector card.
Internal cells are half open; the final upper edge is included. Sparse pad
histories retain zero-count intervals through their complete saved axes.

The relative-time direction is DIRECTION. The first arrival edge is t=0;
the sign and its source-z coordinate are saved with each history. Both signs are
supported. With 2 mm slices and 0.1 mm/ms drift, each slice is 20 ms. Absolute
event arrival offsets are a separate convention from these saved relative times.
Drift, diffusion and grid RNG streams use SHA-256 keys built from run seed,
sample, event ID and stage; event processing order preserves each stream.

\section*{Validation}
Tests cover Petit occupancy parity away from its excluded final boundaries,
exact threshold equality, disconnected-event rejection, full retained-charge
collection, grid extension, integer-pitch phase shifts, final edges, empty
histories, zero slices, both time directions and input ownership. An end-to-end
synthetic signal passes and a two-cloud background fails. Real-event replay is
identical; saved HDF5 histories are read back to check charge conservation.

FIGURES
\clearpage
\section*{Reproduction and next increment}
\begin{verbatim}
make validate-event-step3 PYTHON=/path/to/python3.11-or-newer
# Presentation only:
make report-event-step3 PYTHON=/path/to/python3.11-or-newer
\end{verbatim}
Numerical tables, histories, diffusion tables, axes, seeds and source hashes are in:
\par\noindent\path{results/event_analysis/step3_validation.h5}.
\par The next increment establishes the production output schema and common
ensemble axes/denominators, then runs parameter scans and the full samples.
\section*{Configuration used}
\footnotesize\begin{verbatim}
CONFIG
\end{verbatim}
\end{document}
"""
        direction = "increasing source z gives later arrival" if summary["drift_direction"] == 1 else "decreasing source z gives later arrival"
        for marker, value in {"TESTS": str(summary["passed"]), "ROWS": "\n".join(rows),
                              "DIRECTION": direction, "FIGURES": "\n".join(figures), "CONFIG": config_text}.items():
            text = text.replace(marker, value)
        (OUT / "report.tex").write_text(text)
    print(OUT / "report.tex")


if __name__ == "__main__":
    render()
