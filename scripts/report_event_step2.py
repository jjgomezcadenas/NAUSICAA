"""Present the Julia ion/diffusion validation HDF5 as PNGs and standalone LaTeX.

All physical calculations, moments, losses and scan summaries come from Julia.
This stage formats the saved data and draws the final figures.
"""
from pathlib import Path
import json
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "latex_generated/event_analysis_step2"


def render():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    with h5py.File(ROOT / "results/event_analysis/step2_validation.h5") as f:
        example = f["gaussian_example"]
        scan = f["synthetic_scan"]
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained")
        axes[0].step(example["x_centres_mm"][:], example["observed_density_per_mm"][:],
                     where="mid", label="Sampled histogram", color="#2166ac")
        axes[0].plot(example["x_centres_mm"][:], example["expected_density_per_mm"][:],
                     label="Truncated Gaussian", color="#d6604d", ls="--")
        axes[0].set(xlabel="x [mm]", ylabel="Density [mm$^{-1}$]", title="Retained ions: one-point input")
        axes[0].legend(frameon=False)
        for padding in np.unique(scan["nsigma"][:]):
            choose = scan["nsigma"][:] == padding
            bins = scan["nbins"][:][choose]
            values = scan["variance_x_mm2"][:][choose]
            line, = axes[1].plot(bins, values, "o-", label=f"{padding:g} sigma padding")
            axes[1].axhline(scan["expected_variance_mm2"][:][choose][0], ls=":", lw=1,
                           color=line.get_color())
        axes[1].set(xlabel="Bins per axis", ylabel="Retained x variance [mm$^2$]", title="Histogram resolution scan")
        axes[1].legend(frameon=False, fontsize=8)
        fig.savefig(OUT / "diffusion_validation.png", dpi=180)
        plt.close(fig)
        summary = {"passed": int(f["tests/passed"][()]), "samples": []}
        rows = []
        runtime_rows = []
        for name, g in f["samples"].items():
            s = g["summary"]
            label = "Double beta" if name == "0nubb" else "Electron background"
            item = {"sample": label, **{key: float(value[()]) for key, value in s.items()}}
            summary["samples"].append(item)
            rows.append(
                f"{label} & {int(s['events'][()])} & {s['mean_rounding_difference'][()]:.2f} & "
                f"{100*s['max_absolute_relative_rounding'][()]:.3g} & "
                f"{s['max_absolute_relative_expectation_error'][()]:.3g} \\\\"
            )
            runtime_rows.append(
                f"{label} & {len(g['diffusion_cases/L_cm'])} & {1000*s['mean_seconds'][()]:.2f} & "
                f"{100*s['max_outside_fraction'][()]:.3g} \\\\"
            )
        # These metadata also supply the measured checkpoint table in the HTML.
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        config = f["configuration/card"].asstr()[()]
        info = {"tests": summary["passed"], "w": f["configuration/ionisation_energy_eV"][()],
                "temperature": f["configuration/temperature_K"][()], "field": f["configuration/field_V_cm"][()]}
        scan_rows = []
        baseline = scan["nbins"][:] == f["configuration/nbins"][()]
        for i in np.flatnonzero(baseline):
            scan_rows.append(f"{scan['nsigma'][i]:g} & {int(scan['sampled'][i])} & "
                             f"{int(scan['outside'][i])} & {100*scan['expected_outside_fraction'][i]:.5g} \\\\")
        text = r"""\documentclass[11pt]{article}
\usepackage[margin=20mm]{geometry}
\usepackage{booktabs,graphicx,hyperref}
\begin{document}
\title{Event analysis: ion production and diffusion}
\author{NAUSICAA}\date{}\maketitle
This checkpoint implements Petit's per-hit rounding and per-ion Gaussian Monte
Carlo histogram, with explicit RNG inputs, Float32 physical arrays, Int64 counts,
saved edges and a charge ledger. The cumulative Julia suite passes TESTS assertions.
The pinned replay fixture agrees with Petit's original Float64 diffusion routine
away from bin boundaries. Geometry, conservation, stochastic moments, finite-window
losses, exact edges and reproducibility have separate tests.

\section*{Model and numerical conventions}
The mean pair-creation energy is WVALUE eV. Counts are rounded independently per
hit, with ties to even. The model uses a common drift distance per event, sampled
uniformly from the configured interval. The thermal-width formula retains Petit's
$k_B/e=8.6173\times10^{-5}$ V/K, at TEMPERATURE K and FIELD V/cm. The baseline uses
the same width on all axes. The zero-longitudinal-width comparison is supported.
Coordinates remain in the source frame.

Histogram bounds span the hit extent plus the configured sigma padding. Internal
bins are half open; the final upper edge is included. A coincident axis with zero
diffusion uses one bin centred at that coordinate and adjacent representable
Float32 edges. Every ion consumes three Gaussian draws. Outside counts remain in
the ledger; retained counts are never renormalised. Int64 counts satisfy
$N_{\rm sampled}=N_{\rm retained}+N_{\rm outside}$ exactly.

\section*{Full-sample ion production}
\begin{center}\small\begin{tabular}{lrrrr}\toprule
Sample & Events & Mean rounding [ions] & Max rounding [\%] & Max arithmetic error \\\midrule
ROWS
\bottomrule\end{tabular}\end{center}
Rounding means sampled minus expected ions. The maximum rounding column uses
the absolute event-relative difference. Arithmetic error compares Float32 ion
expectations with an independent Float64 E/W sum; it is a relative error.

\section*{Diffusion checks}
\includegraphics[width=\linewidth]{diffusion_validation.png}
The point-source example has sigma=1 mm. The right panel compares retained
histogram-centre variance with the analytic truncated-Gaussian variance (dotted
lines). The same Gaussian draws are replayed at every bin count and padding.
Increasing padding recovers the Gaussian tails; changing bin count preserves
the retained count for fixed bounds. Finite-sample fluctuations remain after
the histogram resolution converges.

\begin{center}\begin{tabular}{rrrr}\toprule
Padding [sigma] & Sampled ions & Outside ions & Analytic outside [\%] \\\midrule
SCANROWS
\bottomrule\end{tabular}\end{center}

\begin{center}\begin{tabular}{lrrr}\toprule
Sample & Cases & Mean diffusion [ms] & Max outside [\%] \\\midrule
RUNTIMEROWS
\bottomrule\end{tabular}\end{center}
The real-event checks cover events spread across each file, at the configured
fixed drift distances and one sampled distance per event. Timings include histogram
allocation and output construction after compilation; measured values depend on
the machine. The HDF5 also records allocated bytes and seeds. These checks
establish the diffusion stage. Selection efficiencies and pad histories belong
to the next integration checkpoint.

\section*{Reproduction}
\begin{verbatim}
make validate-event-step2 PYTHON=/path/to/python3.11-or-newer
# Regenerate presentation from saved HDF5:
make report-event-step2 PYTHON=/path/to/python3.11-or-newer
\end{verbatim}
Source hashes, resolved inputs and measured tables are saved in:
\par\noindent\path{results/event_analysis/step2_validation.h5}.
\par
The HTML review contains the implemented signatures and source docstrings.
\section*{Configuration used}
\small\begin{verbatim}
CONFIG
\end{verbatim}
\end{document}
"""
        for marker, value in {
            "TESTS": str(info["tests"]), "WVALUE": f"{info['w']:.4g}",
            "TEMPERATURE": f"{info['temperature']:.6g}", "FIELD": f"{info['field']:.6g}",
            "RUNTIMEROWS": "\n".join(runtime_rows), "SCANROWS": "\n".join(scan_rows),
            "ROWS": "\n".join(rows), "CONFIG": config,
        }.items():
            text = text.replace(marker, value)
        (OUT / "report.tex").write_text(text)
    print(OUT / "report.tex")


if __name__ == "__main__":
    render()
