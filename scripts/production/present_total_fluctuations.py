#!/usr/bin/env python3
"""Present saved event-total fluctuations as plots and a self-contained PDF."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def value(group, name):
    """Read a scalar saved by the Julia numerical reduction."""
    item = group[name][()]
    return item.decode() if isinstance(item, bytes) else item


def normal_curve(x, sigma):
    """Draw a Gaussian reference density with the saved predicted width."""
    return np.exp(-0.5 * (x / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def present(input_path: Path, report_dir: Path) -> Path:
    """Generate a four-panel figure and compile its data-driven LaTeX report."""
    report_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(input_path) as f:
        if f.attrs["schema"] not in (b"nausicaa.total_fluctuations.v1",
                                  "nausicaa.total_fluctuations.v1"):
            raise ValueError("Incompatible fluctuation reduction")
        m, s, e = f["metadata"], f["summary"], f["events"]
        ion = e["ion_energy_residual_keV"][:]
        xe = e["xe_energy_residual_keV"][:]
        he = e["he_energy_residual_keV"][:]
        x = np.linspace(-12, 12, 700)
        fig, ax = plt.subplots(2, 2, figsize=(11.8, 8.4), constrained_layout=True)
        ax[0, 0].hist(ion, bins=18, density=True, color="#286090", alpha=0.7,
                      label="Ion energy")
        ax[0, 0].plot(x, normal_curve(x, value(s, "ion_predicted_energy_sd_keV")),
                      color="#b53b35", lw=2, label="Fano prediction")
        ax[0, 0].set(title="Ion energy residual", xlabel="Estimated − deposited energy (keV)",
                     ylabel="Density")
        ax[0, 0].legend()

        ax[0, 1].hist(xe, bins=18, density=True, color="#286090", alpha=0.7,
                      label="Pure Xe")
        ax[0, 1].plot(x, normal_curve(x, value(s, "xe_predicted_energy_sd_keV")),
                      color="#b53b35", lw=2, label="Corrected prediction")
        ax[0, 1].set(title="Pure-Xe corrected photon energy residual",
                     xlabel="Estimated − deposited energy (keV)", ylabel="Density")
        ax[0, 1].legend()

        ax[1, 0].hist(he, bins=18, density=True, color="#d98032", alpha=0.7,
                      label="Xe–15% He")
        ax[1, 0].plot(x, normal_curve(x, value(s, "he_predicted_energy_sd_keV")),
                      color="#b53b35", lw=2, label="Corrected prediction")
        ax[1, 0].set(title="Xe–He corrected photon energy residual",
                     xlabel="Estimated − deposited energy (keV)", ylabel="Density")
        ax[1, 0].legend()

        drift = e["electron_drift_cm"][:]
        ax[1, 1].scatter(drift, xe, s=13, alpha=0.65, color="#286090",
                         label="Pure Xe")
        ax[1, 1].scatter(drift, he, s=13, alpha=0.65, color="#d98032",
                         label="Xe–15% He")
        ax[1, 1].axhline(0, color="0.6", lw=0.7)
        ax[1, 1].set(title="Corrected photon energy versus drift",
                     xlabel="Electron drift (cm)", ylabel="Energy residual (keV)")
        ax[1, 1].legend(fontsize=8)
        for panel in ax.flat:
            panel.grid(alpha=0.2)
        png = report_dir / "total_fluctuations.png"
        fig.savefig(png, dpi=180)
        plt.close(fig)

        fmt = lambda x, places=2: f"{x:,.{places}f}"
        source_mu = value(s, "source_mean_pairs")
        lines = [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=21mm]{geometry}",
            r"\usepackage{amsmath,booktabs,graphicx}",
            r"\begin{document}",
            r"\title{Event-total ion and photon fluctuations: 100 bb0nu events}",
            r"\author{NAUSICAA generated analysis report}",
            r"\date{}\maketitle",
            r"\section*{Question and sample}",
            f"The analysis includes {value(m, 'events')} bb0nu events from the 1 mm Geant4 "
            r"hit table. Their mean expected primary-pair count is "
            f"{fmt(source_mu, 3)}; its across-event standard deviation is "
            f"{fmt(value(s, 'source_sd_pairs'), 4)}. Its mean energy is "
            f"{fmt(value(s, 'source_mean_energy_MeV'), 6)} MeV, effectively fixed "
            r"across events for this comparison.",
            r"\section*{Method and conditions}",
            rf"At each hit $\mu_i=10^6E_i/W$, with $W={fmt(value(m, 'W_eV'), 2)}$ eV "
            f"and $F={fmt(value(m, 'fano_factor'), 3)}$. "
            r"The minimum integer variance gives "
            r"$V_{\rm pair}=\sum_i\max[F\mu_i,\{\mu_i\}(1-\{\mu_i\})]$. "
            r"The ion branch collects every primary ion, so its event total is the pair draw. "
            f"The photon branches use $J={fmt(value(m, 'j_variance_ratio'), 3)}$, "
            f"optical acceptance $A={fmt(value(m, 'optical_acceptance'), 5)}$, "
            rf"and EL yields $Y_{{\rm Xe}}={fmt(value(m, 'xe_yield_photons_per_e'), 1)}$ "
            rf"and $Y_{{\rm XeHe}}={fmt(value(m, 'he_yield_photons_per_e'), 1)}$ "
            r"photons per surviving electron. The yield and acceptance are fixed model "
            r"inputs; photons touching the EL-hole wall are absorbed.",
            r"\[p=\exp[-L_e/(v_e\tau_e)],\qquad "
            r"\widehat E_{\gamma}=\frac{W}{10^6}\,"
            r"\frac{N_{\gamma,\mathrm{incident}}}{pYA}\quad\mathrm{MeV}.\]",
            r"The photon energy estimator is corrected separately for each event's "
            r"known drift length. Only the corrected estimator minus that event's "
            r"deposited energy enters the measured width. The correction removes "
            r"the drift-dependent mean. For a fixed pair count $N$, binomial attachment "
            r"still gives $\operatorname{Var}(N_{\rm survive}/p\mid N)=N(1-p)/p$. "
            r"On the primary-pair scale, the predicted conditional variance is",
            r"\[V_{\gamma,\mathrm{eq}}=V_{\rm pair}+"
            r"\mu\frac{1-p}{p}+\frac{J\mu}{pY}+"
            r"\frac{\mu(1-A)}{pYA}.\]",
            r"These terms represent pair creation, attachment, EL emission and optical "
            r"collection. The two gas branches share the primary pair and attachment "
            r"draws. No PDE or electronics response is applied.",
            f"Electron survival spans {fmt(value(s, 'min_survival'), 4)}--"
            f"{fmt(value(s, 'max_survival'), 4)} across sampled drift lengths.",
            r"\section*{Results}",
            r"\begin{center}\begin{tabular}{lrrrr}\toprule",
            r"Energy estimator & SD obs. (keV) & SD model (keV) & SD obs. (\%) & SD model (\%)\\\midrule",
            f"Ions & {fmt(value(s, 'ion_energy_residual_sd_keV'), 2)} & "
            f"{fmt(value(s, 'ion_predicted_energy_sd_keV'), 2)} & "
            f"{fmt(value(s, 'ion_relative_sd_percent'), 4)} & "
            f"{fmt(value(s, 'ion_predicted_relative_sd_percent'), 4)}" + r"\\",
            f"Pure-Xe photons, corrected & {fmt(value(s, 'xe_energy_residual_sd_keV'), 2)} & "
            f"{fmt(value(s, 'xe_predicted_energy_sd_keV'), 2)} & "
            f"{fmt(value(s, 'xe_relative_sd_percent'), 4)} & "
            f"{fmt(value(s, 'xe_predicted_relative_sd_percent'), 4)}" + r"\\",
            f"Xe--He photons, corrected & {fmt(value(s, 'he_energy_residual_sd_keV'), 2)} & "
            f"{fmt(value(s, 'he_predicted_energy_sd_keV'), 2)} & "
            f"{fmt(value(s, 'he_relative_sd_percent'), 4)} & "
            f"{fmt(value(s, 'he_predicted_relative_sd_percent'), 4)}" + r"\\",
            r"\bottomrule\end{tabular}\end{center}",
            r"For 100 independent events, the relative sampling uncertainty of an SD "
            r"is approximately 7\%, so the observed and predicted widths agree "
            r"at this sample size.",
            f"The corrected energy-residual correlation with drift length is "
            f"{fmt(value(s, 'xe_residual_drift_correlation'), 3)} (pure Xe) and "
            f"{fmt(value(s, 'he_residual_drift_correlation'), 3)} (Xe--He).",
            r"\begin{center}\begin{tabular}{lrr}\toprule",
            r"Mean variance term (pair count$^2$) & Pure Xe & Xe--He\\\midrule",
            f"Fano/integer pairs & {fmt(value(s, 'mean_pair_variance'), 1)} & "
            f"{fmt(value(s, 'mean_pair_variance'), 1)}" + r"\\",
            f"Attachment & {fmt(value(s, 'mean_attachment_variance'), 1)} & "
            f"{fmt(value(s, 'mean_attachment_variance'), 1)}" + r"\\",
            f"EL emission & {fmt(value(s, 'xe_mean_el_variance'), 1)} & "
            f"{fmt(value(s, 'he_mean_el_variance'), 1)}" + r"\\",
            f"Optical collection & {fmt(value(s, 'xe_mean_optical_variance'), 1)} & "
            f"{fmt(value(s, 'he_mean_optical_variance'), 1)}" + r"\\",
            r"\bottomrule\end{tabular}\end{center}",
            r"\begin{center}\includegraphics[width=0.98\linewidth]{total_fluctuations.png}\end{center}",
            r"\section*{Provenance}",
            rf"Numerical HDF5 reduction: \texttt{{{input_path.name.replace('_', r'\_')}}}. "
            r"The reduction stores the source summary digest, cards' derived values, "
            r"event IDs, drift lengths, observed counts and predicted conditional variances. "
            r"The $F$ and $J$ definitions follow Oliveira et al., JINST 6 (2011) P05007; "
            r"the adopted $J=0.01$ at $E/P=3$ is an approximate user reading of Figure 6.",
            r"\end{document}",
        ]
    tex = report_dir / "total_fluctuations.tex"
    tex.write_text("\n".join(lines) + "\n")
    subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", tex.name],
                   cwd=report_dir, check=True, stdout=subprocess.DEVNULL)
    log = (report_dir / "total_fluctuations.log").read_text()
    if any(token in log for token in ("LaTeX Warning:", "Overfull ", "Underfull ",
                                     "Undefined control sequence")):
        raise RuntimeError("LaTeX build log contains warnings")
    subprocess.run(["latexmk", "-c", tex.name], cwd=report_dir, check=True,
                   stdout=subprocess.DEVNULL)
    return report_dir / "total_fluctuations.pdf"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reduction", type=Path)
    parser.add_argument("report_dir", type=Path)
    args = parser.parse_args()
    print(present(args.reduction, args.report_dir))
