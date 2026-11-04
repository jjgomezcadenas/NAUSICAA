#!/usr/bin/env python3
"""Present saved full-chain bb0nu response statistics in HTML and LaTeX/PDF."""
from __future__ import annotations

import argparse
import html
import os
from pathlib import Path
import shutil
import subprocess
import tomllib

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work/matplotlib"))
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def scalar(g, key):
    value = g[key][()]
    return value.decode() if isinstance(value, bytes) else value


def tex(value):
    return (str(value).replace("\\", r"\textbackslash{}")
            .replace("_", r"\_").replace("%", r"\%")
            .replace("&", r"\&").replace("µ", r"$\mu$"))


def table(headers, rows):
    head = "".join(f"<th>{html.escape(str(x))}</th>" for x in headers)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(x))}</td>" for x in row) +
                   "</tr>" for row in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def tex_table(headers, rows):
    spec = "l" * len(headers)
    lines = [r"\begin{center}\resizebox{\linewidth}{!}{\begin{tabular}{" + spec +
             r"}\toprule " + " & ".join(map(tex, headers)) + r"\\\midrule"]
    lines += [" & ".join(map(tex, row)) + r"\\" for row in rows]
    lines += [r"\bottomrule\end{tabular}}\end{center}"]
    return "\n".join(lines)


def plots(f, path):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    ions = f["ions"]
    ax = axes[0, 0]
    ax.stairs(ions["snr_counts"][:], ions["snr_edges"][:], color="#245d7a")
    ax.set(xlabel="Noiseless S/N per positive ion pad–frame cell", ylabel="Cells",
           title="Ion CSA · 3 mm / 30 ms", yscale="log")
    ax.grid(alpha=.2)
    ax = axes[0, 1]
    ax.stairs(ions["noisy_bias_counts"][:], ions["noisy_bias_edges"][:], color="#758a3b")
    ax.set(xlabel="Signed event charge bias after noise", ylabel="Events",
           title="Ion event charge bias")
    ax.grid(alpha=.2)
    colors = {"pure_xe": "#245d7a", "xe_he15": "#d47b36"}
    labels = {"pure_xe": "Pure Xe", "xe_he15": "Xe–15% He"}
    for gas in colors:
        g = f[f"photons/{gas}"]
        axes[1, 0].stairs(g["clipped_counts"][:], g["clipped_edges"][:],
                          color=colors[gas], label=labels[gas])
        axes[1, 1].stairs(g["recovery_counts"][:], g["recovery_edges"][:],
                          color=colors[gas], label=labels[gas])
    axes[1, 0].set(xlabel="Clipped cells / all acquired SiPM–bucket cells per event",
                   ylabel="Events", title="ASIC saturation · 250 ns")
    axes[1, 1].set(xlabel="Sensor charge loss / primary plus dark PE per event",
                   ylabel="Events", title="SiPM recovery · 250 ns")
    for ax in axes[1]:
        ax.grid(alpha=.2)
        ax.legend()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def present(source: Path, result_dir: Path, latex_dir: Path):
    result_dir.mkdir(parents=True, exist_ok=True)
    latex_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(source) as f:
        if f.attrs["schema"] != b"nausicaa.response_case_analysis.v1" or not scalar(f["metadata"], "complete"):
            raise ValueError("Incomplete combined response analysis")
        m = f["metadata"]
        ion_card = tomllib.loads(scalar(m, "ion_card_toml"))
        sipm_card = tomllib.loads(scalar(m, "sipm_card_toml"))
        asic_card = tomllib.loads(scalar(m, "asic_card_toml"))
        ion = f["ions"]
        events = len(ion["event_id"])
        ion_rows = [(f"{events:,}", f"{int(scalar(ion,'input_ion_total')):,}",
                     f"{float(scalar(ion,'enc_ions_rms')):.2f}",
                     f"{float(scalar(ion,'mean_event_snr')):.2f}",
                     f"{100*float(scalar(ion,'mean_noiseless_bias_fraction')):.3f}%",
                     f"{100*float(scalar(ion,'mean_noisy_bias_fraction')):.3f}%",
                     f"{100*float(scalar(ion,'compression50_event_fraction')):.2f}%")]
        ion_headers = ("Events", "Input ions", "Frame ENC (ions RMS)",
                       "Mean event S/N", "Mean noiseless bias", "Mean noisy bias",
                       "Events with ≥50% local loss")
        electron_headers = ("Gas", "Incident photons", "Detected PE", "Dark pulses",
                            "Clipped cells / acquired", "Events with clipping",
                            "Mean recovery loss", "Clipped PE / sensor PE")
        electron_rows = []
        for gas, label in (("pure_xe", "Pure Xe"), ("xe_he15", "Xe–15% He")):
            g = f[f"photons/{gas}"]
            electron_rows.append((label, f"{int(scalar(g,'incident_total')):,}",
                                  f"{int(scalar(g,'detected_total')):,}",
                                  f"{int(scalar(g,'dark_total')):,}",
                                  f"{100*float(scalar(g,'clipped_cell_fraction_total')):.3f}%",
                                  f"{100*float(scalar(g,'events_with_clipping_fraction')):.2f}%",
                                  f"{100*float(scalar(g,'mean_recovery_loss_fraction')):.3f}%",
                                  f"{100*float(scalar(g,'clip_loss_per_sensor_pe_fraction')):.3f}%"))
        figure = result_dir / "response_distributions.png"
        plots(f, figure)
        shutil.copy2(figure, latex_dir / figure.name)
        conditions = (f"The bb0nu sample has {events:,} matched events. Positive ions use "
                      f"{scalar(m,'ion_pitch_mm'):g} mm pitch and {scalar(m,'ion_frame_ms'):g} ms "
                      f"frames; pure-Xe and Xe–15% He electrons use "
                      f"{scalar(m,'electron_bucket_us'):g} µs SiPM buckets. "
                      f"The CSA uses {1000*ion_card['averaging_s']:g} ms voltage averages, "
                      f"a three-frame inverse and stationary quiescent-point noise. "
                      f"SiPM PDE is {sipm_card['pde']:.2f}; the ASIC adds "
                      f"{asic_card['noise_rms_pe']:g} PE RMS bucket noise and clips "
                      f"bucket-average current at {1000*asic_card['current_limit_A']:g} mA. "
                      f"All selected channels' acquisition buckets are retained. "
                      f"The ion weighting field is the central square-pad benchmark; "
                      f"arrival times are uniform within their detector buckets.")
        interpretation = ("The clipping-cell denominator is every acquired SiPM–bucket cell, "
                          "including quiet buckets. Clipped PE divided by sensor PE is a "
                          "separate integrated charge-loss measure. Ion S/N is calculated "
                          "on positive input cells; no S/N selection is applied. Pure Xe "
                          "and Xe–15% He are alternative electron gas models.")
        style = ("body{font:16px/1.55 system-ui;max-width:1180px;margin:32px auto;padding:0 20px;"
                 "color:#243348}img{width:100%}table{border-collapse:collapse;width:100%;"
                 "font-size:13px}th,td{border:1px solid #dce3ec;padding:7px;text-align:right}"
                 "th:first-child,td:first-child{text-align:left}th{background:#edf2f6}"
                 "code{overflow-wrap:anywhere}")
        body = ("<h1>bb0nu detector response: ion CSA and SiPM/ASIC</h1>"
                + f"<p>{html.escape(conditions)}</p>"
                + "<h2>Ion CSA · 3 mm / 30 ms</h2>"
                + table(ion_headers, ion_rows)
                + "<h2>Electron SiPM/ASIC · 250 ns</h2>"
                + table(electron_headers, electron_rows)
                + f"<p>{html.escape(interpretation)}</p>"
                + '<img src="response_distributions.png" alt="Ion and SiPM response distributions">'
                + "<h2>Provenance</h2>"
                + f"<p>Numerical analysis: <code>{html.escape(str(source))}</code>. "
                + f"Ion manifest SHA-256: <code>{html.escape(scalar(m,'ion_manifest_sha256'))}</code>. "
                + f"SiPM manifest SHA-256: <code>{html.escape(scalar(m,'sipm_manifest_sha256'))}</code>. "
                + "The HDF5 records cards, source hashes, event arrays and histogram bins.</p>")
        page = result_dir / "response.html"
        page.write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
                        '<title>bb0nu detector response</title>' + f'<style>{style}</style>{body}</html>')
        lines = [r"\documentclass[11pt]{article}",
                 r"\usepackage[margin=18mm]{geometry}",
                 r"\usepackage{graphicx,booktabs}",
                 r"\begin{document}",
                 r"\title{bb0nu detector response: ion CSA and SiPM/ASIC}",
                 r"\author{NAUSICAA generated report}",
                 r"\date{}\maketitle",r"\section*{Operating conditions}",tex(conditions),
                 r"\section*{Ion CSA: 3 mm / 30 ms}",tex_table(ion_headers,ion_rows),
                 r"\section*{Electron SiPM/ASIC: 250 ns}",tex_table(electron_headers,electron_rows),
                 tex(interpretation),
                 r"\begin{center}\includegraphics[width=\linewidth]{response_distributions.png}\end{center}",
                 r"\section*{Provenance}",
                 r"Numerical results and exact cards are in the saved HDF5 analysis. "
                 + r"The source manifest hashes and event ledgers are retained there.",
                 r"\end{document}"]
        report = latex_dir / "response.tex"
        report.write_text("\n".join(lines) + "\n")
    subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", report.name],
                   cwd=latex_dir, check=True, stdout=subprocess.DEVNULL)
    log = (latex_dir / "response.log").read_text()
    if any(x in log for x in ("LaTeX Warning:", "Overfull ", "Underfull ",
                             "Undefined control sequence")):
        raise RuntimeError("LaTeX report log contains warnings")
    subprocess.run(["latexmk", "-c", report.name], cwd=latex_dir, check=True,
                   stdout=subprocess.DEVNULL)
    return page, latex_dir / "response.pdf"


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("analysis", type=Path)
    p.add_argument("--result-dir", type=Path, default=ROOT / "results/response_bb0nu_3mm")
    p.add_argument("--latex-dir", type=Path, default=ROOT / "latex_generated/response_bb0nu_3mm")
    a = p.parse_args()
    print(present(a.analysis, a.result_dir, a.latex_dir))
