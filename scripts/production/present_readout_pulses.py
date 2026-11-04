#!/usr/bin/env python3
"""Present saved ion pixel–frame and photon SiPM–bucket pulse distributions."""
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


SAMPLES = (("bb0nu", "Double beta"), ("single_electrons", "Single electrons"))
CONFIGS = {
    "ions": (("2", "2 mm / 20 ms", "#245d7a"),
             ("3", "3 mm / 30 ms", "#d47b36"),
             ("4", "4 mm / 40 ms", "#758a3b")),
    "photons": (("pure_xe", "Pure Xe", "#245d7a"),
                ("xe_he15", "Xe–15% He", "#d47b36")),
}


def scalar(group, name):
    """Read one saved HDF5 scalar, decoding fixed strings."""
    item = group[name][()]
    return item.decode() if isinstance(item, bytes) else item


def html_table(headers, rows):
    """Render a compact escaped HTML table from saved numerical values."""
    head = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) +
                   "</tr>" for row in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def tex_escape(value):
    """Escape the plain labels and formatted numbers used in report tables."""
    return (str(value).replace("\\", r"\textbackslash{}")
            .replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")
            .replace("µ", r"$\mu$").replace("–", "--"))


def tex_table(headers, rows):
    """Keep wide comparison tables inside the page margin."""
    spec = "l" * len(headers)
    lines = [r"\begin{center}\resizebox{\linewidth}{!}{\begin{tabular}{" + spec +
             r"}\toprule " + " & ".join(map(tex_escape, headers)) + r"\\\midrule"]
    lines += [" & ".join(map(tex_escape, row)) + r"\\" for row in rows]
    lines += [r"\bottomrule\end{tabular}}\end{center}"]
    return "\n".join(lines)


def rows_for(f, kind):
    """Collect exact event-summary rows already reduced by Julia."""
    rows = []
    for sample, sample_label in SAMPLES:
        for case, setting, _ in CONFIGS[kind]:
            g = f["cases"][kind][sample][case]
            rows.append((sample_label, setting, f"{int(scalar(g, 'events_with_charge')):,}",
                         f"{int(scalar(g, 'total_occupied_cells')):,}",
                         f"{int(scalar(g, 'total_signal')):,}",
                         f"{float(scalar(g, 'mean_of_event_means')):,.2f}",
                         f"{float(scalar(g, 'median_of_event_means')):,.2f}",
                         f"{float(scalar(g, 'mean_of_event_maxima')):,.1f}",
                         f"{float(scalar(g, 'median_of_event_maxima')):,.1f}"))
    return rows


def figure(f, kind, output, photon_bucket_us):
    """Plot saved, common-bin event-mean and event-maximum histograms."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.7), layout="constrained")
    for row, (sample, sample_label) in enumerate(SAMPLES):
        for column, (field, title) in enumerate((("mean", "Event mean"),
                                                  ("max", "Event maximum"))):
            ax = axes[row, column]
            edges = f["histogram_edges"][kind][field][:]
            for case, label, color in CONFIGS[kind]:
                g = f["cases"][kind][sample][case]
                counts = g["hist_" + field + "_counts"][:]
                ax.stairs(counts, edges, color=color, linewidth=1.6, label=label)
                if int(counts.sum()) != int(scalar(g, "events_with_charge")):
                    raise ValueError(f"Histogram entry mismatch: {kind}/{sample}/{case}/{field}")
            unit = ("Ions / occupied pixel–frame cell" if kind == "ions" else
                    f"Incident photons / SiPM–{photon_bucket_us:g} µs bucket")
            ax.set(title=f"{sample_label} · {title}", xlabel=unit, ylabel="Events",
                   xscale="log", yscale="log", xlim=(edges[0], edges[-1]))
            ax.grid(alpha=0.18)
            ax.legend(fontsize=9)
    fig.suptitle("Ion pixel–frame pulses" if kind == "ions" else
                 "Photon SiPM–bucket pulses", fontsize=15)
    fig.savefig(output, dpi=170)
    plt.close(fig)


def present(input_path: Path, output_dir: Path, latex_dir: Path) -> tuple[Path, Path]:
    """Write browser HTML, figures and a generated LaTeX/PDF report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    latex_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(input_path) as f:
        schema = f.attrs["schema"]
        schema = schema.decode() if isinstance(schema, bytes) else schema
        if schema != "nausicaa.readout_pulses.v2" or not scalar(f["metadata"], "complete"):
            raise ValueError("Incomplete or incompatible pulse reduction")
        m = f["metadata"]
        realisation = tomllib.loads(scalar(m, "realisation_card"))
        electron = tomllib.loads(scalar(m, "electron_card"))
        ion_cards = {pitch: tomllib.loads(scalar(m, f"ion_card_{pitch}mm"))
                     for pitch in (2, 3, 4)}
        photon_bucket_us = float(electron["bucket_width_us"])
        if photon_bucket_us != 0.25:
            raise ValueError("Expected 0.25 µs electron buckets")
        for sample, _ in SAMPLES:
            for gas, _, _ in CONFIGS["photons"]:
                g = f[f"cases/photons/{sample}/{gas}"]
                if (float(scalar(g, "bucket_duration")) != photon_bucket_us or
                        scalar(g, "bucket_unit") != "us"):
                    raise ValueError(f"Photon bucket mismatch: {sample}/{gas}")
        for pitch, card in ion_cards.items():
            if (card["pad_pitch_mm"] != pitch or
                    card["bucket_width_ms"] != 10 * pitch):
                raise ValueError(f"Unexpected ion pitch/frame card: {pitch} mm")
        ion_rows = rows_for(f, "ions")
        photon_rows = rows_for(f, "photons")
        for kind in ("ions", "photons"):
            png = figures / f"{kind}_mean_max.png"
            figure(f, kind, png, photon_bucket_us)
            shutil.copy2(png, latex_dir / png.name)
        headers = ("Species", "Setting", "Events", "Occupied cells", "Total signal",
                   "Mean of means", "Median mean", "Mean maximum", "Median maximum")
        ion_intro = ("An ion pulse is one pixel–frame cell with positive ion count. "
                     "The three ion grids use 2, 3 and 4 mm transverse pitch and "
                     "20, 30 and 40 ms frames, respectively, corresponding to "
                     "2, 3 and 4 mm along drift. Every event mean is its total ions "
                     "divided by its number of occupied pixel–frame cells; the event "
                     "maximum is the largest cell count.")
        photon_intro = (f"A photon pulse is one SiPM–{photon_bucket_us:g} µs "
                        "(250 ns) bucket with positive incident "
                        "VUV photon count. Pure Xe and Xe–15% He use the same 10 mm "
                        "hole pitch and electron settings. Every event mean divides its "
                        "total incident photons by its occupied SiPM–bucket count; "
                        "the event maximum is the largest bucket count. These are "
                        "incident photons before SiPM PDE and electronics.")
        common = (f"The 1 mm Geant4 source samples contribute 8,000 events each. "
                  f"The same Fano-sampled primary pairs and event placement feed all "
                  f"branches. W = {realisation['w_eV']:.2f} eV, F = "
                  f"{realisation['fano_factor']:.3f}; electron lifetime = "
                  f"{electron['lifetime_ms']:.1f} ms, E/P = "
                  f"{electron['reduced_field_kV_cm_bar']:.1f} kV/cm/bar, "
                  f"J = {electron['j_variance_ratio']:.3f}. All occupied cells and "
                  f"all events with pulses enter the distributions; no count cut is applied. "
                  f"Carriers were sampled directly into physical readout cells, with "
                  f"no diffusion-histogram truncation. "
                  f"No time-integrated pixel or SiPM charge enters these plots. "
                  f"The ion and photon ledgers are exact against their saved histories.")
        style = ("body{font:16px system-ui;max-width:1250px;margin:30px auto;padding:0 20px;"
                 "color:#243348}p{line-height:1.55}img{width:100%}table{width:100%;"
                 "border-collapse:collapse;font-size:13px;margin:16px 0 26px}th,td{border:1px "
                 "solid #dce3ec;padding:7px;text-align:right}th:first-child,td:first-child,"
                 "th:nth-child(2),td:nth-child(2){text-align:left}th{background:#edf2f6}"
                 "code{overflow-wrap:anywhere}")
        body = ("<h1>Ion and photon pulses per readout cell</h1>"
                + f"<p>{html.escape(common)}</p>"
                + "<h2>Ion pixel–frame pulses</h2>"
                + f"<p>{html.escape(ion_intro)}</p>"
                + html_table(headers, ion_rows)
                + '<img src="figures/ions_mean_max.png" alt="Ion event mean and maximum distributions">'
                + "<h2>Photon SiPM–bucket pulses</h2>"
                + f"<p>{html.escape(photon_intro)}</p>"
                + html_table(headers, photon_rows)
                + '<img src="figures/photons_mean_max.png" alt="Photon event mean and maximum distributions">'
                + "<h2>Provenance</h2>"
                + f"<p>Numerical HDF5 reduction: <code>{html.escape(str(input_path))}</code>. "
                + f"Ion grid-scan index SHA-256: <code>{html.escape(scalar(m, 'scan_index_sha256'))}</code>. "
                + f"Electron 250 ns index SHA-256: <code>{html.escape(scalar(m, 'electron_index_sha256'))}</code>. "
                + f"Ion run ID: <code>{html.escape(scalar(m, 'run_id'))}</code>; "
                + f"electron run ID: <code>{html.escape(scalar(m, 'electron_run_id'))}</code>. "
                + "The HDF5 reduction records each branch file, source and code identity, "
                  "per-event pulse totals and shared histogram edges.</p>")
        html_path = output_dir / "readout_pulses.html"
        html_path.write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
                             '<title>Ion and photon pulses</title><style>' + style +
                             '</style>' + body + '</html>')

        tex_lines = [r"\documentclass[11pt]{article}",
                     r"\usepackage[margin=18mm]{geometry}",
                     r"\usepackage{graphicx,booktabs}",
                     r"\begin{document}",
                     r"\title{Ion and photon pulses per readout cell}",
                     r"\author{NAUSICAA generated report}",
                     r"\date{}\maketitle",
                     r"\section*{Sample and assumptions}",
                     tex_escape(common),
                     r"\section*{Ion pixel--frame pulses}",
                     tex_escape(ion_intro),
                     tex_table(headers, ion_rows),
                     r"\begin{center}\includegraphics[width=\linewidth]{ions_mean_max.png}\end{center}",
                     r"\clearpage\section*{Photon SiPM--bucket pulses}",
                     tex_escape(photon_intro),
                     tex_table(headers, photon_rows),
                     r"\begin{center}\includegraphics[width=\linewidth]{photons_mean_max.png}\end{center}",
                     r"\section*{Provenance}",
                     rf"Numerical reduction: \texttt{{{input_path.name.replace('_', r'\_')}}}. "
                     r"\par\noindent "
                     rf"Ion run ID: \texttt{{{scalar(m, 'run_id')}}}. "
                     r"\par\noindent "
                     rf"Electron run ID: \texttt{{{scalar(m, 'electron_run_id')}}}. "
                     r"\par\noindent "
                     r"The HDF5 file preserves per-event values, histogram bins, and "
                     r"resolved cards.\par\noindent It also records source/code hashes "
                     r"and charge/light ledgers.",
                     r"\end{document}"]
        tex = latex_dir / "readout_pulses.tex"
        tex.write_text("\n".join(tex_lines) + "\n")
    subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", tex.name],
                   cwd=latex_dir, check=True, stdout=subprocess.DEVNULL)
    log = (latex_dir / "readout_pulses.log").read_text()
    if any(token in log for token in ("LaTeX Warning:", "Overfull ", "Underfull ",
                                     "Undefined control sequence")):
        raise RuntimeError("LaTeX build log contains warnings")
    subprocess.run(["latexmk", "-c", tex.name], cwd=latex_dir, check=True,
                   stdout=subprocess.DEVNULL)
    return html_path, latex_dir / "readout_pulses.pdf"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reduction", type=Path)
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "results/readout_pulses_8000")
    parser.add_argument("--latex-dir", type=Path,
                        default=ROOT / "latex_generated/readout_pulses_8000")
    args = parser.parse_args()
    html_path, pdf_path = present(args.reduction, args.output_dir, args.latex_dir)
    print("HTML:", html_path)
    print("PDF:", pdf_path)
