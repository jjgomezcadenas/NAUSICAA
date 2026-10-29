"""Figures and a compilable LaTeX report from Julia's saved ion-track analysis."""
from __future__ import annotations

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]


import argparse
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work/matplotlib-cache"))
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def scalar(dataset):
    value = dataset[()]
    return value.decode() if isinstance(value, bytes) else value


def tex(text):
    replacements = {"\\": r"\textbackslash{}", "_": r"\_", "%": r"\%",
                    "&": r"\&", "#": r"\#", "{": r"\{", "}": r"\}",
                    "$": r"\$", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(replacements.get(c, c) for c in str(text))


def load_analysis(path):
    with h5py.File(path, "r") as f:
        schema = f.attrs["schema"]
        if isinstance(schema, bytes):
            schema = schema.decode()
        if schema != "nausicaa.ion_tracks.v1":
            raise ValueError(f"Unknown schema: {schema}")
        result = {group: {name: scalar(ds) for name, ds in f[group].items()
                          if isinstance(ds, h5py.Dataset)}
                  for group in ("metadata", "parameters", "summary", "validation", "bins", "mean", "example")}
        cfg = f.get("metadata/source_configuration")
        result["source_configuration"] = list(zip(cfg["param_key"].asstr()[:], cfg["param_value"].asstr()[:])) if cfg else []
        for group in ("mean", "example"):
            expected = (len(result["bins"]["y_edges_mm"])-1, len(result["bins"]["x_edges_mm"])-1)
            if result[group]["xy"].shape != expected:
                raise ValueError("XY shape must match HDF5 disk order [y,x]")
            axes = f[f"{group}/xy"].attrs["axis_order_on_disk"]
            if axes not in ("y,x", b"y,x"):
                raise ValueError("Unknown matrix axis order")
            if result[group]["z"].shape != (len(result["bins"]["z_edges_mm"])-1,):
                raise ValueError("Z shape does not match edges")
        if "selection" in f:
            result["selection"] = {name: scalar(ds) for name,ds in f["selection"].items()}
        if "diagnostics" in f:
            result["diagnostics"] = {name: scalar(ds) for name,ds in f["diagnostics"].items()}
        result["path"] = Path(path)
        return result


def plot_sample(data):
    bins, p, meta = data["bins"], data["parameters"], data["metadata"]
    output = data["path"].parent
    files = {}
    for group, suffix in (("mean", "mean"), ("example", f'event_{data["example"]["event_id"]}')):
        label = "Mean per event" if group == "mean" else f'Event {data["example"]["event_id"]}'
        units = "Expected ions / bin / event" if group == "mean" else "Expected ions / bin"
        values = data[group]["xy"]
        positive = values[values > 0]
        norm = LogNorm(vmin=float(positive.min()), vmax=float(positive.max())) if positive.size and positive.max()>positive.min() else None
        fig, ax = plt.subplots(figsize=(6.4, 4.8), layout="constrained")
        artist = ax.imshow(np.ma.masked_less_equal(values, 0), origin="lower", interpolation="nearest",
                           extent=(bins["x_edges_mm"][0], bins["x_edges_mm"][-1],
                                   bins["y_edges_mm"][0], bins["y_edges_mm"][-1]),
                           cmap="viridis", norm=norm, aspect="equal")
        fig.colorbar(artist, ax=ax, label=units)
        if not positive.size:
            ax.text(.5,.5,"No retained hits",transform=ax.transAxes,ha="center")
        ax.set(xlabel="x relative to generation vertex (mm)", ylabel="y relative to generation vertex (mm)", title=f'{meta["sample_label"]}: {label}')
        files[f"xy_{group}"] = output / f"xy_{suffix}.png"
        fig.savefig(files[f"xy_{group}"], dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.4, 4.8), layout="constrained")
        ax.stairs(data[group]["z"], bins["z_edges_mm"], color="#245e8b", linewidth=1.2)
        ax.set(xlabel="z relative to generation vertex (mm)", ylabel=units, title=f'{meta["sample_label"]}: {label}',
               xlim=(bins["z_edges_mm"][0],bins["z_edges_mm"][-1]))
        ax.set_yscale("symlog", linthresh=1)
        ax.grid(alpha=.18)
        # Both coordinate systems were computed and saved by Julia.
        top = ax.twiny()
        top.set_xlim(bins["relative_time_edges_ms"][[0,-1]])
        top.set_xlabel("Relative longitudinal time (ms)")
        files[f"z_{group}_full"] = output / f"z_{suffix}_full.png"
        fig.savefig(files[f"z_{group}_full"], dpi=180)
        halfwidth = data.get("selection",{}).get("xy_halfwidth_mm",max(abs(bins["z_edges_mm"][[0,-1]])))
        ax.set_xlim(-halfwidth,halfwidth)
        top.set_xlim(np.interp([-halfwidth,halfwidth], bins["z_edges_mm"], bins["relative_time_edges_ms"]))
        ax.set_title(f'{meta["sample_label"]}: {label}, central detail')
        files[f"z_{group}"] = output / f"z_{suffix}.png"
        fig.savefig(files[f"z_{group}"], dpi=180)
        plt.close(fig)
    if "diagnostics" in data:
        d = data["diagnostics"]
        halfwidth = data["selection"]["xy_halfwidth_mm"]
        for quantity, units in (("ions","Expected ions / bin / event"),("hits","Hits / bin / event")):
            fig, axes = plt.subplots(2,2,figsize=(9,6),layout="constrained")
            for col, axis in enumerate(("x","y")):
                for row in (0,1):
                    ax = axes[row,col]
                    ax.stairs(d[f"uncut_{axis}_{quantity}"],d[f"{axis}_edges_mm"],color="#245e8b")
                    ax.set_yscale("symlog",linthresh=.01)
                    ax.axvline(-halfwidth,color="#b46618",ls="--",lw=1)
                    ax.axvline(halfwidth,color="#b46618",ls="--",lw=1)
                    if row:
                        ax.set_xlim(-3*halfwidth,3*halfwidth)
                    ax.set(xlabel=f"{axis} relative to vertex (mm)",ylabel=units,
                           title=f"{axis.upper()}: {'central detail' if row else 'full distribution'}")
                    ax.grid(alpha=.15)
            fig.suptitle(f'{meta["sample_label"]}: before spatial cut')
            files[f"tails_{quantity}"] = output/f"xy_tails_{quantity}.png"
            fig.savefig(files[f"tails_{quantity}"],dpi=180)
            plt.close(fig)
    return files


def number(value, digits=6):
    return r"\num{" + format(float(value), f".{digits}g") + "}"


def make_report(datasets, figures, report_dir):
    report_dir.mkdir(parents=True, exist_ok=True)
    p = datasets[0]["parameters"]
    macros = {"IonWorkEv": p["w_eV"], "IonPitchMm": p["pitch_mm"], "IonFrameMs": p["frame_ms"],
              "IonZBinMm": p["z_bin_mm"], "IonSpeedMmMs": p["drift_speed_mm_ms"],
              "IonMarginFraction": p["margin_fraction"], "IonTolerance": p["tolerance"]}
    (report_dir / "values.tex").write_text("% Generated from analysis.h5; do not edit.\n" + "".join(
        rf"\newcommand{{\{key}}}{{{value:.12g}}}"+"\n" for key,value in macros.items()))
    sections = []
    for index, (data, plots) in enumerate(zip(datasets, figures)):
        meta, summary = data["metadata"], data["summary"]
        sid = meta["sample_id"]
        rows = []
        for label, key in ((r"Energy (MeV)", "energy_MeV"), (r"Expected ions/event", "expected_ions"),
                           (r"Mean ions/hit in an event", "mean_ions_per_hit"), (r"Hits/event", "n_hits_per_event")):
            rows.append(tex(label)+" & "+" & ".join(number(v) for v in summary[key])+r" \\")
        table = "\n".join([r"\begin{tabular}{lrrrr}\toprule",r"Quantity & Mean & Population SD & Minimum & Maximum \\\midrule",
                            *rows, r"\bottomrule\end{tabular}"])
        (report_dir / f"{sid}_tables.tex").write_text(table+"\n")
        prefix = "IonSample"+chr(ord("A")+index)
        vals = {"Events": summary["n_events"], "Hits": summary["n_hits"],
                "TotalIons": summary["total_expected_ions"], "MeanIons": summary["expected_ions"][0],
                "MeanEnergyMeV": summary["energy_MeV"][0], "PooledIonsPerHit": summary["pooled_mean_ions_per_hit"]}
        (report_dir / f"{sid}_values.tex").write_text("% Generated from analysis.h5\n"+"".join(
            rf"\newcommand{{\{prefix}{key}}}{{{float(value):.12g}}}"+"\n" for key,value in vals.items()))
        block = [r"\clearpage", r"\section{"+tex(meta["sample_label"])+"}",
                 r"\input{"+f"{sid}_values.tex"+"}",
                 "Saved events: "+number(summary["n_events"])+"; hits: "+number(summary["n_hits"])+".",
                 "Total expected ions: "+number(summary["total_expected_ions"])+". Pooled ions/hit: "+number(summary["pooled_mean_ions_per_hit"])+".",
                 r"\begin{center}\small\input{"+f"{sid}_tables.tex"+r"}\end{center}"]
        if "selection" in data:
            sel=data["selection"]
            block.append("The table above uses all original hits. After the vertex-relative XY cut, the mean retained energy is "+number(sel["energy_MeV"][0])+" MeV and the mean retained ion count is "+number(sel["mean_expected_ions"])+" per event.")
            block.append("Discarded hits: "+number(sel["n_discarded_hits"])+" ("+number(100*sel["discarded_hit_fraction"])+r"\%). Discarded ions: "+number(sel["total_discarded_ions"])+" ("+number(100*sel["discarded_ion_fraction"])+r"\%). Empty events: "+number(sel["empty_events"])+".")
            block.append("Maximum discarded ion fraction in an individual event: "+number(100*sel["maximum_event_discarded_fraction"])+r"\%. Discarded energy is never restored by rescaling.")
        for kind in ("mean", "example"):
            block += [r"\begin{figure}[htbp]\centering"]
            for axis in ("xy", "z"):
                relative = os.path.relpath(plots[f"{axis}_{kind}"], report_dir)
                block.append(r"\includegraphics[width=.49\linewidth]{\detokenize{"+relative+"}}")
            caption = "Selected ensemble mean, with zero contributions included in the denominator." if kind == "mean" else "Selected hits of example event "+str(data["example"]["event_id"])+"."
            block += [r"\caption{"+tex(caption)+" XY is integrated over z; Z is integrated over x and y. The right panel shows the central longitudinal detail; it applies no additional Z cut.}", r"\end{figure}"]
        if "diagnostics" in data:
            block += [r"\clearpage\subsection*{Tail diagnostics}"]
            for key in ("tails_hits","tails_ions","z_mean_full"):
                relative=os.path.relpath(plots[key],report_dir)
                block += [r"\begin{figure}[htbp]\centering",
                          r"\includegraphics[width=.82\linewidth]{\detokenize{"+relative+"}}",
                          r"\caption{"+("Full longitudinal distribution of retained hits." if key=="z_mean_full" else "Uncut vertex-relative distributions. Dashed lines mark the XY cut; bottom panels show a detail of the same histograms.")+"}",r"\end{figure}"]
        sections.append("\n".join(block))
    checks = []
    sources = []
    for data in datasets:
        v, meta, bins = data["validation"], data["metadata"], data["bins"]
        checks.append(tex(meta["sample_label"])+" & "+" & ".join(number(v[k],3) for k in
                      ("event_max_relative_error","xy_mean_relative_error","z_mean_relative_error"))+r" \\")
        cfg = dict(data["source_configuration"])
        source_rows = [("Input",meta["input_file"]),("Original",meta["original_file"]),
                       ("Generator",cfg.get("/nexus/RegisterGenerator","not recorded")),
                       ("Generated events",cfg.get("num_events","not recorded")),
                       ("Source lower energy cut",cfg.get("/Actions/DefaultEventAction/min_energy","not recorded")),
                       ("Source upper energy cut",cfg.get("/Actions/DefaultEventAction/max_energy","not recorded")),
                       ("Source radius",cfg.get("/Geometry/NextPrecdr/radius","not recorded")),
                       ("Source height",cfg.get("/Geometry/NextPrecdr/height","not recorded")),
                       ("Source pressure",cfg.get("/Geometry/NextPrecdr/pressure","not recorded"))]
        sources += [r"\subsection*{"+tex(meta["sample_label"])+"}",r"\begin{tabular}{lp{.65\linewidth}}\toprule",
                    *[tex(k)+" & "+tex(val)+r" \\" for k,val in source_rows],r"\bottomrule\end{tabular}",
                    r"\par\small Input SHA-256: \texttt{"+tex(meta["input_sha256"])+"}.",
                    r"\par Analysis code SHA-256: \texttt{"+tex(meta["analysis_code_sha256"])+"}.\normalsize"]
    bins = datasets[0]["bins"]
    grid_text = (f"Common XY grid: {len(bins['x_edges_mm'])-1} by {len(bins['y_edges_mm'])-1} bins, "
                 f"from {bins['x_edges_mm'][0]:g} to {bins['x_edges_mm'][-1]:g} mm on each axis. "
                 f"Z: {len(bins['z_edges_mm'])-1} bins, from {bins['z_edges_mm'][0]:g} to {bins['z_edges_mm'][-1]:g} mm.")
    cut_text = ""
    if "selection" in datasets[0]:
        cut_text = (r"\par The spatial selection keeps $|x'|\leq " + str(datasets[0]["selection"]["xy_halfwidth_mm"]) +
                    r"$ mm and the same bound on $|y'|$. It is applied after translating to the generation vertex. "
                    r"No Z cut or subsequent recentering is applied. Raw energy totals include every original hit; "
                    r"the mean maps use the selected hits and the original number of events. The full uncut marginal "
                    r"histograms and per-event discarded energy are retained. Common grids enclose every selected hit "
                    r"with the configured margin; bins are left-closed/right-open, with the final upper edge included.")
    preamble = r"""\documentclass[10pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern,microtype,amsmath,booktabs,graphicx,siunitx}
\usepackage[margin=20mm]{geometry}
\usepackage[colorlinks=true,urlcolor=blue,citecolor=blue]{hyperref}
\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}
\setlength{\emergencystretch}{2em}
\input{values}
\title{Initial ionization from NEXUS tracks}
\author{NAUSICAA analysis}\date{}
\begin{document}\maketitle
\section{Calculation and interpretation}
For each deposited-energy hit, the expected number of positive ions is
\[
 \overline n_i=10^6 E_i[\mathrm{MeV}]/W[\mathrm{eV}],\qquad
 E_{\rm evt}=\sum_i E_i,\qquad \overline N_{\rm evt}=\sum_i\overline n_i.
\]
We adopt $W=\IonWorkEv$ eV as a configurable nominal mean energy per electron--ion
pair in gaseous xenon. It describes the partition of deposited energy between
ionization and other processes, rather than the atomic ionization threshold;
see Ref.~\cite{Serra:2014zda}. This is an adopted input, not a measurement from these samples.
Counts are continuous expectations, not rounded integer realizations. Each
molecular $\mathrm{Xe}_2^+$ ion carries charge $+e$.

The per-event mean ions per hit is $\overline N_{\rm evt}/N_{\rm hits}$.
The table reports the mean and population dispersion of this quantity across
events. The pooled value instead divides total ions by total hits. Hit-level
quantities depend on the Geant4 deposition subdivision.

Each event is translated to the common generation vertex of its primary
particles, copied from the original NEXUS particle table:
$\mathbf r'_i=\mathbf r_i-\mathbf r_{\rm vertex}$. Orientations are preserved;
no energy barycenter is used to define the origin.
XY bins have pitch \IonPitchMm~mm. With the configured uniform speed
\IonSpeedMmMs~mm/ms, a \IonFrameMs~ms interval corresponds to a Z bin of
\IonZBinMm~mm. The time axis is relative to the generation vertex, not an absolute
arrival-time prediction. Projections are averaged over all saved events,
including zero contributions, without normalizing individual event energies.

This stage describes initial ion production. It includes no ionization
fluctuations, recombination, diffusion, focusing or collection losses.
The two marginal projections do not determine joint pad-and-frame occupancy.
"""
    report = preamble + "\n" + cut_text + "\n" + tex(grid_text) + "\n" + "\n".join(sections)
    report += r"""
\clearpage\section{Validation and provenance}
The Julia calculation checks that each event's XY and Z projections sum to its
retained expected ion count. The summed mean projections must equal the retained mean ion count
per event. No renormalization is applied. Out-of-range hits cause an error.
The following relative errors must be below \num{\IonTolerance}.
\begin{center}\begin{tabular}{lrrr}\toprule
Sample & Worst event & Mean XY & Mean Z\\\midrule
""" + "\n".join(checks) + r"\bottomrule\end{tabular}\end{center}" + "\n" + "\n".join(sources)
    report += r"""
\par Numerical quantities, plots and tables are read from the Julia HDF5 result.
The HDF5 files also retain configuration and software hashes, per-event
quantities, bin edges and per-hit ion expectations in source row order.
\bibliographystyle{unsrt}\bibliography{references}
\end{document}
"""
    (report_dir / "report.tex").write_text(report)
    bibliography = (ROOT / "notes/physics/references.bib").read_text()
    entry = re.search(r"(?ms)^@article\{Serra:2014zda,.*?(?=^@|\Z)", bibliography)
    if entry is None:
        raise ValueError("Canonical bibliography is missing Serra:2014zda")
    (report_dir / "references.bib").write_text(entry.group().strip()+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path)
    parser.add_argument("--report-dir", type=Path, default=ROOT / "latex_generated/ion_tracks")
    args = parser.parse_args()
    paths = args.inputs or sorted((ROOT / "results/ion_tracks").glob("*/analysis.h5"))
    if not paths:
        parser.error("No analysis.h5 files; run the Julia analysis first")
    data = [load_analysis(path) for path in paths]
    for item in data[1:]:
        if item["metadata"]["config_sha256"] != data[0]["metadata"]["config_sha256"]:
            raise ValueError("Mixed configurations: regenerate all sample analyses together")
        for axis in ("x_edges_mm", "y_edges_mm", "z_edges_mm"):
            if not np.array_equal(item["bins"][axis], data[0]["bins"][axis]):
                raise ValueError("Samples do not share common bin edges")
    with plt.rc_context({"font.size": 10, "axes.spines.right": False}):
        figures = [plot_sample(item) for item in data]
    make_report(data, figures, args.report_dir)
    print(f"Generated {sum(len(f) for f in figures)} PNG figures and {args.report_dir / 'report.tex'}")


if __name__ == "__main__":
    main()
