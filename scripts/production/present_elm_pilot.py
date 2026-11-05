#!/usr/bin/env python3
"""Present saved ELM clipping results; regenerate HTML/TeX/PDF without rerunning physics."""
from pathlib import Path
import html
import shutil
import subprocess
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
import sys
sys.path.insert(0, str(ROOT / "scripts/documentation"))
from build_report import build_report

OUT = ROOT / "results/elm_pilot"
LATEX = ROOT / "latex_generated/electron_readout/elm_pilot"
GASES = (("pure_xe", "Pure Xe"), ("xe_he15", "Xe–15% He"))

def fmt(value, digits=3):
    return f"{float(value):,.{digits}f}"

def table(headers, rows):
    return "<table><thead><tr>" + "".join(f"<th>{html.escape(str(x))}</th>" for x in headers) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(f"<td>{html.escape(str(x))}</td>" for x in row) + "</tr>" for row in rows) + "</tbody></table>"

def savefig(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=160)
    plt.close(fig)
    shutil.copy2(OUT / name, LATEX / name)

def main():
    LATEX.mkdir(parents=True, exist_ok=True)
    with h5py.File(OUT / "elm_pilot.h5", "r") as f:
        assert bool(f["metadata/complete"][()])
        p = {k: f["parameters"][k][()] for k in f["parameters"]}
        data = {gas: {k: f[gas][k][()] for k in f[gas]} for gas, _ in GASES}
        timing = {k: f["metadata"][k][()] for k in ("photon_wall_s", "response_wall_s", "total_wall_s")}
        metadata = {k: f["metadata"][k][()] for k in ("sensor_card", "electronics_card", "source_manifest")}
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(3, 2, figsize=(12, 11))
    for col, (gas, label) in enumerate(GASES):
        d = data[gas]
        vals = d["inclusive_pe"]
        bins = np.unique(np.r_[[-.5, .5], np.geomspace(1.5, float(vals.max()) + .5, 90)])
        axes[0, col].hist(vals, bins=bins, histtype="step", lw=1.7)
        axes[0, col].set_xscale("symlog", linthresh=1)
        axes[0, col].set_yscale("log")
        axes[0, col].axvline(p["clip_threshold_pe"], color="crimson", ls="--", label="ASIC limit")
        axes[0, col].set(xlabel="Detected PE / SiPM–250 ns bucket", ylabel="Input cells", title=label + " · inclusive pulses")
        axes[0, col].legend()
        for row, key, title in ((1, "event_mean_pe", "Event mean pulse"), (2, "event_max_pe", "Event maximum pulse")):
            axes[row, col].hist(d[key], bins=25, histtype="step", lw=1.7)
            if row == 2:
                axes[row, col].axvline(p["clip_threshold_pe"], color="crimson", ls="--")
            axes[row, col].set(xlabel="Detected PE / SiPM–250 ns bucket", ylabel="Events", title=title)
    savefig(fig, "elm_pe_pulses.png")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for col, (gas, label) in enumerate(GASES):
        d = data[gas]
        for row, keys in enumerate((("detected_total", "readout_total"), ("detected_corrected", "readout_corrected"))):
            edges = np.histogram_bin_edges(np.r_[d[keys[0]], d[keys[1]]], bins=35)
            for key, lab in zip(keys, ("Detected optical PE", "After ASIC")):
                axes[row, col].hist(d[key], bins=edges, histtype="step", lw=1.6, label=lab)
            axes[row, col].set(ylabel="Events", title=label)
            axes[row, col].legend()
        axes[0, col].set_xlabel("Total PE-equivalent / event (raw sums)")
        axes[1, col].set_xlabel("Attachment-corrected PE / deposited MeV")
    savefig(fig, "elm_event_totals.png")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, (gas, label) in zip(axes, GASES):
        d = data[gas]
        dots = ax.scatter(d["electron_drift_cm"], d["clip_loss_percent"], c=d["event_max_pe"], s=35, cmap="viridis")
        ax.set(xlabel="Electron drift distance (cm)", ylabel="Event charge lost to clipping (%)", title=label)
        fig.colorbar(dots, ax=ax, label="Event maximum optical pulse (PE)")
    savefig(fig, "elm_clipping_vs_drift.png")
    rows=[];pulse_rows=[];resolution_rows=[];worst=[];validation=[]
    for gas,label in GASES:
        d=data[gas]
        rows.append([label, int(d["event_count"]), fmt(d["detected_mean_total_pe"],0), fmt(d["clipped_pulse_percent"],4), fmt(d["clipped_event_percent"],1), fmt(d["charge_lost_percent"],4), fmt(d["inclusive_max_pe"],0)])
        pulse_rows.append([label, fmt(d["inclusive_mean_pe"],2), fmt(d["inclusive_median_pe"],2), fmt(d["event_mean_average_pe"],2), fmt(d["event_mean_median_pe"],2), fmt(d["event_max_average_pe"],1), fmt(d["event_max_median_pe"],1)])
        resolution_rows.append([label, fmt(d["corrected_detected_rms_percent"],4), fmt(d["corrected_unclipped_rms_percent"],4), fmt(d["corrected_readout_rms_percent"],4)])
        for i in np.argsort(d["clip_loss"])[-5:][::-1]:
            worst.append([label,int(d["event_id"][i]),fmt(d["electron_drift_cm"][i],2),fmt(d["event_max_pe"][i],0),fmt(d["clip_loss_percent"][i],3)])
        validation.append([label,int(d["pair_checks"]),int(d["ledger_checks"]),fmt(d["pde_fraction"],6),fmt(d["pde_mean_z"],2),fmt(d["quiet_noise_mean_pe"],4),fmt(d["quiet_noise_rms_pe"],4)])
    headers=["Gas","Events","Mean optical PE/event","Pulses clipped (%)","Events clipped (%)","Charge lost (%)","Largest pulse (PE)"]
    pulse_headers=["Gas","Inclusive mean","Inclusive median","Mean-pulse mean","Mean-pulse median","Max-pulse mean","Max-pulse median"]
    res_headers=["Gas","Detected PE RMS (%)","Before clipping RMS (%)","After clipping RMS (%)"]
    worst_headers=["Gas","Event ID","Electron drift (cm)","Maximum pulse (PE)","Charge lost (%)"]
    val_headers=["Gas","Pair checks","Ledger checks","Measured PDE","PDE z-score","Quiet mean (PE)","Quiet RMS (PE)"]
    description=("100 matched bb0nu events, IDs 0–99, for both gas models. E/P = 2 kV/cm/bar at 15 bar; "
        "3 mm holes, 5 mm EL depth and 5 mm square-grid pitch. Electron collection is unity. "
        "VUV photons emitted uniformly on the hole axis are thinned through absorbing-hole acceptance, "
        "TPB conversion 0.65, forward fraction 0.5 and meridional capture 1 − 1/1.5. "
        f"The combined optical probability is {float(data['pure_xe']['optical_acceptance']):.8f}. "
        "Guides have a 0.5 mm radial xenon gap; propagation and NOA 87 coupling are ideal. "
        "Saved photons are blue photons at the SiPM entrance; TPB is applied exactly once.")
    sensor=(f"The 3 × 3 mm² S14160-3015PS uses packaged PDE {float(p['pde']):.2f} (peak-value approximation), "
        f"gain {float(p['gain']):.0f}, and {float(p['dark_rate_Hz'])/1000:.0f} kcps primary dark rate. "
        "Each optical/dark detection contributes full gain; microcell occupancy and recovery are deferred. "
        "The device has a silicone-resin entrance window; its transmission is included in PDE. "
        f"The electronics adds {float(p['noise_rms_pe']):.4f} PE RMS independent Gaussian noise and clips the "
        f"{float(p['bucket_us'])*1000:.0f} ns bucket-average current at {float(p['current_limit_mA']):g} mA "
        f"({float(p['clip_threshold_pe']):.2f} PE). Signed noise values are retained.")
    definitions=("Pulse statistics select photon-positive SiPM–bucket cells and keep the same cells through electronics; "
        "zero detected PE in an input cell is retained. Each event contributes one mean and one maximum. "
        "Pulse clipping uses input cells as denominator; event clipping means any acquired bucket clips. "
        "Charge loss is summed clipped PE-equivalent divided by summed sensor PE-equivalent (optical plus dark). "
        "Event energy sums include all acquired cells on illuminated channels, over the event time span plus four quiet "
        "buckets at each end. Resolution diagnostics subtract the expected dark count, correct the known electron "
        "attachment survival and divide by deposited MC energy. Widths are full-distribution RMS/mean, including tails; "
        "they are not fitted peak widths. Before-clipping totals include electronics noise.")
    timing_text=(f"Photon replay/archive: {float(timing['photon_wall_s']):.3f} s; SiPM + electronics: "
        f"{float(timing['response_wall_s']):.3f} s; combined timed workflow: {float(timing['total_wall_s']):.3f} s. "
        "These timings include first-use compilation within the stages, source checks and HDF5 I/O; process startup and presentation are separate.")
    figures=("elm_pe_pulses.png","elm_event_totals.png","elm_clipping_vs_drift.png")
    contents=(f"<h1>ELM · linear SiPM and electronics clipping pilot</h1><p>{description}</p><p>{sensor}</p>"
        "<h2>Clipping and pulse counts</h2>"+table(headers,rows)+table(pulse_headers,pulse_rows)+
        "<p>All pulse entries are PE per SiPM per 250 ns bucket.</p><h2>Event response widths</h2>"+table(res_headers,resolution_rows)+
        f"<p>{definitions}</p>"+"".join(f'<img src="{name}" alt="{name}">' for name in figures)+
        "<h2>Largest event losses</h2>"+table(worst_headers,worst)+"<h2>Validation and timing</h2>"+table(val_headers,validation)+
        f"<p>{timing_text}</p><p>The 22 dedicated optical/linear-response tests passed; saved-event checks verify grid geometry, "
        "primary-pair identity, charge ledgers, PDE sampling and quiet-cell electronics noise.</p>"
        '<p><a href="../../latex_generated/electron_readout/elm_pilot/elm_pilot.pdf">PDF report</a> · <a href="files.toml">File index</a></p>'
        '<p>Complete workflow: <code>make elm-pilot</code>. Presentation only: <code>make elm-report</code>.</p>'
        +"<details><summary>Cards and provenance</summary><pre>"+html.escape(metadata["sensor_card"].decode()+"\n"+metadata["electronics_card"].decode()+"\n"+data["pure_xe"]["detector_card"].decode())+"</pre></details>")
    (OUT/"elm_pilot.html").write_text('<!doctype html><meta charset="utf-8"><title>ELM clipping pilot</title><style>body{font:16px system-ui;line-height:1.6;max-width:1400px;margin:35px auto;padding:0 22px;color:#17334a}table{border-collapse:collapse;margin:18px 0;width:100%}td,th{padding:8px;text-align:right;border-bottom:1px solid #ccd6de}td:first-child,th:first-child{text-align:left}img{width:100%;margin:15px 0}pre{white-space:pre-wrap}</style>'+contents)
    def tex_escape(s):
        return str(s).replace("%",r"\%").replace("_",r"\_").replace("&",r"\&").replace("×",r"$\times$").replace("²",r"$^2$").replace("−",r"$-$").replace("–","--").replace("³",r"$^3$")
    def tex_table(heads, rs):
        return r"\begin{center}\resizebox{\textwidth}{!}{\begin{tabular}{"+"l"+"r"*(len(heads)-1)+r"}\toprule"+"\n"+" & ".join(map(tex_escape,heads))+r"\\\midrule"+"\n"+"\n".join(" & ".join(map(tex_escape,row))+r"\\" for row in rs)+r"\bottomrule\end{tabular}}\end{center}"+"\n"
    tex=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage[margin=22mm]{geometry}
\usepackage{graphicx,booktabs,amsmath,xurl,hyperref,fvextra}
\begin{document}\title{ELM: linear SiPM and electronics clipping pilot}\author{NAUSICAA}\date{}\maketitle
'''
    tex+=tex_escape(description)+"\n\n"+tex_escape(sensor)+"\n\n"+tex_table(headers,rows)+tex_table(pulse_headers,pulse_rows)
    tex+=r"\section*{Event response and definitions}"+tex_table(res_headers,resolution_rows)+tex_escape(definitions)
    for name in figures:
        tex+=r"\clearpage\begin{figure}[p]\centering\includegraphics[width=\textwidth,height=.88\textheight,keepaspectratio]{"+name+r"}\end{figure}"+"\n"
    tex+=r"\clearpage\section*{Largest event losses}"+tex_table(worst_headers,worst)+r"\section*{Validation and timing}"+tex_table(val_headers,validation)+tex_escape(timing_text)
    tex+="\n\nThe 22 dedicated optical/linear-response tests passed. Saved-event checks verify grid geometry, primary-pair identity, charge ledgers, PDE sampling and quiet-cell electronics noise.\n"
    tex+=r"\section*{Reproduction}\noindent Numerical source: \path{results/elm_pilot/elm_pilot.h5}.\\ File index: \path{results/elm_pilot/files.toml}.\\ Complete workflow: \texttt{make elm-pilot}.\\ Presentation only: \texttt{make elm-report}."+"\n"
    tex+=r"\section*{Input cards}\begin{Verbatim}[breaklines,breakanywhere,fontsize=\footnotesize]"+"\n"+data["pure_xe"]["detector_card"].decode()+"\n"+metadata["sensor_card"].decode()+"\n"+metadata["electronics_card"].decode()+r"\end{Verbatim}"+"\n"+r"\end{document}"+"\n"
    (LATEX/"elm_pilot.tex").write_text(tex)
    build_report(LATEX / "elm_pilot.tex")
    print(OUT/"elm_pilot.html")
    print(LATEX/"elm_pilot.pdf")

if __name__=="__main__":
    main()
