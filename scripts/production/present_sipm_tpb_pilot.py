#!/usr/bin/env python3
"""Generate TPB pilot figures, HTML and self-contained LaTeX from Julia results."""
from pathlib import Path
import html
import shutil
import subprocess
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
import sys
sys.path.insert(0, str(ROOT / "scripts/documentation"))
from build_report import build_report

OUT = ROOT / "results/sipm_tpb_pilot"
TEX = ROOT / "latex_generated/electron_readout/sipm_tpb_pilot"
TEX.mkdir(parents=True, exist_ok=True)
text = lambda x: x.decode() if isinstance(x, bytes) else str(x)

def table(headers, rows):
    return "<table><thead><tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + "</tr></thead><tbody>" + "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(x))}</td>" for x in r) + "</tr>" for r in rows) + "</tbody></table>"

def tex_table(headers, rows):
    esc = lambda s: str(s).replace("%", r"\%").replace("_", r"\_")
    return (r"\begin{center}\small\begin{tabular}{" + "l" * len(headers) + r"}\toprule" + "\n"
            + " & ".join(map(esc, headers)) + r"\\\midrule" + "\n"
            + "\n".join(" & ".join(map(esc, r)) + r"\\" for r in rows)
            + "\n" + r"\bottomrule\end{tabular}\end{center}")

def main():
    summary = []
    widths = []
    pulses = []
    sensitivity = []
    converged = []
    with h5py.File(OUT / "analysis.h5") as f:
        assert f["metadata/complete"][()]
        nominal = "tau12p5_n1024" if "tau12p5_n1024" in f else "tau12p5_n512"
        for point in ("ep2", "ep2p5"):
            for gas, label in (("pure_xe", "Xe"), ("xe_he15", "Xe/He")):
                g = f[f"{nominal}/{point}/{gas}"]
                ep = "2.0" if point == "ep2" else "2.5"
                val = lambda k: float(g[k][()])
                summary.append([ep, label, f'{val("detected_total_distribution/mean")/1e6:.3f}',
                    f'{val("sensor_charge_retained_percent"):.2f}', f'{val("clip_loss_percent"):.3f}',
                    f'{val("events_clipped_percent"):.0f}', f'{val("positive_cells_clipped_percent"):.3f}'])
                widths.append([ep, label] + [f'{val(k+"_corrected_distribution/relative_rms_percent"):.3f}'
                    for k in ("incident", "detected", "sensor", "asic")])
                pulses.append([ep,label]+[f'{val(k+"_distribution/"+stat):.1f}'
                    for k in ("event_mean_sensor", "event_max_sensor") for stat in ("mean","median")])
                for setting in ("tau12p5_n512", "tau12p5_n1024", "tau25p0_n512"):
                    if setting not in f:
                        continue
                    s = f[f"{setting}/{point}/{gas}"]
                    sensitivity.append([setting.replace("tau", "").replace("p", ".").replace("_n", " / "),ep,label,
                        f'{s["clip_loss_percent"][()]:.3f}',
                        f'{s["sensor_corrected_distribution/relative_rms_percent"][()]:.3f}',
                        f'{s["asic_corrected_distribution/relative_rms_percent"][()]:.3f}'])
                if "tau12p5_n1024" in f:
                    c = f[f"convergence/{point}/{gas}"]
                    converged.append([ep,label,
                        f'{c["mean_sensor_change_percent"][()]:+.4f}',
                        f'{c["clip_loss_change_percentage_points"][()]:+.4f}'])

        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        for row, (gas, label) in enumerate((("pure_xe", "Xe"), ("xe_he15", "Xe/He"))):
            for point, color, ep in (("ep2", "#2764a5", "2.0"), ("ep2p5", "#c45a21", "2.5")):
                g = f[f"{nominal}/{point}/{gas}"]
                p = g["sensor_pulses"]
                axes[row,0].stairs(p["log_counts"][:], p["log_edges"][:], color=color, label=f"E/P = {ep}")
                for stage, style in (("sensor", "--"), ("asic", "-")):
                    d = g[f"{stage}_corrected_distribution"]
                    axes[row,1].stairs(d["counts"][:], d["edges"][:]/d["mean"][()],
                        color=color, linestyle=style,label=f"{ep}, {stage.upper()}")
            axes[row,0].axvline(float(g["bucket_limit_pe"][()]), color="black",ls=":",label="1 mA limit")
            axes[row,0].set(xscale="log",yscale="log",xlabel="Sensor charge before clipping [PE / SiPM / 250 ns]",
                ylabel="Positive-input pulses",title=label)
            axes[row,1].set(xlabel="Attachment-corrected response / MeV, normalized to its mean",
                ylabel="Events",title=label)
            for ax in axes[row]:
                ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.savefig(OUT / "tpb_response.png", dpi=160)
        plt.close(fig)
        fig, axes = plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
        for row, (gas,label) in enumerate((("pure_xe","Xe"),("xe_he15","Xe/He"))):
            for point,color,ep in (("ep2","#2764a5","2.0"),("ep2p5","#c45a21","2.5")):
                g=f[f"{nominal}/{point}/{gas}"]
                for col,stage in enumerate(("detected","asic")):
                    d=g[f"{stage}_total_distribution"]
                    axes[row,col].stairs(d["counts"][:],d["edges"][:],color=color,label=f"E/P = {ep}")
                    axes[row,col].set(xlabel=f"Total {stage} PE per event",ylabel="Events",title=label)
            for ax in axes[row]:ax.legend();ax.grid(alpha=.2)
        fig.savefig(OUT/"tpb_totals.png",dpi=160);plt.close(fig)
        settings_note = f"Nominal results use {nominal.split('_n')[1]} statistical substeps per 250 ns bucket."

    assumptions = (
        "100 matched bb0nu events (IDs 0–99), 1 mm MC stepping, each gas and E/P point; 15 bar, 15% He. "
        "Centred 1.3 × 1.3 mm S14160-1315PS at zero gap, 5 mm hole diameter/depth, 10 mm pitch. "
        "Uniform isotropic axial EL emission; absorbing surroundings. TPB coats the sensor; complete VUV absorption, "
        "intrinsic conversion 0.62, forward fraction 0.5, blue PDE approximated by its 0.32 peak: effective PDE 0.0992. "
        "Prompt TPB stays in the parent bucket. 7,284 cells, gain 360,000, dark rate 120 kHz. "
        "Provisional recharge 12.5 ns with approved 25 ns sensitivity; uniform microcell illumination and constant PDE during recovery. "
        "ASIC: 1 mA, 250 ns, 3.5417 PE RMS (same input charge noise as the previous sensor). "
        "J = 0.030 at E/P = 2.0 and 0.015 at 2.5, for both gases. " + settings_note)
    correction = ("Widths are full-distribution RMS/mean. Known attachment is corrected event by event and totals are divided by "
        "the deposited MC energy; sensor/ASIC sums have expected stationary dark baseline subtracted. "
        "Recovery and clipping remain in the forward response. The detected-PE width isolates optical counting statistics; "
        "the sensor and ASIC widths expose nonlinear response. These are not resolutions after a fitted nonlinearity correction. "
        "The raw total distributions also include energy deposition and drift variation. Each pulse is one SiPM–250 ns cell. "
        "The logarithmic pulse plot displays positive sensor output; zero outputs remain in all event sums and pulse means.")
    sources = [
        ("Benson et al. (2018), intrinsic TPB QE, Figure 19", "https://doi.org/10.1140/epjc/s10052-018-5807-z"),
        ("Hamamatsu S14160 datasheet", "https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s14160-1310ps_etc_kapd1070e.pdf"),
        ("Klest and Zurek, measured pulse decay, slide 16", "https://indico.bnl.gov/event/20473/contributions/84622/attachments/51960/88885/BIC%2520SiPM%2520Choice%2520Intro.pdf")]
    headers = ["E/P", "Gas", "Mean detected PE [million]", "Sensor / detected [%]", "ASIC charge lost [%]", "Events clipping [%]", "Positive cells clipping [%]"]
    wh = ["E/P", "Gas", "VUV RMS [%]", "Detected PE RMS [%]", "Sensor RMS [%]", "ASIC RMS [%]"]
    ph = ["E/P","Gas","Event-mean PE: mean","Event-mean PE: median","Event-max PE: mean","Event-max PE: median"]
    sh = ["Recharge ns / steps","E/P","Gas","ASIC loss [%]","Sensor RMS [%]","ASIC RMS [%]"]
    ch = ["E/P","Gas","Mean sensor change [%]","Clip loss change [percentage points]"]
    source_html = "<ul>"+"".join(f'<li><a href="{u}">{t}</a></li>' for t,u in sources)+"</ul>"
    content = (f"<h1>S14160-1315PS + TPB: E/P pilot</h1><p>{assumptions}</p>"
        + "<h2>Charge collection and clipping</h2>" + table(headers,summary)
        + "<h2>Event-total widths</h2><p>"+correction+"</p>"+table(wh,widths)
        + '<img src="tpb_response.png" alt="Pulse charge and corrected response distributions">'
        + '<img src="tpb_totals.png" alt="Raw total detected and readout PE distributions">'
        + "<h2>Pulse means and maxima before clipping</h2>"+table(ph,pulses)
        + "<h2>Recharge sensitivity and time convergence</h2>"+table(sh,sensitivity)+table(ch,converged)
        + "<h2>TPB reference</h2><p>Figure 19 gives approximately 0.65 at 170 nm and 0.50 at 190 nm; "
          "0.62 is a rounded interpolation near 175 nm. The typical quoted QE uncertainty is about 8% relative. "
          "This intrinsic QE is separate from film-dependent external wavelength-shifting efficiency.</p>"+source_html
        + '<p><a href="../../latex_generated/electron_readout/sipm_tpb_pilot/sipm_tpb_pilot.pdf">PDF report</a></p>')
    (OUT/"sipm_tpb_pilot.html").write_text('<!doctype html><meta charset="utf-8"><title>SiPM TPB pilot</title>'
        '<style>body{max-width:1250px;margin:40px auto;padding:0 20px;font:16px system-ui;line-height:1.5}'
        'table{border-collapse:collapse;width:100%}th,td{padding:8px;border-bottom:1px solid #ddd;text-align:right}'
        'img{width:100%}h2{margin-top:2em}</style>'+content)
    for name in ("tpb_response.png","tpb_totals.png"):shutil.copy2(OUT/name,TEX/name)
    # Explicit substitutions keep the generated source portable to pdfLaTeX.
    latex_text=lambda s:s.replace("×",r"$\times$").replace("–","--").replace("%",r"\%").replace("_",r"\_")
    tex=r"""\documentclass{article}
\usepackage[a4paper,landscape,margin=16mm]{geometry}
\usepackage{graphicx,booktabs,hyperref}
\begin{document}
\section*{S14160-1315PS + TPB: E/P pilot}
"""+latex_text(assumptions)+"\n"+tex_table(headers,summary)+"\n"+latex_text(correction)+"\n"+tex_table(wh,widths)
    tex+=r"\clearpage\includegraphics[width=\textwidth,height=.88\textheight,keepaspectratio]{tpb_response.png}"
    tex+=r"\clearpage\includegraphics[width=\textwidth,height=.88\textheight,keepaspectratio]{tpb_totals.png}"
    tex+=r"\clearpage\section*{Pulse and sensitivity summaries}"+tex_table(ph,pulses)+tex_table(sh,sensitivity)+tex_table(ch,converged)
    tex+=r"\section*{Physical references}"+"\n"+r"\begin{itemize}"+"\n"
    tex+="\n".join(r"\item \href{"+u+"}{"+t+"}" for t,u in sources)
    tex+=r"\end{itemize}\end{document}"
    (TEX/"sipm_tpb_pilot.tex").write_text(tex)
    build_report(TEX / "sipm_tpb_pilot.tex")

if __name__ == "__main__":main()
