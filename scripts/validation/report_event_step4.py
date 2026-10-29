"""Present paired scan and ensemble results calculated in Julia and saved in HDF5."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import json
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"latex_generated/event_analysis_step4"


def render():
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    with h5py.File(ROOT/"results/event_analysis/step4_scans.h5") as f:
        status=f.attrs["status"]
        if isinstance(status,bytes):status=status.decode()
        if status!="complete":raise ValueError("Scan output is incomplete")
        ids=f["variants/id"].asstr()[:].tolist()
        labels=[s.replace("_"," ") for s in ids]
        summary={"passed":int(f["tests/passed"][()]),"variants":ids,"samples":[]}
        fig,axes=plt.subplots(2,1,figsize=(11,7),layout="constrained",sharex=True)
        mapfig,mapaxes=plt.subplots(2,2,figsize=(10,9),layout="constrained")
        profilefig,profileaxes=plt.subplots(2,1,figsize=(10,6),layout="constrained",sharex=True)
        table_sections=[]
        baseline_rows=[]
        colours=["#2166ac","#b35806"]
        peak=max(float(np.nanmax(s["baseline/ensemble/xy_mean_accepted"][:])) for s in f["samples"].values())
        norm=LogNorm(vmin=0.01,vmax=max(1,peak))
        for index,(name,sg) in enumerate(f["samples"].items()):
            title="Double beta" if name=="0nubb" else "Electron background"
            s=sg["summary"]
            record={"sample":title,"rows":[]}
            for i,id in enumerate(ids):
                row={"variant":id}
                for field in s.keys():
                    if field=="variant":continue
                    if field=="resolution_reference":row[field]=s[field].asstr()[i]
                    else:
                        value=s[field][i].item()
                        row[field]=None if isinstance(value,float) and not np.isfinite(value) else value
                record["rows"].append(row)
            summary["samples"].append(record)
            x=np.arange(len(ids))
            axes[0].plot(x,100*s["acceptance_fraction"][:],"o-",color=colours[index],label=title)
            axes[1].plot(x,100*s["mean_shared_pad_distance"][:],"o--",alpha=0.65,color=colours[index],label=title+", full pipeline")
            axes[1].plot(x,100*s["mean_fixed_grid_pad_distance"][:],"s-",color=colours[index],label=title+", fixed grid")
            e=sg["baseline/ensemble"]
            for col,field in enumerate(("xy_mean_accepted","xy_mean_processed")):
                z=e[field][:]
                m=mapaxes[index,col].pcolormesh(e["x_edges_mm"][:],e["y_edges_mm"][:],np.ma.masked_where(z<=0,z),
                    shading="flat",cmap="viridis",norm=norm,rasterized=True)
                suffix="per accepted event" if col==0 else "per processed event"
                mapaxes[index,col].set(xlabel="Relative pad x [mm]",ylabel="Relative pad y [mm]",title=title+"\n"+suffix)
                mapaxes[index,col].set_aspect("equal")
            profileaxes[0].stairs(e["arrival_mean_accepted"][:],e["time_edges_ms"][:],color=colours[index],label=title)
            profileaxes[1].stairs(e["arrival_mean_processed"][:],e["time_edges_ms"][:],color=colours[index],label=title)
            base=record["rows"][0]
            baseline_rows.append(f"{title} & {base['accepted']}/{base['processed']} & {base['mean_ions_accepted']:.1f} & {base['mean_selected_ions_processed']:.1f} \\\\")
            rows=[]
            for row in record["rows"]:
                texid=row["variant"].replace("_",r"\_")
                number=lambda v:"--" if v is None else f"{100*v:.2f}"
                rows.append(f"{texid} & {row['accepted']} & {row['lost']} & {row['gained']} & {row['shared']} & "
                    f"{number(row['mean_shared_pad_distance'])} & {number(row['mean_fixed_grid_pad_distance'])} \\\\")
            table_sections.append(r"\subsection*{"+title+"}\n"+r"\begin{center}\small\begin{tabular}{lrrrrrr}\toprule"+"\n"+
                r"Variant & Pass & Lost & Gained & Shared & Full [\%] & Fixed [\%] \\\midrule"+"\n"+
                "\n".join(rows)+"\n"+r"\bottomrule\end{tabular}\end{center}")
        axes[0].set(ylabel="Accepted events [%]",title="Paired selection scan: same events and ion draws")
        axes[0].legend(frameon=False)
        axes[1].set(ylabel="Mean paired pad/slice distance [%]",xticks=np.arange(len(ids)))
        axes[1].set_xticklabels(labels,rotation=40,ha="right")
        axes[1].legend(frameon=False,ncol=2,fontsize=9)
        fig.savefig(OUT/"paired_scans.png",dpi=170);plt.close(fig)
        mapfig.colorbar(m,ax=mapaxes.ravel().tolist(),label="Mean ions per 2 mm pad",shrink=0.7)
        mapfig.savefig(OUT/"ensemble_xy.png",dpi=170);plt.close(mapfig)
        profileaxes[0].set(ylabel="Mean ions / slice / accepted event",title="First-arrival alignment; 20 ms slices")
        profileaxes[1].set(ylabel="Selected ions / slice / processed event",xlabel="Relative time [ms]")
        for ax in profileaxes:ax.legend(frameon=False)
        profilefig.savefig(OUT/"ensemble_time.png",dpi=170);plt.close(profilefig)
        (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
        resolution_rows=[]
        for sample in summary["samples"]:
            for row in sample["rows"]:
                if row["variant"] not in ("bins_200","bins_400"):continue
                value=row["mean_resolution_fixed_grid_distance"]
                distance="--" if value is None else f"{100*value:.2f}"
                reference_tex=row['resolution_reference'].replace('_',r'\_')
                variant_tex=row['variant'].replace('_',r'\_')
                resolution_rows.append(f"{sample['sample']} & {reference_tex} & "
                    f"{variant_tex} & {row['resolution_shared']} & {distance} \\\\")
        text=r"""\documentclass[11pt]{article}
\usepackage[margin=18mm]{geometry}
\usepackage{booktabs,graphicx,hyperref}
\begin{document}
\title{Event analysis: paired scans and ensemble means}
\author{NAUSICAA}\date{}\maketitle
This checkpoint retains the approved baseline and evaluates numerical/selection
changes on the same events. Julia computes all counts, decisions, paired metrics
and ensemble reductions. The cumulative core suite passes TESTS assertions.

\section*{Ensemble convention and denominators}
The pad containing the selected track's energy-weighted x/y reference becomes
the central pad, with centre at (0,0). This relabels existing physical pads and
preserves integer counts. The relative arrival axis starts at each history's
first saved arrival edge. The configured direction is retained. Every accepted
event contributes its full retained diffusion charge; all cells and zero-count
intervals are represented by the saved axes.

The same exact integer numerator is divided by two explicit denominators:
accepted events, and all processed events (with rejected events contributing
zero selected charge). The means use Float32. With no accepted events, the
conditional mean is undefined (NaN) and the processed mean is zero. The HDF5
stores numerator, both denominators, both means, axes and alignment metadata.
\begin{center}\small\begin{tabular}{lrrr}\toprule
Sample & Accepted/processed & Ions/accepted event & Selected ions/processed event \\\midrule
BASELINE
\bottomrule\end{tabular}\end{center}
\includegraphics[width=\linewidth]{ensemble_xy.png}
\includegraphics[width=\linewidth]{ensemble_time.png}

\section*{Paired scan design}
Each variation changes one setting family from the baseline: diffusion histogram
bin count, Gaussian padding, kNN neighbour count, ion threshold or a common
multiplier of both the $3\sigma_T$ bin width and $6\sigma_T$ edge cutoff. Drift,
ion draws and grid offsets are paired by saved stage seeds. Graph-only variants
reuse the same diffusion table. Changing histogram bin count at fixed bounds
preserves every retained ion; this is checked for each event.

Lost and gained refer to the same baseline event decisions, so equal acceptance
counts cannot conceal exchanged events. Shared counts give the denominator for
the paired pad metric. This subset is a validation sample, not a full-sample
efficiency determination. Baseline cuts remain unchanged.

\section*{Pad-response comparison}
For paired count tables A and B, use $D=\sum|A-B|/(\sum A+\sum B)$: zero means
identical charge allocation; one means disjoint allocation. The Full comparison
uses each variant's reconstructed grid reference and first-arrival alignment.
The Fixed comparison collects both diffusion tables on the baseline accepted
track's physical lattice, using signed cell indices to eliminate coverage and
time-origin relabelling. Graph-only changes therefore have zero Fixed distance.
These distances measure redistribution, not lost total charge.
\includegraphics[width=\linewidth]{paired_scans.png}
TABLES

\section*{Successive histogram resolution checks}
The following comparisons use the same baseline physical lattice and events
accepted in the baseline and both compared settings. The 200-to-400 check was
added after the first scan showed a measurable pad response to resolution.
\begin{center}\small\begin{tabular}{lllrr}\toprule
Sample & Reference & Variant & Shared & Fixed distance [\%] \\\midrule
RESOLUTION
\bottomrule\end{tabular}\end{center}
Selection stability and pad-response convergence are evaluated separately.
The results and the required accuracy for pad histories should be reviewed
before selecting production settings and running the complete samples.

\section*{Validation and reproduction}
Tests cover translation-invariant central-pad alignment, both denominators,
missing-event cells as zero, complete axis coverage, integer conservation,
undefined conditional means for zero acceptance, Float32 reduction precision,
HDF5 round trips and corruption detection. Every written ensemble is read back;
its exact numerator matches the persisted aligned histories.
\begin{verbatim}
make validate-event-step4 PYTHON=/path/to/python3.11-or-newer
# Presentation from the saved numerical output:
make report-event-step4 PYTHON=/path/to/python3.11-or-newer
\end{verbatim}
Numerical data, event decisions, aligned histories, paired metrics and provenance:
\par\noindent\path{results/event_analysis/step4_scans.h5}.
\par The ensemble subgroups use schema \texttt{nausicaa.event.ensemble.v1}.
\section*{Configuration used}
\footnotesize\begin{verbatim}
CONFIG
\end{verbatim}
\end{document}
"""
        for marker,value in {"TESTS":str(summary["passed"]),"BASELINE":"\n".join(baseline_rows),
            "TABLES":"\n".join(table_sections),"RESOLUTION":"\n".join(resolution_rows),
            "CONFIG":f["configuration/card"].asstr()[()]}.items():text=text.replace(marker,value)
        (OUT/"report.tex").write_text(text)
    print(OUT/"report.tex")


if __name__=="__main__":render()
