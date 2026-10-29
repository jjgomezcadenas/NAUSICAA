#!/usr/bin/env python3
"""Present saved Julia ion/PE grid readout summaries as figure and LaTeX PDF."""
from pathlib import Path
import argparse
import subprocess
import tomllib
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
STREAMS=[(sample,carrier) for sample in ("bb0nu","electrons")
         for carrier in ("ion","pure_xe","xe_he")]


def value(obj):
    item=obj[()]
    return item.decode() if isinstance(item,bytes) else item


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,default=ROOT/"results/grid_readout_direct_8000/grid_readout_summary.h5")
    parser.add_argument("--output",type=Path,default=ROOT/"latex_generated/grid_readout_direct")
    args=parser.parse_args()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    with h5py.File(args.input) as f:
        assert bool(f["summary_complete"][()])
        ec=tomllib.loads(value(f["electron_card_toml"]))
        ic=tomllib.loads(value(f["ion_card_toml"]))
        dc=tomllib.loads(value(f["detector_card_toml"]))
        n=int(f["events_per_stream"][()]);wall=float(f["wall_seconds"][()])
        workers=int(f["workers"][()])
        rows=[];ratios={};window_losses={}
        for sample,carrier in STREAMS:
            g=f[f"stream_event_totals/{sample}_{carrier}"]
            inputs=g["input_carriers"][()]
            signals=g["signal_count"][()]
            assert len(inputs)==len(signals)==n
            rows.append((sample,carrier,float(np.mean(inputs)),float(np.mean(signals)),
                         float(g["mean_signal_per_input_carrier"][()]),
                         float(np.mean(g["occupied_bucket_count"][()])),
                         float(g["pitch_mm"][()]),float(g["bucket_mm"][()]),
                         float(g["el_photons_per_electron"][()]),
                         float(g["axial_acceptance"][()])))
            window_losses[sample,carrier]=float(g["sum_outside_carriers"][()])/float(g["sum_sampled_carriers"][()])
            if carrier!="ion":ratios[sample,carrier]=signals/inputs
    fig,axes=plt.subplots(1,2,figsize=(9.2,3.2),sharey=True)
    for ax,sample in zip(axes,("bb0nu","electrons")):
        for carrier,label,color in (("pure_xe","Pure Xe","#3d65af"),
                                    ("xe_he","Xe + 15% He","#bd5a3b")):
            values=ratios[sample,carrier]
            ax.hist(values,bins=36,histtype="step",linewidth=1.8,label=label,color=color)
        ax.set_title("Double beta" if sample=="bb0nu" else "Single electron")
        ax.set_xlabel("Detected PE / incoming electron")
        ax.grid(alpha=.18)
    axes[0].set_ylabel("Events / bin")
    axes[1].legend(frameon=False,fontsize=8)
    fig.tight_layout()
    fig.savefig(output/"pe_per_electron.png",dpi=180)
    plt.close(fig)
    table=[]
    names={"ion":"ions","pure_xe":"pure Xe e$^{-}$","xe_he":"Xe--He e$^{-}$"}
    for sample,carrier,mean_input,mean_signal,ratio,pulses,pitch,bucket,yield_photons,acceptance in rows:
        table.append(f"{sample} & {names[carrier]} & {mean_input:,.0f} & "
                     f"{mean_signal:,.0f} & {ratio:.2f} & {pulses:,.0f} " + r"\\")
    el=ec["elcc"]
    pure_yield=next(row[8] for row in rows if row[0]=="bb0nu" and row[1]=="pure_xe")
    he_yield=next(row[8] for row in rows if row[0]=="bb0nu" and row[1]=="xe_he")
    direct_acceptance=next(row[9] for row in rows if row[0]=="bb0nu" and row[1]=="xe_he")
    ion_bucket_mm=float(ic["pads"]["z_bin_width_mm"])
    ion_speed_m_s=float(dc["gas"]["bulk_speed_m_s"])
    electron_bucket_mm=float(ec["readout"]["bucket_width_mm"])
    electron_speed_m_s=float(ec["readout"]["electron_drift_speed_m_s"])
    tex=r"""\documentclass{article}
\usepackage[margin=21mm]{geometry}
\usepackage{booktabs,graphicx}
\title{Direct ion pad and electron SiPM readout}
\author{NAUSICAA analysis}
\date{29 October 2026}
\begin{document}
\maketitle
"""+f"{workers} Julia processes sampled and collected {n:,} events from each MC sample in {wall:.2f}~s. " + r"""
Each retained sampled carrier was assigned directly to its physical pad or EL hole
and longitudinal bucket. Pure-Xe and Xe--15\%~He electron arrivals were
reduced to detected photoelectrons (PE) per SiPM and time bucket.
\begin{center}\small
\begin{tabular}{llrrrr}
\toprule
Sample & Readout & Input/event & Signal/event & Signal/input & Buckets/event\\
\midrule
"""+"\n".join(table)+r"""
\bottomrule
\end{tabular}\end{center}
\begin{center}\includegraphics[width=.94\linewidth]{pe_per_electron.png}\end{center}
\paragraph{Grid and timing inputs.}
Ion pads have """+f"{rows[0][6]:g}~mm pitch and {ion_bucket_mm:g}~mm " + r"""longitudinal buckets
"""+f"({ion_bucket_mm/ion_speed_m_s:g}~ms at {ion_speed_m_s:g}~m/s drift speed). " + r"""
Electron holes have """+f"{el['pitch_mm']:g}~mm pitch and {el['hole_diameter_mm']:g}~mm diameter, " + \
          f"with {el['hole_depth_mm']:g}~mm EL depth. " + \
          f"PE use {electron_bucket_mm:g}~mm buckets ({1000*electron_bucket_mm/electron_speed_m_s:g}~microseconds) " + r"""at the configured electron drift speed. The ion and electron
arrival orders are opposite. An independent uniform phase is used for each
physical grid; the electron phase is shared by the two gas models of an event.
"""+f"The EL reduced field is {el['reduced_field_kV_cm_bar']:g}~kV/cm/bar at " + \
          f"{ec['gas']['pressure_bar']:g}~bar, with gate/anode potentials " + \
          f"{el['entry_gate_voltage_V']:g}/{el['exit_anode_voltage_V']:g}~V " + r"""in a relative voltage gauge.
\paragraph{Optical model.}
The fitted mean EL yields from the Xe--He study are """+ \
          f"{pure_yield:.1f} and {he_yield:.1f} " + r"""photons per electron
for pure Xe and Xe--15\%~He. Photons are isotropic along the hole axis;
VUV rays touching the walls are absorbed. The analytic direct-light
acceptance is """+f"{direct_acceptance:.6f}"+r""" and the configured SiPM PDE is """ + \
          f"{100*ec['readout']['sipm_pde']:.0f}\\%. " + r"""The fitted yield is held at its mean;
PE counts are sampled with an exact binomial draw from the combined optical
acceptance and PDE. PE are placed in the parent electron's arrival bucket;
in-hole transit delay and electronics noise are outside this readout stage.
The Xe--He EL fits were measured near 1.2~bar; their use at 15~bar applies
reduced-field scaling.
\paragraph{Validation.}
The ion signal equals the retained input ion count in every event. Saved
grid indices and time buckets account for each signal; pure-Xe and Xe--He
readouts share a transverse grid phase per event. The finite diffusion window
is the source-hit extent plus """+f"{ic['transport']['nsigma']:g}" + r""" diffusion widths
on each side. The fractions outside that window in the two ion streams are
"""+f"{100*window_losses['bb0nu','ion']:.4f}\\% and {100*window_losses['electrons','ion']:.4f}\\%, respectively. " + r"""The complete event-wise validation results accompany
the HDF5 run summary.
\end{document}
"""
    (output/"grid_readout.tex").write_text(tex)
    with (output/"latexmk.stdout").open("w") as log:
        subprocess.run(["latexmk","-pdf","-interaction=nonstopmode","-halt-on-error",
                        "grid_readout.tex"],cwd=output,check=True,stdout=log,stderr=subprocess.STDOUT)
    assert "! LaTeX Error" not in (output/"grid_readout.log").read_text()


if __name__=="__main__":main()
