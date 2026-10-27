#!/usr/bin/env python3
"""Present the saved Julia pixel-area calculation for the manuscript.

Run from any directory with Python containing h5py and matplotlib.
No physical calculation is performed here; lengths/areas come from analysis.h5.
"""
from pathlib import Path
import hashlib
import json
import os

import h5py

HERE = Path(__file__).resolve()
NOTE = HERE.parents[1]
ROOT = HERE.parents[3]
os.environ.setdefault("MPLCONFIGDIR", str(NOTE / "work/matplotlib-cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch


def main():
    source = ROOT / "results/pixel_area/analysis.h5"
    out = NOTE / "generated"
    with h5py.File(source) as h:
        # Reject stale results if any input to the physical calculation changed.
        for path, digest in zip(h["source_paths"].asstr()[:], h["source_sha256"].asstr()[:]):
            assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest, path
        a = dict(h.attrs)
        g = h["wide_margins"]
        b = dict(g.attrs)
        r = h["regions"]
        areas = dict(zip(r["names"].asstr()[:], r["area_m2"][:]))
        devices = dict(zip(h["devices/names"].asstr()[:],
                           zip(h["devices/W_m"][:], h["devices/L_m"][:])))
        colors = {"Pad and clearance": "#dce4ec", "Shared services": "#dbccec",
                  "Buffer and selection": "#d4e7d1", "Control and routing": "#f3dfb7"}
        labels = {"Pad and clearance": "Pad clearance", "Shared services": "Shared circuitry allocation",
                  "Buffer and selection": "Buffer and selection reserve", "Control and routing": "Control and routing reserve"}
        with plt.rc_context({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}):
            fig, (ax, detail) = plt.subplots(1, 2, figsize=(9, 6),
                                           gridspec_kw={"width_ratios": [2.5, 1]}, layout="constrained")
            pitch = a["pitch_m"] * 1e3
            side = a["envelope_side_m"] * 1e3
            offset = (pitch - side) / 2
            ax.add_patch(Rectangle((0, 0), pitch, pitch, fill=False, edgecolor="0.2"))
            for name, x, y, w, height in zip(r["names"].asstr()[:], r["x_m"][:], r["y_m"][:], r["w_m"][:], r["h_m"][:]):
                if name == "CSA strip":
                    continue
                ax.add_patch(Rectangle((offset + x * 1e3, offset + y * 1e3), w * 1e3, height * 1e3,
                                       facecolor=colors[name], edgecolor="white", linewidth=1))
            pad = a["pad_width_m"] * 1e3
            ax.add_patch(Rectangle(((pitch-pad)/2, (pitch-pad)/2), pad, pad,
                                   facecolor="#626d76", edgecolor="black"))
            ax.text(pitch/2, pitch/2, f"Pad\n{pad:g} mm", ha="center", va="center", color="white")
            cw, ch = b["csa_width_m"]*1e3, b["csa_height_m"]*1e3
            ax.add_patch(Rectangle((offset+side-cw, offset+(side-ch)/2), cw, ch,
                                   facecolor="#35779b", edgecolor="white"))
            ax.text(offset+side-cw/2, pitch/2, "CSA", rotation=90, ha="center", va="center", color="white")
            ax.add_patch(Rectangle((offset, offset), side, side, fill=False, ls="--", edgecolor="black"))
            ax.text(pitch/2, offset-.12, f"{a['envelope_area_m2']*1e6:g} mm² envelope", ha="center", va="top")
            ax.set(xlim=(-.05, pitch+.05), ylim=(-.05, pitch+.05), xlabel="x (mm)", ylabel="y (mm)", title="(a) Pixel cell")
            ax.set_aspect("equal")
            ax.legend(handles=[Patch(facecolor=colors[n], label=labels[n]) for n in colors],
                      loc="upper center", bbox_to_anchor=(.5, -.17), frameon=False, ncol=2, fontsize=8)
            margin = b["device_margin_m"]
            for name, x, y, w, height in zip(g["names"].asstr()[:], g["x_m"][:], g["y_m"][:], g["w_m"][:], g["h_m"][:]):
                detail.add_patch(Rectangle((x*1e6, y*1e6), w*1e6, height*1e6,
                                           facecolor="#dce4ec", edgecolor="#35779b"))
                if name in devices:
                    dw, dl = devices[name]
                    detail.add_patch(Rectangle(((x+margin)*1e6, (y+margin)*1e6), min(dw,dl)*1e6,
                                               max(dw,dl)*1e6, facecolor="#b46b19"))
                label = name.replace("_", ".").replace("Cf reserve", "Cf\nreserve")
                detail.text((x+w/2)*1e6, (y+height/2)*1e6, label, ha="center", va="center", fontsize=7)
            detail.set(xlim=(0, b["csa_width_m"]*1e6), ylim=(0, b["csa_height_m"]*1e6),
                       xlabel="x (µm)", ylabel="y (µm)", title="(b) CSA packing")
            detail.set_aspect("equal")
            fig.savefig(out/"pixel_area.pdf", metadata={"CreationDate": None, "ModDate": None})
            plt.close(fig)
        rows = [("Pad and clearance", areas["Pad and clearance"]),
                ("CSA and feedback", b["csa_area_m2"]),
                ("Buffer and selection reserve", areas["Buffer and selection"]),
                ("Control and routing reserve", areas["Control and routing"]),
                ("Shared circuitry allocation", areas["Shared services"])]
        table = [r"\begin{tabular}{lr}", r"\toprule", r"Block & Area (mm$^2$) \\", r"\midrule"]
        table += [f"{name} & {value*1e6:.2f}" + r" \\" for name, value in rows]
        table += [r"\midrule", f"Total & {b['budget_used_m2']*1e6:.2f}" + r" \\", r"\bottomrule", r"\end{tabular}"]
        (out/"pixel_area_table.tex").write_text("\n".join(table)+"\n")
        values = {"PixelCellArea": a["cell_area_m2"]*1e6,
                  "PixelBudgetArea": a["envelope_area_m2"]*1e6,
                  "PixelClearanceUm": a["pad_clearance_m"]*1e6,
                  "PixelDeviceMarginUm": margin*1e6,
                  "PixelOuterMarginUm": b["outer_margin_m"]*1e6,
                  "PixelRowGapUm": b["row_gap_m"]*1e6,
                  "PixelCsaWidthUm": b["csa_width_m"]*1e6,
                  "PixelCsaHeightUm": b["csa_height_m"]*1e6}
        (out/"pixel_area_values.tex").write_text("% Generated by figures/build_pixel_area.py\n" +
              "".join(f"\\newcommand{{\\{name}}}{{{value:g}}}\n" for name,value in values.items()))
    provenance = {"source": str(source.relative_to(ROOT)), "scenario": "wide_margins",
                  "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                  "generator": str(HERE.relative_to(ROOT)),
                  "generator_sha256": hashlib.sha256(HERE.read_bytes()).hexdigest()}
    (out/"pixel_area_provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")


if __name__ == "__main__":
    main()
