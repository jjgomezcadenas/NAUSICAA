#!/usr/bin/env python3
"""Present the Julia optical-geometry estimate saved in HDF5."""
from pathlib import Path
import html
import h5py

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sipm_geometry_estimate"
TEX = ROOT / "latex_generated/sipm_geometry_estimate"
TEX.mkdir(parents=True, exist_ok=True)
with h5py.File(OUT / "estimate.h5", "r") as f:
    assert f["metadata/complete"][()]
    rows = []
    for gap in (0, 10):
        g = f[f"gap_{gap}"]
        rows.append([str(gap), "Circle, diameter 5 mm",
                     f'{g["circle_acceptance_percent"][()]:.4f}', "100.00", "1.00"])
        for side in (3, 2, 1):
            r = g[f"square_{side}"]
            rows.append([str(gap), f"Square, {side} x {side} mm",
                         f'{r["acceptance_percent"][()]:.4f}',
                         f'{r["relative_to_circle_percent"][()]:.2f}',
                         f'{r["average_occupancy_ratio"][()]:.2f}'])

title = "SiPM size and distance: direct optical acceptance"
assumptions = ("Uniform isotropic emission along the axis of a 5 mm deep hole, "
               "5 mm hole diameter, absorbing walls, and a centred sensor parallel "
               "to the exit plane. Values describe light from one hole reaching its "
               "aligned sensor. Array totals additionally include neighbouring sensors. "
               "PDE and microcell pitch are held fixed.")
interpretation = ("The occupancy ratio is incident light per microcell relative to the "
                  "circular sensor at the same gap. At zero gap, reducing sensor size "
                  "reduces total light but increases illumination per microcell. At "
                  "10 mm gap, light approximately scales with sensor area. These "
                  "optical ratios precede SiPM recovery and electronics clipping.")
headers = ["Gap (mm)", "Sensor", "Acceptance (%)", "Relative light (%)", "Occupancy ratio"]
table = "<table><tr>" + "".join(f"<th>{x}</th>" for x in headers) + "</tr>"
table += "".join("<tr>" + "".join(f"<td>{html.escape(x)}</td>" for x in r) + "</tr>" for r in rows)
table += "</table>"
(OUT / "estimate.html").write_text(
    '<!doctype html><meta charset="utf-8"><title>' + title + '</title>'
    '<style>body{max-width:1050px;margin:40px auto;font:17px system-ui;line-height:1.5}'
    'td,th{padding:9px;border-bottom:1px solid #ccc;text-align:left}</style>'
    f"<h1>{title}</h1><p>{assumptions}</p>{table}<p>{interpretation}</p>"
    '<p>Julia depth quadrature was checked at 50,000 and 100,000 intervals, '
    'with absolute acceptance changes below 10⁻⁹.</p>')

tex_rows = "\n".join(" & ".join(r) + r" \\" for r in rows)
tex = r"""\documentclass{article}
\usepackage[margin=18mm]{geometry}
\usepackage{amsmath,amssymb,booktabs}
\begin{document}
\section*{SiPM size and distance: direct optical acceptance}
""" + assumptions + r"""

For a square of half-side $a$, gap $g$, and hole depth $h$, the accepted
fraction of all emitted photons is
\[
A_\square=\frac{1}{\pi h}\int_0^h
\arctan\!\left[\frac{a^2}{(g+u)\sqrt{(g+u)^2+2a^2}}\right]du.
\]
The square corners lie within the hole aperture for all sizes tabulated.
For the circular reference with radius $r$ equal to the hole radius,
\[
A_\circ=\frac{h+\sqrt{g^2+r^2}-\sqrt{(g+h)^2+r^2}}{2h}.
\]
\begin{center}\small
\begin{tabular}{rlrrr}\toprule
Gap (mm)&Sensor&Acceptance (\%)&Relative light (\%)&Occupancy ratio\\\midrule
""" + tex_rows + r"""
\bottomrule\end{tabular}\end{center}
""" + interpretation + r"""

Depth quadrature at 50,000 and 100,000 intervals agreed within $10^{-9}$
in absolute acceptance. Numerical results and card/source provenance are
saved in \texttt{results/sipm\_geometry\_estimate/estimate.h5}.
\end{document}
"""
(TEX / "estimate.tex").write_text(tex)
