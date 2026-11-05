#!/usr/bin/env python3
"""Render the saved 250 ns electron/photon run as a self-contained TeX/PDF."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tomllib

import h5py


ROOT = Path(__file__).resolve().parents[2]


def value(group: h5py.Group, key: str):
    """Read a scalar dataset, decoding an HDF5 string if present."""
    result = group[key][()]
    return result.decode() if isinstance(result, bytes) else result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    summary = args.summary.resolve()
    with h5py.File(summary) as source:
        if source.attrs["schema"] != b"nausicaa.electron_rebin_summary.v1":
            raise ValueError("Unexpected electron-rebin summary schema")
        if not bool(value(source["metadata"], "complete")):
            raise ValueError("Electron-rebin summary is incomplete")
        metadata = source["metadata"]
        run_id = value(metadata, "run_id")
        electron = tomllib.loads(value(metadata, "electron_card"))
        realisation = tomllib.loads(value(metadata, "realisation_card"))
        events = int(value(metadata, "events_per_sample"))
        bucket_us = float(value(metadata, "bucket_width_us"))
        bucket_mm = float(value(metadata, "bucket_width_mm"))
        rows = []
        for sample, label in (("bb0nu", r"$0\nu\beta\beta$"),
                              ("single_electrons", "Single electrons")):
            group = source[sample]
            if int(value(group, "events")) != events:
                raise ValueError(f"Event count mismatch for {sample}")
            pairs = int(value(group, "initial_pairs"))
            survivors = int(value(group, "surviving_electrons"))
            if survivors > pairs:
                raise ValueError(f"Electron ledger mismatch for {sample}")
            for gas, title in (("pure_xe", "Pure Xe"), ("xe_he15", "Xe--15\\% He")):
                emitted = int(value(group, f"{gas}_emitted"))
                incident = int(value(group, f"{gas}_incident"))
                if incident > emitted:
                    raise ValueError(f"Photon ledger mismatch for {sample}/{gas}")
                rows.append(
                    f"{label} & {title} & {events:,} & {pairs:,} & "
                    f"{survivors:,} & {emitted:,} & {incident:,} " + r"\\"
                )
    out = (args.output or ROOT / "latex_generated/electron_readout/electron_250ns" / run_id).resolve()
    out.mkdir(parents=True, exist_ok=True)
    tex = r"""\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm]{geometry}
\usepackage{booktabs,longtable,hyperref}
\title{Electron photon histories with 250 ns buckets}
\author{NAUSICAA production analysis}
\date{\today}
\begin{document}\maketitle
\section*{Run and detector settings}
The detector-stage electron branch was replayed from the saved common
primary-pair realisation. The two samples use the same source events,
placement, pair draws, and attachment streams as the earlier 2~$\mu$s run.
The ion histories were retained. Pure Xe and Xe--15\% He are alternative gas
models for each event. The HDF5 inputs to the SiPM response are VUV photons
incident on each SiPM per time bucket, before PDE or electronics.

\begin{center}\begin{tabular}{lr}\toprule
Quantity & Value\\\midrule
Electron time bucket & """ + f"{bucket_us:g}" + r"""~$\mu$s\\
Longitudinal bucket width & """ + f"{bucket_mm:g}" + r"""~mm\\
Electron drift speed & """ + f"{electron['drift_speed_mm_us']:g}" + r"""~mm/$\mu$s\\
Pressure & """ + f"{electron['pressure_bar']:g}" + r"""~bar\\
Helium fraction in mixture & """ + f"{100*electron['helium_fraction']:g}" + r"""\%\\
EL hole pitch, diameter, depth & """ + (
        f"{electron['hole_pitch_mm']:g}, {electron['hole_diameter_mm']:g}, "
        f"{electron['hole_depth_mm']:g}"
    ) + r"""~mm\\
Reduced EL field & """ + f"{electron['reduced_field_kV_cm_bar']:g}" + r"""~kV/cm/bar\\
Electron lifetime & """ + f"{electron['lifetime_ms']:g}" + r"""~ms\\
Primary $W$, Fano $F$ & """ + f"{realisation['w_eV']:g}~eV, {realisation['fano_factor']:g}" + r"""\\
EL variance ratio $J$ & """ + f"{electron['j_variance_ratio']:g}" + r"""\\
\bottomrule\end{tabular}\end{center}

\section*{Production ledgers}
\small
\begin{longtable}{llrrrrr}\toprule
Sample & Gas & Events & Initial pairs & Surviving $e^-$ & Emitted photons & Incident photons\\\midrule
""" + "\n".join(rows) + r"""
\bottomrule\end{longtable}
\normalsize
Every event replays and checks its primary-pair checksum against the saved
realisation. The electron/photon ledger requires surviving electrons no
greater than initial pairs and incident photons no greater than emitted
photons. Photon groups retain event, gas, SiPM and signed bucket indices.

\section*{Timing and optical assumptions}
Electrons are diffused directly into the physical 10~mm SiPM grid and
""" + f"{bucket_mm:g}" + r"""~mm longitudinal buckets. EL yield and its $J$ fluctuation,
electron attachment and analytic direct optical acceptance follow the
resolved detector card and \texttt{docs/production\_formulas.tex}.
Photon time is assigned to the parent electron's arrival bucket at the EL
hole. The model does not resolve electron flight or photon emission time
within the 5~mm hole. The finer 250~ns bucket therefore resolves the saved
arrival-coordinate model; it does not add a measured EL time profile.

The numerical source and full provenance are stored in the HDF5 summary
under the production pilot results directory. Its filename is
\texttt{electron\_250ns\_summary\_} followed by the run identity
\texttt{""" + run_id + r"""} and \texttt{.h5}. Large branch histories remain
local beside that summary.
\end{document}
"""
    tex_path = out / "electron_250ns.tex"
    tex_path.write_text(tex)
    with (out / "build.log").open("w") as log:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode",
                        "-halt-on-error", tex_path.name], cwd=out, check=True,
                       stdout=log, stderr=subprocess.STDOUT)
    print(f"TeX: {tex_path}")
    print(f"PDF: {out / 'electron_250ns.pdf'}")


if __name__ == "__main__":
    main()
