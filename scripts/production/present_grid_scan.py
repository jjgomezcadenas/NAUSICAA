#!/usr/bin/env python3
"""Generate the 2/3/4 mm detector-production index as a LaTeX/PDF report."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tomllib

import h5py


def scalar(group, name):
    """Read one saved scalar with HDF5 byte strings decoded."""
    item = group[name][()]
    return item.decode() if isinstance(item, bytes) else item


def present(summary_path: Path, report_dir: Path) -> Path:
    """Format the Julia scan index; compile and check its LaTeX build log."""
    report_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(summary_path) as f:
        schema = f.attrs["schema"]
        schema = schema.decode() if isinstance(schema, bytes) else schema
        if schema != "nausicaa.detector_grid_scan.v1" or not scalar(f["metadata"], "complete"):
            raise ValueError("Incomplete or incompatible grid-scan summary")
        m = f["metadata"]
        r = tomllib.loads(scalar(m, "realisation_card"))
        e = tomllib.loads(scalar(m, "electron_card"))
        ion_cards = {pitch: tomllib.loads(scalar(m, f"ion_card_{pitch}mm"))
                     for pitch in (2, 3, 4)}
        lines = [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=20mm]{geometry}",
            r"\usepackage{booktabs}",
            r"\begin{document}",
            r"\title{Detector production: 8,000-event ion-grid scan}",
            r"\author{NAUSICAA generated report}",
            r"\date{}\maketitle",
            r"\section*{Configuration}",
            r"Each of the two 1 mm Geant4 samples uses one shared event placement "
            r"and one Fano-sampled primary-pair realisation. Its ions are collected "
            r"on three pad grids. The electron branch uses the same realization "
            r"and is propagated once per sample, producing separate pure-Xe and "
            r"Xe--15\% He incident-photon histories. The electron settings are "
            r"identical for all ion-grid choices.",
            r"\begin{center}\begin{tabular}{lr}\toprule Parameter & Value\\\midrule",
            f"TPC length & {r['chamber_length_cm']:.1f} cm" + r"\\",
            f"Pair-creation energy $W$ & {r['w_eV']:.2f} eV" + r"\\",
            f"Fano factor $F$ & {r['fano_factor']:.3f}" + r"\\",
            f"Electron lifetime & {e['lifetime_ms']:.1f} ms" + r"\\",
            f"Electron hole pitch / bucket & {e['hole_pitch_mm']:.1f} mm / "
            f"{e['bucket_width_us']:.1f} $\\mu$s" + r"\\",
            f"EL $E/P$ and $J$ & {e['reduced_field_kV_cm_bar']:.1f} kV/cm/bar; "
            f"{e['j_variance_ratio']:.3f}" + r"\\",
            r"\bottomrule\end{tabular}\end{center}",
            r"The $J=0.01$ operating-point input is the user's approximate "
            r"reading of Oliveira et al., Figure 6. The saved histories are "
            r"ions per pad/time frame and photons incident on each SiPM/time bucket. "
            r"SiPM PDE and electronics are subsequent stages.",
            r"\section*{Event and charge ledgers}",
            r"\begin{center}\small\begin{tabular}{lrrrrr}\toprule",
            r"Sample & Events & Pairs/ions & Surviving $e^-$ & Xe photons & Xe--He photons\\\midrule",
        ]
        for sample, label in (("bb0nu", r"bb0nu"),
                              ("single_electrons", r"single electrons")):
            g = f[sample]
            ions = [int(scalar(g["ions"][str(p)], "ions")) for p in (2, 3, 4)]
            pairs = int(scalar(g, "pairs"))
            if ions != [pairs] * 3:
                raise ValueError(f"Ion charge mismatch in {sample}")
            xe = int(scalar(g["photons/pure_xe"], "incident_photons"))
            he = int(scalar(g["photons/xe_he15"], "incident_photons"))
            lines.append(f"{label} & {int(scalar(g, 'events')):,} & {pairs:,} & "
                         f"{int(scalar(g, 'survivors')):,} & {xe:,} & {he:,}" + r"\\")
        lines += [
            r"\bottomrule\end{tabular}\end{center}",
            r"Every ion-grid file has exactly one collected ion per sampled pair. "
            r"Each ion-only replay verified every event's pair checksum against "
            r"the shared realization before collecting charge.",
            r"\section*{Ion-grid outputs}",
            r"\begin{center}\begin{tabular}{llrrrr}\toprule",
            r"Sample & Grid pitch & Frame & Ions & File size & Ion-only time\\\midrule",
        ]
        for sample, label in (("bb0nu", "bb0nu"),
                              ("single_electrons", "single electrons")):
            for pitch in (2, 3, 4):
                g = f[sample]["ions"][str(pitch)]
                frame = ion_cards[pitch]["bucket_width_ms"]
                if abs(float(scalar(g, "bucket_width_ms")) - frame) > 1e-8:
                    raise ValueError("Ion frame differs from its card")
                seconds = float(scalar(g, "elapsed_s"))
                runtime = "joint run" if pitch == 3 else f"{seconds:.1f} s"
                lines.append(f"{label} & {pitch} mm & {frame:.0f} ms & "
                             f"{int(scalar(g, 'ions')):,} & "
                             f"{float(scalar(g, 'file_bytes'))/1e6:.1f} MB & "
                             f"{runtime}" + r"\\")
        lines += [
            r"\bottomrule\end{tabular}\end{center}",
            r"The 3 mm histories were produced in the joint ion/electron run; "
            r"their separate ion-only runtime was not measured. The index HDF5 "
            r"stores full paths and source/configuration identities for every branch.",
            r"\section*{Provenance}",
            rf"Run ID: \texttt{{{scalar(m, 'run_id')}}}. "
            rf"Numerical index: \texttt{{{summary_path.name.replace('_', r'\_')}}}. "
            r"The index records the SHA-256 digest of the joint-run summary and "
            r"the resolved realisation, ion and electron TOML cards.",
            r"\end{document}",
        ]
    tex = report_dir / "grid_scan.tex"
    tex.write_text("\n".join(lines) + "\n")
    subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", tex.name],
                   cwd=report_dir, check=True, stdout=subprocess.DEVNULL)
    log = (report_dir / "grid_scan.log").read_text()
    if any(token in log for token in ("LaTeX Warning:", "Overfull ", "Underfull ",
                                     "Undefined control sequence")):
        raise RuntimeError("LaTeX build log contains warnings")
    subprocess.run(["latexmk", "-c", tex.name], cwd=report_dir, check=True,
                   stdout=subprocess.DEVNULL)
    return report_dir / "grid_scan.pdf"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("report_dir", type=Path)
    args = parser.parse_args()
    print(present(args.summary, args.report_dir))
