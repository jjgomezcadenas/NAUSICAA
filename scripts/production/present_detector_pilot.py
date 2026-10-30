#!/usr/bin/env python3
"""Present saved detector-pilot reductions in a generated LaTeX/PDF report."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tomllib

import h5py


def scalar(group: h5py.Group, name: str):
    value = group[name][()]
    value = value.item() if hasattr(value, "item") else value
    return value.decode() if isinstance(value, bytes) else value


def num(value: float, decimals: int = 3) -> str:
    return f"{value:,.{decimals}f}"


def make_report(summary_path: Path, report_dir: Path) -> Path:
    report_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(summary_path, "r") as f:
        schema = f.attrs["schema"]
        schema = schema.decode() if isinstance(schema, bytes) else str(schema)
        if schema != "nausicaa.detector_pilot_summary.v1" or not scalar(f["metadata"], "complete"):
            raise ValueError("Incomplete or incompatible pilot summary")
        r = tomllib.loads(scalar(f["metadata"], "realisation_card"))
        i = tomllib.loads(scalar(f["metadata"], "ion_card"))
        e = tomllib.loads(scalar(f["metadata"], "electron_card"))
        d = f["derived"]
        rows = []
        for sample in ("bb0nu", "single_electrons"):
            if sample not in f:
                continue
            g = f[sample]
            if not scalar(g, "ion_charge_conserved") or not scalar(g, "photon_ledgers_valid"):
                raise ValueError(f"Failed detector ledger: {sample}")
            rows.append((sample, *(int(scalar(g, key)) for key in (
                "events", "pairs", "ions", "survivors", "pure_xe_incident", "xe_he15_incident"))))
        run_id = scalar(f["metadata"], "run_id")
        code_hash = scalar(f["metadata"], "code_hash")
        peak_rss = int(scalar(f["metadata"], "peak_rss_bytes"))
        lines = [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=20mm]{geometry}",
            r"\usepackage{amsmath,booktabs,longtable}",
            r"\begin{document}",
            r"\title{Detector production pilot: formulas and results}",
            r"\author{NAUSICAA generated report}",
            r"\date{}",
            r"\maketitle",
            r"\noindent This calculation uses the recorded parameter cards. "
            r"The approximate $J$ values at $E/P=1.5, 2, 3$ are user readings of "
            r"Oliveira et al., Figure 6; the form requires an explicit $J$ at intermediate fields.",
            r"\section*{Resolved inputs}",
            r"\begin{tabular}{ll}\toprule Parameter & Value\\\midrule",
            f"TPC length & {num(r['chamber_length_cm'], 1)} cm\\\\",
            f"Pressure & {num(e['pressure_bar'], 1)} bar\\\\",
            f"Ion drift field & {num(i['drift_field_V_cm'], 1)} V/cm\\\\",
            f"Ion pitch, frame & {num(i['pad_pitch_mm'], 1)} mm, {num(i['bucket_width_ms'], 1)} ms\\\\",
            f"Electron hole pitch, bucket & {num(e['hole_pitch_mm'], 1)} mm, {num(e['bucket_width_us'], 1)} $\\mu$s\\\\",
            f"Electron lifetime & {num(e['lifetime_ms'], 1)} ms\\\\",
            f"EL $E/P$ & {num(e['reduced_field_kV_cm_bar'], 2)} kV/cm/bar\\\\",
            f"Pair-creation energy $W$, Fano factor $F$ & {num(r['w_eV'], 2)} eV, {num(r['fano_factor'], 3)}\\\\",
            f"EL variance ratio $J$ & {num(e['j_variance_ratio'], 3)}\\\\",
            f"$J$ source & {e['j_status'].replace('_', r'\_')}\\\\",
            r"\bottomrule\end{tabular}",
            r"\section*{Calculations and assumptions}",
            r"\begin{align*}",
            r"\mu_i &= 10^6 E_i/W, &\operatorname{Var}(N_i)&\simeq F\mu_i,\\",
            r"L_e &= L_{\max}-L_i, & L_i&\sim U(L_{\min},L_{\max}),\\",
            r"\sigma_{i,T}=\sigma_{i,L}&=10\sqrt{2(k_B/e)TL_i/E_d}\ {\rm mm},\\",
            r"\sigma_{e,T/L}&=d_{T/L}\sqrt{L_e/P}, &p_{\rm survive}&=e^{-L_e/(v_e\tau_e)},\\",
            r"Y_g&=a_g\max(E/P-x_{0,g},0)Ph, &\operatorname{Var}(N_{\rm EL}|n_e)&=n_e JY_g,\\",
            r"A&=(h+r-\sqrt{h^2+r^2})/(2h), &N_{\rm incident}|N_{\rm EL}&\sim\mathrm{Binomial}(N_{\rm EL},A).",
            r"\end{align*}",
            r"At each source hit, one integer pair draw creates both initial carrier populations. "
            r"The minimum integer variance limits the Fano target at low occupancy. "
            r"Electron attachment is sampled once per pair realisation and shared by the two gas alternatives. "
            r"Diffused carriers enter physical grid cells directly; there is no diffusion-window charge loss. "
            r"VUV photons touching a hole wall are absorbed; a circular SiPM covers its exit.",
            r"\par\medskip",
            r"\begin{tabular}{lrr}\toprule Derived quantity & Pure Xe & Xe--15\% He\\\midrule",
            f"Max transverse electron $\\sigma$ (mm) & {num(scalar(d['pure_xe'], 'sigma_t_max_mm'))} & {num(scalar(d['xe_he15'], 'sigma_t_max_mm'))}\\\\",
            f"Max longitudinal electron $\\sigma$ (mm) & {num(scalar(d['pure_xe'], 'sigma_l_max_mm'))} & {num(scalar(d['xe_he15'], 'sigma_l_max_mm'))}\\\\",
            f"Mean EL yield (photons/electron) & {num(scalar(d['pure_xe'], 'yield_photons_per_e'), 1)} & {num(scalar(d['xe_he15'], 'yield_photons_per_e'), 1)}\\\\",
            r"\bottomrule\end{tabular}",
            f"\\par Ion width at maximum drift: {num(scalar(d, 'ion_sigma_max_mm'))} mm. "
            f"Ion frame: {num(scalar(d, 'ion_bucket_mm'), 1)} mm. "
            f"Electron bucket: {num(scalar(d, 'electron_bucket_mm'), 1)} mm. "
            f"Minimum electron survival: {num(scalar(d, 'electron_min_survival'), 4)}. "
            f"Direct optical acceptance: {num(scalar(d, 'optical_acceptance'), 5)}. "
            f"Gate/anode potentials relative to cathode: {num(scalar(d, 'gate_voltage_V') / 1000, 1)} / "
            f"{num(scalar(d, 'anode_voltage_V') / 1000, 1)} kV.",
            r"\section*{Saved event totals and validation}",
            r"\begin{center}\small\begin{tabular}{lrrrrrr}\toprule Sample & Events & Pairs & Ions & Surviving $e^-$ & Xe photons & Xe--He photons\\\midrule",
        ]
        for sample, events, pairs, ions, survivors, xe, he in rows:
            label = sample.replace("_", r"\_")
            lines.append(f"{label} & {events:,} & {pairs:,} & {ions:,} & {survivors:,} & {xe:,} & {he:,}\\\\")
        lines += [
            r"\bottomrule\end{tabular}\end{center}",
            r"\par The ion ledger is exact for both samples: ions in pad/time cells equal initial pairs. "
            r"For each gas, incident photons do not exceed emitted photons. "
            r"The three branch files store matching per-event pair checksums.",
            r"\section*{Provenance and scope}",
            f"Run ID: \\texttt{{{run_id}}}. Code SHA-256 prefix: \\texttt{{{code_hash[:16]}}}. "
            r"The complete digest is in the HDF5 summary.",
            f"Peak process RSS: {peak_rss / 1e6:.1f} MB. Source hit tables are the 1 mm Geant4 samples.",
            r"Ion thermal diffusion follows ITACA Revisited, Eq. (2.1) (arXiv:2604.02443). "
            r"The pure-Xe transverse coefficient is its Table 1 input; other diffusion coefficients "
            r"are recorded model inputs. EL yield coefficients follow Fernandes et al., JHEP 04 (2020) 034, Table 1. "
            r"The definition of $J$ follows Oliveira et al., JINST 6 (2011) P05007, Eq. (4.6); "
            r"the recorded field-point values are approximate user readings of its Figure 6.",
            r"\end{document}",
        ]
    output = report_dir / "formulas.tex"
    output.write_text("\n".join(lines) + "\n")
    subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", output.name],
                   cwd=report_dir, check=True, stdout=subprocess.DEVNULL)
    log = (report_dir / "formulas.log").read_text()
    if any(token in log for token in ("LaTeX Warning:", "Overfull ", "Underfull ", "Undefined control sequence")):
        raise RuntimeError("LaTeX build log contains warnings")
    subprocess.run(["latexmk", "-c", output.name], cwd=report_dir, check=True,
                   stdout=subprocess.DEVNULL)
    return report_dir / "formulas.pdf"


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("summary", type=Path)
    p.add_argument("report_dir", type=Path)
    args = p.parse_args()
    print(make_report(args.summary, args.report_dir))
