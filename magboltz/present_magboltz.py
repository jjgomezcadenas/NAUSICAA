#!/usr/bin/env python3
"""Make plots and an HTML report from saved Magboltz HDF5 results."""

from pathlib import Path
import csv
import html
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def strings(dataset):
    return np.array([value.decode() if isinstance(value, bytes) else str(value)
                     for value in dataset[:]])


def present(paths):
    paths = [Path(path).resolve() for path in paths]
    if not paths or len({path.parent for path in paths}) != 1:
        raise ValueError("Provide result files from one directory")
    out = paths[0].parent
    plots = out / "plots"
    plots.mkdir(exist_ok=True)
    columns = {name: [] for name in ("gas", "density_mode", "temperature_K",
                                      "pressure_bar", "collision_target",
                                      "field_V_cm", "reduced_field_Td",
                                      "results/velocity_mm_us",
                                      "results/velocity_error_mm_us",
                                      "results/drift_time_150cm_ms",
                                      "results/transverse_card_mm_sqrtbar_sqrtcm",
                                      "results/longitudinal_card_mm_sqrtbar_sqrtcm",
                                      "results/transverse_card_error_mm_sqrtbar_sqrtcm",
                                      "results/longitudinal_card_error_mm_sqrtbar_sqrtcm",
                                      "elapsed_s")}
    source_ids = set()
    for path in paths:
        with h5py.File(path) as f:
            source_ids.add((f["source_sha256"][()].decode(),
                            f["source_filename_version"][()].decode(),
                            f["source_header_version"][()].decode()))
            for name, parts in columns.items():
                parts.append(strings(f[name]) if name in ("gas", "density_mode")
                             else f[name][:])
    if len(source_ids) != 1:
        raise ValueError("Magboltz sources differ across result files")
    data = {name: np.concatenate(parts) for name, parts in columns.items()}
    gas = data["gas"]
    mode = data["density_mode"]
    pressure = data["pressure_bar"]
    temperature = data["temperature_K"]
    collisions = data["collision_target"]
    field = data["field_V_cm"]
    reduced = data["reduced_field_Td"]
    velocity = data["results/velocity_mm_us"]
    velocity_error = data["results/velocity_error_mm_us"]
    drift_time = data["results/drift_time_150cm_ms"]
    transverse = data["results/transverse_card_mm_sqrtbar_sqrtcm"]
    longitudinal = data["results/longitudinal_card_mm_sqrtbar_sqrtcm"]
    transverse_error = data["results/transverse_card_error_mm_sqrtbar_sqrtcm"]
    longitudinal_error = data["results/longitudinal_card_error_mm_sqrtbar_sqrtcm"]
    elapsed = data["elapsed_s"]
    source_sha, source_file_version, source_header_version = source_ids.pop()

    literature = Path(__file__).resolve().parents[1] / "docs/handoffs/data/xe_he_electron_drift_velocity_literature.csv"
    with literature.open(newline="") as stream:
        points = list(csv.DictReader(stream))

    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    for label, marker, color in (("Xe/He 85/15", "o", "#247a9e"),):
        measured = [r for r in points if r["gas"] == label and r["kind"] == "measured"]
        measured = [r for r in measured if 14 <= float(r["reduced_drift_field_V_cm_bar"]) <= 26]
        ax.scatter([float(r["reduced_drift_field_V_cm_bar"]) for r in measured],
                   [float(r["drift_velocity_mm_us"]) for r in measured],
                   marker=marker, s=26, alpha=.45, color=color, label=f"Measured {label}")
    # Literature points use E/P; Magboltz points use the selected physical or
    # model-equivalent pressure, labelled separately.
    for g in np.unique(gas):
        for m in np.unique(mode):
            for p in np.unique(pressure[(gas == g) & (mode == m)]):
                chosen = (gas == g) & (mode == m) & (pressure == p)
                order = np.argsort(field[chosen])
                x = (field[chosen] / pressure[chosen])[order]
                y = velocity[chosen][order]
                err = velocity_error[chosen][order]
                ax.errorbar(x, y, yerr=err, fmt="o-" if m == "geant4_material" else "s--",
                            capsize=2,
                            label=f"{'Density-corrected' if m == 'geant4_material' else 'Plain'} {p:g} bar")
    ax.set(xlabel="Physical E/P (V/cm/bar)", ylabel="Electron drift velocity (mm/µs)",
           title="Electron drift · literature and Magboltz")
    ax.set_xlim(14, 26)
    ax.set_ylim(.2, 1.3)
    ax.grid(alpha=.2)
    ax.legend(fontsize=8, loc="upper left")
    operating = (gas == "Xe_85_He_15") & (mode == "geant4_material") & np.isclose(pressure, 17.5)
    if np.any(operating):
        summary = ", ".join(f"{e:g} V/cm: {v:.3f} mm/µs"
                            for e, v in sorted(zip(field[operating], velocity[operating])))
        ax.text(.98, .04, f"17.5 bar corrected · {summary}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8,
                bbox=dict(facecolor="white", alpha=.85, edgecolor="#ccd6dc"))
    fig.savefig(plots / "velocity_reduced.png", dpi=160)
    plt.close(fig)

    chosen = (gas == "Xe_85_He_15") & np.isclose(pressure, 17.5)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
    for m in np.unique(mode[chosen]):
        mask = chosen & (mode == m)
        order = np.argsort(field[mask])
        x = field[mask][order]
        axes[0].errorbar(x, velocity[mask][order], yerr=velocity_error[mask][order],
                         fmt="o-", capsize=2,
                         label="Density-corrected" if m == "geant4_material" else "Plain")
        axes[1].plot(x, drift_time[mask][order], "o-",
                     label="Density-corrected" if m == "geant4_material" else "Plain")
    axes[0].axhline(.6, color="#777", linestyle="--", linewidth=1,
                    label="Existing 15 bar card: 0.6")
    axes[0].set(ylabel="Velocity (mm/µs)", title="Xe–15% He at 17.5 bar")
    axes[1].set(ylabel="150 cm drift time (ms)", title="From saved Magboltz velocity")
    for ax in axes:
        ax.set_xlabel("Drift field (V/cm)")
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.savefig(plots / "velocity_field.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
    for m in np.unique(mode[chosen]):
        mask = chosen & (mode == m)
        order = np.argsort(field[mask])
        x = field[mask][order]
        for ax, y, err in ((axes[0], transverse, transverse_error),
                           (axes[1], longitudinal, longitudinal_error)):
            ax.errorbar(x, y[mask][order], yerr=err[mask][order], fmt="o-", capsize=2,
                        label="Density-corrected" if m == "geant4_material" else "Plain")
    for ax, title, reference in ((axes[0], "Transverse", 1.6), (axes[1], "Longitudinal", .75)):
        ax.axhline(reference, color="#777", linestyle="--", linewidth=1,
                   label=f"Existing 15 bar card: {reference:g}")
        ax.set(xlabel="Drift field (V/cm)",
               ylabel="Card coefficient (mm·√bar/√cm)", title=title)
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.savefig(plots / "diffusion_field.png", dpi=160)
    plt.close(fig)

    rows = "".join(
        "<tr>" + "".join(f"<td>{html.escape(value)}</td>" for value in (
            gas[i], "Density-corrected" if mode[i] == "geant4_material" else "Plain",
            f"{pressure[i]:.3f}", f"{field[i]:.1f}",
            f"{reduced[i]:.5f}", f"{velocity[i]:.5f} ± {velocity_error[i]:.5f}",
            f"{transverse[i]:.4f} ± {transverse_error[i]:.4f}",
            f"{longitudinal[i]:.4f} ± {longitudinal_error[i]:.4f}",
            f"{elapsed[i]:.1f}")) + "</tr>"
        for i in range(len(gas)))
    result_links = " · ".join(f'<a href="{html.escape(path.name)}">{html.escape(path.stem)}</a>'
                              for path in paths)
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Xe–He electron transport · Magboltz</title>
<style>body{{font:16px/1.5 system-ui;max-width:1200px;margin:30px auto;padding:0 20px;color:#203343}}
img{{max-width:100%}}table{{border-collapse:collapse;width:100%}}th,td{{padding:7px;border-bottom:1px solid #ccd6dc;text-align:right}}
th:first-child,td:first-child{{text-align:left}}a{{color:#185b83}}</style></head><body>
<h1>Xe–He electron transport · Magboltz</h1>
<p>Magboltz source filename: {source_file_version}; internal header: {source_header_version}.
SHA-256: <code>{source_sha}</code>. Gas composition is by mole fraction. The model-density
setting follows the Geant4 material density at the physical pressure.
Temperatures: {', '.join(f'{t:g} K' for t in np.unique(temperature))};
collision targets: {', '.join(f'{int(n):,}' for n in np.unique(collisions))}.
Measured Xe/He points are from McDonald et al.; their physical E/P convention
is described in <a href="../../../docs/handoffs/magboltz_xe_he.md">the handoff</a>.</p>
<h2>Drift velocity</h2><img src="plots/velocity_reduced.png" alt="Velocity and literature points">
<img src="plots/velocity_field.png" alt="Velocity and drift time against field">
<h2>Diffusion</h2><img src="plots/diffusion_field.png" alt="Diffusion coefficients against field">
<h2>Numerical results</h2><table><thead><tr><th>Gas</th><th>Density</th><th>P (bar)</th>
<th>E (V/cm)</th><th>E/N (Td)</th><th>v (mm/µs)</th>
<th>Dₜ (mm·√bar/√cm)</th><th>Dₗ (mm·√bar/√cm)</th><th>Time (s)</th></tr></thead>
<tbody>{rows}</tbody></table>
<p>HDF5 data, inputs and raw output: {result_links} ·
<a href="https://magboltz.web.cern.ch/magboltz/">CERN Magboltz source</a></p>
</body></html>"""
    (out / "index.html").write_text(page)
    print(out / "index.html")


if __name__ == "__main__":
    present(sys.argv[1:])
