#!/usr/bin/env python3
"""Present saved Julia HDF5 parameters and validation/pilot results as LaTeX."""
from pathlib import Path
import argparse
import json
import tomllib
import h5py

ROOT = Path(__file__).resolve().parents[2]

def escape(value):
    return str(value).replace("_", r"\_").replace("%", r"\%")

def render(summary, output, validation):
    output.mkdir(parents=True, exist_ok=True)
    with h5py.File(summary) as h:
        p = h["parameters"]
        rows = [
            ("PDE", "pde", ""), ("Effective microcells", "microcells", ""),
            ("Gain", "gain", r"electrons/PE"), ("Recovery time", "recovery_s", "s"),
            ("Sensor diameter", "diameter_m", "m"), ("Primary dark rate", "dark_rate_Hz", "Hz"),
            ("Charge per fully recovered PE", "charge_per_pe_C", "C"),
            ("Bucket duration", "bucket_width_s", "s"),
            ("ASIC noise RMS", "noise_rms_pe", "PE"), ("ASIC noise RMS", "noise_rms_A", "A"),
            ("Upper average-current limit", "current_limit_A", "A"),
            ("Equivalent upper range", "current_limit_pe", "PE/bucket"),
            ("Mean primary dark pulses", "dark_mean_per_bucket", "per bucket"),
        ]
        table = [r"\begin{tabular}{lrl}\toprule Parameter & Value & Unit \\\midrule"]
        for label, key, unit in rows:
            table.append(f"{label} & {p[key][()]:.7g} & {unit} " + r"\\")
        table.append(r"\bottomrule\end{tabular}")
        block = "\n".join(table)
        (output / "parameters.tex").write_text(r"\renewcommand{\SiPMGeneratedParameters}{" + "\n" + block + "\n}\n")
        pilot = [r"\begin{tabular}{llrrrr}\toprule Sample & Gas & Event & Detected PE & Recovery loss & Clipped buckets \\\midrule"]
        details = []
        sources = {"summary": str(summary), "validation": str(validation), "inputs": []}
        for sample in ("bb0nu", "single_electrons"):
            for gas in ("pure_xe", "xe_he15"):
                g = h[f"{sample}/{gas}"]
                sources["inputs"].append({k: g[k][()].decode() for k in ("source_path", "source_sha256", "sensor_path", "electronics_path")})
                for key in sorted(g.keys()):
                    if not isinstance(g[key], h5py.Group):
                        continue
                    e = g[key]
                    pilot.append(f"{escape(sample)} & {escape(gas)} & {e['event_id'][()]} & "
                                 f"{e['detected_pe'][()]:,} & {e['recovery_loss_fraction'][()]:.5f} & "
                                 f"{e['saturated_buckets'][()]:,} / {e['acquired_channel_buckets'][()]:,} " + r"\\")
                    details.append(f"{escape(sample)}, {escape(gas)}: maximum sensor charge "
                                   f"{e['max_sensor_pe'][()]:.3f} PE; dark pulses {e['dark_pulses'][()]:,}; "
                                   f"ASIC noise RMS {e['noise_rms_pe'][()]:.5f} PE.")
        pilot.append(r"\bottomrule\end{tabular}")
        pilot_block = "\n".join(pilot)
    validation_lines = [r"\begin{tabular}{lr}\toprule Validation quantity & Saved value \\\midrule"]
    with h5py.File(validation) as h:
        if not h["all_tests_passed"][()]:
            raise ValueError("Validation did not pass")
        for key in sorted(h.keys()):
            if key != "all_tests_passed" and isinstance(h[key], h5py.Dataset):
                validation_lines.append(f"{escape(key)} & {h[key][()]:.7g} " + r"\\")
    validation_lines.append(r"\bottomrule\end{tabular}")
    validation_block = "\n".join(validation_lines)
    fragment = r"""\section{Saved pilot and validation results}
The pilot uses event 0 in each of the four existing 250 ns photon histories.
Channels illuminated anywhere in the event are acquired through the global
photon-bucket span with four quiet buckets at either end. Zero-photon buckets
contain dark pulses and ASIC noise. Recovery loss includes optical and dark
charge. Clipping counts include every acquired channel--bucket cell.
These four events establish operation of the pipeline; they do not estimate
the full-sample saturation fraction.

\begin{center}\small
""" + pilot_block + "\n\\end{center}\n\n" + "\n\n".join(details) + r"""

\subsection{Independent numerical checks}
The two-trigger expectation is independent midpoint quadrature of the
uniform-arrival separation density. The stationary test uses independent
Poisson photon counts, 40 microcells and 80 detected triggers per bucket;
the transient test begins with recovered cells. The inverse means demonstrate
finite-bin bias and the difference between stationary and transient illumination.
Dark-count tests use the 5 mm effective sensor; binomial tests use 100 incident
photons at PDE 0.30. ASIC tests include variance, lag-one covariance, charge/current
units and exact range boundaries. Interface checks include signed gaps, empty
events, HDF5 round trips, seed replay and reuse of the sensor output.

\begin{center}\small
""" + validation_block + "\n\\end{center}\n"
    (output / "pilot_results.tex").write_text(fragment)
    report = r"""\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm]{geometry}
\usepackage{booktabs,amsmath,hyperref}
\providecommand{\SiPMGeneratedParameters}{}
\input{parameters.tex}
\title{SiPM and ASIC response: first pilot}
\author{NAUSICAA}\date{}
\begin{document}\maketitle
The Julia response starts from saved incident photons, applies binomial PDE,
uniform microcell illumination and exponential recovery, and adds Poisson
primary dark pulses. The independent ASIC stage adds Gaussian PE noise and
clips the bucket-average current at the configured upper limit. All avalanche
charge is assigned to the trigger bucket. Crosstalk, afterpulsing and gain
spread are deferred. Dark initial ages are stationary; the PDE is constant
during recovery. Float64 charge accumulators preserve the sampled sensor charge.
The physical model and source discussion are in \texttt{docs/sipm\_electronics\_model.tex}.
\section{Saved parameters}
\begin{center}\SiPMGeneratedParameters\end{center}
\input{pilot_results.tex}
\end{document}
"""
    (output / "report.tex").write_text(report)
    (output / "sources.json").write_text(json.dumps(sources, indent=2) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "latex_generated/sipm_response")
    parser.add_argument("--validation", type=Path, default=ROOT / "results/sipm_response_validation/validation.h5")
    args = parser.parse_args()
    summary = args.summary
    if summary is None:
        summary = Path(tomllib.loads((ROOT / "results/sipm_response_pilot/latest.toml").read_text())["summary"])
    render(summary, args.output, args.validation)
    print(args.output / "report.tex")
