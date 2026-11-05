#!/usr/bin/env python3
"""Format saved ELM calculations as LaTeX macros, without repeating physics."""
from pathlib import Path
import h5py

ROOT = Path(__file__).resolve().parents[2]
source = ROOT / "results/elm/elm_example.h5"
target = ROOT / "latex_generated/detector/elm/elm_values.tex"
target.parent.mkdir(parents=True, exist_ok=True)
with h5py.File(source, "r") as data:
    lines = ["% Generated from results/elm/elm_example.h5; do not edit values."]
    for key in sorted(data["values"]):
        value = float(data["values"][key][()])
        lines.append(r"\newcommand{\elm" + key + "}{" + f"{value:.8f}" + "}")
    for label, prefix in (("four", "Four"), ("five", "Five")):
        pilot = data[f"pilots/{label}_tau_12p5"]
        for key, suffix in (("event_count", "Events"), ("retained_events", "Retained"),
                            ("asic_charge_lost_percent", "ASICLoss"),
                            ("corrected_readout_rms_percent", "FullRMS"),
                            ("retained_readout_rms_percent", "RetainedRMS"),
                            ("corrected_unclipped_rms_percent", "UnclippedRMS"),
                            ("event_mean_average_pe", "MeanPE"),
                            ("event_max_average_pe", "MaxPE")):
            lines.append(r"\newcommand{\elmPilot" + prefix + suffix + "}{"
                         + f"{float(pilot[key][()]):.8f}" + "}")
    rows = [r"\begin{tabular}{rrrrrrr}", r"\toprule",
            r"Pitch & $d$ & $\tau_r$ & SiPM loss & ASIC loss & Retained & Full / retained RMS\\",
            r"(mm) & (mm) & (ns) & (\%) & (\%) & events & (\%)\\", r"\midrule"]
    for label, pitch, diameter in (("five", 5, 3), ("four", 4, 2.5)):
        for tag, recovery in (("tau_12p5", 12.5), ("tau_25", 25)):
            pilot = data[f"pilots/{label}_{tag}"]
            v = lambda k: float(pilot[k][()])
            rows.append(f"{pitch} & {diameter:g} & {recovery:g} & "
                        f"{v('sipm_charge_lost_percent'):.4f} & {v('asic_charge_lost_percent'):.4f} & "
                        f"{int(v('retained_events'))}/{int(v('event_count'))} & "
                        f"{v('corrected_readout_rms_percent'):.4f} / {v('retained_readout_rms_percent'):.4f}"
                        + r"\\")
    rows += [r"\bottomrule", r"\end{tabular}"]
    (target.parent / "elm_pilot_table.tex").write_text("\n".join(rows) + "\n")
target.write_text("\n".join(lines) + "\n")
print(target)
