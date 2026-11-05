#!/usr/bin/env python3
"""Terminal parameter form, reproducible cards and complete detector pilot run."""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tomllib


import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from ion_signal.project_paths import resolve_config_paths, resolve_recorded_path
ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "config/specs"
FIGURE6_J_READINGS = {1.5: 0.1, 2.0: 0.05, 3.0: 0.01}


def figure6_j_reading(reduced_field: float) -> float | None:
    """Return a user-supplied approximate Oliveira Fig. 6 reading at an anchor."""
    return next((j for field, j in FIGURE6_J_READINGS.items()
                 if abs(reduced_field - field) < 1e-9), None)


def toml_value(value):
    """Serialize the scalar and string-array values used by pilot cards."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, list) and all(isinstance(item, str) for item in value):
        return "[" + ", ".join(toml_value(item) for item in value) + "]"
    raise TypeError(f"Unsupported pilot TOML value: {type(value).__name__}")


def toml_dump(card: dict) -> str:
    """Write the known one-level pilot card schema with stable key ordering."""
    lines = [f"{key} = {toml_value(card[key])}" for key in sorted(card)
             if not isinstance(card[key], dict)]
    for table in sorted(key for key in card if isinstance(card[key], dict)):
        lines += ["", f"[{table}]"]
        lines += [f"{key} = {toml_value(value)}" for key, value in sorted(card[table].items())]
    return "\n".join(lines) + "\n"


def templates() -> tuple[dict, dict, dict, dict]:
    """Read checked-in realisation, ion, electron and run defaults."""
    paths = (
        TEMPLATE / "realisations/pilot_1mm.toml",
        TEMPLATE / "detector/ions/pitch3mm_bucket30ms.toml",
        TEMPLATE / "detector/electrons/elcc_pilot.toml",
        TEMPLATE / "runs/production_pilot.toml",
    )
    return tuple(resolve_config_paths(tomllib.loads(path.read_text())) for path in paths)


def ask(label: str, default, cast):
    """Prompt once with a visible default and a clear conversion error."""
    answer = input(f"{label} [{default}]: ").strip()
    if not answer:
        return default
    try:
        return cast(answer)
    except ValueError as exc:
        raise ValueError(f"Invalid value for {label}: {answer!r}") from exc


def collect_values(defaults: tuple[dict, dict, dict, dict], interactive: bool) -> dict:
    """Collect the physics controls exposed by this pilot parameter form."""
    r, i, e, run = defaults
    if not interactive:
        return {
            "samples": list(run["selected_samples"]),
            "max_events_per_sample": r["max_events_per_sample"],
            "seed": r["seed"],
            "chamber_length_cm": r["chamber_length_cm"],
            "min_ion_drift_cm": r["min_ion_drift_cm"],
            "w_eV": r["w_eV"],
            "fano_factor": r["fano_factor"],
            "drift_field_V_cm": i["drift_field_V_cm"],
            "pad_pitch_mm": i["pad_pitch_mm"],
            "ion_bucket_ms": i["bucket_width_ms"],
            "electron_lifetime_ms": e["lifetime_ms"],
            "hole_pitch_mm": e["hole_pitch_mm"],
            "electron_bucket_us": e["bucket_width_us"],
            "reduced_field_kV_cm_bar": e["reduced_field_kV_cm_bar"],
            "j_variance_ratio": e["j_variance_ratio"],
            "j_status": e["j_status"],
        }
    print("Detector pilot parameter form. Press Enter to keep each default.")
    print("Approximate Oliveira Fig. 6 J readings: E/P 1.5 → 0.1; 2 → 0.05; 3 → 0.01.")
    choice = ask("Sample (both, bb0nu, single_electrons)", "both", str)
    if choice not in ("both", "bb0nu", "single_electrons"):
        raise ValueError("Sample must be both, bb0nu or single_electrons")
    values = collect_values(defaults, False)
    values["samples"] = ["bb0nu", "single_electrons"] if choice == "both" else [choice]
    prompts = (
        ("max_events_per_sample", "Events per sample", int),
        ("seed", "Random seed", int),
        ("chamber_length_cm", "Chamber length (cm)", float),
        ("min_ion_drift_cm", "Minimum ion drift (cm)", float),
        ("w_eV", "Pair-creation energy W (eV)", float),
        ("fano_factor", "Fano factor F", float),
        ("drift_field_V_cm", "Ion drift field (V/cm)", float),
        ("pad_pitch_mm", "Ion pad pitch (mm)", float),
        ("ion_bucket_ms", "Ion frame duration (ms)", float),
        ("electron_lifetime_ms", "Electron lifetime (ms)", float),
        ("hole_pitch_mm", "EL hole pitch (mm)", float),
        ("electron_bucket_us", "Electron time bucket (us)", float),
        ("reduced_field_kV_cm_bar", "EL E/P (kV/cm/bar, max 3)", float),
    )
    for key, label, cast in prompts:
        values[key] = ask(label, values[key], cast)
    anchor = figure6_j_reading(values["reduced_field_kV_cm_bar"])
    if anchor is None:
        raw = input("EL variance ratio J (required at this E/P): ").strip()
        if not raw:
            raise ValueError("Enter J explicitly when E/P is between Figure 6 anchors")
        try:
            values["j_variance_ratio"] = float(raw)
        except ValueError as exc:
            raise ValueError(f"Invalid EL variance ratio J: {raw!r}") from exc
        values["j_status"] = "user_specified"
    else:
        values["j_variance_ratio"] = ask("EL variance ratio J", anchor, float)
        values["j_status"] = ("oliveira_fig6_user_reading"
                              if abs(values["j_variance_ratio"] - anchor) < 1e-9
                              else "user_specified")
    return values


def validate_values(v: dict) -> None:
    """Reject settings outside the implemented pilot sampler's domain."""
    if not v["samples"] or set(v["samples"]) - {"bb0nu", "single_electrons"}:
        raise ValueError("Unknown sample selection")
    if v["max_events_per_sample"] <= 0 or v["seed"] < 0:
        raise ValueError("Event count must be positive and seed nonnegative")
    if not 0 <= v["min_ion_drift_cm"] < v["chamber_length_cm"]:
        raise ValueError("Ion drift bounds must lie within the chamber")
    if v["w_eV"] <= 0 or not 0 <= v["fano_factor"] < 1:
        raise ValueError("W must be positive and F in [0,1)")
    if any(v[key] <= 0 for key in ("drift_field_V_cm", "pad_pitch_mm", "ion_bucket_ms",
                                    "electron_lifetime_ms", "hole_pitch_mm", "electron_bucket_us")):
        raise ValueError("Field, pitches, bucket widths and lifetime must be positive")
    if not 1.5 <= v["reduced_field_kV_cm_bar"] <= 3:
        raise ValueError("Pilot E/P must be in [1.5,3] kV/cm/bar")
    if not 0 <= v["j_variance_ratio"] < 1:
        raise ValueError("Pilot J must be in [0,1)")
    if v["j_status"] not in ("oliveira_fig6_user_reading", "user_specified"):
        raise ValueError("Unknown J source")
    if v["j_status"] == "oliveira_fig6_user_reading":
        anchor = figure6_j_reading(v["reduced_field_kV_cm_bar"])
        if anchor is None or abs(v["j_variance_ratio"] - anchor) >= 1e-9:
            raise ValueError("Figure 6 J reading requires its matching E/P anchor")


def write_cards(values: dict, output_root: Path = ROOT / "config/specs/generated") -> Path:
    """Save five resolved cards; the content digest gives a stable directory name."""
    validate_values(values)
    realisation, ions, electrons, run = map(deepcopy, templates())
    realisation.update({key: values[key] for key in (
        "max_events_per_sample", "seed", "chamber_length_cm", "min_ion_drift_cm",
        "w_eV", "fano_factor")})
    ions.update({
        "drift_field_V_cm": values["drift_field_V_cm"],
        "pad_pitch_mm": values["pad_pitch_mm"],
        "bucket_width_ms": values["ion_bucket_ms"],
    })
    electrons.update({
        "lifetime_ms": values["electron_lifetime_ms"],
        "hole_pitch_mm": values["hole_pitch_mm"],
        "bucket_width_us": values["electron_bucket_us"],
        "reduced_field_kV_cm_bar": values["reduced_field_kV_cm_bar"],
        "j_variance_ratio": values["j_variance_ratio"],
        "j_status": values["j_status"],
    })
    run["selected_samples"] = list(values["samples"])
    digest = sha256(json.dumps(values, sort_keys=True).encode()).hexdigest()[:12]
    directory = output_root / digest
    directory.mkdir(parents=True, exist_ok=True)
    paths = {name: directory / (name + ".toml") for name in
             ("realisation", "ions", "electrons", "detector", "run")}
    def reference(path: Path) -> str:
        try:
            return str(path.relative_to(ROOT))
        except ValueError:
            return str(path)
    detector = {
        "schema": "nausicaa.detector_configuration.v1",
        "realisation": reference(paths["realisation"]),
        "ions": reference(paths["ions"]),
        "electrons": reference(paths["electrons"]),
    }
    run["detector_configuration"] = reference(paths["detector"])
    for name, card in (("realisation", realisation), ("ions", ions),
                       ("electrons", electrons), ("detector", detector), ("run", run)):
        paths[name].write_text(toml_dump(card))
    return paths["run"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--defaults", action="store_true", help="use all checked-in defaults")
    parser.add_argument("--dry-run", action="store_true", help="write cards and stop before Julia")
    args = parser.parse_args()
    values = collect_values(templates(), interactive=not args.defaults)
    card = write_cards(values)
    print(f"Resolved cards: {card}")
    if args.dry_run:
        return
    env = os.environ.copy()
    env["PYTHON"] = sys.executable
    subprocess.run([env.get("JULIA", "julia"), "--project=.",
                    "scripts/production/run_detector_pilot.jl", str(card)],
                   cwd=ROOT, env=env, check=True)


if __name__ == "__main__":
    main()
