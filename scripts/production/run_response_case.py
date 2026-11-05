#!/usr/bin/env python3
"""Run or reuse the first bb0nu detector response and generate its report."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]


def manifest_from_run(command: list[str], prefix: str) -> Path:
    result = subprocess.run(command, cwd=ROOT, check=True, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout, end="")
    lines = [line.removeprefix(prefix).strip() for line in result.stdout.splitlines()
             if line.startswith(prefix)]
    if len(lines) != 1:
        raise RuntimeError(f"Expected one {prefix!r} line from production runner")
    return Path(lines[0])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ion-manifest", type=Path)
    p.add_argument("--sipm-manifest", type=Path)
    p.add_argument("--present-only", action="store_true")
    p.add_argument("--result-dir", type=Path, default=ROOT / "results/response_bb0nu_3mm")
    p.add_argument("--latex-dir", type=Path, default=ROOT / "latex_generated/production/response_bb0nu_3mm")
    a = p.parse_args()
    julia = os.environ.get("JULIA", "julia")
    a.result_dir.mkdir(parents=True, exist_ok=True)
    analysis = a.result_dir / "analysis.h5"
    if not a.present_only:
        ion = a.ion_manifest or manifest_from_run(
            [julia, "--startup-file=no", "--project=.",
             str(ROOT / "scripts/production/run_ion_csa.jl"),
             "config/specs/runs/ion_csa_bb0nu_3mm_production.toml"],
            "Ion response manifest: ")
        sipm = a.sipm_manifest or manifest_from_run(
            [julia, "--startup-file=no", "--project=.",
             str(ROOT / "scripts/production/run_sipm_production.jl"),
             "config/specs/runs/sipm_bb0nu_production.toml"],
            "SiPM production manifest: ")
        subprocess.run([julia, "--startup-file=no", "--project=.",
                        str(ROOT / "scripts/production/reduce_response_case.jl"),
                        str(ion), str(sipm), str(analysis)], cwd=ROOT, check=True)
    if not analysis.is_file():
        raise FileNotFoundError(analysis)
    subprocess.run([sys.executable, str(ROOT / "scripts/production/present_response_case.py"),
                    str(analysis), "--result-dir", str(a.result_dir),
                    "--latex-dir", str(a.latex_dir)], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
