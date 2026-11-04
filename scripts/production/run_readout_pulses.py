#!/usr/bin/env python3
"""Reduce saved readout histories and regenerate the HTML/LaTeX/PDF comparison."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-index", type=Path,
                        default=ROOT / "results/simulations/production_pilot/"
                                       "grid_scan_summary_238ddc5a75e0fd9b.h5")
    parser.add_argument("--electron-index", type=Path,
                        default=ROOT / "results/simulations/production_pilot/"
                                       "electron_250ns_summary_464425eac32a2c34.h5")
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "results/readout_pulses_8000")
    parser.add_argument("--latex-dir", type=Path,
                        default=ROOT / "latex_generated/readout_pulses_8000")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reduction = args.output_dir / "readout_pulses.h5"
    subprocess.run([os.environ.get("JULIA", "julia"), "--project=.",
                    str(ROOT / "scripts/production/reduce_readout_pulses.jl"),
                    str(args.scan_index), str(args.electron_index), str(reduction)],
                   cwd=ROOT, check=True)
    subprocess.run([sys.executable,
                    str(ROOT / "scripts/production/present_readout_pulses.py"),
                    str(reduction), "--output-dir", str(args.output_dir),
                    "--latex-dir", str(args.latex_dir)], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
