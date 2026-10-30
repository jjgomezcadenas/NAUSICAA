#!/usr/bin/env python3
"""Run the full 2/3/4 mm ion scan with one electron pass per species."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parents[2]


def julia_stage(script: str, args: list[Path], log: Path) -> str:
    """Run one Julia stage, retain its log and display sparse progress lines."""
    log.parent.mkdir(parents=True, exist_ok=True)
    command = [os.environ.get("JULIA", "julia"), "--project=.",
               str(ROOT / "scripts/production" / script), *(str(path) for path in args)]
    last = ""
    with log.open("w") as output:
        with subprocess.Popen(command, cwd=ROOT, env=os.environ.copy(), text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as process:
            assert process.stdout is not None
            for line in process.stdout:
                output.write(line)
                last = line.strip()
                match = re.search(r"(?:event |: )(\d+)(?:/8000)?(?: events)?", line)
                if ("Summary:" in line or "events in" in line or
                        "reused completed" in line or
                        match and int(match.group(1)) % 1000 == 0):
                    print(last, flush=True)
            if process.wait() != 0:
                raise RuntimeError(f"{script} failed; see {log}")
    return log.read_text()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p2", type=Path,
                        default=ROOT / "config/specs/generated/4725adea70b9/run.toml")
    parser.add_argument("--p3", type=Path,
                        default=ROOT / "config/specs/generated/df1414afcdc2/run.toml")
    parser.add_argument("--p4", type=Path,
                        default=ROOT / "config/specs/generated/91f0e47f63fc/run.toml")
    parser.add_argument("--full-summary", type=Path,
                        help="reuse a completed 3 mm joint run")
    args = parser.parse_args()
    for card in (args.p2, args.p3, args.p4):
        if not card.is_file():
            raise FileNotFoundError(card)
    if args.full_summary is None:
        print("Running one joint ion/electron pass per species at 3 mm", flush=True)
        log = julia_stage("run_detector_pilot.jl", [args.p3],
                          ROOT / "work/grid_scan_joint.log")
        found = re.findall(r"^Summary: (.+\.h5)$", log, flags=re.MULTILINE)
        if len(found) != 1:
            raise RuntimeError("Joint run did not announce exactly one HDF5 summary")
        full_summary = Path(found[0])
    else:
        full_summary = args.full_summary.resolve()
        if not full_summary.is_file():
            raise FileNotFoundError(full_summary)
        print(f"Reusing completed joint summary: {full_summary}", flush=True)
    for pitch, card in ((2, args.p2), (4, args.p4)):
        print(f"Collecting ions at {pitch} mm", flush=True)
        julia_stage("run_ion_readout.jl", [card],
                    ROOT / f"work/grid_scan_ions_{pitch}mm.log")
    index_log = julia_stage("summarize_grid_scan.jl",
                            [args.p2, args.p3, args.p4, full_summary],
                            ROOT / "work/grid_scan_summary.log")
    index = Path(index_log.strip().splitlines()[-1])
    if not index.is_file():
        raise RuntimeError("Numerical scan index was not written")
    run_id = index.stem.removeprefix("grid_scan_summary_")
    report_dir = ROOT / "latex_generated/production_grid_scan" / run_id
    subprocess.run([os.environ.get("PYTHON", "python3.13"),
                    str(ROOT / "scripts/production/present_grid_scan.py"),
                    str(index), str(report_dir)], cwd=ROOT, check=True)
    print(f"Numerical index: {index}")
    print(f"Report: {report_dir / 'grid_scan.pdf'}")


if __name__ == "__main__":
    main()
