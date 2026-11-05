#!/usr/bin/env python3
"""Compile a report with all LaTeX intermediates in work/latex.

The source working directory is retained for relative figures and inputs. Only
the completed PDF is copied into the report directory. Authored notes can set
the repository as their working directory and a separate PDF destination.
"""
from pathlib import Path
import argparse
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def build_report(source, *, cwd=None, destination=None):
    """Compile one TeX source; return its delivered PDF and keep the build log."""
    source = Path(source).resolve()
    cwd = Path(cwd).resolve() if cwd else source.parent
    destination = Path(destination).resolve() if destination else source.with_suffix(".pdf")
    try:
        identity = source.relative_to(ROOT).with_suffix("")
    except ValueError:
        raise ValueError("Report sources must be inside the repository") from None
    build = ROOT / "work/latex" / identity
    build.mkdir(parents=True, exist_ok=True)
    with (build / "latexmk.stdout").open("w") as stream:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
                        f"-outdir={build}", str(source)], cwd=cwd, check=True,
                       stdout=stream, stderr=subprocess.STDOUT)
    log = (build / (source.stem + ".log")).read_text(errors="replace")
    errors = ("Undefined control sequence", "There were undefined references",
              "There were undefined citations", "LaTeX Error:")
    if any(token in log for token in errors):
        raise RuntimeError(f"Unresolved LaTeX diagnostics: {build}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(build / source.with_suffix(".pdf").name, destination)
    warnings = [line for line in log.splitlines()
                if "Warning:" in line or "Overfull" in line or "Underfull" in line]
    if warnings:
        print(f"{len(warnings)} LaTeX layout/package warnings; log: {build / (source.stem + '.log')}")
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--cwd", type=Path)
    parser.add_argument("--destination", type=Path)
    args = parser.parse_args()
    print(build_report(args.source, cwd=args.cwd, destination=args.destination))
