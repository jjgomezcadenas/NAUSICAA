"""Locate a pinned LTspice installation and describe its executable identity.

Selection is explicit argument, LTSPICE environment variable, then the local
tooling card. Detector and study cards retain physical settings and scientific
provenance; this module handles machine-specific application paths.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import plistlib
import tomllib
from typing import Mapping

ROOT = Path(__file__).resolve().parents[2]


def resolve_simulator(explicit: str | Path | None = None, *, root: Path = ROOT,
                      environment: Mapping[str, str] | None = None,
                      required: bool = True) -> Path | None:
    """Return the configured executable, checking its presence and execute bit.

    ``explicit`` and ``LTSPICE`` may name either the application bundle or its
    executable. A local card at ``config/local/tools.toml`` uses
    ``[ltspice] executable = "..."`` and may pin ``expected_sha256``. Relative
    paths resolve from the repository root and ``~`` expands to the user's home.
    The pin is checked for the local-card selection; explicit overrides retain
    the caller's existing calibration-hash validation.

    With ``required=False``, an entirely unconfigured simulator returns None so
    the viewer can display saved/linear results. Invalid configured paths and
    mismatched local-card hashes raise an error at selection time.
    """
    environment = os.environ if environment is None else environment
    root = Path(root).resolve()
    selected = explicit or environment.get("LTSPICE")
    pin = None
    if not selected:
        card = root / "config/local/tools.toml"
        if card.exists():
            settings = tomllib.loads(card.read_text()).get("ltspice", {})
            selected = settings.get("executable")
            pin = settings.get("expected_sha256")
    if not selected:
        if not required:
            return None
        raise FileNotFoundError("Configure LTspice with --simulator, LTSPICE, or "
                                "config/local/tools.toml; see docs/ltspice_setup.md")
    path = Path(selected).expanduser()
    if not path.is_absolute():
        path = root / path
    if path.suffix == ".app":
        path = path / "Contents/MacOS/LTspice"
    path = path.resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise FileNotFoundError(f"LTspice executable is missing or not executable: {path}")
    if pin and executable_sha256(path) != pin:
        raise ValueError(f"LTspice executable differs from the pinned local installation: {path}")
    return path


def executable_sha256(path: str | Path) -> str:
    """Return the SHA-256 used by existing acquisition/calibration provenance."""
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def simulator_identity(path: str | Path) -> dict[str, str]:
    """Describe an executable without launching it; read app version when present."""
    path = Path(path).resolve()
    version = "unrecorded"
    info = path.parent.parent / "Info.plist"
    if info.exists():
        metadata = plistlib.loads(info.read_bytes())
        version = str(metadata.get("CFBundleShortVersionString", metadata.get("CFBundleVersion", version)))
    return {"executable": str(path), "version": version,
            "sha256": executable_sha256(path)}


def source_path(path: str | Path, *, root: Path = ROOT) -> str:
    """Record repository-relative inputs and absolute external runtime paths."""
    path = Path(path)
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--simulator", type=Path)
    parser.add_argument("--path", action="store_true", help="Print only the resolved executable path")
    args = parser.parse_args()
    simulator = resolve_simulator(args.simulator)
    print(simulator if args.path else json.dumps(simulator_identity(simulator), indent=2))
