"""Resolve recorded repository paths through the explicit relocation maps.

Scientific cards and manifests retain the paths and hashes recorded by their
original runs. Callers resolve those paths when opening inputs or writing a
presentation; directory roots remain free of compatibility aliases.
"""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]


def resolve_recorded_path(path, *, root=ROOT):
    """Return the current location of a relative or absolute recorded path."""
    root = Path(root).absolute()
    value = Path(path).expanduser()
    if value.is_absolute():
        try:
            relative = value.relative_to(root).as_posix()
        except ValueError:
            return value
    else:
        relative = value.as_posix()
    moves = []
    mapping = root / "docs/ltspice_paths.json"
    if mapping.exists():
        moves.extend((r["old"], r["new"]) for r in json.loads(mapping.read_text())["moves"])
    reports = root / "docs/generated_reports.json"
    if reports.exists():
        moves.extend((f"latex_generated/{name}", f"latex_generated/{spec['group']}/{name}")
                     for name, spec in json.loads(reports.read_text())["reports"].items())
    moves.sort(key=lambda pair: len(pair[0]), reverse=True)
    for _ in range(8):
        replacement = next((new + relative[len(old):] for old, new in moves
                            if relative == old or relative.startswith(old + "/")), relative)
        if replacement == relative:
            return root / relative
        relative = replacement
    raise ValueError(f"Cyclic relocation mapping for {path}")


def resolve_config_paths(value, *, root=ROOT):
    """Map path values in an in-memory card, leaving its saved bytes unchanged."""
    if isinstance(value, dict):
        return {key: resolve_config_paths(item, root=root) for key, item in value.items()}
    if isinstance(value, list):
        return [resolve_config_paths(item, root=root) for item in value]
    if isinstance(value, str) and value.startswith(("LTSpice/", "latex_generated/")):
        return resolve_recorded_path(value, root=root).relative_to(Path(root).absolute()).as_posix()
    return value
