"""Present the calibrated Julia/HDF5 averaging scan."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from scripts.paper.build_spice_paper import plot_averaging

if __name__ == "__main__":
    plot_averaging()
