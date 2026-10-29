#!/usr/bin/env python3
"""Launch the six-panel browser viewer for completed diffusion clouds."""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from cloud_viewer.server import serve


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--clouds", type=Path, default=ROOT / "results/three_clouds_8000")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8767)
    parser.add_argument("--open-browser", action="store_true")
    args = parser.parse_args()
    serve(args.clouds, args.host, args.port, args.open_browser)


if __name__ == "__main__":
    main()
