#!/usr/bin/env python3
"""Launch the MC event viewer: python scripts/viewers/view_events.py [--no-browser]."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=ROOT/'config/event_viewer.toml')
    parser.add_argument('--port',type=int)
    parser.add_argument('--work-dir',type=Path)
    parser.add_argument('--no-browser',action='store_true')
    args=parser.parse_args()
    from event_viewer.server import run
    run(args.config.resolve(),port=args.port,open_browser=not args.no_browser,work_dir=args.work_dir)

if __name__=='__main__':main()
