#!/usr/bin/env python3
"""Open the local three-panel ion pad-history viewer."""
from pathlib import Path
import argparse,sys,importlib.util,os
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
def main():
 fallback=ROOT/'work/event-viewer-venv/bin/python'
 if importlib.util.find_spec("plotly") is None and fallback.exists() and Path(sys.executable).absolute()!=fallback.absolute():
  os.execv(str(fallback),[str(fallback),str(Path(__file__)),*sys.argv[1:]])
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--config',type=Path,default=ROOT/'config/pad_history_viewer.toml')
 p.add_argument('--port',type=int);p.add_argument('--no-browser',action='store_true')
 p.add_argument('--work-dir',type=Path);p.add_argument('--simulator',type=Path);p.add_argument('--julia')
 a=p.parse_args()
 from pad_history_viewer.server import run
 run(a.config,a.port,not a.no_browser,a.work_dir,a.simulator,a.julia)
if __name__=='__main__':main()
