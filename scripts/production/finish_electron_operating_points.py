#!/usr/bin/env python3
"""Wait for complete operating-point manifests, then build the comparison report.

This read-only waiter checks completion every 30 seconds. Production owns the
photon and response files; this step streams them through the Julia reducer and
Python presenter. Its log records failures and completion.
"""
from pathlib import Path
import os,subprocess,time,tomllib
ROOT=Path(__file__).resolve().parents[2]

def finish():
    """Run presentation once both new operating-point productions are complete."""
    while True:
        ready=[]
        for tag in ('ep2p5','ep3'):
            p=ROOT/f'results/electron_operating_points/{tag}.toml'
            if p.exists() and tomllib.loads(p.read_text()).get('complete'):ready.append(tag)
        if len(ready)==2:break
        print('Waiting for completed productions; ready: '+(', '.join(ready) or 'none'),flush=True)
        time.sleep(30)
    env=os.environ.copy()
    env['JULIA_DEPOT_PATH']=f'{ROOT}/work/julia-depot:{ROOT}/notes/physics/work/julia-depot:{Path.home()}/.julia'
    env['MPLCONFIGDIR']=str(ROOT/'work/matplotlib')
    subprocess.run([str(Path.home()/'.juliaup/bin/julia'),'--project=.','scripts/production/reduce_electron_operating_points.jl'],cwd=ROOT,env=env,check=True)
    subprocess.run([str(ROOT/'work/event-viewer-venv/bin/python'),'scripts/production/present_electron_operating_points.py'],cwd=ROOT,env=env,check=True)
    print('Complete: results/electron_operating_points/electron_operating_points.html',flush=True)
if __name__=='__main__':finish()
