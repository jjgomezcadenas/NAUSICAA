#!/usr/bin/env python3
"""Reproduce both preserved retention scans and their combined report source."""
from pathlib import Path
import argparse,os,subprocess,sys,tomllib
ROOT=Path(__file__).resolve().parents[3]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--simulator',required=True,type=Path)
    p.add_argument('--julia',default=str(Path.home()/'.juliaup/bin/julia'))
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--config',type=Path,default=ROOT/'config/retention_comparison.toml')
    a=p.parse_args();cp=a.config.resolve();c=tomllib.loads(cp.read_text())
    env=os.environ.copy()
    env.setdefault('JULIA_DEPOT_PATH',str(ROOT/'notes/physics/work/julia-depot')+':'+str(Path.home()/'.julia'))
    env['PYTHONPATH']=str(ROOT/'src')
    def run(*parts):subprocess.run(list(map(str,parts)),cwd=ROOT,env=env,check=True)
    for card in c['study_cards']:
        run(sys.executable,ROOT/'scripts/electronics/feedback/run_retention_scan.py',
            '--config',ROOT/card,'--simulator',a.simulator.resolve(),'--julia',a.julia,'--workers',a.workers)
    run(a.julia,'--startup-file=no','--project=.',ROOT/'scripts/electronics/feedback/combine_retention_scans.jl',cp)
    run(sys.executable,ROOT/'scripts/electronics/feedback/present_retention_scan.py','--config',cp)

if __name__=='__main__':main()
