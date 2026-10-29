#!/usr/bin/env python3
"""Reproduce the feedback-device scan, operating-point tuning and report."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse,os,subprocess,sys,tomllib
import h5py
ROOT=Path(__file__).resolve().parents[3]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--simulator',required=True,type=Path);p.add_argument('--julia',default='julia');a=p.parse_args();sim=a.simulator.resolve()
 env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'src');env.setdefault('JULIA_DEPOT_PATH',str(ROOT/'notes/physics/work/julia-depot')+':'+str(Path.home()/'.julia'))
 def run(*args):subprocess.run(list(map(str,args)),cwd=ROOT,env=env,check=True)
 c=tomllib.loads((ROOT/'config/feedback_optimization.toml').read_text());out=ROOT/c['output']
 for v in c['variants']:
  path=out/(v['name']+'.h5')
  if path.exists():
   with h5py.File(path) as h:
    if 'acquisition_failure' in h.attrs:continue
  result=subprocess.run([sys.executable,'scripts/electronics/feedback/run_feedback_optimization.py','--simulator',str(sim),'--names',v['name']],cwd=ROOT,env=env)
  if result.returncode:
   with h5py.File(path) as h:
    if 'acquisition_failure' not in h.attrs:raise RuntimeError('Unrecorded acquisition failure')
 for name in ('tune_feedback_optimization','tune_feedback_combined'):
  run(a.julia,'--project=.',f'scripts/electronics/feedback/{name}.jl',sys.executable,sim)
 run(sys.executable,'scripts/electronics/feedback/run_feedback_candidate.py','--simulator',sim,'--variant','length_400um_gate0','--m4-gate','0.6164')
 run(a.julia,'--project=.','scripts/electronics/feedback/analyze_feedback_optimization.jl')
 run(sys.executable,'scripts/electronics/feedback/present_feedback_optimization.py')
if __name__=='__main__':main()
