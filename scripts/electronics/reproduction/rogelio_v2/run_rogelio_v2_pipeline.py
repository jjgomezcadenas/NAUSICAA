#!/usr/bin/env python3
"""Acquire v2 results, verify the S/N threshold, and generate the study report."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[4]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse,subprocess,sys,tomllib
import h5py
ROOT=Path(__file__).resolve().parents[4]
def main():
 p=argparse.ArgumentParser();p.add_argument('--simulator',type=Path,required=True);p.add_argument('--julia',default='julia');a=p.parse_args()
 def run(*args):subprocess.run(list(map(str,args)),cwd=ROOT,check=True)
 acquisition=[sys.executable,'scripts/electronics/reproduction/rogelio_v2/run_rogelio_v2.py','--simulator',a.simulator.resolve()]
 run(*acquisition)
 run(a.julia,'--project=.','scripts/electronics/reproduction/rogelio_v2/analyze_rogelio_v2.jl')
 c=tomllib.loads((ROOT/'config/rogelio_v2.toml').read_text())
 with h5py.File(ROOT/c['output']/'analysis.h5') as h:
  cases=[(v['name'],h[v['name']+'/suggested_counts'][:]) for v in c['variants']]
 for name,counts in cases:run(*acquisition,'--variants',name,'--counts',*counts)
 run(a.julia,'--project=.','scripts/electronics/reproduction/rogelio_v2/analyze_rogelio_v2.jl')
 run(a.julia,'--project=.','scripts/electronics/reproduction/rogelio_v2/retune_rogelio_v2.jl',sys.executable,a.simulator.resolve())
 tc='config/rogelio_v2_retuned.toml'
 run(*acquisition,'--config',tc)
 run(a.julia,'--project=.','scripts/electronics/reproduction/rogelio_v2/analyze_rogelio_v2.jl',tc)
 with h5py.File(ROOT/'results/ltspice/rogelio_v2_retuned/analysis.h5') as h:counts=h['v2_retuned/suggested_counts'][:]
 run(*acquisition,'--config',tc,'--counts',*counts)
 run(sys.executable,'scripts/electronics/reproduction/rogelio_v2/compact_v2_raw.py')
 run(a.julia,'--project=.','scripts/electronics/reproduction/rogelio_v2/analyze_rogelio_v2.jl')
 run(a.julia,'--project=.','scripts/electronics/reproduction/rogelio_v2/analyze_rogelio_v2.jl',tc)
 run(sys.executable,'scripts/electronics/reproduction/rogelio_v2/diagnose_rogelio_v2.py','--simulator',a.simulator.resolve())
 run(a.julia,'--project=.','scripts/electronics/reproduction/rogelio_v2/analyze_rogelio_v2.jl','--diagnostics')
 run(sys.executable,'scripts/electronics/reproduction/rogelio_v2/present_rogelio_v2.py')
if __name__=='__main__':main()
