#!/usr/bin/env python3
"""Acquire v2 results, verify the S/N threshold, and generate the study report."""
from pathlib import Path
import argparse,subprocess,sys,tomllib
import h5py
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--simulator',type=Path,required=True);p.add_argument('--julia',default='julia');a=p.parse_args()
 def run(*args):subprocess.run(list(map(str,args)),cwd=ROOT,check=True)
 acquisition=[sys.executable,'scripts/run_rogelio_v2.py','--simulator',a.simulator.resolve()]
 run(*acquisition)
 run(a.julia,'--project=.','scripts/analyze_rogelio_v2.jl')
 c=tomllib.loads((ROOT/'config/rogelio_v2.toml').read_text())
 with h5py.File(ROOT/c['output']/'analysis.h5') as h:
  cases=[(v['name'],h[v['name']+'/suggested_counts'][:]) for v in c['variants']]
 for name,counts in cases:run(*acquisition,'--variants',name,'--counts',*counts)
 run(a.julia,'--project=.','scripts/analyze_rogelio_v2.jl')
 run(a.julia,'--project=.','scripts/retune_rogelio_v2.jl',sys.executable,a.simulator.resolve())
 tc='config/rogelio_v2_retuned.toml'
 run(*acquisition,'--config',tc)
 run(a.julia,'--project=.','scripts/analyze_rogelio_v2.jl',tc)
 with h5py.File(ROOT/'results/ltspice/rogelio_v2_retuned/analysis.h5') as h:counts=h['v2_retuned/suggested_counts'][:]
 run(*acquisition,'--config',tc,'--counts',*counts)
 run(sys.executable,'scripts/compact_v2_raw.py')
 run(a.julia,'--project=.','scripts/analyze_rogelio_v2.jl')
 run(a.julia,'--project=.','scripts/analyze_rogelio_v2.jl',tc)
 run(sys.executable,'scripts/diagnose_rogelio_v2.py','--simulator',a.simulator.resolve())
 run(a.julia,'--project=.','scripts/analyze_rogelio_v2.jl','--diagnostics')
 run(sys.executable,'scripts/present_rogelio_v2.py')
if __name__=='__main__':main()
