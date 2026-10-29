#!/usr/bin/env python3
"""Run the bias sweep, convergence checks, Julia reduction and report."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse,subprocess,sys
ROOT=Path(__file__).resolve().parents[3]
REFINED=['paper_rounded_3','paper_rounded_4','paper_rounded_5','v2_rounded_3','v2_rounded_4','paper_m4_6','v2_m4_6',
         'paper_m3_4','paper_m3_5','paper_m3_7','paper_m3_8','paper_m4_3','paper_m4_9','v2_m4_3','v2_m4_9']
def main():
 p=argparse.ArgumentParser();p.add_argument('--simulator',type=Path,required=True);p.add_argument('--julia',default='julia');a=p.parse_args()
 def run(*args):subprocess.run(list(map(str,args)),cwd=ROOT,check=True)
 if not (ROOT/'results/ltspice/bias_sensitivity/points.h5').exists():run(a.julia,'--project=.','scripts/electronics/bias/build_bias_sensitivity.jl')
 acquire=[sys.executable,'scripts/electronics/bias/run_bias_sensitivity.py','--simulator',a.simulator.resolve()]
 run(*acquire,'--stage','op');run(*acquire,'--stage','transient')
 run(*acquire,'--stage','transient','--refine','--names',*REFINED)
 run(a.julia,'--project=.','scripts/electronics/bias/analyze_bias_sensitivity.jl')
 run(a.julia,'--project=.','tests/test_bias_sensitivity.jl')
 run(sys.executable,'scripts/electronics/bias/present_bias_sensitivity.py')
if __name__=='__main__':main()
