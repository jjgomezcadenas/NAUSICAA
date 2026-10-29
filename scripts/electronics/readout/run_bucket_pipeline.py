#!/usr/bin/env python3
"""Run inputs, LTspice, Julia reduction, integer-threshold check and presentation."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse, subprocess, sys, tomllib
import h5py
ROOT=Path(__file__).resolve().parents[3]
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--simulator',required=True,type=Path)
    p.add_argument('--julia',default='julia')
    a=p.parse_args();cfg=tomllib.loads((ROOT/'config/bucket_study.toml').read_text())
    def run(*args):subprocess.run(list(map(str,args)),cwd=ROOT,check=True)
    def julia(script):run(a.julia,'--project=.',f'scripts/electronics/readout/{script}.jl')
    paths=[ROOT/cfg['output']/f'{round(f*1000)}ms/input.h5' for f in cfg['bucket_sizes_s']]
    if not any(p.exists() for p in paths):julia('build_bucket_inputs')
    if not all(p.exists() for p in paths):raise RuntimeError('Incomplete input set: inspect preserved directories before resuming')
    acquisition=[sys.executable,'scripts/electronics/readout/run_bucket_study.py','--simulator',a.simulator.resolve()]
    run(*acquisition);julia('analyze_bucket_study')
    for frame in cfg['bucket_sizes_s']:
        ms=round(frame*1000)
        with h5py.File(ROOT/cfg['output']/f'{ms}ms/analysis.h5') as h:
            n=int(h.attrs['linear_N_snr5'])
        run(*acquisition,'--bucket-ms',ms,'--counts',n-1,n,n+1)
    julia('analyze_bucket_study')
    run(sys.executable,'scripts/electronics/readout/present_bucket_study.py')
if __name__=='__main__':main()
