#!/usr/bin/env python3
"""Run the discharge calibration, threshold validation and report generation."""
from pathlib import Path
import argparse,subprocess,sys,tomllib
import h5py
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--simulator',required=True,type=Path);p.add_argument('--julia',default='julia');a=p.parse_args()
    cfg=tomllib.loads((ROOT/'config/discharge_calibration.toml').read_text());out=ROOT/cfg['output']
    def run(*args):subprocess.run(list(map(str,args)),cwd=ROOT,check=True)
    def julia(name):run(a.julia,'--project=.',f'scripts/{name}.jl')
    if not (out/'input.h5').exists():julia('build_discharge_inputs')
    acquisition=[sys.executable,'scripts/run_discharge_calibration.py','--simulator',a.simulator.resolve()]
    run(*acquisition);julia('analyze_discharge_calibration')
    with h5py.File(out/'analysis.h5') as h:
        counts=set()
        for kind in ('old','calibrated'):
            n=int(h[kind].attrs['linear_N_snr5']);counts.update((n-1,n,n+1))
    run(*acquisition,'--counts',*sorted(counts));julia('analyze_discharge_calibration')
    run(sys.executable,'scripts/present_discharge_calibration.py')
if __name__=='__main__':main()
