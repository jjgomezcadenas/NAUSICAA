#!/usr/bin/env python3
"""Correct the saved CSA pilot and regenerate diagnostics, using the same noise."""
from pathlib import Path
import argparse,os,subprocess,sys,tomllib
ROOT=Path(__file__).resolve().parents[2]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=ROOT/'config/specs/runs/ion_csa_v3_corrected_pilot.toml');p.add_argument('--julia',default=str(Path.home()/'.juliaup/bin/julia'));a=p.parse_args()
 c=tomllib.loads(a.config.read_text());env=os.environ.copy();env['JULIA_DEPOT_PATH']=':'.join(map(str,(ROOT/'work/julia-depot',ROOT/'notes/physics/work/julia-depot',Path.home()/'.julia')))
 subprocess.run([a.julia,'--startup-file=no','--project=.',str(ROOT/'scripts/production/ion_csa_correction.jl'),str(a.config.resolve())],cwd=ROOT,env=env,check=True)
 directory=max((ROOT/c['results_directory']).glob('*/summary.h5'),key=lambda p:p.stat().st_mtime).parent
 subprocess.run([sys.executable,str(ROOT/'scripts/production/present_ion_csa_pilot.py'),'--results',str(directory)],cwd=ROOT,env=env,check=True)
 print('Corrected pilot:',directory)
if __name__=='__main__':main()
