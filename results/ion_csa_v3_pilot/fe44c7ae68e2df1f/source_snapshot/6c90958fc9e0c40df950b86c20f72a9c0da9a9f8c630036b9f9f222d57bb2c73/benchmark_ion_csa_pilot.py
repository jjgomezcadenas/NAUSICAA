#!/usr/bin/env python3
"""Measure fresh LTspice acquisition time for a spread of saved pilot histories."""
from pathlib import Path
import argparse,concurrent.futures,json,shutil,sys,tempfile,time,tomllib
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from ion_signal.spice_acquisition import run

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);p.add_argument('--cases',type=int,default=100);p.add_argument('--workers',type=int,default=4)
 p.add_argument('--simulator',type=Path,default=ROOT/'notes/physics/work/ltspice-review/legacy/LTspice.pkg/Payload/Applications/LTspice.app/Contents/MacOS/LTspice')
 a=p.parse_args();jobs=tomllib.loads((a.results/'jobs.toml').read_text())['jobs'];n=min(a.cases,len(jobs));chosen=[jobs[round(i*(len(jobs)-1)/max(1,n-1))] for i in range(n)]
 with tempfile.TemporaryDirectory(dir=ROOT/'notes/physics/work',prefix='csa_benchmark_') as folder:
  folder=Path(folder)
  for i,job in enumerate(chosen):
   dest=folder/str(i);dest.mkdir();src=Path(job['folder'])
   for name in ('current.pwl','corrected.inc','history.cir'):shutil.copyfile(src/name,dest/name)
  started=time.monotonic()
  def acquire(i):
   run(a.simulator,folder/str(i)/'history.cir',timeout_s=180)
  with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:list(pool.map(acquire,range(n)))
  elapsed=time.monotonic()-started
 result={'fresh_cases':n,'workers':a.workers,'elapsed_s':elapsed,'wall_seconds_per_history':elapsed/n,'estimated_acquisition_s_for_pilot_selected_histories':elapsed/n*len(jobs),'scope':'Fresh acquisitions only; current generation, reduction, noise sampling and presentation are separate.'}
 (a.results/'acquisition_benchmark.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
