#!/usr/bin/env python3
"""Run independent calibrated retention history studies with private Julia workers."""
from pathlib import Path
import argparse,concurrent.futures,json,os,subprocess,sys,time,tomllib
ROOT=Path(__file__).resolve().parents[3]

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--simulator',type=Path,required=True)
 p.add_argument('--config',type=Path,default=ROOT/'config/retention_comparison.toml')
 p.add_argument('--parallel',type=int,default=2)
 p.add_argument('--workers',type=int,default=2)
 p.add_argument('--only',nargs='+')
 a=p.parse_args();c=tomllib.loads(a.config.read_text());jobs=[]
 directory=ROOT/c['output'];directory.mkdir(parents=True,exist_ok=True)
 logs=directory/'batch_logs';logs.mkdir(exist_ok=True)
 for card in c['study_cards']:
  study=tomllib.loads((ROOT/card).read_text())
  for variant in study['variants']:
   name=variant['name']
   if a.only and name not in a.only:continue
   qualified=ROOT/study['output']/'cards'/(name+'_validated_worker_viewer.toml')
   viewer=qualified if qualified.exists() else ROOT/study['output']/'cards'/(name+'_viewer.toml')
   jobs.append((name,viewer))
 started=time.monotonic()
 def execute(job):
  name,card=job;log=logs/(name+'.stdout')
  env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'src')
  env.setdefault('JULIA_DEPOT_PATH',str(ROOT/'notes/physics/work/julia-depot')+':'+str(Path.home()/'.julia'))
  env.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/retention_matplotlib'))
  command=[sys.executable,str(ROOT/'scripts/electronics/histories/run_pad_compression.py'),
           '--config',str(card),'--simulator',str(a.simulator.resolve()),'--workers',str(a.workers),
           '--work-dir',str(ROOT/'notes/physics/work'/('retention_parallel_'+name))]
  with log.open('w') as stream:
   result=subprocess.run(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT)
  print(name,'complete' if result.returncode==0 else 'FAILED',f'{time.monotonic()-started:.1f} s',flush=True)
  return dict(name=name,returncode=result.returncode,log=str(log.relative_to(ROOT)))
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.parallel) as pool:
  results=list(pool.map(execute,jobs))
 (directory/('batch_runs_'+str(time.time_ns())+'.json')).write_text(json.dumps(dict(results=results,elapsed_s=time.monotonic()-started),indent=2)+'\n')
 if any(r['returncode'] for r in results):raise RuntimeError('Retention batch failed; inspect the saved logs')

if __name__=='__main__':main()
