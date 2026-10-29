#!/usr/bin/env python3
"""Run the screened nonlinear CSA study and fixed-seed control histories."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse,concurrent.futures,json,os,sys,time,subprocess
import h5py
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src'))
from pad_history_viewer.engine import Engine

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--config',type=Path,default=ROOT/'config/pad_history_viewer_direct.toml')
 p.add_argument('--simulator',type=Path);p.add_argument('--julia');p.add_argument('--workers',type=int,default=4)
 p.add_argument('--limit',type=int);p.add_argument('--refine',type=int,default=24);p.add_argument('--retry-failed',action='store_true')
 args=p.parse_args()
 engine=Engine(args.config,ROOT/'notes/physics/work/pad_compression_batch',args.simulator,args.julia)
 failures=[];started=time.monotonic()
 try:
  plan=engine.plan();tasks=[]
  with h5py.File(plan) as h:
   for name,g in h.items():
    for event,px,py,stratum in zip(g['event_id'][:],g['pixel_x'][:],g['pixel_y'][:],g['stratum'][:]):
     tasks.append((name,int(event),int(px),int(py),int(stratum)))
  # Distribute configurations throughout the batch so progress statistics cover all samples.
  tasks.sort(key=lambda t:(t[4],t[1],t[0],t[2],t[3]))
  alltasks=tasks.copy()
  if args.retry_failed:
   prior=json.loads((engine.directory/'batch_status.json').read_text())
   failed={tuple(item['history']) for item in prior['failed']}
   tasks=[t for t in tasks if t in failed]
  if args.limit is not None:tasks=tasks[:args.limit]
  print('Output:',engine.directory,'planned:',len(tasks),flush=True)
  def calculate(task):
   engine.simulate(*task[:4]);return task
  with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
   futures={pool.submit(calculate,t):t for t in tasks}
   for i,future in enumerate(concurrent.futures.as_completed(futures),1):
    try:future.result()
    except Exception as error:
     task=futures[future];failures.append({'history':task,'error':str(error)});print('FAILED',task,error,flush=True)
    if i%250==0 or i==len(tasks):
     engine.aggregate();print(f'{i}/{len(tasks)}, {time.monotonic()-started:.1f} s, failures={len(failures)}',flush=True)
  # A reproducible spread through completed cases checks the finer current and solver grids.
  controls=[next((t for t in alltasks if t[0]==name and t[4]>0),None) for name in engine.samples]
  controls=[t for t in controls if t is not None]
  chosen=alltasks[::max(1,len(alltasks)//max(1,args.refine))][:max(0,args.refine-len(controls))]+controls
  for task in chosen:
   try:engine.simulate(*task[:4],refine=True)
   except Exception as error:failures.append({'history':task,'refinement_error':str(error)})
  engine.aggregate()
  (engine.directory/'batch_status.json').write_text(json.dumps({'planned':len(alltasks),'scheduled':len(tasks),'failed':failures,'elapsed_s':time.monotonic()-started},indent=2)+'\n')
  print('Complete:',engine.directory,flush=True)
  if failures:raise RuntimeError(f'{len(failures)} simulations failed; see batch_status.json')
  subprocess.run([engine.worker.process.args[0],'--project=.',str(ROOT/'scripts/electronics/histories/report_pad_compression.jl'),str(engine.directory)],cwd=ROOT,env=os.environ.copy()|{'JULIA_DEPOT_PATH':str(ROOT/'notes/physics/work/julia-depot')+':'+str(Path.home()/'.julia')},check=True)
  subprocess.run([sys.executable,str(ROOT/'scripts/electronics/histories/present_pad_compression.py'),'--results',str(engine.directory),'--config',str(engine.card)],cwd=ROOT,check=True)
 finally:engine.close()
if __name__=='__main__':main()
