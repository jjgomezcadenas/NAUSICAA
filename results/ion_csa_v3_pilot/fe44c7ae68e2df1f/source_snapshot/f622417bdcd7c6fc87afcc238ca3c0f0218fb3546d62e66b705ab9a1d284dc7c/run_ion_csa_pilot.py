#!/usr/bin/env python3
"""Acquire nonlinear pad controls and run the Julia ion CSA response pilot."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,os,shutil,subprocess,sys,time,tomllib
import h5py,numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from ion_signal.spice_acquisition import run as acquire
LOADED_SHA=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--config',type=Path,default=ROOT/'config/specs/runs/ion_csa_v3_pilot.toml')
 p.add_argument('--simulator',type=Path,default=ROOT/'notes/physics/work/ltspice-review/legacy/LTspice.pkg/Payload/Applications/LTspice.app/Contents/MacOS/LTspice')
 p.add_argument('--workers',type=int,default=4);p.add_argument('--julia',default=str(Path.home()/'.juliaup/bin/julia'))
 p.add_argument('--reuse-cases',type=Path,help='Reuse acquisitions only after current and model content checks')
 p.add_argument('--prepared',type=Path);p.add_argument('--present-only',action='store_true')
 a=p.parse_args();a.config=a.config.resolve();run=tomllib.loads(a.config.read_text());c=tomllib.loads((ROOT/run['response_card']).read_text())
 env=os.environ.copy();env['JULIA_DEPOT_PATH']=':'.join(map(str,(ROOT/'work/julia-depot',ROOT/'notes/physics/work/julia-depot',Path.home()/'.julia')))
 def julia(stage,directory=None):
  cmd=[a.julia,'--startup-file=no','--project=.',str(ROOT/'scripts/production/ion_csa_pilot.jl'),stage,str(a.config)]
  if directory:cmd.append(str(directory))
  print('Julia stage',stage,flush=True);subprocess.run(cmd,cwd=ROOT,env=env,check=True)
 if not a.prepared:julia('prepare')
 directory=a.prepared.resolve() if a.prepared else max((ROOT/run['results_directory']).glob('*/plan.h5'),key=lambda p:p.stat().st_mtime).parent
 started=time.monotonic();summary={'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'runner_sha256':LOADED_SHA,'iterations':[]}
 model=ROOT/c['calibration_directory']/'decks'/c['variant']/'corrected.inc'
 source=ROOT/c['calibration_directory']/'decks'/c['variant']/'30ms_single_N100_r2.cir'
 body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('is1 ','.tran','.save','.end')))+'\n'
 def transient(job,refined=False):
  folder=Path(job['folder']);suffix='refined' if refined else 'history';target=folder/('refined_raw.h5' if refined else 'raw.h5')
  if target.exists():return False
  if a.reuse_cases:
   prior=a.reuse_cases/'spice_cases'/job['key'];cached=prior/target.name
   if cached.exists() and (prior/'current.pwl').exists() and sha(prior/'current.pwl')==sha(folder/'current.pwl') and sha(prior/'corrected.inc')==sha(model):
    for filename in (target.name,suffix+'.cir',suffix+'.run.txt','corrected.inc'):
     old=prior/filename
     if old.exists():shutil.copyfile(old,folder/filename)
    return False
  config=tomllib.loads((folder/'input.toml').read_text());shutil.copyfile(model,folder/'corrected.inc')
  step=c['solver_step_s']/(2 if refined else 1)
  circuit=body
  if refined and (folder/'raw.h5').exists():
   with h5py.File(folder/'raw.h5') as h:
    if h.attrs.get('solver_method')=='gear':circuit=circuit.replace('.options reltol=','.options method=gear reltol=')
  deck=folder/(suffix+'.cir')
  deck.write_text(circuit+f'Is1 0 in PWL FILE="current.pwl"\n.tran 0 {config["stoptime_s"]:.17e} 0 {step:.17e}\n.save V(in) V(out) I(Is1)\n.end\n')
  try:wave,_=acquire(a.simulator,deck,timeout_s=120)
  except subprocess.TimeoutExpired:
   deck.write_text(deck.read_text().replace('.options reltol=','.options method=gear reltol='))
   wave,_=acquire(a.simulator,deck,timeout_s=240)
  with h5py.File(str(target)+'.tmp','w') as h:
   for key,value in wave.items():h.create_dataset(key,data=np.asarray(value).real,compression='gzip',shuffle=True)
   h.attrs.update(deck_sha256=sha(deck),model_sha256=sha(model),simulator_sha256=sha(a.simulator),acquisition_sha256=LOADED_SHA,solver_step_s=step,solver_method='gear' if 'method=gear' in deck.read_text() else 'modified_trap')
  os.replace(str(target)+'.tmp',target)
  for ext in ('.raw','.log','.op.raw'):deck.with_suffix(ext).unlink(missing_ok=True)
  return True
 if not a.present_only:
  for iteration in range(8):
   jobs=tomllib.loads((directory/'jobs.toml').read_text())['jobs'];new=0
   print('Acquiring',len(jobs),'candidate/control histories',flush=True)
   with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
    futures={pool.submit(transient,j):j for j in jobs}
    for i,future in enumerate(concurrent.futures.as_completed(futures),1):
     try:new+=bool(future.result())
     except Exception as error:
      (directory/'acquisition_failure.json').write_text(json.dumps({'job':futures[future],'error':str(error)},indent=2));raise
     if i%50==0 or i==len(jobs):print(f'{i}/{len(jobs)} histories, {time.monotonic()-started:.1f} s',flush=True)
   julia('reduce',directory);julia('assess',directory)
   v=tomllib.loads((directory/'fast_path_validation.toml').read_text());summary['iterations'].append({'new_transients':new,**v})
   (directory/'acquisition_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
   if v['passed']:break
  else:raise RuntimeError('Fast-path validation did not converge')
  # Independent solver refinement of controls spread through the selected set.
  jobs=tomllib.loads((directory/'jobs.toml').read_text())['jobs']
  selected=jobs[::max(1,len(jobs)//12)][:12]
  with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:list(pool.map(lambda j:transient(j,True),selected))
  julia('validate',directory)
  julia('finish',directory);julia('summary',directory)
  summary['elapsed_s']=time.monotonic()-started;summary['direct_histories']=len(jobs);summary['refined_histories']=len(selected)
  (directory/'acquisition_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 subprocess.run([sys.executable,str(ROOT/'scripts/production/present_ion_csa_pilot.py'),'--results',str(directory)],cwd=ROOT,env=env,check=True)
 print('Study:',directory,flush=True)
if __name__=='__main__':main()
