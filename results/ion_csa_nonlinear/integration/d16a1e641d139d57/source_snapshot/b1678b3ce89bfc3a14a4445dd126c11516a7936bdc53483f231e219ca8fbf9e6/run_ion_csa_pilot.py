#!/usr/bin/env python3
"""Acquire nonlinear pad controls and run the Julia ion CSA response pilot."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,os,shutil,subprocess,sys,time,tomllib
import h5py,numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from ion_signal.spice_acquisition import run as acquire
from ion_signal.spice_runtime import resolve_simulator, simulator_identity
LOADED_SHA=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def nonlinear_case(job):
 """Identify an accepted calibrated response; controls retain their SPICE response."""
 path=Path(job['folder'])/'response.h5'
 if not path.exists():return False
 with h5py.File(path) as h:
  value=h.attrs.get('response_model')
  if isinstance(value,bytes):value=value.decode('utf-8')
  return value=='calibrated_nonlinear'

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--config',type=Path,default=ROOT/'config/specs/runs/ion_csa_v3_pilot.toml')
 p.add_argument('--simulator',type=Path,default=None)
 p.add_argument('--workers',type=int,default=4);p.add_argument('--julia',default=str(Path.home()/'.juliaup/bin/julia'))
 p.add_argument('--reuse-cases',type=Path,help='Reuse acquisitions only after current and model content checks')
 p.add_argument('--prepared',type=Path);p.add_argument('--present-only',action='store_true')
 a=p.parse_args();a.simulator=resolve_simulator(a.simulator, root=ROOT);a.config=a.config.resolve();run=tomllib.loads(a.config.read_text());c=tomllib.loads((ROOT/run['response_card']).read_text())
 env=os.environ.copy();env['JULIA_DEPOT_PATH']=':'.join(map(str,(ROOT/'work/julia-depot',Path.home()/'.julia')))
 def julia(stage,directory=None):
  cmd=[a.julia,'--startup-file=no','--project=.',str(ROOT/'scripts/production/ion_csa_pilot.jl'),stage,str(a.config)]
  if directory:cmd.append(str(directory))
  print('Julia stage',stage,flush=True);subprocess.run(cmd,cwd=ROOT,env=env,check=True)
 calibration_dir=ROOT/c['calibration_directory']
 calibration_file=calibration_dir/'analysis.h5'
 if calibration_file.exists():
  with h5py.File(calibration_file) as saved:
   recorded=saved[c['variant']].attrs['simulator_sha256']
   if isinstance(recorded,bytes):recorded=recorded.decode()
   if recorded!=sha(a.simulator):
    raise RuntimeError('CSA calibration simulator changed; regenerate calibration')
 if not calibration_file.exists():
  calibration_dir.mkdir(parents=True,exist_ok=True)
  design=tomllib.loads((ROOT/c['design_card']).read_text())
  with h5py.File(ROOT/design['study_directory']/'operating_points.h5') as src,h5py.File(calibration_dir/'operating_points.h5','a') as dst:
   if c['variant'] not in dst:src.copy(src[c['variant']],dst,name=c['variant'])
  subprocess.run([sys.executable,str(ROOT/'scripts/electronics/feedback/run_retention_scan.py'),'--config',str(ROOT/c['calibration_card']),'--simulator',str(a.simulator),'--julia',a.julia,'--stage','calibrate'],cwd=ROOT,env=env,check=True)
 if not a.prepared:julia('prepare')
 directory=a.prepared.resolve() if a.prepared else max((ROOT/run['results_directory']).glob('*/plan.h5'),key=lambda p:p.stat().st_mtime).parent
 started=time.monotonic();summary={'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'runner_sha256':LOADED_SHA,'simulator':simulator_identity(a.simulator),'iterations':[],'charge_resolution_retries':[]}
 model=ROOT/c['calibration_directory']/'decks'/c['variant']/'corrected.inc'
 source=ROOT/c['calibration_directory']/'decks'/c['variant']/'30ms_single_N100_r2.cir'
 model_sha=sha(model);simulator_sha=sha(a.simulator)
 body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('is1 ','.tran','.save','.end')))+'\n'
 def transient(job,refined=False,*,step_override=None,force=False):
  if nonlinear_case(job):return False
  folder=Path(job['folder']);suffix='refined' if refined else 'history';target=folder/('refined_raw.h5' if refined else 'raw.h5')
  if target.exists() and not force:return False
  if target.exists() and force:
   archive=folder/(target.stem+'_rejected_'+str(time.time_ns())+'.h5');target.rename(archive)
  if a.reuse_cases and not force:
   prior=a.reuse_cases/'spice_cases'/job['key'];cached=prior/target.name
   same_inputs=cached.exists() and (prior/'current.pwl').exists() and (prior/'corrected.inc').exists() and sha(prior/'current.pwl')==sha(folder/'current.pwl') and sha(prior/'corrected.inc')==model_sha
   if same_inputs:
    with h5py.File(cached) as saved:
     same_inputs=(saved.attrs['model_sha256']==model_sha and saved.attrs['simulator_sha256']==simulator_sha and saved.attrs['solver_step_s']==c['solver_step_s']/(2 if refined else 1))
   if same_inputs:
    for filename in (target.name,suffix+'.cir',suffix+'.run.txt','corrected.inc'):
     old=prior/filename
     if old.exists():shutil.copyfile(old,folder/filename)
    return False
  config=tomllib.loads((folder/'input.toml').read_text());shutil.copyfile(model,folder/'corrected.inc')
  step=step_override if step_override is not None else c['solver_step_s']/(2 if refined else 1)
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
   h.attrs.update(deck_sha256=sha(deck),model_sha256=model_sha,simulator_sha256=simulator_sha,acquisition_sha256=LOADED_SHA,solver_step_s=step,solver_method='gear' if 'method=gear' in deck.read_text() else 'modified_trap')
  os.replace(str(target)+'.tmp',target)
  for ext in ('.raw','.log','.op.raw'):deck.with_suffix(ext).unlink(missing_ok=True)
  return True
 def reduce_with_resolution(stage):
  attempts={};marker=directory/'charge_retry.toml'
  while True:
   marker.unlink(missing_ok=True)
   try:julia(stage,directory);return
   except subprocess.CalledProcessError:
    if not marker.exists():raise
    failed=tomllib.loads(marker.read_text());folder=Path(failed['folder']);refined=failed['refined'];key=(str(folder),refined)
    attempts[key]=attempts.get(key,0)+1
    if attempts[key]>3:raise RuntimeError('Current quadrature failed after three finer acquisitions: '+str(folder))
    raw=folder/('refined_raw.h5' if refined else 'raw.h5')
    with h5py.File(raw) as h:oldstep=float(h.attrs['solver_step_s'])
    newstep=min(oldstep/4,c['current_step_s']/2**attempts[key])
    print('Refining current quadrature:',folder.name,'step',newstep,flush=True)
    transient({'folder':str(folder),'row':failed['row'],'key':folder.name},refined,step_override=newstep,force=True)
    summary['charge_resolution_retries'].append({**failed,'old_step_s':oldstep,'new_step_s':newstep})
    (directory/'acquisition_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 if not a.present_only:
  for iteration in range(8):
   julia('nonlinear',directory)
   jobs=tomllib.loads((directory/'jobs.toml').read_text())['jobs'];new=0
   nonlinear_summary=directory/'nonlinear_summary.toml'
   if nonlinear_summary.exists():
    expected=tomllib.loads(nonlinear_summary.read_text())['nonlinear_histories']
    if sum(nonlinear_case(j) for j in jobs)!=expected:
     raise RuntimeError('Julia and Python disagree on accepted nonlinear histories')
   print('Acquiring',len(jobs),'candidate/control histories',flush=True)
   with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
    futures={pool.submit(transient,j):j for j in jobs}
    for i,future in enumerate(concurrent.futures.as_completed(futures),1):
     try:new+=bool(future.result())
     except Exception as error:
      (directory/'acquisition_failure.json').write_text(json.dumps({'job':futures[future],'error':str(error)},indent=2));raise
     if i%50==0 or i==len(jobs):print(f'{i}/{len(jobs)} histories, {time.monotonic()-started:.1f} s',flush=True)
   reduce_with_resolution('reduce');julia('assess',directory)
   v=tomllib.loads((directory/'fast_path_validation.toml').read_text());summary['iterations'].append({'new_transients':new,**v})
   (directory/'acquisition_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
   if v['passed']:break
  else:raise RuntimeError('Fast-path validation did not converge')
  # Independent solver refinement of controls spread through the selected set.
  jobs=tomllib.loads((directory/'jobs.toml').read_text())['jobs']
  spice_jobs=[j for j in jobs if not nonlinear_case(j)]
  selected=spice_jobs[::max(1,len(spice_jobs)//12)][:12]
  with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:list(pool.map(lambda j:transient(j,True),selected))
  reduce_with_resolution('validate')
  julia('finish',directory);julia('summary',directory)
  summary['elapsed_s']=time.monotonic()-started;summary['direct_histories']=len(spice_jobs);summary['nonlinear_histories']=len(jobs)-len(spice_jobs);summary['refined_histories']=len(selected)
  (directory/'acquisition_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 subprocess.run([sys.executable,str(ROOT/'scripts/production/present_ion_csa_pilot.py'),'--results',str(directory)],cwd=ROOT,env=env,check=True)
 print('Study:',directory,flush=True)
if __name__=='__main__':main()
