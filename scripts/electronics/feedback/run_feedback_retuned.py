#!/usr/bin/env python3
"""Acquire a feedback variant at a specified M4 bias, retaining all run inputs."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse,hashlib,re,tomllib
import h5py,numpy as np
from ion_signal.spice_acquisition import run
from scripts.electronics.feedback.run_feedback_optimization import ROOT,sha,preserve

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--simulator',required=True,type=Path);p.add_argument('--variant',required=True);p.add_argument('--m4-gate',required=True,type=float);p.add_argument('--op-only',action='store_true');a=p.parse_args()
 cp=ROOT/'config/feedback_optimization.toml';c=tomllib.loads(cp.read_text());v=next(v for v in c['variants'] if v['name']==a.variant);out=ROOT/c['output'];sim=a.simulator.resolve()
 if a.op_only:
  tag=hashlib.sha256(f'{a.variant}/{a.m4_gate:.17e}'.encode()).hexdigest()[:14];name=f'{a.variant}_{tag}';folder=out/'decks/tuning'/name;raw=out/'tuning'/(name+'.h5')
 else:name=a.variant+'_retuned';folder=out/'decks'/name;raw=out/(name+'.h5')
 folder.mkdir(parents=True,exist_ok=True);raw.parent.mkdir(parents=True,exist_ok=True)
 source=ROOT/c['circuit'];model=ROOT/c['model'];preserve(folder/'corrected.inc',model.read_bytes())
 body='\n'.join(x for x in source.read_text().splitlines() if x.strip().lower() != '.op' and not x.lower().startswith(('is1 ','.noise','.tran','.save','.end')))+'\n'
 for pattern,value in [(r'^Vf1 Vdis 0 .*$',f"Vf1 Vdis 0 {v['gate_V']:.17e}"),(r'^\.parameter WMf=.*$',f".parameter WMf={v['width_m']:.17e}"),(r'^\.parameter LMf=.*$',f".parameter LMf={v['length_m']:.17e}"),(r'^Vp1 Vncasc 0 .*$',f'Vp1 Vncasc 0 {a.m4_gate:.17e}')]:body,n=re.subn(pattern,value,body,flags=re.M);assert n==1
 with h5py.File(raw,'a') as h:
  attrs={**v,'name':name,'m4_gate_V':a.m4_gate,'config_sha256':sha(cp),'model_sha256':sha(model),'circuit_sha256':sha(source),'acquisition_sha256':sha(Path(__file__)),'simulator_sha256':sha(sim)}
  for k,val in attrs.items():
   if k in h.attrs:assert h.attrs[k]==val,(name,k)
   else:h.attrs[k]=val
  def acquire(key,tail):
   deck=folder/(key.replace('/','_')+'.cir');preserve(deck,body+tail+'\n.end\n')
   if key in h:assert h[key].attrs['deck_sha256']==sha(deck);return
   d,_=run(sim,deck);log=deck.with_suffix('.run.txt').read_text();assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I)
   g=h.create_group(key);g.attrs.update(log=log,deck_sha256=sha(deck))
   for k,val in d.items():g.create_dataset(k,data=np.asarray(val).real,compression='gzip',shuffle=True)
   h.flush()
  acquire('op','Is1 0 in 0\n.op');h.attrs['baseline_V']=float(h['op/V(out)'][0])
  if a.op_only:print(raw.relative_to(ROOT),flush=True);return
  assert .2<h.attrs['baseline_V']<1.6
  for density in (200,400):acquire('noise_'+str(density),f'Is1 0 in 0 AC 1\n.noise V(out) Is1 dec {density} .1 1Meg')
  with h5py.File(out/'inputs.h5') as inp:
   h.attrs['inputs_sha256']=sha(out/'inputs.h5')
   for frame in c['frames_s']:
    key=f'{round(frame*1000)}ms'
    for r,step in enumerate(c['max_steps_s'],1):
     stop=inp[f'{key}/r{r}/time_s'][-1]
     for n in c['counts']:
      pwl=f'../inputs/{key}_N{n}_r{r}.pwl'
      acquire(f'{key}/single_N{n}_r{r}',f'Is1 0 in PWL FILE="{pwl}"\n.tran 0 {stop:.17e} 0 {step:.17e}\n.save V(in) V(out) I(Is1)')
    print(name,key,'complete',flush=True)
if __name__=='__main__':main()
