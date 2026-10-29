#!/usr/bin/env python3
"""Acquire the configured feedback-device variants using the saved physical input."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse,hashlib,re,tomllib
import h5py,numpy as np
from ion_signal.spice_acquisition import run
ROOT=Path(__file__).resolve().parents[3]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def preserve(p,data):
 if isinstance(data,str):data=data.encode()
 if p.exists() and p.read_bytes()!=data:raise ValueError(f'Preserving {p}')
 p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--simulator',required=True,type=Path);p.add_argument('--names',nargs='+');a=p.parse_args()
 cfgfile=ROOT/'config/feedback_optimization.toml';c=tomllib.loads(cfgfile.read_text());out=ROOT/c['output'];out.mkdir(parents=True,exist_ok=True)
 source=ROOT/c['circuit'];model=ROOT/c['model'];sim=a.simulator.resolve()
 preserve(out/'configuration.toml',cfgfile.read_bytes())
 preserve(out/'.gitignore','decks/**/*.raw\ndecks/**/*.log\n')
 with h5py.File(out/'inputs.h5','a') as inp:
  for frame,inputfile in zip(c['frames_s'],c['input_sources']):
   key=f'{round(frame*1000)}ms';path=ROOT/inputfile
   if key not in inp:
    with h5py.File(path) as original:
     g=inp.create_group(key);g.attrs['source']=inputfile;g.attrs['source_sha256']=sha(path)
     for k,v in original.attrs.items():g.attrs[k]=v
     for r in (1,2):original.copy(original[f'r{r}'],g,name=f'r{r}')
   else:assert inp[key].attrs['source_sha256']==sha(path)
   for r in (1,2):
    t=inp[f'{key}/r{r}/time_s'][:];unit=inp[f'{key}/r{r}/single_A_per_ion'][:]
    for n in c['counts']:
     preserve(out/'decks/inputs'/f'{key}_N{n}_r{r}.pwl',''.join(f'{x:.17e} {y:.17e}\n' for x,y in zip(t,n*unit)))
 with h5py.File(out/'inputs.h5') as inp:
  for variant in c['variants']:
   name=variant['name']
   if a.names and name not in a.names:continue
   folder=out/'decks'/name;folder.mkdir(parents=True,exist_ok=True);preserve(folder/'corrected.inc',model.read_bytes())
   body='\n'.join(x for x in source.read_text().splitlines() if x.strip().lower() != '.op' and not x.lower().startswith(('is1 ','.noise','.tran','.save','.end')))+'\n'
   for pattern,value in [(r'^Vf1 Vdis 0 .*$',f"Vf1 Vdis 0 {variant['gate_V']:.17e}"),(r'^\.parameter WMf=.*$',f".parameter WMf={variant['width_m']:.17e}"),(r'^\.parameter LMf=.*$',f".parameter LMf={variant['length_m']:.17e}")]:body,n=re.subn(pattern,value,body,flags=re.M);assert n==1
   with h5py.File(out/(name+'.h5'),'a') as h:
    attrs={**variant,'config_sha256':sha(cfgfile),'circuit_sha256':sha(source),'model_sha256':sha(model),'inputs_sha256':sha(out/'inputs.h5'),'simulator_sha256':sha(sim),'acquisition_sha256':sha(Path(__file__))}
    for k,v in attrs.items():
     if k in h.attrs:assert h.attrs[k]==v,(name,k)
     else:h.attrs[k]=v
    def acquire(tag,tail):
     deck=folder/(tag.replace('/','_')+'.cir');preserve(deck,body+tail+'\n.end\n')
     if tag in h:assert h[tag].attrs['deck_sha256']==sha(deck);return
     try:d,head=run(sim,deck)
     except Exception as e:
      h.attrs['acquisition_failure']=str(e);raise
     log=deck.with_suffix('.run.txt').read_text();assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I)
     g=h.create_group(tag);g.attrs.update(deck_sha256=sha(deck),log=log)
     for k,v in d.items():g.create_dataset(k,data=np.asarray(v).real,compression='gzip',shuffle=True)
     h.flush()
    acquire('op','Is1 0 in 0\n.op')
    baseline=float(h['op/V(out)'][0]);h.attrs['baseline_V']=baseline
    if not c['minimum_output_V']<baseline<c['maximum_output_V']:
     h.attrs['rejected_operating_point']=True;print(name,'rejected DC',baseline,flush=True);continue
    for density in (200,400):acquire('noise_'+str(density),f'Is1 0 in 0 AC 1\n.noise V(out) Is1 dec {density} .1 1Meg')
    for frame in c['frames_s']:
     key=f'{round(frame*1000)}ms'
     for r,step in enumerate(c['max_steps_s'],1):
      stop=inp[f'{key}/r{r}/time_s'][-1]
      for n in c['counts']:
       pwl=f'../inputs/{key}_N{n}_r{r}.pwl'
       acquire(f'{key}/single_N{n}_r{r}',f'Is1 0 in PWL FILE="{pwl}"\n.tran 0 {stop:.17e} 0 {step:.17e}\n.save V(in) V(out) I(Is1)')
     print(name,key,'complete',flush=True)
if __name__=='__main__':main()
