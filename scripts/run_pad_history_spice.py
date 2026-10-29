#!/usr/bin/env python3
"""Validate complete demanding pad histories with the optimized transistor circuit."""
from pathlib import Path
import argparse,hashlib,re,tomllib
import h5py,numpy as np
from build_ltspice_study import run
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def preserve(p,data):
 if isinstance(data,str):data=data.encode()
 if p.exists() and p.read_bytes()!=data:raise ValueError(f'Preserving {p}')
 p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--simulator',required=True,type=Path);p.add_argument('--cases',nargs='+');p.add_argument('--config',type=Path,default=ROOT/'config/pad_history_readout_direct.toml');a=p.parse_args()
 cp=a.config.resolve();c=tomllib.loads(cp.read_text());out=ROOT/c['output'];sim=a.simulator.resolve()
 preserve(out/'configuration.toml',cp.read_bytes());preserve(out/'.gitignore','decks/**/*.raw\ndecks/**/*.log\n')
 with h5py.File(out/'spice_inputs.h5') as inp:
  for name in inp:
   if a.cases and name not in a.cases:continue
   case=inp[name];folder=out/'decks'/name;folder.mkdir(parents=True,exist_ok=True)
   model=ROOT/c['source_decks']/'corrected.inc';preserve(folder/'corrected.inc',model.read_bytes())
   source=ROOT/c['source_decks']/'40ms_single_N400_r2.cir'
   body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('is1 ','.tran','.save','.end')))+'\n'
   result=out/'spice'/(name+'.h5');result.parent.mkdir(parents=True,exist_ok=True)
   with h5py.File(result,'a') as h:
    attrs={**dict(case.attrs),'config_sha256':sha(cp),'inputs_sha256':sha(out/'spice_inputs.h5'),'circuit_sha256':sha(source),'model_sha256':sha(model),'simulator_sha256':sha(sim),'acquisition_sha256':sha(Path(__file__))}
    for k,v in attrs.items():
     if k in h.attrs:
      if k!='acquisition_sha256':assert h.attrs[k]==v,(name,k)
     else:h.attrs[k]=v
    for r,maxstep in enumerate(c['tail_max_steps_s'],1):
     t=case[f'r{r}/time_s'][:];current=case[f'r{r}/current_A'][:]
     for baseline in (True,False):
      key=('baseline' if baseline else 'history')+f'_r{r}'
      pwl=folder/(key+'.pwl');preserve(pwl,''.join(f'{x:.17e} {y:.17e}\n' for x,y in zip(t,np.zeros_like(current) if baseline else current)))
      deck=folder/(key+'.cir');preserve(deck,body+f'Is1 0 in PWL FILE="{pwl.name}"\n.tran 0 {t[-1]:.17e} 0 {maxstep:.17e}\n.save V(in) V(out) I(Is1)\n.end\n')
      if key in h:assert h[key].attrs['deck_sha256']==sha(deck);continue
      try:d,_=run(sim,deck)
      except Exception as e:
       h.attrs['acquisition_failure']=str(e);print(name,key,'FAILED',e,flush=True);break
      log=deck.with_suffix('.run.txt').read_text();assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I)
      g=h.create_group(key);g.attrs.update(deck_sha256=sha(deck),pwl_sha256=sha(pwl),acquisition_sha256=sha(Path(__file__)),log=log)
      for k,v in d.items():g.create_dataset(k,data=np.asarray(v).real,compression='gzip',shuffle=True)
      h.flush();print(name,key,'complete',flush=True)
     if 'acquisition_failure' in h.attrs:break
if __name__=='__main__':main()
