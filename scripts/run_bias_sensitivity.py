#!/usr/bin/env python3
"""Acquire the configured bias sweep; reuse the established physical input."""
from pathlib import Path
import argparse,re,tomllib,hashlib
import h5py,numpy as np
from run_rogelio_v2 import ROOT,preserve,save
from build_ltspice_study import run
def main():
 p=argparse.ArgumentParser();p.add_argument('--simulator',type=Path,required=True)
 p.add_argument('--stage',choices=['op','transient'],default='op');p.add_argument('--names',nargs='+');p.add_argument('--refine',action='store_true');a=p.parse_args()
 cp=ROOT/'config/bias_sensitivity.toml';c=tomllib.loads(cp.read_text());out=ROOT/c['output']
 with h5py.File(out/'points.h5') as points,h5py.File(ROOT/c['input']) as inp:
  names=a.names or sorted(points.keys(),key=lambda n:(not ('rounded' in n),n))
  for name in names:
   attr=dict(points[name].attrs);attr={k:v.decode() if isinstance(v,bytes) else v for k,v in attr.items()}
   folder=out/'decks'/name;folder.mkdir(parents=True,exist_ok=True)
   preserve(folder/'corrected.inc',(ROOT/c['model']).read_text())
   original=(ROOT/c[attr['family']]['source']).read_text()
   body='\n'.join(x for x in original.splitlines() if not x.lower().startswith(('.noise','.op','.tran','.save','.end','.options','is1 ')))+'\n'
   body=re.sub(r'^Vnload Vnload 0 .*$',f"Vnload Vnload 0 {attr['m3_gate_V']:.17e}",body,flags=re.M)
   body=re.sub(r'^Vp1 Vncasc 0 .*$',f"Vp1 Vncasc 0 {attr['m4_gate_V']:.17e}",body,flags=re.M)
   raw=out/(name+'.h5')
   with h5py.File(raw,'a') as h:
    h.attrs.update(attr);h.attrs['config_sha256']=hashlib.sha256(cp.read_bytes()).hexdigest()
    h.attrs['input_sha256']=hashlib.sha256((ROOT/c['input']).read_bytes()).hexdigest()
    h.attrs['acquisition_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    h.attrs['single_ion_charge_C']=inp.attrs['single_ion_charge_C']
    h.attrs['assumptions']=inp.attrs['assumptions']
    def acquire(tag,tail,opt):
     deck=folder/(tag+'.cir');preserve(deck,body+opt+'\n'+tail+'\n.end\n')
     if tag in h:return
     data,head=run(a.simulator.resolve(),deck);log=deck.with_suffix('.run.txt').read_text()
     assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I),log
     save(h,tag,data,header=head,log=log,deck_sha256=hashlib.sha256(deck.read_bytes()).hexdigest())
    if a.stage=='op':
     precise='.options reltol=1e-6 abstol=1e-19 vntol=1e-10 numdgt=15 plotwinsize=0'
     acquire('op','Is1 0 in 0\n.op',precise)
     # Preserve complex AC values separately; save() handles real acquisitions.
     deck=folder/'ac.cir';preserve(deck,body+precise+'\nIs1 0 in 0 AC 1\n.ac dec 30 1 1Meg\n.save V(out) V(in)\n.end\n')
     if 'ac' not in h:
      d,head=run(a.simulator.resolve(),deck);g=h.create_group('ac')
      for k,v in d.items():g.create_dataset(k,data=v)
      g.attrs['header']=head
     print(name,'Vout',h['op/V(out)'][0],flush=True)
    else:
     for density in ([200,400] if a.refine else [400]):
      acquire('noise_'+str(density),f'Is1 0 in 0 AC 1\n.noise V(out) Is1 dec {density} .1 1Meg','.options numdgt=15 reltol=1e-6 abstol=1e-19 vntol=1e-10')
     for r in ([1,2] if a.refine else [1]):
      t=inp[f'r{r}/time_s'][:];i=inp[f'r{r}/single_A_per_ion'][:]
      opt='.options reltol=1e-5 abstol=1e-17 vntol=1e-8 numdgt=15 plotwinsize=0' if r==1 else '.options reltol=5e-6 abstol=5e-18 vntol=5e-9 numdgt=15 plotwinsize=0'
      for n in c['counts']:
       tag=f'single_N{n}_r{r}';pwl=folder/(tag+'.pwl')
       preserve(pwl,''.join(f'{x:.17e} {y:.17e}\n' for x,y in zip(t,n*i)))
       acquire(tag,f'Is1 0 in PWL FILE="{pwl.name}"\n.tran 0 {t[-1]:.17e} 0 {5e-6/r:.17e}\n.save V(out) V(in) I(Is1)',opt)
     print(name,'transients complete',flush=True)
if __name__=='__main__':main()
