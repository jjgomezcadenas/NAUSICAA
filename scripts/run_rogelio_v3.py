#!/usr/bin/env python3
"""Reproduce the supplied v3 EXP pulse and compare the recorded ion input."""
from pathlib import Path
import argparse,re,hashlib
import h5py,numpy as np
from run_rogelio_v2 import ROOT,preserve,save
from build_ltspice_study import run
from ion_signal.spice import schematic_netlist
def main():
 p=argparse.ArgumentParser();p.add_argument('--simulator',required=True,type=Path);a=p.parse_args()
 source=ROOT/'LTSpice/Spice2/v3/CSAfolded180nm-2_Astronio.asc';model=source.parent/'180nm_professional.inc'
 out=ROOT/'results/ltspice/rogelio_v3';folder=out/'decks';folder.mkdir(parents=True,exist_ok=True)
 preserve(folder/model.name,model.read_text())
 elements,directives=schematic_netlist(source)
 body='* Rogelio v3; connectivity from supplied schematic\n'+'\n'.join(elements)+'\n'
 body+='\n'.join(x for x in directives if not x.lower().startswith('.tran'))+'\n'
 expected={'M1':['s','in','med','Vdda'],'M2':['Vdda','Vnmos','s','0'],'M3':['med','Vnload','0','0'],'M4':['med','Vncasc','out','0'],'M5':['vd6','Vp5','out','Vdda'],'M6':['Vdda','Vp6','vd6','Vdda'],'Mf':['out','Vdis','in','0']}
 for line in elements:
  parts=line.split()
  if parts[0] in expected:assert parts[1:5]==expected[parts[0]],line
 with h5py.File(out/'raw.h5','a') as h:
  h.attrs['source_sha256']=hashlib.sha256(source.read_bytes()).hexdigest();h.attrs['model_sha256']=hashlib.sha256(model.read_bytes()).hexdigest()
  h.attrs['acquisition_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
  def acquire(name,text):
   deck=folder/(name+'.cir');preserve(deck,text)
   if name in h:return
   data,header=run(a.simulator.resolve(),deck);log=deck.with_suffix('.run.txt').read_text()
   assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I),log
   save(h,name,data,header=header,log=log,deck_sha256=hashlib.sha256(deck.read_bytes()).hexdigest())
   print(name,'complete',flush=True)
  acquire('original',body+'.tran 40\n.end\n')
  acquire('original_uncompressed',body+'.options plotwinsize=0\n.tran 40\n.end\n')
  acquire('op',body+'.op\n.end\n')
  for r in (1,2):
   step=5e-6/r
   options='.options reltol=1e-6 abstol=1e-19 vntol=1e-10 numdgt=15 plotwinsize=0\n'
   acquire(f'exp_r{r}',body+options+f'.tran 0 1 0 {step:.17e}\n.save V(out) V(in) I(Is1) V(Qmip_C)\n.end\n')
   with h5py.File(ROOT/'results/ltspice/discharge_calibration/20ms/input.h5') as inp:
    t=inp[f'r{r}/time_s'][:];i=1000*inp[f'r{r}/single_A_per_ion'][:]
    pwl=folder/f'uniform_1000_r{r}.pwl';preserve(pwl,''.join(f'{x:.17e} {y:.17e}\n' for x,y in zip(t,i)))
   uniform=re.sub(r'^Is1 (\S+) (\S+) .*$',rf'Is1 \1 \2 PWL FILE="{pwl.name}"',body,flags=re.M)
   acquire(f'uniform_1000_r{r}',uniform+options+f'.tran 0 {t[-1]:.17e} 0 {step:.17e}\n.save V(out) V(in) I(Is1) V(Qmip_C)\n.end\n')
if __name__=='__main__':main()
