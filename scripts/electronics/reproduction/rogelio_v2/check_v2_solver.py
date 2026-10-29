#!/usr/bin/env python3
"""Short zero-input transient checks for the high-output-bias v2 circuit."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[4]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import subprocess, argparse, re
from ion_signal.spice_storage import preserve
ROOT = _REPO_ROOT
from ion_signal.spice import read_raw
import h5py
from ion_signal.spice_storage import save
def main():
 p=argparse.ArgumentParser();p.add_argument('--simulator',type=Path,required=True);a=p.parse_args()
 folder=ROOT/'results/ltspice/rogelio_v2/decks/solver';folder.mkdir(parents=True,exist_ok=True)
 source=ROOT/'results/ltspice/rogelio_v2/decks/v2_native_corrected/op.cir'
 body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('.op','.options','.end','is1 ')))+'\n'
 preserve(folder/'corrected.inc',(source.parent/'corrected.inc').read_text())
 with h5py.File(folder.parent.parent/'solver.h5','a') as h:
  for name,opt in [('relaxed','.options reltol=1e-5 abstol=1e-17 vntol=1e-8'),('relaxed_gear','.options method=gear reltol=1e-5 abstol=1e-17 vntol=1e-8'),('default','')]:
   deck=folder/f'{name}.cir';preserve(deck,body+opt+'\n.options numdgt=15 plotwinsize=0\nIs1 0 in 0\n.tran 0 0.0001 0 1e-6\n.save V(in) V(out)\n.end\n')
   if name in h:continue
   try:
    result=subprocess.run([str(a.simulator.resolve()),'-b',str(deck)],cwd=folder,capture_output=True,timeout=30)
    if result.returncode:raise RuntimeError(result.returncode)
    data,header=read_raw(deck.with_suffix('.raw'));save(h,name,data,header=header)
    print(name,'completed',flush=True)
   except Exception as e:
    h.attrs[name+'_failure']=str(e);print(name,e,flush=True)
if __name__=='__main__':main()
