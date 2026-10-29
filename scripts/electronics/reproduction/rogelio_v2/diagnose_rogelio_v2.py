#!/usr/bin/env python3
"""Record M4 small-signal parameters and amplifier voltage gain at fixed bias."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[4]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse, hashlib, re
import h5py
from ion_signal.spice_acquisition import run
from ion_signal.spice_storage import preserve, save
ROOT = _REPO_ROOT

def main():
 p=argparse.ArgumentParser();p.add_argument('--simulator',type=Path,required=True);a=p.parse_args()
 out=ROOT/'results/ltspice/rogelio_v2';folder=out/'decks/diagnostics';folder.mkdir(parents=True,exist_ok=True)
 preserve(folder/'corrected.inc',(ROOT/'LTSpice/variants/m4_x4_m3_x4/corrected.inc').read_text())
 cases={'paper':ROOT/'LTSpice/variants/m4_x4_m3_x4/noise.cir',
        'v2_native':out/'decks/v2_native_corrected/op.cir',
        'v2_matched':out/'decks/v2_matched/op.cir',
        'v2_retuned':ROOT/'results/ltspice/rogelio_v2_retuned/decks/v2_retuned/op.cir'}
 with h5py.File(out/'diagnostics.h5','a') as h:
  for name,source in cases.items():
   if not source.exists():continue
   body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('.noise','.op','.save','.end','is1 ')))+'\n'
   deck=folder/f'{name}_device.cir'
   preserve(deck,body+'Is1 0 in 0\n.op\n.save V(out) V(in) V(med) V(vd6) @m4[gds] @m4[gm] @m4[gmbs] @m4[id] @mf[id]\n.end\n')
   if name+'_device' not in h:
    d,head=run(a.simulator.resolve(),deck);save(h,name+'_device',d,header=head,deck_sha256=hashlib.sha256(deck.read_bytes()).hexdigest())
   d={k:h[name+'_device'][k][:] for k in h[name+'_device']}
   print(name,{k:float(v[0]) for k,v in d.items()},flush=True)
   # The actual output/input voltage ratio includes the physical feedback loading.
   # Current excitation preserves the operating point without forcing a rounded bias.
   deck=folder/f'{name}_ac.cir'
   preserve(deck,body+'Is1 0 in 0 AC 1\n.ac dec 200 1 1Meg\n.save V(out) V(in)\n.end\n')
   if name+'_ac' not in h:
    data,head=run(a.simulator.resolve(),deck);g=h.create_group(name+'_ac');g.attrs['header']=head
    for k,v in data.items():g.create_dataset(k,data=v)
    h.flush()
if __name__=='__main__':main()
