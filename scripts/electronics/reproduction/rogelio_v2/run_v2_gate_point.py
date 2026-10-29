#!/usr/bin/env python3
"""Acquire one unchanged-geometry v2 operating point for Julia bias tuning."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[4]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

import argparse,hashlib
from pathlib import Path
import h5py
from ion_signal.spice_storage import preserve, save
ROOT = _REPO_ROOT
from ion_signal.spice_acquisition import run
def main():
 p=argparse.ArgumentParser();p.add_argument('--gate',type=float,required=True);p.add_argument('--simulator',type=Path,required=True);a=p.parse_args()
 folder=ROOT/'results/ltspice/rogelio_v2_retuned/decks/tuning';folder.mkdir(parents=True,exist_ok=True)
 source=ROOT/'results/ltspice/rogelio_v2/decks/v2_matched/op.cir'
 preserve(folder/'corrected.inc',(source.parent/'corrected.inc').read_text())
 name=f'gate_{a.gate:.14f}'
 deck=folder/f'{name}.cir';preserve(deck,source.read_text().replace('Vp1 Vncasc 0 0.64',f'Vp1 Vncasc 0 {a.gate:.17e}'))
 with h5py.File(folder.parent.parent/'tuning_raw.h5','a') as h:
  if name in h:return
  data,header=run(a.simulator.resolve(),deck)
  save(h,name,data,header=header,gate_V=a.gate,deck_sha256=hashlib.sha256(deck.read_bytes()).hexdigest())
if __name__=='__main__':main()
