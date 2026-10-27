#!/usr/bin/env python3
"""Acquire one unchanged-geometry v2 operating point for Julia bias tuning."""
import argparse,hashlib
from pathlib import Path
import h5py
from run_rogelio_v2 import ROOT,preserve,save
from build_ltspice_study import run
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
