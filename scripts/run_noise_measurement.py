"""Snapshot Rogelio's updated v3 and reproduce its integrated noise measurement."""
from pathlib import Path
import argparse, hashlib
import h5py
from ion_signal.spice import schematic_netlist
from run_rogelio_v2 import ROOT, preserve, save
from build_ltspice_study import run

def main():
    p=argparse.ArgumentParser();p.add_argument('--simulator',type=Path,required=True);a=p.parse_args()
    source=ROOT/'LTSpice/Spice2/v3/CSAfolded180nm-2_Astronio.asc'
    model=source.parent/'180nm_professional.inc'
    folder=ROOT/'results/ltspice/noise_measurement';decks=folder/'decks';decks.mkdir(parents=True,exist_ok=True)
    preserve(decks/source.name,source.read_text(encoding='latin1'))
    preserve(decks/model.name,model.read_text())
    elements,directives=schematic_netlist(source)
    body='* v3 integrated-noise comparison\n'+'\n'.join(elements)+'\n'
    body+='\n'.join(x for x in directives if not x.lower().startswith(('.tran','.noise','.meas')))+'\n'
    with h5py.File(folder/'raw.h5','a') as h:
        for n in (200,400):
            name=f'noise_{n}';deck=decks/(name+'.cir')
            preserve(deck,body+f'.options numdgt=15\n.noise V(out) Is1 dec {n} .1 1Meg\n.meas NOISE int_sigma INTEG V(onoise) FROM .1 TO 1Meg\n.end\n')
            if name not in h:
                data,header=run(a.simulator.resolve(),deck)
                save(h,name,data,header=header,log=deck.with_suffix('.run.txt').read_text(),deck_sha256=hashlib.sha256(deck.read_bytes()).hexdigest())
                print(name, 'complete', flush=True)
        h.attrs['source_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
        h.attrs['model_sha256']=hashlib.sha256(model.read_bytes()).hexdigest()

if __name__=='__main__':main()
