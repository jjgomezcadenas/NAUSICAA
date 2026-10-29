#!/usr/bin/env python3
"""Acquire Rogelio's updated noise runs; derived copies fix syntax only.
Physical simulation is performed by LTspice. Julia reduces the saved HDF5.
"""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[4]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse, hashlib, re
import h5py
import numpy as np
from ion_signal.spice_acquisition import run
from ion_signal.spice import schematic_netlist
ROOT=Path(__file__).resolve().parents[4]

def corrected_model(text):
    text='\n'.join(line.rstrip() for line in text.splitlines()).strip()+'\n'
    assert len(re.findall(r'\bNOIBE\b',text))==1
    text=re.sub(r'\bNOIBE\b','NOIB',text)
    # Both closing parentheses were left on commented alternate-model lines.
    assert sum(line.count(')') for line in text.splitlines() if not line.startswith('*'))==0
    text,count=re.subn(r'(?m)^(\+EF[^\n]*)$',r'\1\n+)',text)
    assert count==2
    active='\n'.join(line for line in text.splitlines() if not line.startswith('*'))
    assert active.count('(')==active.count(')')==2
    return text

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--simulator',type=Path,required=True)
    args=ap.parse_args(); simulator=args.simulator.resolve()
    source=ROOT/'LTSpice/Spice2'; out=ROOT/'results/ltspice/spice2_reproduction'
    decks=out/'decks';decks.mkdir(parents=True,exist_ok=True)
    asc=source/'CSAfolded180nm.asc'; model=source/'180nm_professional.inc'
    elements,dirs=schematic_netlist(asc)
    old_elements,_=schematic_netlist(ROOT/'LTSpice/CSAfolded180nm.asc')
    assert elements==old_elements, 'Received circuit changed beyond analysis directives'
    assert 'Cd1 in 0 250f' in elements
    dirs=[d for d in dirs if not d.lower().startswith(('.lib','.noise','.tran','.meas','.op'))]
    revised=corrected_model(model.read_text())
    (decks/'corrected.inc').write_text(revised)
    (decks/'previous.inc').write_text((ROOT/'LTSpice/180nm_professional.inc').read_text())
    with h5py.File(out/'raw.h5','w') as h:
        h.attrs['simulator']=str(simulator)
        h.attrs['corrections']='NMOS NOIBE -> NOIB; close each .MODEL after EF/EM. All numerical coefficients unchanged.'
        h.attrs['conditions']='Original circuit, Cd1=250 fF, Cf=3 fF, original biases and transistor sizes; default 27 C; 1000 frequency points/decade.'
        files=[asc,model,source/'CSAfolded180nm.plt',source/'ltspice noise models.docx',source/'tsmc180 nm noise models spice 3f5.docx',ROOT/'LTSpice/180nm_professional.inc',Path(__file__).resolve(),ROOT/'src/ion_signal/spice.py',simulator]
        h['source_paths']=np.array([str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p) for p in files],dtype=h5py.string_dtype())
        h['source_sha256']=np.array([hashlib.sha256(p.read_bytes()).hexdigest() for p in files],dtype=h5py.string_dtype())
        cases=[('updated_1Hz','corrected','.noise V(out) Is1 dec 1000 1 1G'),
               ('updated_0p1Hz','corrected','.noise V(out) Is1 dec 1000 .1 1G'),
               ('previous_0p1Hz','previous','.noise V(out) Is1 dec 1000 .1 1G'),
               ('updated_op','corrected','.op'),('previous_op','previous','.op')]
        for name,library,analysis in cases:
            deck=decks/(name+'.cir')
            deck.write_text('* Reproduce supplied Spice2 circuit; syntax corrections in derived model only\n'+'\n'.join(elements+dirs)+f'\n.lib {library}.inc\n.options numdgt=15 plotwinsize=0\n{analysis}\n.end\n')
            data,header=run(simulator,deck)
            log=deck.with_suffix('.run.txt').read_text()
            assert not re.search(r'unknown|unrecognized|unrecognised|fatal|syntax error',log,re.I), log
            group=h.create_group(name)
            for col,values in data.items():group[col]=np.asarray(values).real
            group.attrs['log']=log;group.attrs['analysis']=analysis
            print(name, len(next(iter(data.values()))),'points',flush=True)
    print(out/'raw.h5')
if __name__=='__main__':main()
