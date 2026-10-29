#!/usr/bin/env python3
"""Acquire nonlinear LTspice waveforms for the specified ten-bucket example."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse, hashlib, os, subprocess, sys, tomllib
import h5py
import numpy as np
from ion_signal.spice_acquisition import run
ROOT=Path(__file__).resolve().parents[3]
def sha(p): return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
def preserve(p,content):
    if isinstance(content,str):content=content.encode()
    if p.exists() and p.read_bytes()!=content:raise ValueError(f'Preserving {p}')
    p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(content)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--simulator',type=Path,required=True)
    parser.add_argument('--julia',default='julia')
    args=parser.parse_args();sim=args.simulator.resolve()
    cp=ROOT/'config/pad_sequence_example.toml';c=tomllib.loads(cp.read_text());out=ROOT/c['output']
    env=os.environ.copy();env.setdefault('JULIA_DEPOT_PATH',str(ROOT/'notes/physics/work/julia-depot')+':'+str(Path.home()/'.julia'))
    if not (out/'inputs.h5').exists():
        subprocess.run([args.julia,'--project=.','scripts/electronics/histories/pad_sequence_example.jl','--build'],cwd=ROOT,env=env,check=True)
    preserve(out/'configuration.toml',cp.read_bytes())
    preserve(out/'.gitignore','decks/**/*.raw\ndecks/**/*.log\n')
    source=ROOT/c['source_decks']/'40ms_single_N400_r2.cir';model=source.parent/'corrected.inc'
    body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('is1 ','.tran','.save','.end')))+'\n'
    with h5py.File(out/'inputs.h5') as inp:
        assert inp.attrs['config_sha256'].decode()==sha(cp)
        assert inp.attrs['generator_sha256'].decode()==sha(ROOT/'src/ion_signal/pad_sequence_example.jl')
        for name,case in inp.items():
            folder=out/'decks'/name;preserve(folder/'corrected.inc',model.read_bytes())
            with h5py.File(out/(name+'.h5'),'a') as h:
                for key,val in {'config_sha256':sha(cp),'inputs_sha256':sha(out/'inputs.h5'),'circuit_sha256':sha(source),'model_sha256':sha(model),'simulator_sha256':sha(sim),'acquisition_sha256':sha(Path(__file__))}.items():
                    if key in h.attrs:assert h.attrs[key]==val
                    else:h.attrs[key]=val
                for r,step in enumerate(c['solver_max_steps_s'],1):
                    t=case[f'r{r}/time_s'][:];current=case[f'r{r}/current_A'][:]
                    for baseline in (True,False):
                        key=('baseline' if baseline else 'history')+f'_r{r}'
                        pwl=folder/(key+'.pwl')
                        preserve(pwl,''.join(f'{x:.17e} {y:.17e}\n' for x,y in zip(t,np.zeros_like(current) if baseline else current)))
                        deck=folder/(key+'.cir')
                        preserve(deck,body+f'Is1 0 in PWL FILE="{pwl.name}"\n.tran 0 {t[-1]:.17e} 0 {step:.17e}\n.save V(in) V(out) I(Is1)\n.end\n')
                        if key in h:
                            assert h[key].attrs['deck_sha256']==sha(deck)
                            continue
                        wave,_=run(sim,deck)
                        group=h.create_group(key);group.attrs.update(deck_sha256=sha(deck),pwl_sha256=sha(pwl),log=deck.with_suffix('.run.txt').read_text())
                        for var,val in wave.items():group.create_dataset(var,data=np.asarray(val).real,compression='gzip',shuffle=True)
                        h.flush();print(name,key,'complete',flush=True)
    subprocess.run([args.julia,'--project=.','scripts/electronics/histories/pad_sequence_example.jl'],cwd=ROOT,env=env,check=True)
if __name__=='__main__':main()
