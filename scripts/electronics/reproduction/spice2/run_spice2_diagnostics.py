#!/usr/bin/env python3
"""Acquire capacitance and short-injection controls with the revised noise card."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[4]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse, hashlib, re
import h5py
import numpy as np
from ion_signal.spice_acquisition import run
ROOT=Path(__file__).resolve().parents[4]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--simulator',required=True,type=Path)
    a=p.parse_args(); sim=a.simulator.resolve()
    out=ROOT/'results/ltspice/spice2_diagnostics'; decks=out/'decks'
    decks.mkdir(parents=True,exist_ok=True)
    base=ROOT/'results/ltspice/spice2_enc/decks'
    source=base/'noise_400.cir'
    body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('is1 ','.noise','.end')))+'\n'
    model=base/'corrected.inc'; (decks/'corrected.inc').write_bytes(model.read_bytes())
    pwl=ROOT/'results/ltspice/decks/uniform_r2.pwl'
    data=np.loadtxt(pwl)
    raw=out/'raw.h5'
    with h5py.File(raw,'w') as h:
        h.attrs.update(frame_s=.01,averaging_s=300e-6,correction_tau_s=.3,temperature_K=293.15,
          frequency_low_Hz=.1,frequency_high_Hz=1e6,simulator=str(sim),
          assumptions='External capacitance scan; intrinsic MOS capacitances retained. Original device dimensions and biases. Stationary NOIMOD=2 noise. Uniform finite-gap ion current, no renormalization. Short injection 1 us plateau, 10 ns edges, one elementary charge. Trapezoid rise T, flat G=T, no pole-zero correction.')
        files=[source,model,pwl,Path(__file__).resolve(),sim]
        h['source_paths']=np.array([str(x.relative_to(ROOT)) for x in files],dtype=h5py.string_dtype())
        h['source_sha256']=np.array([hashlib.sha256(x.read_bytes()).hexdigest() for x in files],dtype=h5py.string_dtype())
    def acquire(cap,name,tail):
        tag=f'cap_{cap:g}pF'; path=decks/f'{tag}_{name}.cir'
        circuit=re.sub(r'^Cd1 .*$',f'Cd1 in 0 {cap:g}p',body,flags=re.M)
        path.write_text(circuit+tail+'\n.end\n')
        d,_=run(sim,path); log=path.with_suffix('.run.txt').read_text()
        assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I)
        with h5py.File(raw,'a') as h:
            group=h.require_group(tag);group.attrs['external_capacitance_F']=cap*1e-12
            g=group.create_group(name);g.attrs['log']=log
            for key,value in d.items():g[key]=np.asarray(value).real
        print(tag,name,flush=True)
    for cap in (.75,1.25,3.):
        for density in (200,400):
            acquire(cap,f'noise_{density}',f'Is1 0 in 0 AC 1\n.noise V(out) Is1 dec {density} .1 1Meg')
        for r in (1,2):
            for n in (0,1):
                filename=f'uniform_N{n}.pwl'
                np.savetxt(decks/filename,np.c_[data[:,0],n*data[:,1]],fmt='%.17e')
                acquire(cap,f'ions_N{n}_r{r}',f'Is1 0 in PWL FILE="{filename}"\n.tran 0 .032 0 {5e-6/r:.17e}\n.save V(in) V(out) I(Is1)')
        if cap==1.25:
            for r in (1,2):
                for n in (0,1):
                    amp=n*1.602176634e-19/1.01e-6
                    acquire(cap,f'injection_N{n}_r{r}',f'Is1 0 in PULSE(0 {amp:.17e} .005 10n 10n 1u 1)\n.tran 0 .014 0 {1e-6/r:.17e}\n.save V(in) V(out) I(Is1)')
if __name__=='__main__':main()
