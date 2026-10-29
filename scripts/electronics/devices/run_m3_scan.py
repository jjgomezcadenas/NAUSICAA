#!/usr/bin/env python3
"""Run derived M3 geometries at fixed W/L and unchanged bias voltages."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse, hashlib, re
import h5py
import numpy as np
from ion_signal.spice_acquisition import run
ROOT=Path(__file__).resolve().parents[3]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--simulator',required=True,type=Path)
    p.add_argument('--scales',nargs='+',type=int,default=[1])
    p.add_argument('--gate',type=float,default=None)
    p.add_argument('--op-only',action='store_true')
    p.add_argument('--tight-injection',action='store_true')
    p.add_argument('--injection-only',action='store_true')
    p.add_argument('--injection-counts',nargs='+',type=int,default=[0,1,10,100])
    p.add_argument('--refinements',nargs='+',type=int,default=[1,2])
    p.add_argument('--counts',nargs='+',type=int,default=[0,1,1000])
    a=p.parse_args();sim=a.simulator.resolve()
    out=ROOT/'results/ltspice/m3_scan';decks=out/'decks';decks.mkdir(parents=True,exist_ok=True)
    source=ROOT/'results/ltspice/m4_scan/decks/scale_4_gate_0.5568034172_noise_400.cir'
    model=source.parent/'corrected.inc'
    body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('is1 ','.noise','.end')))+'\n'
    body=body.replace('* M4 geometry scan; all other device dimensions and bias voltages unchanged','* M3 geometry scan; M4 x4 and its tuned bias retained')
    (decks/'corrected.inc').write_bytes(model.read_bytes())
    pwl=ROOT/'results/ltspice/decks/uniform_r2.pwl';data=np.loadtxt(pwl)
    raw=out/'raw.h5'
    with h5py.File(raw,'a') as h:
        h.attrs.update(external_capacitance_F=1.25e-12,feedback_F=3e-15,frame_s=.01,averaging_s=300e-6,correction_tau_s=.3,temperature_K=293.15,
          frequency_low_Hz=.1,frequency_high_Hz=1e6,simulator=str(sim),
          assumptions='Scale M3 W and L together; retain W/L=100 and multiplier=1; M4 dimensions x4 and gate 0.556803417206 V. Optional M3 gate adjustment is stored per variant. All other dimensions, biases and noise coefficients retained. External 1.25 pF plus intrinsic MOS capacitances. Uniform finite-gap current over 10 ms; no renormalization. Stationary NOIMOD=2 noise. Short injection: 1 us plateau and 10 ns edges; one elementary charge.')
        files=[source,model,pwl,Path(__file__).resolve(),ROOT/'scripts/electronics/reproduction/initial/build_ltspice_study.py',ROOT/'src/ion_signal/spice.py',sim]
        for name,values in [('source_paths',[str(x.relative_to(ROOT)) for x in files]),('source_sha256',[hashlib.sha256(x.read_bytes()).hexdigest() for x in files])]:
            if name in h:del h[name]
            h[name]=np.array(values,dtype=h5py.string_dtype())
    def acquire(scale,name,tail):
        tag=f'scale_{scale}' + (f'_gate_{a.gate:.10f}' if a.gate is not None else '')
        with h5py.File(raw) as h:
            if f'{tag}/{name}' in h:return
        circuit=body.replace('.param Lnload=0.5u',f'.param Lnload={.5*scale:g}u').replace('.param Wnload=50u',f'.param Wnload={50*scale:g}u')
        if a.gate is not None:circuit=circuit.replace('Vnload Vnload 0 0.6',f'Vnload Vnload 0 {a.gate:.12g}')
        if name.startswith('injection_precise'):
            circuit=circuit.replace('reltol=1e-6 abstol=1e-19 vntol=1e-10','reltol=1e-7 abstol=1e-21 vntol=1e-11')
        path=decks/f'{tag}_{name}.cir';path.write_text(circuit+tail+'\n.end\n')
        d,_=run(sim,path);log=path.with_suffix('.run.txt').read_text()
        assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I)
        with h5py.File(raw,'a') as h:
            group=h.require_group(tag);group.attrs.update(linear_scale=scale,W_m=50e-6*scale,L_m=.5e-6*scale,area_scale=scale**2,gate_V=.6 if a.gate is None else a.gate)
            g=group.create_group(name);g.attrs['log']=log;g.attrs['deck_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
            for key,value in d.items():g.create_dataset(key,data=np.asarray(value).real,compression="gzip",shuffle=True)
        print(tag,name,flush=True)
    for scale in a.scales:
        assert scale>=1
        acquire(scale,'op','Is1 0 in 0\n.op')
        if a.op_only:continue
        for density in (200,400):
            acquire(scale,f'noise_{density}',f'Is1 0 in 0 AC 1\n.noise V(out) Is1 dec {density} .1 1Meg')
        for r in a.refinements:
            for n in ([] if a.injection_only else a.counts):
                filename=f'uniform_N{n}.pwl'
                np.savetxt(decks/filename,np.c_[data[:,0],n*data[:,1]],fmt='%.17e')
                acquire(scale,f'ions_N{n}_r{r}',f'Is1 0 in PWL FILE="{filename}"\n.tran 0 .032 0 {5e-6/r:.17e}\n.save V(in) V(out) I(Is1)')
            for n in a.injection_counts:
                amp=n*1.602176634e-19/1.01e-6
                kind='injection_precise' if a.tight_injection else 'injection'
                acquire(scale,f'{kind}_N{n}_r{r}',f'Is1 0 in PULSE(0 {amp:.17e} .005 10n 10n 1u 1)\n.tran 0 .014 0 {1e-6/r:.17e}\n.save V(in) V(out) I(Is1)')
if __name__=='__main__':main()
