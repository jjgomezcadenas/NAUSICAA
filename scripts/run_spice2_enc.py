#!/usr/bin/env python3
"""Run 1.25 pF ion-frame SPICE study; Julia controls the threshold search.
Reuse the validated finite-gap PWL, preserve all transistor dimensions/biases.
"""
from pathlib import Path
import argparse, hashlib, re, subprocess
import h5py
import numpy as np
from build_ltspice_study import run
from reproduce_spice2 import corrected_model
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--simulator',required=True,type=Path);p.add_argument('--julia',default='julia')
    args=p.parse_args();sim=args.simulator.resolve()
    out=ROOT/'results/ltspice/spice2_enc';decks=out/'decks';decks.mkdir(parents=True,exist_ok=True)
    source=ROOT/'results/ltspice/capacitance/pad_1.25pF/decks/ions_N1_r2.cir'
    pwl=ROOT/'results/ltspice/decks/uniform_r2.pwl'
    model=ROOT/'LTSpice/Spice2/180nm_professional.inc'
    body='\n'.join(s for s in source.read_text().splitlines() if not s.lower().startswith(('.tran','.save','.end','is1 ')))+'\n'
    body=body.replace('.lib 180nm_professional.inc','.lib corrected.inc')
    assert '.temp 20' in body and 'Cd1 in 0 1.24999999999999992e-12' in body
    (decks/'corrected.inc').write_text(corrected_model(model.read_text()))
    rawpwl=np.loadtxt(pwl)
    raw=out/'raw.h5'
    with h5py.File(raw,'w') as h:
        h.attrs.update(external_pad_F=1.25e-12,feedback_F=3e-15,frame_s=.01,averaging_s=300e-6,
            correction_tau_s=.3,temperature_K=293.15,frequency_low_Hz=.1,frequency_high_Hz=1e6,
            simulator=str(sim),assumptions='Uniform mean ion arrivals; 0.5 mm square pad in grounded plane; finite 5 mm approach; no charge renormalization; intrinsic MOS capacitances added to external pad load; stationary NOIMOD=2 noise; continuous centered averages; no ADC noise.')
        h['input_time_s']=rawpwl[:,0];h['input_current_per_ion_A']=rawpwl[:,1]
        files=[source,pwl,model,Path(__file__).resolve(),ROOT/'scripts/reproduce_spice2.py',ROOT/'src/ion_signal/spice.py',sim]
        h['source_paths']=np.array([str(x.relative_to(ROOT)) for x in files],dtype=h5py.string_dtype())
        h['source_sha256']=np.array([hashlib.sha256(x.read_bytes()).hexdigest() for x in files],dtype=h5py.string_dtype())
    def acquire(name,tail):
        deck=decks/(name+'.cir');deck.write_text(body+tail+'\n.end\n')
        data,header=run(sim,deck);log=deck.with_suffix('.run.txt').read_text()
        assert not re.search(r'unknown|unrecognized|unrecognised|fatal|syntax error',log,re.I)
        with h5py.File(raw,'a') as h:
            group=h.create_group(name)
            for col,values in data.items():group[col]=np.asarray(values).real
            group.attrs['log']=log
        print(name,flush=True)
    for density in (200,400):
        acquire(f'noise_{density}',f'Is1 0 in 0 AC 1\n.noise V(out) Is1 dec {density} .1 1Meg')
    def transient(n,r):
        name=f'ions_N{n}_r{r}'
        with h5py.File(raw) as h:
            if name in h:return
        filename=name+'.pwl'
        np.savetxt(decks/filename,np.c_[rawpwl[:,0],n*rawpwl[:,1]],fmt='%.17e')
        acquire(name,f'Is1 0 in PWL FILE="{filename}"\n.tran 0 .032 0 {5e-6/r:.17e}\n.save V(in) V(out) I(Is1)')
    for n in (0,1):
        for r in (1,2):transient(n,r)
    for iteration in range(12):
        subprocess.run([args.julia,'--project=.',str(ROOT/'scripts/analyze_spice2_enc.jl')],cwd=ROOT,check=True)
        with h5py.File(out/'analysis.h5') as h:pending=h['pending_counts'][:]
        if len(pending)==0:break
        for n in pending:
            for r in (1,2):transient(int(n),r)
    else:raise RuntimeError('Threshold search did not converge')
if __name__=='__main__':main()
