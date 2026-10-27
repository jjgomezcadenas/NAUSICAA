#!/usr/bin/env python3
"""Acquire extended LTspice decay tails and a specified sequence of ion packets."""
from pathlib import Path
import argparse,hashlib,re,tomllib
import h5py
import numpy as np
from build_ltspice_study import run
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def string(x):return x.decode() if isinstance(x,bytes) else x
def preserve(p,data):
    if p.exists():
        if p.read_bytes()!=data:raise ValueError(f'Refusing overwrite: {p}')
    else:p.write_bytes(data)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--simulator',required=True,type=Path);p.add_argument('--counts',nargs='+',type=int)
    args=p.parse_args();sim=args.simulator.resolve()
    cp=ROOT/'config/discharge_calibration.toml';cfg=tomllib.loads(cp.read_text())
    out=ROOT/cfg['output'];decks=out/'decks';decks.mkdir(parents=True,exist_ok=True)
    source=ROOT/cfg['circuit'];model=ROOT/cfg['model'];history=ROOT/cfg['historical']/'raw.h5'
    body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('is1 ','.tran','.save','.end')))+'\n'
    preserve(out/'configuration.toml',cp.read_bytes());preserve(decks/'corrected.inc',model.read_bytes())
    with h5py.File(out/'input.h5') as inp,h5py.File(out/'raw.h5','a') as h,h5py.File(history) as old:
        assert string(inp.attrs['config_sha256'])==sha(cp)
        assert old.attrs['circuit_sha256']==sha(source) and old.attrs['model_sha256']==sha(model)
        provenance={'config_sha256':sha(cp),'input_sha256':sha(out/'input.h5'),'circuit_sha256':sha(source),'model_sha256':sha(model),'noise_source_sha256':sha(history),'simulator_sha256':sha(sim),'acquisition_sha256':sha(Path(__file__).resolve()),'parser_sha256':sha(ROOT/'src/ion_signal/spice.py')}
        for k,v in provenance.items():
            if k in h.attrs and string(h.attrs[k])!=v:raise ValueError(f'Provenance mismatch {k}')
            h.attrs[k]=v
        for k,v in inp.attrs.items():
            if k not in h.attrs:h.attrs[k]=v
        for name in ('noise_200','noise_400'):
            if name not in h:old.copy(old[name],h,name=name)
        for r,step in enumerate(cfg['max_steps_s'],1):
            cases=[('single',n) for n in (cfg['counts'] if args.counts is None else args.counts)]
            if args.counts is None:cases.append(('train',1))
            for kind,n in cases:
                assert n>=0
                name=f'single_N{n}_r{r}' if kind=='single' else f'train_r{r}'
                t=inp[f'r{r}/time_s'][:]
                current=n*inp[f'r{r}/single_A_per_ion'][:] if kind=='single' else inp[f'r{r}/train_A'][:]
                pwl=decks/f'{name}.pwl';preserve(pwl,''.join(f'{x:.17e} {y:.17e}\n' for x,y in zip(t,current)).encode())
                deck=decks/f'{name}.cir'
                preserve(deck,(body+f'Is1 0 in PWL FILE="{pwl.name}"\n.tran 0 {t[-1]:.17e} 0 {step:.17e}\n.save V(in) V(out) I(Is1)\n.end\n').encode())
                if name in h:
                    assert h[name].attrs['deck_sha256']==sha(deck) and h[name].attrs['pwl_sha256']==sha(pwl)
                    continue
                d,_=run(sim,deck);log=deck.with_suffix('.run.txt').read_text()
                assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I)
                g=h.create_group(name);g.attrs.update(deck_sha256=sha(deck),pwl_sha256=sha(pwl),log=log)
                for key,val in d.items():g.create_dataset(key,data=np.asarray(val).real,compression='gzip',shuffle=True)
                h.flush();print(name,flush=True)
if __name__=='__main__':main()
