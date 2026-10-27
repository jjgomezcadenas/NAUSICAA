#!/usr/bin/env python3
"""Acquire LTspice transients for configured buckets; append only missing cases.

Julia supplies physical input currents. Raw stationary spectra are copied from
an identical circuit's saved noise run. Existing records are checked by hash.
"""
from pathlib import Path
import argparse, hashlib, re, tomllib
import h5py
import numpy as np
from build_ltspice_study import run
ROOT=Path(__file__).resolve().parents[1]
def string(value): return value.decode() if isinstance(value, bytes) else value
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def preserve(path, data):
    if path.exists():
        if path.read_bytes()!=data: raise ValueError(f'Refusing to overwrite {path}')
    else: path.write_bytes(data)
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--simulator',required=True,type=Path)
    p.add_argument('--bucket-ms',type=int,nargs='+')
    p.add_argument('--counts',type=int,nargs='+')
    args=p.parse_args(); sim=args.simulator.resolve()
    cfgpath=ROOT/'config/bucket_study.toml';cfg=tomllib.loads(cfgpath.read_text())
    source=ROOT/cfg['circuit'];model=ROOT/cfg['model'];noise=ROOT/cfg['noise_source']
    body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('is1 ','.tran','.save','.end')))+'\n'
    # Verify that reused noise actually belongs to the same transistor circuit.
    original_noise=ROOT/'results/ltspice/m3_scan/decks'/f"{cfg['noise_group']}_noise_400.cir"
    noisebody='\n'.join(x for x in original_noise.read_text().splitlines() if not x.lower().startswith(('is1 ','.noise','.end')))+'\n'
    assert body==noisebody, 'Noise and transient circuits differ'
    assert sha(model)==sha(original_noise.parent/'corrected.inc')
    for frame in cfg['bucket_sizes_s']:
        ms=round(frame*1000)
        if args.bucket_ms and ms not in args.bucket_ms: continue
        out=ROOT/cfg['output']/f'{ms}ms';decks=out/'decks';decks.mkdir(parents=True,exist_ok=True)
        inputs=out/'input.h5';raw=out/'raw.h5'
        preserve(out/'configuration.toml',cfgpath.read_bytes())
        preserve(decks/'corrected.inc',model.read_bytes())
        with h5py.File(inputs) as inp, h5py.File(raw,'a') as h:
            assert string(inp.attrs['config_sha256'])==sha(cfgpath)
            assert string(inp.attrs['detector_sha256'])==sha(ROOT/cfg['detector'])
            provenance={'config_sha256':sha(cfgpath),'input_sha256':sha(inputs),'circuit_sha256':sha(source),'model_sha256':sha(model),'noise_source_sha256':sha(noise),'simulator_sha256':sha(sim),'acquisition_sha256':sha(Path(__file__).resolve()),'parser_sha256':sha(ROOT/'src/ion_signal/spice.py')}
            for k,v in provenance.items():
                if k in h.attrs and h.attrs[k]!=v: raise ValueError(f'{raw}: provenance mismatch {k}')
                h.attrs[k]=v
            for k,v in inp.attrs.items():
                if k not in h.attrs:h.attrs[k]=v
            h.attrs.update(external_capacitance_F=1.25e-12,feedback_F=3e-15,temperature_K=293.15,frequency_low_Hz=.1,frequency_high_Hz=1e6,noise_source=cfg['noise_source'],noise_group=cfg['noise_group'],simulator=str(sim))
            for density in (200,400):
                name=f'noise_{density}'
                if name not in h:
                    with h5py.File(noise) as old:old.copy(old[cfg['noise_group']][name],h,name=name)
            for r,step in enumerate(cfg['max_steps_s'],1):
                cases=[('ions',n) for n in (args.counts if args.counts is not None else cfg['counts'])]
                if args.counts is None:cases += [('train',cfg['train_count'])]
                for kind,n in cases:
                    if n<0:raise ValueError('Negative ion count')
                    name=f'{kind}_N{n}_r{r}'
                    pwl=decks/f'{name}.pwl'
                    t=inp[f'{kind}_r{r}/time_s'][:];i=inp[f'{kind}_r{r}/current_A_per_ion'][:]*n
                    # Formatting/acquisition only; physical waveform was calculated in Julia.
                    data=''.join(f'{x:.17e} {y:.17e}\n' for x,y in zip(t,i)).encode()
                    preserve(pwl,data)
                    path=decks/f'{name}.cir'
                    tail=f'Is1 0 in PWL FILE="{pwl.name}"\n.tran 0 {t[-1]:.17e} 0 {step:.17e}\n.save V(in) V(out) I(Is1)\n.end\n'
                    preserve(path,(body+tail).encode())
                    if name in h:
                        assert h[name].attrs['deck_sha256']==sha(path)
                        assert h[name].attrs['pwl_sha256']==sha(pwl)
                        continue
                    d,_=run(sim,path);log=path.with_suffix('.run.txt').read_text()
                    assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I)
                    group=h.create_group(name)
                    group.attrs.update(deck_sha256=sha(path),pwl_sha256=sha(pwl),log=log)
                    for k,v in d.items():group.create_dataset(k,data=np.asarray(v).real,compression='gzip',shuffle=True)
                    h.flush();print(f'{ms} ms {name}',flush=True)
if __name__=='__main__':main()
