#!/usr/bin/env python3
"""Acquire separate v2 reproduction and matched ion/noise simulations."""
from pathlib import Path
import argparse, hashlib, re, tomllib
import h5py
import numpy as np
from build_ltspice_study import run
from ion_signal.spice import read_raw

ROOT = Path(__file__).resolve().parents[1]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def preserve(p, text):
    if p.exists():
        assert p.read_text() == text, f'Refusing overwrite {p}'
    else: p.write_text(text)
def save(h, name, data, **attrs):
    g=h.create_group(name);g.attrs.update(attrs)
    for k,v in data.items():g.create_dataset(k,data=np.asarray(v).real,compression='gzip',compression_opts=9,shuffle=True)
    h.file.flush()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--simulator',type=Path,required=True)
    ap.add_argument('--variants',nargs='+');ap.add_argument('--counts',nargs='+',type=int)
    ap.add_argument('--noise-only',action='store_true')
    ap.add_argument('--config',default='config/rogelio_v2.toml')
    a=ap.parse_args();sim=a.simulator.resolve()
    cp=ROOT/a.config;c=tomllib.loads(cp.read_text())
    out=ROOT/c['output'];out.mkdir(parents=True,exist_ok=True)
    preserve(out/'configuration.toml',cp.read_text())
    original=(ROOT/c['source']).read_text()
    body='* Derived v2 circuit; first line is the SPICE title\n'+'\n'.join(x for x in original.splitlines() if x and not x.startswith('*') and
                  not x.lower().startswith(('.lib','.model','.noise','.backanno','.end','.meas','is1 ','vion1 ','b1 ')))+'\n'
    options='.options reltol=1e-6 abstol=1e-19 vntol=1e-10 numdgt=15 plotwinsize=0\n'
    with h5py.File(out/'raw.h5','a') as h, h5py.File(ROOT/c['input']) as inp:
        for key,path in [('config',cp),('input',ROOT/c['input']),('source',ROOT/c['source']),
                         ('source_model',ROOT/c['source_model']),('corrected_model',ROOT/c['corrected_model']),('simulator',sim)]:
            value=sha(path)
            if key+'_sha256' in h.attrs:assert h.attrs[key+'_sha256']==value
            h.attrs[key+'_sha256']=value
        h.attrs['acquisition_sha256']=sha(Path(__file__))
        h.attrs['assumptions']=inp.attrs['assumptions']
        h.attrs['single_ion_charge_C']=inp.attrs['single_ion_charge_C']
        def acquire(parent,name,deck,allow_warning=False):
            if name in parent:
                assert parent[name].attrs['deck_sha256']==sha(deck)
                return
            data,header=run(sim,deck);log=deck.with_suffix('.run.txt').read_text()
            if not allow_warning:assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I),log
            save(parent,name,data,deck_sha256=sha(deck),log=log,header=header)
            print(parent.name,name,flush=True)
        # Matched to Rogelio's own 1 Hz--1 GHz supplied spectrum.
        if not a.variants and not c.get('retuned',False):
            rdir=out/'decks'/'reproduction';rdir.mkdir(parents=True,exist_ok=True)
            preserve(rdir/'original.inc',(ROOT/c['source_model']).read_text())
            preserve(rdir/'corrected.inc',(ROOT/c['corrected_model']).read_text())
            group=h.require_group('reproduction')
            for k in ('original','corrected'):
                deck=rdir/f'{k}.cir'
                preserve(deck,body+f'.lib {k}.inc\nIs1 0 in 0 AC 1\n.noise V(out) Is1 dec 1000 1 1G\n.end\n')
                acquire(group,k,deck,allow_warning=k=='original')
            for name,filename in [('supplied_noise','CSAfolded180nm-2.raw'),('supplied_op','CSAfolded180nm-2.op.raw')]:
                source=(ROOT/c['source']).parent/filename
                if name not in group:
                    data,header=read_raw(source);save(group,name,data,source_sha256=sha(source),header=header)
        for v in c['variants']:
            name=v['name']
            if a.variants and name not in a.variants:continue
            g=h.require_group(name);g.attrs.update(v)
            folder=out/'decks'/name;folder.mkdir(parents=True,exist_ok=True)
            preserve(folder/'corrected.inc',(ROOT/c['corrected_model']).read_text())
            circuit=re.sub(r'^Cd1 .*$',f"Cd1 in 0 {v['capacitance_F']:.17e}",body,flags=re.M)
            if 'm4_gate_V' in v:circuit=circuit.replace('Vp1 Vncasc 0 0.64',f"Vp1 Vncasc 0 {v['m4_gate_V']:.17e}")
            if v.get('m3_large'):
                circuit=circuit.replace('.param Lnload=0.5u','.param Lnload=2u').replace('.param Wnload=50u','.param Wnload=200u')
                circuit=circuit.replace('Vnload Vnload 0 0.6','Vnload Vnload 0 0.538477478921')
            circuit+=f".lib corrected.inc\n.temp {v['temperature_C']}\n"+options
            deck=folder/'op.cir';preserve(deck,circuit+'Is1 0 in 0\n.op\n.end\n');acquire(g,'op',deck)
            for density in (200,400):
                deck=folder/f'noise_{density}.cir'
                preserve(deck,circuit+f'Is1 0 in 0 AC 1\n.noise V(out) Is1 dec {density} .1 1Meg\n.end\n')
                acquire(g,f'noise_{density}',deck)
            if a.noise_only:continue
            for r,step in enumerate(c['max_steps_s'],1):
                t=inp[f'r{r}/time_s'][:];current=inp[f'r{r}/single_A_per_ion'][:]
                for n in a.counts if a.counts is not None else c['counts']:
                    case=f'single_N{n}_r{r}';pwl=folder/f'{case}.pwl'
                    preserve(pwl,''.join(f'{x:.17e} {y:.17e}\n' for x,y in zip(t,n*current)))
                    deck=folder/f'{case}_tolerance2.cir'
                    tol='.options reltol=1e-5 abstol=1e-17 vntol=1e-8' if r==1 else '.options reltol=5e-6 abstol=5e-18 vntol=5e-9'
                    transient=circuit.replace(options,tol+' numdgt=15 plotwinsize=0\n')
                    preserve(deck,transient+f'Is1 0 in PWL FILE="{pwl.name}"\n.tran 0 {t[-1]:.17e} 0 {step:.17e}\n.save V(in) V(out) V(med) V(vd6) I(Is1)\n.end\n')
                    acquire(g,case,deck)
if __name__=='__main__':main()
