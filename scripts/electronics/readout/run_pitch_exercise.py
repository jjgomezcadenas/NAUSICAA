#!/usr/bin/env python3
"""Run the 4 mm pitch / 40 ms CSA exercise without changing previous studies."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse,hashlib,os,re,subprocess,sys,tomllib
import h5py
import numpy as np
from ion_signal.spice_acquisition import run
ROOT=Path(__file__).resolve().parents[3]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def preserve(p,data):
    if p.exists() and p.read_bytes()!=data:raise ValueError(f'Preserving {p}')
    p.write_bytes(data)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--simulator',required=True,type=Path);p.add_argument('--julia',default='julia');a=p.parse_args()
    cfgfile=ROOT/'config/pitch_4mm_40ms.toml';c=tomllib.loads(cfgfile.read_text());out=ROOT/c['output'];decks=out/'decks';decks.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env.setdefault('JULIA_DEPOT_PATH',str(ROOT/'notes/physics/work/julia-depot')+':'+str(Path.home()/'.julia'));env.setdefault('PYTHONPATH',str(ROOT/'src'))
    def julia(stage):subprocess.run([a.julia,'--project=.','scripts/electronics/readout/pitch_exercise.jl',stage],cwd=ROOT,env=env,check=True)
    if not (out/'input.h5').exists():julia('inputs')
    source=ROOT/c['circuit'];model=ROOT/c['model'];noise=ROOT/c['noise_source'];sim=a.simulator.resolve()
    body='\n'.join(x for x in source.read_text().splitlines() if not x.lower().startswith(('is1 ','.tran','.save','.end')))+'\n'
    preserve(decks/'corrected.inc',model.read_bytes());preserve(out/'configuration.toml',cfgfile.read_bytes())
    (out/'.gitignore').write_text('decks/*.raw\ndecks/*.log\n')
    with h5py.File(out/'input.h5') as inp,h5py.File(out/'raw.h5','a') as h,h5py.File(noise) as old:
        assert old.attrs['circuit_sha256']==sha(source) and old.attrs['model_sha256']==sha(model)
        assert inp.attrs['config_sha256'].decode() == sha(cfgfile) if isinstance(inp.attrs['config_sha256'], bytes) else inp.attrs['config_sha256'] == sha(cfgfile)
        attrs={'config_sha256':sha(cfgfile),'input_sha256':sha(out/'input.h5'),'circuit_sha256':sha(source),'model_sha256':sha(model),'noise_source_sha256':sha(noise),'simulator_sha256':sha(sim),'acquisition_sha256':sha(Path(__file__))}
        for k,v in attrs.items():
            if k in h.attrs and h.attrs[k]!=v:raise ValueError(f'Provenance mismatch {k}')
            h.attrs[k]=v
        for k,v in inp.attrs.items():
            if k not in h.attrs:h.attrs[k]=v
        for name in ('noise_200','noise_400'):
            if name not in h:old.copy(old[name],h,name=name)
        for r,step in enumerate(c['max_steps_s'],1):
            t=inp[f'r{r}/time_s'][:];unit=inp[f'r{r}/single_A_per_ion'][:]
            for n in c['counts']:
                name=f'single_N{n}_r{r}';pwl=decks/(name+'.pwl');deck=decks/(name+'.cir')
                preserve(pwl,''.join(f'{x:.17e} {y:.17e}\n' for x,y in zip(t,n*unit)).encode())
                preserve(deck,(body+f'Is1 0 in PWL FILE="{pwl.name}"\n.tran 0 {t[-1]:.17e} 0 {step:.17e}\n.save V(in) V(out) I(Is1)\n.end\n').encode())
                if name in h:
                    assert h[name].attrs['deck_sha256']==sha(deck) and h[name].attrs['pwl_sha256']==sha(pwl)
                    continue
                d,_=run(sim,deck);log=deck.with_suffix('.run.txt').read_text()
                assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I)
                g=h.create_group(name);g.attrs.update(deck_sha256=sha(deck),pwl_sha256=sha(pwl),log=log)
                for k,v in d.items():g.create_dataset(k,data=np.asarray(v).real,compression='gzip',shuffle=True)
                h.flush();print(name,flush=True)
    julia('analyze');julia('checks')
    subprocess.run([sys.executable,'scripts/electronics/readout/present_pitch_exercise.py'],cwd=ROOT,env=env,check=True)
if __name__=='__main__':main()
