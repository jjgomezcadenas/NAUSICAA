#!/usr/bin/env python3
"""Acquire and reproduce the isolated longer-retention CSA study."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
import argparse,hashlib,json,os,re,shutil,subprocess,time,tomllib
import h5py,numpy as np
from ion_signal.spice_acquisition import run

def sha(path):
    with Path(path).open('rb') as source:return hashlib.file_digest(source,'sha256').hexdigest()

def preserve(path,data):
    path=Path(path);data=data.encode() if isinstance(data,str) else data
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists() and path.read_bytes()!=data:raise ValueError(f'Preserving existing {path}')
    if not path.exists():path.write_bytes(data)

def body(c,variant,gate):
    text='\n'.join(line for line in (ROOT/c['source_deck']).read_text().splitlines()
                   if not line.lower().startswith(('is1 ','.tran','.noise','.save','.end')))+'\n'
    replacements=((r'^Vp1 Vncasc 0 .*$',f'Vp1 Vncasc 0 {gate:.12g}'),
                  (r'^Vf1 Vdis 0 .*$',f'Vf1 Vdis 0 {c["feedback_gate_V"]:.12g}'),
                  (r'^\.parameter WMf=.*$',f'.parameter WMf={c["feedback_width_m"]:.17e}'),
                  (r'^\.parameter LMf=.*$',f'.parameter LMf={variant["length_m"]:.17e}'))
    for pattern,value in replacements:
        text,n=re.subn(pattern,value,text,flags=re.M);assert n==1,(pattern,n)
    return text

def archive_sources(c):
    sources=[Path(__file__),ROOT/'src/ion_signal/retention_scan.jl',
             ROOT/'src/ion_signal/spice_acquisition.py',ROOT/'src/ion_signal/bucket_study.jl',
             ROOT/'src/ion_signal/discharge_calibration.jl',ROOT/'src/ion_signal/averaging_scan.jl']
    out=ROOT/c['output']
    for source in sources:
        preserve(out/'source_snapshot'/sha(source)/source.name,source.read_bytes())

def acquire(c,cp,variant,gate,simulator,op_only=False):
    out=ROOT/c['output'];name=variant['name']
    if op_only:
        token=hashlib.sha256(f'{name}/{gate:.17e}'.encode()).hexdigest()[:14]
        folder=out/'tuning'/token;path=folder/'op.h5'
    else:folder=out/'decks'/name;path=out/'variants'/(name+'.h5')
    folder.mkdir(parents=True,exist_ok=True);path.parent.mkdir(parents=True,exist_ok=True)
    preserve(folder/'corrected.inc',(ROOT/c['model']).read_bytes())
    circuit=body(c,variant,gate)
    with h5py.File(path,'a') as h:
        attributes=dict(name=name,length_m=variant['length_m'],width_m=c['feedback_width_m'],
                        gate_V=c['feedback_gate_V'],m4_gate_V=gate,config_sha256=sha(cp),
                        model_sha256=sha(ROOT/c['model']),source_sha256=sha(ROOT/c['source_deck']),
                        simulator_sha256=sha(simulator))
        for key,value in attributes.items():
            if key in h.attrs:assert h.attrs[key]==value,(name,key)
            else:h.attrs[key]=value
        def case(key,tail):
            deck=folder/(key.replace('/','_')+'.cir');preserve(deck,circuit+tail+'\n.end\n')
            if key in h:
                assert h[key].attrs['deck_sha256']==sha(deck);return
            data,_=run(simulator,deck);log=deck.with_suffix('.run.txt').read_text()
            assert not re.search(r'unknown|unrecognized|fatal|syntax error',log,re.I),log
            group=h.create_group(key)
            group.attrs.update(deck_sha256=sha(deck),acquisition_sha256=sha(Path(__file__)),
                               runner_sha256=sha(ROOT/'src/ion_signal/spice_acquisition.py'))
            for label,values in data.items():group.create_dataset(label,data=np.asarray(values).real,compression='gzip',shuffle=True)
            h.flush()
            for extension in ('.raw','.log','.op.raw'):deck.with_suffix(extension).unlink(missing_ok=True)
        case('op','Is1 0 in 0\n.op')
        h.attrs['baseline_V']=float(h['op/V(out)'][0])
        if op_only:return path
        assert c['minimum_output_V']<h.attrs['baseline_V']<c['maximum_output_V']
        for density in c['noise_densities']:
            case('noise_'+str(density),f'Is1 0 in 0 AC 1\n.noise V(out) Is1 dec {density} {c["noise_low_Hz"]:.12g} {c["noise_high_Hz"]:.12g}')
        with h5py.File(out/'inputs.h5') as inp:
            h.attrs['inputs_sha256']=sha(out/'inputs.h5')
            for frame in c['frames_s']:
                key=f'{round(frame*1000)}ms'
                for r,step in enumerate(c['max_steps_s'],1):
                    t=inp[f'{key}/fast_r{r}/time_s'][:];unit=inp[f'{key}/fast_r{r}/single_A_per_ion'][:]
                    for n in c['counts']:
                        pwl=out/'inputs'/f'{key}_N{n}_r{r}.pwl'
                        if not pwl.exists():
                            if n==0:preserve(pwl,f'{t[0]:.17e} 0\n{t[-1]:.17e} 0\n')
                            else:preserve(pwl,''.join(f'{x:.17e} {n*y:.17e}\n' for x,y in zip(t,unit)))
                        case(f'{key}/single_N{n}_r{r}',f'Is1 0 in PWL FILE="../../inputs/{pwl.name}"\n.tran 0 {t[-1]:.17e} 0 {step:.17e}\n.save V(in) V(out) I(Is1)')
                    tailt=inp[f'{key}/tail_r{r}/time_s'][:];taili=inp[f'{key}/tail_r{r}/single_A_per_ion'][:]
                    for n in (0,10,100):
                        pwl=out/'inputs'/f'{key}_tail_N{n}_r{r}.pwl'
                        if not pwl.exists():
                            if n==0:preserve(pwl,f'{tailt[0]:.17e} 0\n{tailt[-1]:.17e} 0\n')
                            else:preserve(pwl,''.join(f'{x:.17e} {n*y:.17e}\n' for x,y in zip(tailt,taili)))
                        case(f'{key}/tail/single_N{n}_r{r}',f'Is1 0 in PWL FILE="../../inputs/{pwl.name}"\n.tran 0 {tailt[-1]:.17e} 0 {c["tail_max_steps_s"][r-1]:.17e}\n.save V(in) V(out) I(Is1)')
                print(name,key,'acquired',flush=True)
    return path

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=ROOT/'config/retention_scan.toml')
    parser.add_argument('--simulator',required=True,type=Path);parser.add_argument('--julia',default=str(Path.home()/'.juliaup/bin/julia'))
    parser.add_argument('--stage',choices=['all','calibrate','inputs','tune','op','acquire','analyze','histories','summary','present'],default='all')
    parser.add_argument('--variant');parser.add_argument('--m4-gate',type=float);parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();cp=args.config.resolve();c=tomllib.loads(cp.read_text());sim=args.simulator.resolve()
    out=ROOT/c['output'];out.mkdir(parents=True,exist_ok=True)
    preserve(out/'configuration.toml',cp.read_bytes())
    preserve(out/'.gitignore','*.raw\n*.log\nvariants/*.h5\ntuning/\ninputs/*.pwl\nscreens/**/*.h5\ncompression/**/cases/\n')
    if args.stage=='op':
        variant=next(v for v in c['variants'] if v['name']==args.variant)
        print(acquire(c,cp,variant,args.m4_gate,sim,True).relative_to(ROOT));return
    archive_sources(c)
    env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'src')
    env.setdefault('JULIA_DEPOT_PATH',str(ROOT/'notes/physics/work/julia-depot')+':'+str(Path.home()/'.julia'))
    env.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/retention_matplotlib'))
    def command(*parts):subprocess.run(list(map(str,parts)),cwd=ROOT,env=env,check=True)
    def julia(stage):
        command(args.julia,'--startup-file=no','--project=.',ROOT/'scripts/electronics/feedback/retention_scan.jl',stage,cp,sys.executable,sim)
    if args.stage in ('all','calibrate','inputs'):julia('inputs')
    if args.stage in ('all','calibrate','tune'):julia('tune')
    if args.stage in ('all','calibrate','acquire'):
        with h5py.File(out/'operating_points.h5') as points:
            for variant in c['variants']:
                if args.variant and variant['name']!=args.variant:continue
                acquire(c,cp,variant,float(points[variant['name']].attrs['m4_gate_V']),sim)
    if args.stage in ('all','calibrate','analyze'):julia('analyze')
    if args.stage in ('all','histories'):
        julia('screens')
        for variant in c['variants']:
            command(sys.executable,ROOT/'scripts/electronics/histories/run_pad_compression.py',
                    '--config',out/'cards'/(variant['name']+'_viewer.toml'),'--simulator',sim,'--workers',args.workers,
                    '--work-dir',ROOT/'notes/physics/work'/('retention_'+out.name+'_'+variant['name']))
    if args.stage in ('all','summary'):julia('summary')
    if args.stage in ('all','present'):
        command(sys.executable,ROOT/'scripts/electronics/feedback/present_retention_scan.py','--config',cp)

if __name__=='__main__':main()
