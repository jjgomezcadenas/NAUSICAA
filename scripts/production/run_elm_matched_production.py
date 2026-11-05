#!/usr/bin/env python3
"""Run coherent detector placement, both electronics arms and the combined graph."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,os,subprocess,time,tomllib
ROOT=Path(__file__).resolve().parents[2]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=['smoke','full'],default='full');p.add_argument('--julia',default=str(Path.home()/'.juliaup/bin/julia'));a=p.parse_args()
    out=ROOT/'results'/('elm_matched_smoke' if a.kind=='smoke' else 'elm_matched_production');out.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env['JULIA_DEPOT_PATH']=':'.join(map(str,[ROOT/'work/julia-depot',Path.home()/'.julia']))
    def command(script,*args):return [a.julia,'--startup-file=no','--project=.',str(ROOT/'scripts/production'/script),*args]
    started=time.monotonic();print(f'Matched production PID {os.getpid()}, kind={a.kind}',flush=True)
    subprocess.run(command('run_matched_detector.jl',a.kind),cwd=ROOT,env=env,check=True)
    detector_s=time.monotonic()-started
    def stage(name,cmd):
        begin=time.monotonic();log=out/(name+'.log')
        with log.open('w') as stream:
            process=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT)
            print(f'{name} PID {process.pid}; log={log}',flush=True)
            if process.wait():raise RuntimeError(f'{name} failed; inspect {log}')
        elapsed=time.monotonic()-begin;print(f'{name} complete: {elapsed:.2f} s',flush=True);return elapsed
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        sipm=pool.submit(stage,'sipm',command('run_matched_sipm.jl',a.kind))
        csa=pool.submit(stage,'csa',command('run_ion_csa.jl',f'config/specs/runs/elm_matched_csa_{a.kind}.toml'))
        sipm_s=sipm.result();csa_s=csa.result()
    sm=tomllib.loads((out/'sipm_files.toml').read_text())['manifest']
    subprocess.run(command('assemble_matched_response.jl',a.kind,sm),cwd=ROOT,env=env,check=True)
    (out/'run_timings.toml').write_text(f'complete = true\ndetector_wall_s = {detector_s:.9f}\nsipm_wall_s = {sipm_s:.9f}\ncsa_wall_s = {csa_s:.9f}\ntotal_wall_s = {time.monotonic()-started:.9f}\nrunner_sha256 = "{hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}"\n')
    print(f'Matched production complete: {out / "manifest.toml"}',flush=True)
if __name__=='__main__':main()
