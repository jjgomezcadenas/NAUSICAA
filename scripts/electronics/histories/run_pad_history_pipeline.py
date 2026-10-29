#!/usr/bin/env python3
"""Run or resume the complete stored-pad-history electronics study."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tomllib
import h5py

ROOT = Path(__file__).resolve().parents[3]

def sha(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()

def value(attr):
    return attr.decode() if isinstance(attr, bytes) else str(attr)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--simulator', required=True, type=Path)
    parser.add_argument('--julia', default='julia')
    parser.add_argument('--config',type=Path,default=ROOT/'config/pad_history_readout_direct.toml')
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = tomllib.loads(config_path.read_text())
    out = ROOT/config['output']
    out.mkdir(parents=True, exist_ok=True)
    code = ROOT/'src/ion_signal/pad_history_readout.jl'
    sources = [config_path, code, ROOT/config['circuit_results'], ROOT/config['circuit_analysis'],
               ROOT/config['source_decks']/'corrected.inc',
               ROOT/config['source_decks']/'40ms_single_N400_r2.cir']
    sources += [ROOT/s['path'] for s in config['samples']]
    manifest = {str(p.relative_to(ROOT)):sha(p) for p in sources}
    manifest_path = out/'pipeline_inputs.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise RuntimeError('Inputs changed; select a new output directory to preserve this study.')
    env = os.environ.copy()
    env['PYTHONPATH'] = str(ROOT/'src')
    env['NAUSICAA_PAD_HISTORY_CARD'] = str(config_path)
    env.setdefault('JULIA_DEPOT_PATH', str(ROOT/'notes/physics/work/julia-depot')+':'+str(Path.home()/'.julia'))
    env.setdefault('MPLCONFIGDIR', str(ROOT/'notes/physics/work/matplotlib'))
    def run(*command):
        subprocess.run(list(map(str,command)),cwd=ROOT,env=env,check=True)
    screens = [out/(s['name']+'.h5') for s in config['samples']]
    if any(p.exists() for p in screens) and not all(p.exists() for p in screens):
        raise RuntimeError('Partial linear screen; use a separate output directory for a clean run.')
    if not all(p.exists() for p in screens):
        run(args.julia,'--project=.','scripts/electronics/histories/screen_pad_histories.jl')
    for sample, path in zip(config['samples'], screens):
        with h5py.File(path) as h:
            for key, expected in {'config_sha256':sha(config_path), 'code_sha256':sha(code),
                                  'source_sha256':manifest[sample['path']]}.items():
                if value(h.attrs[key]) != expected:
                    raise RuntimeError(f'{path}: mismatched {key}')
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    inputs = out/'spice_inputs.h5'
    if not inputs.exists():
        run(args.julia,'--project=.','scripts/electronics/histories/build_pad_history_spice_inputs.jl')
    with h5py.File(inputs) as h:
        if value(h.attrs['generator_sha256']) != sha(code):
            raise RuntimeError('PWL generator changed.')
        for group in h.values():
            sample = value(group.attrs['sample'])
            if value(group.attrs['source_sha256']) != sha(out/(sample+'.h5')):
                raise RuntimeError('Linear screen changed since PWL generation.')
    run(sys.executable,'scripts/electronics/histories/run_pad_history_spice.py','--simulator',args.simulator.resolve(),'--config',config_path)
    with h5py.File(inputs) as h:
        for name in h:
            with h5py.File(out/'spice'/(name+'.h5')) as raw:
                required = {'baseline_r1','history_r1','baseline_r2','history_r2'}
                if not required.issubset(raw) or 'acquisition_failure' in raw.attrs:
                    raise RuntimeError(f'Incomplete SPICE acquisition: {name}')
    run(args.julia,'--project=.','scripts/electronics/histories/analyze_pad_history_spice.jl')
    run(sys.executable,'scripts/electronics/histories/present_pad_histories.py','--config',config_path)

if __name__ == '__main__':
    main()
