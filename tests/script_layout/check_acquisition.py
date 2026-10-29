"""Compare original and extracted LTspice acquisition on identical frozen inputs."""
from pathlib import Path
import ast
import json
import shutil
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from ion_signal.spice import read_raw
from ion_signal.spice_acquisition import run

manifest = json.loads((ROOT/'docs/script_paths.json').read_text())
source_path = next(r['old'] for r in manifest['moves'] if r['new'].endswith('/initial/build_ltspice_study.py'))
source = subprocess.check_output(['git','show',manifest['source_commit']+':'+source_path],cwd=ROOT).decode()
node = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == 'run')
scope = dict(subprocess=subprocess, read_raw=read_raw)
exec(compile(ast.Module(body=[node], type_ignores=[]),'<original runner>','exec'),scope)
original = scope['run']
cases = ROOT/'results/ltspice/pad_compression_direct/8b0d8b859f56/cases'
case = next(p for p in sorted(cases.glob('p2_bb0nu_*')) if (p/'history.cir').exists())
work = ROOT/'work/reorg_acquisition'
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--simulator',type=Path,required=True)
sim=parser.parse_args().simulator.resolve()
answers=[]
for label,runner in (('original',original),('shared',run)):
    folder=work/label;folder.mkdir(parents=True,exist_ok=True)
    for name in ('history.cir','history.pwl','corrected.inc'):
        shutil.copy2(case/name,folder/name)
    wave,_=runner(sim,folder/'history.cir')
    assert (folder/'history.run.txt').exists()
    answers.append(wave)
assert answers[0].keys()==answers[1].keys()
for key in answers[0]:
    np.testing.assert_array_equal(answers[0][key],answers[1][key])
print('Fresh original/shared LTspice acquisitions give identical waveform arrays:',case.name)
