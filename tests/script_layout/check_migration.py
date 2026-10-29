"""Check moved command paths, shared helper behaviour and CLI imports.

Run from the repository root with the scientific Python environment.
The subprocess checks request help and perform no acquisition or PDF compilation.
"""
from pathlib import Path
import ast
import json
import os
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]


def function(tree, name):
    return next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)


def main():
    manifest = json.loads((ROOT/'docs/script_paths.json').read_text())
    moves = manifest['moves']
    for item in moves:
        assert not (ROOT/item['old']).exists(), item['old']
        assert (ROOT/item['new']).is_file(), item['new']
    assert len({item['new'] for item in moves}) == len(moves)
    assert not list((ROOT/'scripts').glob('*.py'))
    assert not list((ROOT/'scripts').glob('*.jl'))
    print(len(moves), 'relocated paths exist, with no old top-level commands')

    for original, current, names in (
        ('scripts/build_ltspice_study.py', 'src/ion_signal/spice_acquisition.py', ('run',)),
        ('scripts/run_rogelio_v2.py', 'src/ion_signal/spice_storage.py', ('sha', 'preserve', 'save')),
    ):
        old = ast.parse(subprocess.check_output(['git','show',manifest['source_commit']+':'+original],cwd=ROOT).decode())
        new = ast.parse((ROOT/current).read_text())
        for name in names:
            assert ast.dump(function(old,name)) == ast.dump(function(new,name)), name
    print('Shared acquisition and storage functions preserve their original syntax trees')

    commands = list((ROOT/'scripts').rglob('*.py'))
    for path in commands:
        ast.parse(path.read_text(), filename=str(path))
    known = {item['old'] for item in moves}
    files = commands + [ROOT/'Makefile']
    files += list((ROOT/'src').rglob('*.py')) + list((ROOT/'src').rglob('*.jl'))
    files += list((ROOT/'tests').rglob('*.py')) + list((ROOT/'tests').rglob('*.jl'))
    stale = []
    for path in files:
        if path.resolve() == Path(__file__).resolve():
            continue
        for value in re.findall(r'scripts/[A-Za-z0-9_./-]+[.](?:py|jl)', path.read_text()):
            if value in known:
                stale.append((str(path.relative_to(ROOT)),value))
    assert not stale, stale
    print(len(commands), 'Python command files parse; active callers contain no old command paths')

    selected = (
        'viewers/view_events.py', 'viewers/serve_cloud_viewer.py', 'viewers/view_pad_histories.py',
        'electronics/histories/run_pad_history_pipeline.py',
        'electronics/histories/run_pad_compression.py',
        'electronics/histories/run_pad_history_spice.py',
        'electronics/histories/cache_standard_pad_histories.py',
        'electronics/feedback/run_feedback_pipeline.py',
        'electronics/feedback/run_feedback_candidate.py',
        'electronics/bias/run_bias_sensitivity.py',
        'electronics/reproduction/initial/build_ltspice_study.py',
        'electronics/reproduction/spice2/run_spice2_enc.py',
        'electronics/reproduction/rogelio_v3/run_rogelio_v3.py',
        'documentation/build_event_review.py',
    )
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLCONFIGDIR=str(ROOT/'work/reorg-matplotlib'))
    env.pop('PYTHONPATH',None)
    for relative in selected:
        result = subprocess.run([sys.executable,str(ROOT/'scripts'/relative),'--help'],cwd='/private/tmp',env=env,capture_output=True,text=True,timeout=60)
        assert result.returncode == 0, (relative,result.stdout,result.stderr)
        assert 'usage:' in result.stdout.lower(), relative
    print(len(selected), 'CLI entry points resolve imports from another working directory')


if __name__ == '__main__':
    main()
