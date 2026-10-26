#!/usr/bin/env python3
"""Package git-tracked calculation sources and results for collaborator review."""
from pathlib import Path
import subprocess
import zipfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'notes/physics/work/handoff/calculos_CSA_Rogelio.zip'
PREFIXES=('LTSpice/','src/','scripts/','tests/','config/','results/ltspice/')
FILES={'requirements.txt','pyproject.toml','Makefile'}

def main():
    paths=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0')
    chosen=sorted(p for p in paths if p and (p.startswith(PREFIXES) or p in FILES))
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(OUT,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for name in chosen:
            path=ROOT/name
            if path.is_file():archive.write(path,arcname='calculos_CSA/'+name)
    with zipfile.ZipFile(OUT) as archive:
        assert archive.testzip() is None
        assert 'calculos_CSA/LTSpice/calculos_para_rogelio.md' in archive.namelist()
    print(f'{OUT}: {OUT.stat().st_size} bytes, {len(chosen)} files')
    # Smaller review package retains reports, spectra and the inputs needed by
    # the three reproduction commands; regenerated bulk transients stay in full ZIP.
    required={'results/ltspice/study.json', 'results/ltspice/noise_200.csv',
        'results/ltspice/noise_400.csv', 'results/ltspice/decks/uniform_N1_r3.cir',
        'results/ltspice/decks/uniform_r2.pwl', 'results/ltspice/uniform_N1_r1.csv.gz',
        'results/ltspice/uniform_N1_r2.csv.gz'}
    compact=OUT.with_name('calculos_CSA_Rogelio_revision.zip')
    with zipfile.ZipFile(compact,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for name in chosen:
            path=ROOT/name
            keep=(not name.startswith('results/') or name in required
                or name.startswith('results/ltspice/noise_audit/')
                or path.suffix in ('.md','.json','.pdf')
                or path.name in ('shaping_scan.csv','comparison.csv'))
            if keep and path.is_file():archive.write(path,arcname='calculos_CSA/'+name)
    with zipfile.ZipFile(compact) as archive:assert archive.testzip() is None
    print(f'{compact}: {compact.stat().st_size} bytes')

if __name__=='__main__':main()
