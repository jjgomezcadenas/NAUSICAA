#!/usr/bin/env python3
"""Package git-tracked calculation sources and results for collaborator review."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import subprocess
import zipfile

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'work/handoff/calculos_CSA_Rogelio.zip'
PREFIXES=('LTSpice/','src/','scripts/','tests/','config/','results/ltspice/',
          'assets/','docs/event_analysis_review/content/','docs/script_paths.json',
          'docs/ltspice_paths.json','docs/path_migrations.toml','docs/generated_reports.json',
          'latex_generated/ion_readout/rogelio_notes/')
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
        assert 'calculos_CSA/LTSpice/notes/calculos_para_rogelio.md' in archive.namelist()
    print(f'{OUT}: {OUT.stat().st_size} bytes, {len(chosen)} files')
    # Smaller review package retains reports, spectra and the inputs needed by
    # the three reproduction commands; regenerated bulk transients stay in full ZIP.
    required={'results/ltspice/reproduction/initial/study.json', 'results/ltspice/reproduction/initial/noise_200.csv',
        'results/ltspice/reproduction/initial/noise_400.csv', 'results/ltspice/reproduction/initial/decks/uniform_N1_r3.cir',
        'results/ltspice/reproduction/initial/decks/uniform_r2.pwl', 'results/ltspice/reproduction/initial/uniform_N1_r1.csv.gz',
        'results/ltspice/reproduction/initial/uniform_N1_r2.csv.gz'}
    compact=OUT.with_name('calculos_CSA_Rogelio_revision.zip')
    with zipfile.ZipFile(compact,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for name in chosen:
            path=ROOT/name
            keep=(not name.startswith('results/') or name in required
                or name.startswith('results/ltspice/reproduction/noise_audit/')
                or path.suffix in ('.md','.json','.pdf')
                or path.name in ('shaping_scan.csv','comparison.csv'))
            if keep and path.is_file():archive.write(path,arcname='calculos_CSA/'+name)
    with zipfile.ZipFile(compact) as archive:assert archive.testzip() is None
    print(f'{compact}: {compact.stat().st_size} bytes')

if __name__=='__main__':main()
