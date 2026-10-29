#!/usr/bin/env python3
"""Audit the original noise setup and a synthetic flicker positive control.

The control KF is arbitrary and is NOT a technology estimate or revised design.
Preserves the received schematic/models; writes separate derived decks.
"""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[4]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse
import hashlib
import json
import re
import numpy as np
from ion_signal.spice_acquisition import run
from ion_signal.spice import schematic_netlist

ROOT=Path(__file__).resolve().parents[4]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--simulator',required=True,type=Path)
    args=ap.parse_args()
    out=ROOT/'results/ltspice/noise_audit'
    decks=out/'decks';decks.mkdir(parents=True,exist_ok=True)
    model=ROOT/'LTSpice/180nm_professional.inc'
    asc=ROOT/'LTSpice/CSAfolded180nm.asc'
    elements,directives=schematic_netlist(asc)
    directives=[d for d in directives if not d.lower().startswith(('.tran','.meas'))]
    report={'simulator':str(args.simulator.resolve()),'version':'LTspice 17.2.4 native macOS',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (asc,model,Path(__file__).resolve(),ROOT/'LTSpice/Noise Analysis CSA 180 nm.pdf')},
        'conditions':'Original 250 fF external load, original sources and transistor dimensions, default 27 C; PDF noise sweep 1 Hz–1 GHz, 1000 points/decade.',
        'control':'KF=1e-25 added to both models only in positive-control copy; arbitrary diagnostic, not a physical noise calibration.',
        'cases':{}}
    for name,kf in [('original',None),('positive_control_kf',1e-25)]:
        lib=decks/f'{name}.inc'
        text='\n'.join(l.rstrip() for l in model.read_text().splitlines()).strip()+'\n'
        if kf is not None:
            assert text.count(')')==2
            text=text.replace(')',f' KF={kf:.17e} )')
        lib.write_text(text)
        deck=decks/f'{name}.cir'
        dirs=[f'.lib {lib.name}' if d.lower().startswith('.lib') else d for d in directives]
        deck.write_text('* Original schematic noise audit\n'+'\n'.join(elements+dirs)+'\n.options numdgt=15 plotwinsize=0\n.noise V(out) Is1 dec 1000 1 1G\n.end\n')
        d,h=run(args.simulator.resolve(),deck)
        names=list(d);a=np.column_stack([d[n].real for n in names])
        np.savetxt(out/f'{name}.csv.gz',a,delimiter=',',header=','.join(names),comments='')
        flicker={n:float(np.max(d[n].real)) for n in names if '1overf' in n}
        assert len(flicker)==7
        if kf is None:assert all(v==0 for v in flicker.values())
        else:assert all(v>0 for v in flicker.values())
        # Independent check: device total variances sum to output variance.
        totals=[n for n in names if re.fullmatch(r'V\(m\w+\)',n)]
        variance=sum(d[n].real**2 for n in totals)
        closure=float(np.max(np.abs(variance/d['V(onoise)'].real**2-1)))
        assert closure<1e-6
        points=[]
        for f in (1.,10.,100.,1000.,1e6,1e9):
            j=int(np.argmin(abs(a[:,0]-f)))
            points.append({'frequency_Hz':float(a[j,0]),'output_ASD_V_sqrtHz':float(d['V(onoise)'][j].real),
                'input_ASD_A_sqrtHz':float(d['inoise'][j].real),
                'flicker_output_ASD_V_sqrtHz':float(np.sqrt(sum(d[n][j].real**2 for n in flicker)))})
        report['cases'][name]={'flicker_max_ASD_V_sqrtHz':flicker,'variance_closure_relative':closure,'spectral_points':points}
    (out/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['cases'],indent=2))

if __name__=='__main__':main()
