#!/usr/bin/env python3
"""Publish the completed matched-pitch ion CSA response and measured timing."""
from pathlib import Path
import argparse,html,shutil,sys,tomllib
import h5py
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/documentation'))
from build_report import build_report

def main(pitch):
    root=ROOT/f'results/ion_csa_nonlinear/matched_pitch{pitch}'
    latest=tomllib.loads((root/'response/latest.toml').read_text())
    manifest=Path(latest['manifest']);run=tomllib.loads(manifest.read_text())
    timing=tomllib.loads((manifest.parent/'production_timing.toml').read_text())
    collection=tomllib.loads((root/'collection.toml').read_text())
    assert timing['complete'] and timing['fresh_filtered_outputs'] and run['complete']
    node=run['responses'][0]
    with h5py.File(node['path']) as h:
        m=h['metadata']
        stats={k:int(m[k][()]) for k in ('event_count','input_ion_total','pads','direct_spice_pads','nonlinear_pads')}
        assert m['complete'][()] and m['mode'][()].decode()=='production'
    assert stats['event_count']==collection['events']==8000
    assert stats['input_ion_total']==collection['total_ions']
    execution=Path(timing.get('execution_directory',manifest.parent/f'execution/pitch{pitch}'))
    source=next((execution/'corrected_report').glob('*/report.tex'))
    out=ROOT/f'latex_generated/ion_readout/ion_csa_nonlinear_pitch{pitch}'
    out.mkdir(parents=True,exist_ok=True)
    for name in ('diagnostics.pdf','diagnostics.png','sources.json'):
        shutil.copy2(source.parent/name,out/name)
    shutil.copy2(next((execution/'corrected').glob('*/ion_csa.html')),out/'response.html')
    rows=[('Ion diffusion and collection',collection['collection_wall_s']),
          ('Covariance calibration setup',timing['calibration_setup_s'])]
    labels={'prepare_s':'Linear predictions and current preparation','nonlinear_s':'Nonlinear response and storage',
            'response_noise_summary_s':'Noise, filtered response and summary','fallback_workflow_s':'SPICE fallback workflow',
            'correction_summary_s':'Three-frame reconstruction and summary','presentation_s':'Figures and presentation'}
    rows.extend((label,timing[key]) for key,label in labels.items() if key in timing)
    rows.extend([('Event work excluding calibration and final presentation',timing['event_work_s']),
                 ('Event work per event',timing['event_work_seconds_per_event'])])
    if 'complete_wall_s' in timing:
        rows.append(('Complete CSA wall time',timing['complete_wall_s']))
    if 'recovery_wall_s' in timing:
        rows.append(('Publication recovery and registration',timing['recovery_wall_s']))
    latexrows='\n'.join(f'{label} & {value:.4f} '+r'\\' for label,value in rows)
    summary=(f'This production contains {stats["event_count"]:,} bb0nu events at {pitch} mm pitch and {10*pitch} ms frames. '
        'The shared realization places ion drift uniformly between 5 and 145 cm, with complementary electron drift in the 150 cm chamber. '
        'Every event ID, primary-pair checksum and drift distance was checked against the existing matched common realization and 4 mm ion file. '
        'The electron readout was preserved. All events and signed response buckets are retained. '
        f'The ion ledger is {stats["input_ion_total"]:,}, exactly equal to the primary-pair total. '
        f'The CSA used {stats["nonlinear_pads"]:,} nonlinear histories and {stats["direct_spice_pads"]:,} SPICE fallback histories among '
        f'{stats["pads"]:,} pad histories. Production creates no SPICE control histories. '
        'Fresh filtered outputs were asserted before the run; actual circuit-deck, model, kernel-waveform and calibration checksums are saved in preflight.toml.')
    if timing.get('publication_recovered'):
        summary += (' Julia sendfile reported an error after producing a byte-identical publication copy. The saved copy was verified by SHA256, its production metadata completed and its event ledgers validated before registration. Numerical timing covers the completed fresh calculation; recovery and registration are timed separately.')
    block=(r'\section{Matched production and timing}'+'\n'+summary+'\n'+
        r'\begin{center}\small\begin{tabular}{p{11cm}r}\toprule Stage & Seconds\\\midrule'+'\n'+latexrows+'\n'+
        r'\bottomrule\end{tabular}\end{center}'+'\n'+
        'The numerical timer includes screening, nonlinear response, noise, reconstruction, output writing and consistency checks. '
        'Covariance calibration and final presentation are separated. Event work includes fallback orchestration and its internal report generation when SPICE fallbacks occur. Stage and complete times are nested; the per-event value divides measured event work by 8,000.\n')
    (out/'production_timing.tex').write_text(block)
    tex=source.read_text().replace(r'\begin{document}',r'\author{Generated production results}'+'\n'+r'\begin{document}')
    tex=tex.replace(r'\maketitle',r'\maketitle'+'\n'+r'\input{production_timing.tex}')
    (out/'report.tex').write_text(tex)
    table='<table><tr><th>Stage</th><th>Seconds</th></tr>'+''.join(
        f'<tr><td>{html.escape(label)}</td><td>{value:.4f}</td></tr>' for label,value in rows)+'</table>'
    (out/'index.html').write_text(f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>ITACA matched {pitch} mm ion production</title>
<style>body{{font:17px system-ui;line-height:1.5;max-width:1100px;margin:40px auto;padding:0 20px;color:#193c4a}}td,th{{padding:9px;border:1px solid #ccd8df}}table{{border-collapse:collapse}}img{{width:100%}}</style>
<h1>ITACA matched {pitch} mm ion production</h1><p>{summary}</p>{table}
<p><a href="report.pdf">PDF report</a> · <a href="response.html">Interactive response diagnostics</a></p>
<img src="diagnostics.png" alt="Ion CSA response diagnostics"></html>''')
    build_report(out/'report.tex');print(out/'index.html')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('pitch',type=int,choices=(2,3))
    main(parser.parse_args().pitch)
