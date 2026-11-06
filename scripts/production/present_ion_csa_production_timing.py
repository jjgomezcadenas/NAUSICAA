#!/usr/bin/env python3
"""Present the saved production-mode timing and response report without simulation."""
from pathlib import Path
import html
import shutil
import sys
import tomllib
import h5py

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/documentation'))
from build_report import build_report

def main():
    latest=tomllib.loads((ROOT/'results/ion_csa_nonlinear/production_mode/latest.toml').read_text())
    manifest=Path(latest['manifest']);run=tomllib.loads(manifest.read_text())
    total=tomllib.loads((manifest.parent/'timing.toml').read_text())
    base=manifest.parent/'execution/pitch4'
    times=tomllib.loads((base/'timing.toml').read_text())
    assert times['mode']=='production'
    source=next((base/'corrected_report').glob('*/report.tex'))
    out=ROOT/'latex_generated/ion_readout/ion_csa_nonlinear/production_100'
    out.mkdir(parents=True,exist_ok=True)
    node=run['responses'][0]
    with h5py.File(node['path']) as h:
        m=h['metadata'];events=int(m['event_count'][()]);ions=int(m['input_ion_total'][()])
        nonlinear=int(m['nonlinear_pads'][()]);spice=int(m['direct_spice_pads'][()])
    rows=[('Input preparation and linear predictions',times['prepare_s']),
          ('Nonlinear evaluation and storage',times['nonlinear_s']),
          ('Response, noise and summary',times['response_noise_summary_s']),
          ('Three-frame correction and summary',times['correction_summary_s']),
          ('Numerical branch total',times['numerical_s']),
          ('Presentation',times['presentation_s']),
          ('Complete workflow, including module loading and publication',total['entrypoint_elapsed_s'])]
    for name in ('diagnostics.pdf','diagnostics.png','sources.json'):
        shutil.copy2(source.parent/name,out/name)
    response=next((base/'corrected').glob('*/ion_csa.html'))
    shutil.copy2(response,out/'response.html')
    latexrows='\n'.join(f'{label} & {seconds:.2f} '+r'\\' for label,seconds in rows)
    caveat=(('This execution reused a filtered-response file. Its elapsed time is not a fresh-output scaling measurement. '
             if total.get('filtered_response_reused',False) else '')+
            'The separate fresh-output scaling report measures event processing after warm-up and separates covariance calibration and presentation. ')
    timing=(r'\section{Production mode timing}'+'\n'+caveat+'\n\n'+
        f'This run processes {events} bb0nu events at 4 mm pitch and 40 ms frames, '
        f'with {nonlinear} nonlinear pad histories, {spice} SPICE histories and '
        f'{ions:,} input ions. The run card selects '+r'\texttt{mode = "production"}.'+'\n'+
        'The saved calibration is reused. SPICE controls and solver-refinement histories '
        'are generated only in validation mode. Charge, noise, reconstruction, provenance '
        'and storage checks remain active.\n'+
        r'\begin{center}\small\begin{tabular}{p{11cm}r}\toprule Stage & Seconds\\\midrule'+'\n'+
        latexrows+'\n'+r'\bottomrule\end{tabular}\end{center}'+'\n'+
        'The numerical total includes the preceding numerical stages; the complete total '
        'also includes module loading, presentation and publication. These are nested timings. '
        'The entrypoint timer starts before module loading. The run uses one Julia process '
        'for its numerical stages.\n')
    (out/'timing.tex').write_text(timing)
    tex=source.read_text().replace(r'\begin{document}',r'\author{Generated production results}'+'\n'+r'\begin{document}')
    tex=tex.replace(r'\maketitle',r'\maketitle'+'\n'+r'\input{timing.tex}')
    (out/'report.tex').write_text(tex)
    table='<table><tr><th>Stage</th><th>Seconds</th></tr>'+''.join(
        f'<tr><td>{html.escape(label)}</td><td>{seconds:.2f}</td></tr>' for label,seconds in rows)+'</table>'
    page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>ITACA CSA production: 100 events</title>
<style>body{{font:17px system-ui;line-height:1.5;max-width:950px;margin:40px auto;padding:0 20px;color:#193c4a}}td,th{{border:1px solid #cbd8df;padding:9px}}table{{border-collapse:collapse}}img{{width:100%}}</style>
<h1>ITACA CSA production: {events} events</h1>
<p><strong>Timing correction:</strong> {caveat}<a href="../scaling/index.html">Fresh-output scaling report</a>.</p>
<p><strong>Complete workflow: {total['entrypoint_elapsed_s']:.2f} s.</strong> bb0nu, 4 mm pitch, 40 ms frames; mode: production.</p>
<p>{nonlinear} nonlinear pad histories, {spice} SPICE histories, zero repeated SPICE controls. Charge ledger: {ions:,} ions. Numerical consistency checks passed.</p>
{table}<p>Numerical and complete totals include the preceding stages. The complete timer starts before module loading and includes reconstruction, presentation and publication. The stages run in one Julia process.</p>
<p><code>mode = "production"</code> reuses the validated calibration; SPICE serves domain fallback.
<code>mode = "validation"</code> additionally generates SPICE controls and solver refinements.</p>
<p>Card: <code>config/specs/runs/ion_csa_nonlinear_production_100.toml</code>.</p>
<p><a href="report.pdf">PDF report</a> · <a href="response.html">Interactive response diagnostics</a> · <a href="../index.html">Earlier validation report</a></p>
<img src="diagnostics.png" alt="Production response diagnostics"></html>'''
    (out/'index.html').write_text(page)
    build_report(out/'report.tex')
    print(out/'index.html')

if __name__=='__main__':main()
