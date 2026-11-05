#!/usr/bin/env python3
"""Render measured photon-production timing and its file index to HTML and LaTeX/PDF."""
from pathlib import Path
import html, os, subprocess, tomllib
import h5py
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/photon_production_pilot'
LATEX=ROOT/'latex_generated/electron_readout/photon_production_pilot'

def text(group,key):
    """Decode a saved HDF5 scalar string."""
    value=group[key][()]
    return value.decode() if isinstance(value,bytes) else value

def main():
    """Present the saved Julia timings; estimates are read from the numerical summary."""
    LATEX.mkdir(parents=True,exist_ok=True)
    rows=[];links=[];texrows=[];file_rows=[];all_full=True
    with h5py.File(OUT/'timing.h5') as f:
        assert f['metadata/complete'][()]
        warm=f['metadata/warmup_s'][()]
        n=f['metadata/pilot_events'][()]
        for tag in ['ep2p5','ep3']:
            g=f[tag];e=tomllib.loads(text(g,'electron_card_toml'))
            row=(e['reduced_field_kV_cm_bar'],e['j_variance_ratio'],int(g['events'][()]),
                 float(g['elapsed_s'][()]),float(g['estimated_full_s'][()]),
                 float(g['pure_xe/full_recorded_elapsed_s'][()]) if g['pure_xe/matching_full_exists'][()] else None)
            values=[f'{row[0]:g}',f'{row[1]:g}',str(row[2]),f'{row[3]:.2f}',f'{row[4]:.1f}',f'{row[5]:.1f}' if row[5] is not None else 'N/A']
            rows.append('<tr>'+''.join(f'<td>{v}</td>' for v in values)+'</tr>')
            texrows.append(' & '.join(values)+r' \\')
            for gas,label in [('pure_xe','Xe'),('xe_he15','Xe–15% He')]:
                group=g[gas]
                all_full=all_full and bool(group['matching_full_exists'][()])
                choices=[('pilot',text(group,'pilot_path'))]
                if group['matching_full_exists'][()]:choices.append(('full 8,000',text(group,'full_path')))
                anchor=' · '.join(f'<a href="{html.escape(os.path.relpath(path,OUT))}">{title}</a>' for title,path in choices)
                links.append(f'<li>E/P = {row[0]:g}, {label}: {anchor}</li>')
                file_rows.append(f'{row[0]:g} & '+('Xe' if gas=='pure_xe' else r'Xe--15\% He')+f' & {int(group["pilot_bytes"][()]):,} & '+(f'{int(group["full_bytes"][()]):,}' if 'full_bytes' in group else 'N/A')+r' \\')
        e=tomllib.loads(text(f['ep3'],'electron_card_toml'))
        r=tomllib.loads(text(f['ep3'],'realisation_card_toml'))
    description=(f'Fresh pilot: {n:,} bb0nu events per operating point, each propagated through pure Xe and Xe–15% He. '
        'Every reported time covers BOTH gases together, processed sequentially. '
        f'Julia warm-up was measured separately ({warm:.2f} s); times exclude process startup and report generation. '
        'Full-run estimates scale the measured pilot time to 8,000 events. '
        'The pilot checks all event identities, primary-pair replay checksums and photon/electron ledgers.')
    if all_full:description+=' Existing matching 8,000-event files are complete and available.'
    conditions=(f'Pressure {e["pressure_bar"]:g} bar; chamber length {r["chamber_length_cm"]:g} cm; '
        f'EL pitch {e["hole_pitch_mm"]:g} mm; hole diameter/depth {e["hole_diameter_mm"]:g}/{e["hole_depth_mm"]:g} mm; '
        f'time bucket {e["bucket_width_us"]:g} microseconds; electron lifetime {e["lifetime_ms"]:g} ms. '
        'The saved counts are VUV photons incident on each SiPM per bucket. PDE, SiPM recovery and electronics are subsequent stages. '
        'The model assigns photons to the parent electron arrival bucket, uses the recorded gas diffusion and EL laws, '
        'and applies analytical direct optical acceptance with absorbing hole walls.')
    heads=['E/P (kV/cm/bar)','J','Events','Pilot (s)','Estimate: 8,000 (s)','Existing full run (s)']
    page='<!doctype html><meta charset="utf-8"><title>Photon production pilot</title><style>body{max-width:1100px;margin:40px auto;padding:0 20px;font:17px system-ui;color:#17334a;line-height:1.6}td,th{padding:10px;border-bottom:1px solid #ddd;text-align:right}</style><h1>Photon production pilot</h1><p>'+description+'</p><table><tr>'+''.join(f'<th>{v}</th>' for v in heads)+'</tr>'+''.join(rows)+'</table><p>'+conditions+'</p><h2>Photon files</h2><ul>'+''.join(links)+'</ul><p><a href="files.toml">File index (TOML)</a> · <a href="../../latex_generated/electron_readout/photon_production_pilot/photon_production_pilot.pdf">PDF report</a></p>'
    (OUT/'photon_production_pilot.html').write_text(page)
    def latex(s):return s.replace('–','--').replace('%',r'\%')
    tex=r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm]{geometry}
\usepackage{booktabs,hyperref,xurl}
\begin{document}\title{Photon production pilot}\author{NAUSICAA}\date{}\maketitle
'''+latex(description)+r'''
\begin{center}\begin{tabular}{rrrrrr}\toprule
$E/P$ & $J$ & Events & Pilot (s) & 8,000 estimate (s) & Existing full run (s)\\\midrule
'''+ '\n'.join(texrows)+r'''
\bottomrule\end{tabular}\end{center}
$E/P$ is in kV/cm/bar. All times include both gases.

'''+latex(conditions)+r'''
\section*{Saved file sizes}
\begin{center}\begin{tabular}{rlrr}\toprule
$E/P$ & Gas & Pilot (bytes) & Full 8,000 (bytes)\\\midrule
'''+ '\n'.join(file_rows)+r'''
\bottomrule\end{tabular}\end{center}
Numerical timing and provenance: \path{results/photon_production_pilot/timing.h5}.
The photon files are indexed in \path{results/photon_production_pilot/files.toml}.
The HTML report in that directory links directly to the pilot and full photon files.

The production implementation is \path{scripts/production/run_electron_250ns.jl};
the timing driver is \path{scripts/production/benchmark_photon_production.jl}.
The driver rejects cached pilot outputs so a timing measurement always performs
fresh diffusion, collection, EL sampling and HDF5 writing. A subsequent measurement
can select another event limit. Presentation is regenerated from the saved summary
with \path{scripts/production/present_photon_production_pilot.py}.
\end{document}
'''
    name='photon_production_pilot.tex';(LATEX/name).write_text(tex)
    (LATEX/'.gitignore').write_text('*.aux\n*.log\n*.out\n*.fls\n*.fdb_latexmk\n')
    with (LATEX/'build.log').open('w') as log:
        subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error',name],cwd=LATEX,stdout=log,stderr=subprocess.STDOUT,check=True)
    print(OUT/'photon_production_pilot.html')
if __name__=='__main__':main()
