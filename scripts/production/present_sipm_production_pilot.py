#!/usr/bin/env python3
"""Present the saved statistical SiPM plus ASIC pilot and production-time estimate."""
from pathlib import Path
from contextlib import ExitStack
import html,os,subprocess,tomllib
import h5py
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sipm_statistical_production_pilot'
LATEX=ROOT/'latex_generated/electron_readout/sipm_statistical_production_pilot'

def scalar(g,k):
    """Read a scalar value or UTF-8 string from the numerical result."""
    v=g[k][()]
    return v.decode() if isinstance(v,bytes) else v

def main():
    """Render Julia timing estimates and validation metrics without recomputing response."""
    LATEX.mkdir(parents=True,exist_ok=True)
    timing=[];checks=[];totals=[];links=[]
    with ExitStack() as stack:
        f=stack.enter_context(h5py.File(OUT/'timing.h5'))
        assert f['metadata/complete'][()]
        cases=[(tag,f[tag]) for tag in ['ep2p5','ep3']]
        if (OUT/'ep2/timing.h5').exists():
            ep2=stack.enter_context(h5py.File(OUT/'ep2/timing.h5'))
            assert ep2['metadata/complete'][()]
            cases.insert(0,('ep2',ep2['ep2']))
        for tag,g in cases:
            e=tomllib.loads(scalar(g,'electron_card'));sensor=tomllib.loads(scalar(g,'sensor_card'))
            electronics=tomllib.loads(scalar(g,'electronics_card'))
            ep=e['reduced_field_kV_cm_bar']
            totals.append([f'{ep:g}',f'{g["total_s"][()]:.2f}',f'{g["projected_8000_minutes"][()]:.1f}'])
            manifest=scalar(g,'manifest')
            links.append(f'<li>E/P = {ep:g}: <a href="{html.escape(os.path.relpath(manifest,OUT))}">complete pilot manifest</a></li>')
            for gas,label in [('pure_xe','Xe'),('xe_he15','Xe–15% He')]:
                r=g[gas];assert r['all_checks_passed'][()]
                timing.append([f'{ep:g}',label,str(int(r['events'][()])),f'{r["sensor_s"][()]:.2f}',f'{r["asic_s"][()]:.2f}',
                    f'{r["projected_8000_s"][()]/60:.1f}',f'{r["sensor_bytes"][()]/2**20:.1f}',f'{r["asic_bytes"][()]/2**20:.1f}'])
                checks.append([f'{ep:g}',label,f'{r["noise_mean_pe"][()]:.5f}',f'{r["noise_rms_pe"][()]:.5f}',
                    f'{100*r["clipped_cell_fraction"][()]:.2f}',f'{100*r["total_charge_retention"][()]:.2f}',
                    str(int(r['state_repairs'][()])),str(int(r['variance_limits'][()]))])
        description=(f'100 bb0nu events at each operating point, propagated through both gas models. '
            f'The statistical SiPM model uses {sensor["statistical_substeps"]} intervals per {electronics["bucket_width_s"]*1e9:g} ns bucket. '
            f'SiPM PDE is {sensor["pde"]:g}, recovery time {sensor["recovery_s"]*1e9:g} ns, with {sensor["microcell_count"]:,} microcells. '
            f'The ASIC adds {electronics["noise_rms_pe"]:g} PE RMS noise and clips at {electronics["current_limit_A"]*1e3:g} mA. '
            'Julia warm-up was measured separately for each invocation. Full 8,000-event response production has not been started.')
        method=('Each SiPM/ASIC timing includes computation and compressed HDF5 writing. The estimate scales processing and manifest/summary finalization '
            'from 100 to 8,000 events; input resolution is counted once. Each measurement used one Julia process with sequential gas branches. '
            'Warm-up, independent saved-file validation, process startup and report rendering are outside the projection. '
            'The estimate uses the first 100 source events; event complexity and storage throughput can change the full-run rate.')
        validation=('All saved-history checks passed: event IDs and pair checksums, pixel and signed bucket correspondence, parent-file hashes, '
            'detected/dark/sensor/readout charge ledgers and the ASIC current limit. Quiet-cell noise is measured where the sensor input charge is zero. '
            'Clipped fraction uses every acquired cell. Charge retention is the summed ASIC output divided by summed SiPM input charge. '
            'Recovery-state repairs and variance limiting are recorded below.')
    def table(heads,rows):
        return '<table><tr>'+''.join(f'<th>{html.escape(h)}</th>' for h in heads)+'</tr>'+''.join('<tr>'+''.join(f'<td>{html.escape(v)}</td>' for v in row)+'</tr>' for row in rows)+'</table>'
    th=['E/P','Gas','Events','SiPM + write (s)','ASIC + write (s)','8,000 estimate (min)','SiPM MiB','ASIC MiB']
    ch=['E/P','Gas','Quiet mean (PE)','Quiet RMS (PE)','Clipped cells %','Charge retained %','State repairs','Variance limits']
    totalheads=['E/P','Pilot total: both gases (s)','8,000: both gases (min)']
    comparison_link='<p><a href="clipping_comparison.html">E/P comparison: clipping and event-response spread</a></p>' if (OUT/'clipping_comparison.html').exists() else ''
    headline='100-event pilots with fixed SiPM and ASIC settings; production estimates cover both gas models at each E/P.'
    page='<!doctype html><meta charset="utf-8"><title>SiPM and ASIC production pilot</title><style>body{font:16px system-ui;max-width:1250px;margin:35px auto;padding:0 20px;color:#17334a;line-height:1.6}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:right}table{border-collapse:collapse}</style><h1>Statistical SiPM + ASIC production pilot</h1><p><strong>'+headline+'</strong></p><p>'+description+'</p>'+table(totalheads,totals)+'<h2>Timing by gas</h2>'+table(th,timing)+'<p>'+method+'</p><h2>Saved-response checks</h2><p>'+validation+'</p>'+table(ch,checks)+comparison_link+'<h2>Outputs</h2><ul>'+''.join(links)+'</ul><p><a href="files.toml">File index</a> · <a href="../../latex_generated/electron_readout/sipm_statistical_production_pilot/sipm_production_pilot.pdf">PDF report</a></p>'
    (OUT/'sipm_production_pilot.html').write_text(page)
    def esc(v):return v.replace('–','--').replace('%',r'\%')
    def tex_table(heads,rows):
        return '\\begin{center}\\resizebox{\\textwidth}{!}{\\begin{tabular}{'+'l'*len(heads)+'}\\toprule\n'+' & '.join(esc(h) for h in heads)+r' \\ \midrule'+'\n'+'\n'.join(' & '.join(esc(v) for v in row)+r' \\' for row in rows)+'\n\\bottomrule\\end{tabular}}\\end{center}\n'
    tex=r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm]{geometry}\usepackage{booktabs,graphicx,hyperref,xurl}
\begin{document}\title{Statistical SiPM and ASIC production pilot}\author{NAUSICAA}\date{}\maketitle
'''+esc(headline)+'\n\n'+esc(description)+'\n'+tex_table(totalheads,totals)+r'\section*{Timing by gas}'+'\n'+tex_table(th,timing)+esc(method)+r'\section*{Saved-response checks}'+'\n'+esc(validation)+tex_table(ch,checks)+r'''
The numerical source is \path{results/sipm_statistical_production_pilot/timing.h5}.
The independent sensor and ASIC files are linked through the complete manifests
in \path{results/sipm_statistical_production_pilot/files.toml}.
Settings, source hashes and event/channel seeds accompany every response shard.
The photon histories are reused from their completed production.
\end{document}
'''
    name='sipm_production_pilot.tex';(LATEX/name).write_text(tex)
    (LATEX/'.gitignore').write_text('*.aux\n*.log\n*.out\n*.fls\n*.fdb_latexmk\n')
    with (LATEX/'build.log').open('w') as log:
        subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error',name],cwd=LATEX,stdout=log,stderr=subprocess.STDOUT,check=True)
    print(OUT/'sipm_production_pilot.html')
if __name__=='__main__':main()
