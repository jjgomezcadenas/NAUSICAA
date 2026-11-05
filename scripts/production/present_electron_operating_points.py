#!/usr/bin/env python3
"""Present saved pre-electronics pulse and event-energy distributions."""
from pathlib import Path
import shutil,subprocess,html,tomllib
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
import sys
sys.path.insert(0, str(ROOT / "scripts/documentation"))
from build_report import build_report


def present(source=None, report_dir=None):
    """Generate HTML, figures and self-contained LaTeX from Julia reductions."""
    source=Path(source) if source else ROOT/'results/electron_operating_points/electron_operating_points.h5'
    out=source.parent;texdir=Path(report_dir) if report_dir else ROOT/'latex_generated/electron_readout/electron_operating_points';texdir.mkdir(parents=True,exist_ok=True)
    plots=[];rows=[]
    with h5py.File(source) as f:
        if not f['complete'][()]:raise ValueError('Incomplete reduction')
        available=list(f.keys()); available=[x for x in available if isinstance(f[x],h5py.Group)]
        reference=f[available[0]+'/pure_xe']
        sensor=tomllib.loads(reference['sensor_card_toml'][()].decode())
        electron=tomllib.loads(reference['electron_card_toml'][()].decode())
        electronics=tomllib.loads(reference['electronics_card_toml'][()].decode())
        for gas,title in [('pure_xe','Pure Xe'),('xe_he15','Xe + 15% He')]:
            fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
            metrics=[('inclusive_sensor_pe','SiPM output entering electronics [PE / 250 ns bucket]'),
                     ('total_sensor_pe','Total SiPM output per event [PE]'),
                     ('attachment_corrected_pe','Dark-subtracted, attachment-corrected event sum [PE]'),
                     ('corrected_pe_per_MeV','Corrected response / deposited energy [PE / MeV]')]
            for tag in (t for t in ('legacy','ep3','ep2p5') if t in f):
                g=f[f'{tag}/{gas}'];ep=float(g['reduced_field_kV_cm_bar'][()]);j=float(g['j_variance_ratio'][()])
                label=f'E/P {ep:g}, J {j:g}'+(' (previous)' if tag=='legacy' else '')
                for ax,(key,xlabel) in zip(axes.flat,metrics):
                    d=g[key];ax.stairs(d['counts'][()],d['edges'][()],label=label)
                    ax.set_xlabel(xlabel);ax.set_ylabel('Pulses' if key=='inclusive_sensor_pe' else 'Events');ax.grid(alpha=.2)
                    if key=='inclusive_sensor_pe':ax.set_yscale('log')
                c=g['attachment_corrected_pe'];r=g['corrected_pe_per_MeV']
                rows.append((title,ep,j,len(g['event_id']),float(g['inclusive_sensor_pe/mean'][()]),
                    float(c['relative_rms_percent'][()]),float(c['gaussian_equivalent_fwhm_percent'][()]),
                    float(r['gaussian_equivalent_fwhm_percent'][()]),100*float(g['clipped_pulse_fraction'][()]),
                    100*float(g['readout_charge_retention'][()])))
            for ax in axes.flat:ax.legend(fontsize=8)
            fig.suptitle(title+f' · {len(g["event_id"]):,} bb0nu events · 10 mm pitch · 250 ns buckets')
            name=f'{gas}_pre_electronics.png';fig.savefig(out/name,dpi=160);plt.close(fig)
            shutil.copy2(out/name,texdir/name);plots.append(name)
    parameters=(f"Pressure {electron['pressure_bar']:g} bar; holes {electron['hole_diameter_mm']:g} mm in diameter and {electron['hole_depth_mm']:g} mm deep; "
        f"PDE {sensor['pde']:g}; gain {sensor['gain']:g}; {sensor['microcell_count']} microcells; recovery time {sensor['recovery_s']*1e9:g} ns. ")
    completion='Completed productions shown: '+', '.join(available)+'. '
    missing=[x for x in ('ep2p5','ep3') if x not in available]
    if missing:completion+='Awaiting completion: '+', '.join(missing)+'. '
    explanation=completion+parameters+('The signal injected into electronics is the SiPM avalanche charge divided by the fully recovered single-PE charge. '
        f"It includes microcell recovery and dark avalanches, before the {electronics['noise_rms_pe']:g} PE electronic noise and {electronics['current_limit_A']*1e3:g} mA clipping. "
        'The inclusive pulse distribution uses photon-positive SiPM–bucket cells. Event sums use every acquired cell. '
        f"Expected stationary dark charge is subtracted before multiplying by exp(drift time / {electron['lifetime_ms']:g} ms). "
        'The dark baseline uses the stationary dark-only recovery factor; illumination can further reduce dark avalanche charge. '
        'The corrected sum retains SiPM nonlinearity. The tabulated width is the RMS divided by the mean of the complete event distribution. '
        'The response per deposited MeV additionally removes variations in the simulated deposited energy. '
        'Gaussian-equivalent FWHM is 2.355 times RMS; the distribution may be non-Gaussian. '
        'Both gases use the user-specified Oliveira Figure 6 readings: J = 0.015 at E/P = 2.5 and J = 0.025 at 3. '
        'The previous E/P = 3, J = 0.01 production is retained as a reference.')
    headers=['Gas','E/P','J','Events','Mean pulse PE','Corrected RMS %','Equivalent FWHM %','PE/MeV equiv. FWHM %','Clipped pulses %','Output charge retained %']
    table=''.join('<tr>'+''.join('<td>'+html.escape(str(v) if i in (0,3) else f'{v:.3f}')+'</td>' for i,v in enumerate(row))+'</tr>' for row in rows)
    page='<!doctype html><meta charset="utf-8"><title>Electron operating points</title><style>body{font:16px system-ui;color:#17334a;max-width:1350px;margin:35px auto}p{line-height:1.6}table{border-collapse:collapse}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:right}img{width:100%}</style><h1>Electron operating points · pre-electronics response</h1><p>'+explanation+'</p><p><a href="../../latex_generated/electron_readout/electron_operating_points/electron_operating_points.pdf">PDF</a> · <a href="http://127.0.0.1:8769/">Interactive viewer</a> · <a href="../response_overview/response_overview.html">Ion and electron overview</a></p><table><tr>'+''.join(f'<th>{h}</th>' for h in headers)+'</tr>'+table+'</table>'+''.join(f'<img src="{name}">' for name in plots)
    (out/'electron_operating_points.html').write_text(page)
    tex=r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=15mm]{geometry}
\usepackage{graphicx,booktabs,hyperref}
\begin{document}
\title{Electron operating points: signal entering electronics}\author{NAUSICAA}\date{}\maketitle
'''+explanation.replace('%',r'\%')+'\n\\begin{center}\\scriptsize\\resizebox{\\textwidth}{!}{\\begin{tabular}{lrrrrrrrrr}\\toprule\nGas & E/P & J & Events & Mean PE & RMS (\\%) & FWHM (\\%) & PE/MeV FWHM (\\%) & Clipped (\\%) & Retained (\\%) \\\\ \\midrule\n'
    for row in rows:tex+=' & '.join(str(v).replace('%',r'\%') if i in (0,3) else f'{v:.3f}' for i,v in enumerate(row))+r' \\'+'\n'
    tex+='\\bottomrule\\end{tabular}}\\end{center}\n'
    for name in plots:tex+='\\clearpage\\noindent\\includegraphics[width=\\textwidth]{'+name+'}\n'
    tex+='\\end{document}\n';(texdir/'electron_operating_points.tex').write_text(tex)
    build_report(texdir / "electron_operating_points.tex")
    print(out/'electron_operating_points.html')
if __name__=='__main__':present()
