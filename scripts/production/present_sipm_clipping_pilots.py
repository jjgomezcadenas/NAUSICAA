#!/usr/bin/env python3
"""Present the saved E/P=2,2.5,3 SiPM/ASIC clipping comparison."""
from pathlib import Path
import html,shutil,subprocess,tomllib,os
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sipm_statistical_production_pilot'
LATEX=ROOT/'latex_generated/electron_readout/sipm_statistical_production_pilot'

def main():
    """Render the Julia comparison as plots, tables, HTML and a self-contained PDF."""
    LATEX.mkdir(parents=True,exist_ok=True);rows=[];links=[]
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    with h5py.File(OUT/'clipping_comparison.h5') as f:
        assert f['metadata/complete'][()]
        for row,(gas,title) in enumerate([('pure_xe','Pure Xe'),('xe_he15','Xe + 15% He')]):
            for tag in ['ep2','ep2p5','ep3']:
                g=f[tag+'/'+gas];ep=float(g['reduced_field_kV_cm_bar'][()]);j=float(g['J'][()])
                limit_pe=float(g['bucket_limit_pe'][()])
                sensor=tomllib.loads(g['sensor_card'][()].decode());electronics=tomllib.loads(g['electronics_card'][()].decode())
                label=f'E/P {ep:g}, J {j:g}'
                a=g['event_retention'];axes[row,0].stairs(a['counts'][()],100*a['edges'][()],label=label)
                a=g['normalized_asic_response'];axes[row,1].stairs(a['counts'][()],a['edges'][()],label=label)
                rows.append([title,f'{ep:g}',f'{j:g}',str(int(g['event_count'][()])),
                    f'{g["charge_retained_percent"][()]:.1f}',f'{g["clipped_photon_positive_cells_percent"][()]:.1f}',
                    f'{g["events_with_clipping_percent"][()]:.0f}',f'{g["sensor_response_rms_percent"][()]:.1f}',
                    f'{g["asic_response_rms_percent"][()]:.1f}'])
                if row==0:
                    path=g['manifest'][()].decode()
                    links.append(f'<li>E/P = {ep:g}: <a href="{html.escape(os.path.relpath(path,OUT))}">100-event output manifest</a></li>')
            axes[row,0].set_xlabel('Event charge retained [%]')
            axes[row,1].set_xlabel('Corrected ASIC response / mean response')
            for ax in axes[row]:ax.set_title(title);ax.set_ylabel('Events');ax.legend(fontsize=8);ax.grid(alpha=.2)
    figure='clipping_comparison.png';fig.savefig(OUT/figure,dpi=150);plt.close(fig);shutil.copy2(OUT/figure,LATEX/figure)
    description=(f'The same 100 bb0nu source events, primary-pair checksums and drift distances are compared in both gases at E/P = 2, 2.5 and 3 kV/cm/bar. '
        'The user-specified EL variance ratios are J = 0.030, 0.015 and 0.025, respectively. '
        f'The statistical SiPM uses {sensor["statistical_substeps"]} subdivisions per 250 ns bucket, PDE {sensor["pde"]:g}, gain {sensor["gain"]:g} '
        f'and recovery time {sensor["recovery_s"]*1e9:g} ns. Electronics remain fixed at {electronics["noise_rms_pe"]:g} PE RMS noise '
        f'and a {electronics["current_limit_A"]*1e3:g} mA upper current limit ({limit_pe:.1f} PE per bucket). Full production remains stopped.')
    definitions=('Charge retained is the sum of ASIC output divided by the sum of SiPM input charge over the sample. '
        'The clipped-cell percentage uses only SiPM–bucket cells containing incident photons. '
        'The event percentage counts events with at least one clipped bucket. These three quantities have distinct denominators.')
    spread=('Event-total response spread is RMS/mean after subtracting the expected stationary dark baseline, correcting known attachment '
        'and dividing by deposited energy. SiPM recovery nonlinearity is retained in both columns. '
        'The before/after comparison measures the forward response spread under these corrections. '
        'All 100 events enter; the quoted widths are full-distribution RMS values. Each plotted ASIC distribution is divided by its own mean.')
    heads=['Gas','E/P','J','Events','Charge retained %','Photon cells clipped %','Events with clipping %','Before ASIC RMS %','After ASIC RMS %']
    table='<table><tr>'+''.join(f'<th>{h}</th>' for h in heads)+'</tr>'+''.join('<tr>'+''.join(f'<td>{html.escape(v)}</td>' for v in row)+'</tr>' for row in rows)+'</table>'
    page='<!doctype html><meta charset="utf-8"><title>E/P scan: SiPM and ASIC clipping</title><style>body{max-width:1350px;margin:35px auto;padding:0 20px;font:16px system-ui;line-height:1.6;color:#17334a}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:right}img{width:100%}</style><h1>E/P = 2, 2.5, 3: clipping pilot</h1><p>'+description+'</p>'+table+'<p>'+definitions+'</p><p>'+spread+'</p><img src="'+figure+'"><ul>'+''.join(links)+'</ul><p><a href="../../latex_generated/electron_readout/sipm_statistical_production_pilot/clipping_comparison.pdf">PDF</a> · <a href="sipm_production_pilot.html">Timing and saved-history checks</a></p>'
    (OUT/'clipping_comparison.html').write_text(page)
    esc=lambda s:s.replace('%',r'\%').replace('–','--')
    tex=r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=16mm]{geometry}\usepackage{graphicx,booktabs,hyperref,xurl}
\begin{document}\title{E/P = 2, 2.5, 3: SiPM and ASIC clipping pilot}\author{NAUSICAA}\date{}\maketitle
'''+esc(description)+'\n\\begin{center}\\resizebox{\\textwidth}{!}{\\begin{tabular}{llrrrrrrr}\\toprule\n'+' & '.join(esc(h) for h in heads)+r' \\ \midrule'+'\n'+'\n'.join(' & '.join(esc(v) for v in row)+r' \\' for row in rows)+r'''
\bottomrule\end{tabular}}\end{center}
'''+esc(definitions)+'\n\n'+esc(spread)+r'''
\clearpage\noindent\includegraphics[width=\textwidth]{clipping_comparison.png}

Numerical results: \path{results/sipm_statistical_production_pilot/clipping_comparison.h5}.
The HTML counterpart links the complete response manifests for each operating point.
\end{document}
'''
    (LATEX/'clipping_comparison.tex').write_text(tex)
    with (LATEX/'clipping_build.log').open('w') as log:
        subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error','clipping_comparison.tex'],cwd=LATEX,stdout=log,stderr=subprocess.STDOUT,check=True)
    print(OUT/'clipping_comparison.html')
if __name__=='__main__':main()
