#!/usr/bin/env python3
"""Present saved Julia pixel/frame charge reductions as 2×2 figures and HTML/LaTeX/PDF.

Read HDF5 inclusive/mean/maximum/minimum histograms without repeating transport,
collection or statistical reductions. Run from any directory after
analyze_mc_selection.jl --pixel-charge. Numeric values come from saved results.
"""
from pathlib import Path
import html,shutil,subprocess,tomllib,argparse
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--results-dir',type=Path,default=ROOT/'results/pixel_charge')
parser.add_argument('--latex-dir',type=Path,default=ROOT/'latex_generated/pixel_charge')
args=parser.parse_args()
OUT=args.results_dir.resolve();FIG=OUT/'figures';TEX=args.latex_dir.resolve()
FIG.mkdir(parents=True,exist_ok=True);TEX.mkdir(parents=True,exist_ok=True)
samples=('bb0nu','electrons');labels=('Double beta','Single electrons')
files=[h5py.File(OUT/(s+'.h5')) for s in samples]
try:
    limits={}
    for name in ('inclusive','mean','maximum','minimum'):
        occupied=[]
        for f in files:
            g=f['histograms/'+name];ii=np.flatnonzero(g['counts'][:]);edges=g['edges_ions'][:]
            if len(ii):occupied.extend((edges[max(0,ii[0]-1)],edges[min(len(edges)-1,ii[-1]+2)]))
        limits[name]=(min(occupied),max(occupied)) if occupied else (float(files[0]['histograms/'+name+'/edges_ions'][0]),float(files[0]['histograms/'+name+'/edges_ions'][-1]))
    card=tomllib.loads(files[0]['provenance/card_toml'].asstr()[()])
    step=card['run'].get('simulation_step_mm',0.1)
    slice_width=card['pads']['z_bin_width_mm']
    rows=[];timings=[];statistics=[];retention=[];losses=[]
    for sample,label,f in zip(samples,labels,files):
        summary=f['summary'];n_min=int(f['settings/n_min'][()]);frame=float(f['settings/frame_ms'][()]);pitch=float(f['settings/pitch_mm'][()])
        cut_label='no charge cuts' if n_min==0 else f'minimum histogram: min > {n_min}'
        n=int(summary['processed_events'][()]);nwith=int(summary['events_with_charge'][()])
        rows.append([label,str(n),str(nwith),f"{int(summary['pulse_count'][()]):,}",
                     f"{int(summary['collected_ions'][()]):,}",f"{int(summary['outside_ions'][()]):,}"])
        timings.append([label]+[f"{float(f['timing/'+k][()]) if 'timing/'+k in f else 0:.3f}" for k in ('read_seconds','transport_seconds','cache_read_seconds','cache_write_seconds','collection_reduction_seconds','elapsed_seconds')])
        for key,title in (('all_pixels','All pixels'),):
            r=f['ion_retention/'+key]
            retention.append([label,title]+[f'{100*float(r[k][()]):.4f}' for k in ('aggregate_retention','mean_event_retention','median_event_retention','aggregate_loss')])
        loss=f['ion_losses']
        losses.append([label,f"{int(loss['outside_histogram_ions'][()]):,}",f"{100*float(loss['outside_histogram_fraction'][()]):.4f}"])
        for name,title in (('maximum','Maximum'),('minimum','Minimum'),('mean','Event mean')):
            st=f['distribution_statistics/'+name]
            statistics.append([label,title,str(int(st['entries'][()])),f"{float(st['mean_ions'][()]):.2f}" if st['entries'][()] else '—',f"{float(st['median_ions'][()]):.2f}" if st['entries'][()] else '—'])
        fig,axes=plt.subplots(2,2,figsize=(11,8.3),layout='constrained')
        for ax,name,title in zip(axes.flat,('inclusive','mean','maximum','minimum'),
                ('Inclusive: ions per pulse','Mean ions per pulse per event','Maximum ions per pulse per event','Minimum ions per pulse per event' if n_min==0 else f'Event minimum ions (minimum > {n_min})')):
            g=f['histograms/'+name]
            ax.stairs(g['counts'][:],g['edges_ions'][:],color='#215872',lw=1.4)
            ax.set(title=title,xlabel=f'Ions per pixel per {frame:g} ms frame',
                   ylabel='Pixel–frame cells' if name=='inclusive' else 'Events',xlim=limits[name])
            if not np.any(g['counts'][:]):
                ax.set_ylim(.8,1.2)
                ax.text(.5,.5,'No events with minimum > '+str(n_min)+' ions',ha='center',transform=ax.transAxes,fontsize=10)
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.xaxis.set_minor_formatter(NullFormatter())
            if name=='minimum':
                ticks=np.linspace(*limits[name],3)
                ax.set_xticks(ticks,labels=[f'{v:g}' for v in ticks])
            ax.grid(alpha=.18)
            ax.text(.97,.95,f"{int(g['counts'][:].sum()):,} entries",ha='right',va='top',transform=ax.transAxes,fontsize=9)
        fig.suptitle(f'{label} · {step:g} mm MC steps · {n:,} events · {cut_label} · {pitch:g} mm pixels · {frame:g} ms frames',fontsize=12)
        png=FIG/(sample+'_pixel_charge.png');fig.savefig(png,dpi=170);plt.close(fig)
        shutil.copy2(png,TEX/png.name)
    requested=tomllib.loads(files[0]['provenance/card_toml'].asstr()[()])['run']['n_events']
    actual=[int(f['summary/processed_events'][()]) for f in files]
    assumptions=(f'Requested {requested:,} catalogue events per species from the input samples; processed {actual[0]:,} double-beta and {actual[1]:,} electron events, limited by available input events. All events and all MC particles reach collection. MC generation step: {step:g} mm (simulation input). '
        f'Random drift uniform in 10–150 cm; thermal isotropic Gaussian diffusion, 200 V/cm bulk field and 293.15 K; W = 21.9 eV. '
        f'The diffusion histogram has 100 bins per axis and 3σ padding. Out-of-histogram ions are recorded. '
        f'Ideal hole/pixel assignment maps each diffusion-bin centre to its containing {pitch:g} mm pixel. '
        f'Frame width is {frame:g} ms ({slice_width:g} mm slices at 0.1 m/s). Time is relative to the first arrival edge, increasing source z is later. '
        f'Per-hit ion counts are rounded integers. Each event has an independent reproducible drift, diffusion and transverse grid phase. '
        f'A pulse is a pixel–frame cell with positive ion count. Inclusive, event mean and event maximum include every pulse, with no charge threshold. '
        +(f'Every event with pulses enters the minimum histogram, using its true minimum over all positive cells. All four distributions have no charge cut. ' if n_min==0 else f'Compute each event’s true minimum over all its pulses, then include that event in the minimum histogram only if its minimum is greater than {n_min} ions. No individual pulse is discarded when computing the minimum. ')
        +f'The event mean is its total pulse-ion sum divided by its number of pulses. Events without pulses have no distribution entries. '
        f'The statistics table uses the same event populations as each histogram. An empty minimum distribution has unavailable mean/median. All histogram count and ion axes are logarithmic; bins are shared between species. '
        f'Track-reconstruction decisions do not gate this study. No time-integrated pixel charge is used in these plots. '
        f'Ion retention is all collected pixel ions divided by produced event ions, before histogram truncation. The aggregate ratio uses sums over events; event mean and median give every event equal weight. Loss is one minus retention, caused here only by diffusion-histogram truncation. The minimum histogram selection introduces no collection or threshold loss. Produced ions use the saved per-hit rounded counts. '
        f'Diffused ion clouds are saved independently of the collection grid when a cache is enabled. Pixel pitch and frame duration can change without repeating Gaussian draws; changing transport or diffusion-histogram settings selects a different cache. '
        f'Run times include first-event compilation; initial Julia process startup and final histogram/presentation generation are separate.')
    def table(heads,data):
        return '<table><tr>'+''.join('<th>'+html.escape(h)+'</th>' for h in heads)+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in row)+'</tr>' for row in data)+'</table>'
    heads=['Species','Events','Events with charge','Pixel–frame cells with charge','Collected ions','Outside diffusion histogram']
    retentionheads=['Species','Pixels included','Aggregate retention [%]','Mean event retention [%]','Median event retention [%]','Aggregate loss [%]']
    lossheads=['Species','Outside histogram [ions]','Histogram loss [%]']
    statheads=['Species','Distribution','Events','Mean [ions/pulse]','Median [ions/pulse]']
    timeheads=['Species','Source read [s]','Diffusion [s]','Cache read [s]','Cache write [s]','Collection/reduction [s]','Total [s]']
    body='<h1>Pixel charge per electronics frame</h1><p>'+html.escape(assumptions)+'</p>'+table(heads,rows)+'<h2>Ion retention and losses</h2>'+table(retentionheads,retention)+table(lossheads,losses)+'<h2>Distribution means and medians</h2>'+table(statheads,statistics)+'<h2>Run time</h2>'+table(timeheads,timings)
    for sample,label in zip(samples,labels):body+='<h2>'+label+'</h2><img src="figures/'+sample+'_pixel_charge.png" alt="'+label+' pixel-frame charge">'
    body+='<h2>Numerical inputs</h2><ul>'+''.join('<li>'+html.escape(f['provenance/source_path'].asstr()[()])+'<br>SHA-256: '+f['provenance/source_sha256'].asstr()[()]+'</li>' for f in files)+'</ul>'
    (OUT/'pixel_charge.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Pixel charge</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:0 20px;color:#243348}p{line-height:1.5}img{width:100%}table{width:100%;border-collapse:collapse;font-size:13px}th,td{border:1px solid #dce3ec;padding:8px;text-align:left}th{background:#edf2f6}li{overflow-wrap:anywhere}</style>'+body+'</html>')
    def esc(s):
        for a,b in [('_',r'\_'),('%',r'\%'),('&',r'\&'),('μ',r'$\mu$'),('σ',r'$\sigma$'),('±',r'$\pm$'),('–','--'),('—','--')]:s=s.replace(a,b)
        return s
    tex=r'\documentclass{article}\usepackage[margin=18mm]{geometry}\usepackage{graphicx,booktabs}\begin{document}\title{Pixel charge per electronics frame}\author{}\date{\today}\maketitle'+'\n'+esc(assumptions)+'\n'
    for title,h,data in [('Charge accounting',heads,rows),('Ion retention',retentionheads,retention),('Ion loss accounting',lossheads,losses),('Distribution means and medians',statheads,statistics),('Run time',timeheads,timings)]:
        tex+=r'\section*{'+title+r'}\noindent\resizebox{\linewidth}{!}{\begin{tabular}{'+'l'*len(h)+r'}\toprule '+' & '.join(map(esc,h))+r'\\\midrule'+'\n'
        tex+='\n'.join(' & '.join(map(esc,row))+r'\\' for row in data)+r'\bottomrule\end{tabular}}'+'\n'
    for sample,label in zip(samples,labels):
        tex+=r'\clearpage\section*{'+label+r'}\noindent\includegraphics[width=\linewidth]{'+sample+r'_pixel_charge.png}'+'\n'
    tex+=r'\end{document}';(TEX/'pixel_charge.tex').write_text(tex)
    subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error','pixel_charge.tex'],cwd=TEX,check=True)
    print('HTML:',OUT/'pixel_charge.html');print('PDF:',TEX/'pixel_charge.pdf')
    for row in rows+timings:print(' | '.join(row))
finally:
    for f in files:f.close()
