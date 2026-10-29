#!/usr/bin/env python3
"""Generate MC-selection PNGs, HTML and self-contained LaTeX/PDF from Julia HDF5.

Run python scripts/events/present_mc_selection.py from any directory. Plot counts,
efficiencies, uncertainties, validation values and timings are read from the
saved results; physical calculations belong to analyze_mc_selection.jl.
"""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import html
import subprocess
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/mc_selection';FIG=OUT/'figures';TEX=ROOT/'latex_generated/mc_selection'
FIG.mkdir(parents=True,exist_ok=True);TEX.mkdir(parents=True,exist_ok=True)
samples=('bb0nu','electrons');labels=('Double beta (100 μm)','Single electrons (100 μm)')
plots=[('voxel_energy_MeV','Occupied voxel energy','keV','Voxels'),
       ('n_voxels','Occupied voxels per event','Voxels','Events'),
       ('track_length_mm','Longest component backbone','mm','Events'),
       ('voxel_energy_sum_MeV','Total voxel energy per event','MeV','Events'),
       ('n_components','Connected components per event','Components','Events'),
       ('retained_voxel_energy_MeV','Voxel energy after threshold','keV','Voxels'),
       ('n_retained_voxels','Voxels after threshold per event','Voxels','Events'),
       ('retained_energy_sum_MeV','Graph-retained energy per event','MeV','Events'),
       ('component_length_mm','Individual component backbones','mm','Components')]
files=[h5py.File(OUT/(s+'.h5')) for s in samples]
try:
    efficiencies=[];timings=[];validations=[];figures=[]
    for label,f in zip(labels,files):
        s=f['selection'];t=f['timing'];v=f['validation']
        efficiencies.append([label,str(int(s['processed'][()])),str(int(s['accepted'][()])),
            f"{100*s['efficiency'][()]:.1f} ± {100*s['binomial_standard_error'][()]:.1f}",
            f"[{100*s['wilson_68_low'][()]:.1f}, {100*s['wilson_68_high'][()]:.1f}]"])
        timings.append([label]+[f'{t[k][()]:.3f}' for k in ('serial_warmup_seconds','serial_read_seconds',
            'serial_analysis_seconds','serial_batch_seconds','parallel_batch_seconds','parallel_startup_warmup_seconds','parallel_total_seconds')])
        validations.append([label,f"{1000*v['max_energy_closure_MeV'][()]:.6f}",
            f"{1000*v['max_voxel_error_vs_Float64_MeV'][()]:.6f}",str(bool(t['serial_parallel_identical'][()]))])
    for name,title,unit,ylabel in plots:
        energy_scale=1000 if unit=='keV' else 1
        fig,axes=plt.subplots(1,2,figsize=(10,3.4),sharex=True,sharey=True)
        overflow=[]
        for ax,label,f in zip(axes,labels,files):
            g=f['histograms/'+name];edges=g['edges'][:]*energy_scale
            for state,colour in (('before','#567b98'),('after','#cf602b')):
                ax.stairs(g[state+'/counts'][:],edges,label='Before selection' if state=='before' else 'Single component',color=colour,lw=1.5)
            # Restrict the view to occupied bins; saved bin edges remain common across species.
            ax.set(title=label,xlabel=f'{title} [{unit}]' if unit in ('keV','MeV','mm') else title,ylabel=ylabel)
            ax.legend(fontsize=8);ax.grid(alpha=.18)
            overflow.append(f"{label}: before under/overflow {int(g['before/underflow'][()])}/{int(g['before/overflow'][()])}, after {int(g['after/underflow'][()])}/{int(g['after/overflow'][()])}")
        active=[]
        for f in files:
            g=f['histograms/'+name];c=g['before/counts'][:]
            active.extend(i for i,n in enumerate(c) if n)
        if active:
            edges=files[0]['histograms/'+name+'/edges'][:]*energy_scale
            lo=max(0,min(active)-1);hi=min(len(edges)-1,max(active)+2)
            for ax in axes:ax.set_xlim(edges[lo],edges[hi])
        if unit=='keV':
            for ax in axes:ax.set_xlim(0,100)
        fig.suptitle(title,fontsize=12);fig.tight_layout();fig.savefig(FIG/(name+'.png'),dpi=160);plt.close(fig)
        figures.append((name,title,'; '.join(overflow)))
    settings=files[0]['settings']
    setting_text='; '.join(k+' = '+str(settings[k].asstr()[()] if settings[k].dtype.kind in ('O','S','U') else settings[k][()]) for k in settings)
    sample_counts=[int(f['selection/processed'][()]) for f in files]
    requested=int(__import__('tomllib').loads(files[0]['provenance/card_toml'].asstr()[()])['run']['n_events'])
    assumptions=(f'Requested {requested:,} catalogue events per sample; processed {sample_counts[0]:,} double-beta and {sample_counts[1]:,} electron events (limited by available input events). Independent samples, no event matching. Source coordinates, all particles, no diffusion or fiducial cut. '
        'Voxel threshold is inclusive and precedes graph construction. Accept exactly one connected component, including a one-node component. '
        'Event track length is the longest component backbone. Voxel-energy and component-length histograms pool entries over events; all other histograms have one entry per event. '
        'After-selection histograms use the same quantities and bins as before selection. Efficiencies include binomial standard errors and z=1 Wilson intervals. '
        'Float32 is used for physical reductions; saved Float64 sums provide independent precision references. '
        'Timing excludes initial Julia startup from serial warm-up; warm-up uses one event. Parallel uses two separate Julia processes with owned HDF5 readers, warmed behind a start barrier. '
        'Parallel batch time includes dispatch, reading, calculation and shard writing. Parallel total additionally includes process startup and warm-up. '
        'Each sample is benchmarked serial first; OS caches may favour the subsequent parallel run. Numerical output is the serial reference.')
    def table(headers,rows):
        return '<table><thead><tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table>'
    effhead=['Sample','Processed','Accepted','Efficiency [%]','Wilson 68% [%]']
    timehead=['Sample','Serial warm-up [s]','Read [s]','Analysis [s]','Serial batch [s]','2-process batch [s]','Process startup + warm-up [s]','Parallel total [s]']
    valhead=['Sample','Max energy closure [keV]','Max voxel error vs Float64 [keV]','Serial/parallel identical']
    body='<h1>MC selection</h1><p>'+html.escape(assumptions)+'</p><p>'+html.escape(setting_text)+'</p>'
    for title,heads,rows in [('Selection efficiency',effhead,efficiencies),('Timing',timehead,timings),('Energy validation',valhead,validations)]:
        body+='<h2>'+title+'</h2>'+table(heads,rows)
    for name,title,caption in figures:
        body+='<section><h2>'+title+'</h2><img src="figures/'+name+'.png" alt="'+title+'"><p>'+caption+'</p></section>'
    (OUT/'mc_selection.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>MC selection</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:0 20px;color:#243348}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:9px;border:1px solid #dce3ec;text-align:left}th{background:#edf2f6}img{width:100%}p{line-height:1.5}section{margin-top:35px}</style>'+body+'</html>')
    def esc(s):
        for a,b in [('\\',r'\textbackslash{}'),('_',r'\_'),('%',r'\%'),('&',r'\&'),('μ',r'$\mu$'),('±',r'$\pm$')]:s=s.replace(a,b)
        return s
    tex=r'\documentclass{article}\usepackage[margin=20mm]{geometry}\usepackage{graphicx,booktabs}\begin{document}\title{MC selection}\author{}\date{\today}\maketitle'+'\n'+esc(assumptions)+'\n\n'+esc(setting_text).replace('; ',r'\\'+'\n')+'\n'
    for title,heads,rows in [('Selection efficiency',effhead,efficiencies),('Energy validation',valhead,validations)]:
        tex+=r'\section*{'+title+r'}\begin{center}\resizebox{\linewidth}{!}{\begin{tabular}{'+'l'*len(heads)+r'}\toprule '+' & '.join(map(esc,heads))+r'\\\midrule '+'\n'
        tex+='\n'.join(' & '.join(map(esc,row))+r'\\' for row in rows)+r'\bottomrule\end{tabular}}\end{center}'+'\n'
    tex+=r'\section*{Timing}'+'\n'
    for row in timings:
        tex+=esc(row[0])+r'\begin{itemize}'+'\n'+''.join(r'\item '+esc(h)+': '+v+'\n' for h,v in zip(timehead[1:],row[1:]))+r'\end{itemize}'+'\n'
    # Copy figure assets so the generated LaTeX folder is portable and self-contained.
    import shutil
    for name,title,caption in figures:
        shutil.copy2(FIG/(name+'.png'),TEX/(name+'.png'))
        tex+=r'\clearpage\section*{'+esc(title)+r'}\noindent\includegraphics[width=\linewidth]{'+name+r'.png}'+'\n\n'+esc(caption)+'\n'
    tex+=r'\end{document}'
    (TEX/'mc_selection.tex').write_text(tex)
    print('HTML:',OUT/'mc_selection.html');print('LaTeX:',TEX/'mc_selection.tex')
    for row in efficiencies+timings:print(' | '.join(row))
finally:
    for f in files:f.close()
