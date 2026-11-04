#!/usr/bin/env python3
"""Render saved Julia pulse reductions as HTML, PNG and a LaTeX/PDF report."""
from pathlib import Path
import argparse
import html
import shutil
import subprocess
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
LABELS={'pure_xe':'Pure Xe','xe_he15':'Xe + 15% He'}

def render(source):
    """Present saved histograms and summary scalars without repeating reduction."""
    out=source.parent;latex=ROOT/'latex_generated/response_overview';latex.mkdir(parents=True,exist_ok=True)
    plots=[];table=[]
    with h5py.File(source) as f:
        if not f['complete'][()]: raise ValueError('Incomplete reduction')
        events=int(f['validated_paired_events'][()])
        for arm in ('ions','photons','detected','sensor','readout'):
            fig,axes=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
            for ax,metric in zip(axes,('mean','maximum')):
                series=[(f'ions_{pitch}mm_{stage}',f'{pitch} mm / {pitch*10} ms · {stage}',stage=='input')
                        for pitch in (2,3,4) for stage in ('input','corrected')] if arm=='ions' else [
                        (f'{gas}_{arm}',label,False) for gas,label in LABELS.items()]
                for index,(key,label,dashed) in enumerate(series):
                    g=f[key];unit=g['unit'][()].decode()
                    ax.stairs(g[metric+'_histogram'][()],g[metric+'_edges'][()],label=label,linestyle='--' if dashed else '-',color=f'C{index//2 if arm=="ions" else index}')
                ax.set(xlabel=f'Event {metric} [{unit} / channel–bucket]',ylabel='Events',title=f'Event {metric}')
                ax.grid(alpha=.2);ax.legend(fontsize=8)
            fig.suptitle({'ions':'Ions: input and corrected CSA','photons':'Incident photons per SiPM–250 ns bucket',
                          'detected':'Detected PE per SiPM–250 ns bucket','sensor':'SiPM response per 250 ns bucket',
                          'readout':'Electronics PE per SiPM–250 ns bucket'}[arm])
            name=f'{arm}_mean_max.png';fig.savefig(out/name,dpi=160);plt.close(fig)
            shutil.copy2(out/name,latex/name);plots.append(name)
        for key in f:
            if not isinstance(f[key],h5py.Group):continue
            g=f[key];table.append((key.replace('_',' '),len(g['event_id']),g['unit'][()].decode(),
                *(float(g[f'{metric}_ensemble_{stat}'][()]) for metric in ('mean','maximum') for stat in ('mean','median'))))
    explanation=(f'{events:,} paired bb0nu events. Each event contributes one mean and one maximum. '
        'Ions use pixels with positive input ions; electrons use SiPM buckets with positive incident photons. '
        'The same cells enter all response stages, including signed reconstructed values. No amplitude or S/N cut is applied. '
        'Ion pitches 2, 3 and 4 mm use 20, 30 and 40 ms frames. Both electron gases use 10 mm pitch and 250 ns buckets. '
        'CSA results use the corrected reconstruction and the production 0.5 ENC modelling tolerance. '
        'These pulse summaries describe the occupied input support; event energy requires the complete response history.')
    heads=['Stage','Events','Unit','Mean: average','Mean: median','Maximum: average','Maximum: median']
    rows=''.join('<tr>'+''.join(f'<td>{html.escape(str(x) if i<3 else f"{x:.3f}")}</td>' for i,x in enumerate(row))+'</tr>' for row in table)
    document='<!doctype html><meta charset="utf-8"><title>NAUSICAA response overview</title><style>body{font:17px system-ui;max-width:1250px;margin:40px auto;color:#17334a}img{width:100%}td,th{padding:9px;border-bottom:1px solid #ddd;text-align:right}th:first-child,td:first-child{text-align:left}p{line-height:1.6}</style><h1>Ion and electron response · bb0nu</h1><p><a href="../electron_operating_points/electron_operating_points.html">Electron operating points: pre-electronics pulses and event sums</a></p><p>'+html.escape(explanation)+'</p><p><a href="../../latex_generated/response_overview/response_overview.pdf">PDF report</a></p><table><tr>'+''.join(f'<th>{h}</th>' for h in heads)+'</tr>'+rows+'</table>'+''.join(f'<img src="{p}" alt="{p}">' for p in plots)
    (out/'response_overview.html').write_text(document)
    tex=r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=18mm]{geometry}
\usepackage{graphicx,booktabs,hyperref}
\begin{document}
\title{Ion and electron response: bb0nu}\author{NAUSICAA}\date{}\maketitle
'''+explanation.replace('%',r'\%').replace('S/N','signal-to-noise')+'\n\\begin{center}\\scriptsize\n\\begin{tabular}{lrlrrrr}\\toprule\nStage & Events & Unit & Mean avg. & Mean median & Max avg. & Max median \\\\ \\midrule\n'
    for row in table:tex+=' & '.join(str(v) if i<3 else f'{v:.3f}' for i,v in enumerate(row))+r' \\'+'\n'
    tex+='\\bottomrule\\end{tabular}\\end{center}\n'
    for plot in plots:tex+=f'\\begin{{figure}}[htbp]\\centering\\includegraphics[width=\\textwidth]{{{plot}}}\\end{{figure}}\n'
    tex+='\\clearpage\\section*{Reproduction}\nNumerical reduction: \\texttt{scripts/production/reduce\\_response\\_overview.jl}.\\par Presentation: \\texttt{scripts/production/present\\_response\\_overview.py}. The HDF5 reduction retains the source card, SiPM manifest, event IDs and reducer digest.\n\\end{document}\n'
    (latex/'response_overview.tex').write_text(tex)
    subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error','response_overview.tex'],cwd=latex,check=True,stdout=subprocess.DEVNULL)
    print(out/'response_overview.html');print(latex/'response_overview.pdf')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source',nargs='?',type=Path,default=ROOT/'results/response_overview/response_overview.h5');args=parser.parse_args();render(args.source)
