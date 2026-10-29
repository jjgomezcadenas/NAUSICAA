#!/usr/bin/env python3
"""Present saved Julia electron diffusion coefficients and widths in HTML/LaTeX/PDF."""
from pathlib import Path
import subprocess
import h5py
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/electron_diffusion';TEX=ROOT/'latex_generated/electron_diffusion';TEX.mkdir(parents=True,exist_ok=True)
with h5py.File(OUT/'electron_diffusion.h5') as f:
    pressure=float(f['pressure_bar'][()]);helium=100*float(f['helium_fraction'][()])
    coefficients=[];widths=[]
    for gas,label in [('pure_xe','Pure Xe'),('xe_he',f'Xe + {helium:.0f}% He')]:
        coefficients.append([label]+[f"{f[gas+'/d_'+a+'_mm_per_sqrtcm'][()]:.6f}" for a in ('T','L')])
        for i,length in enumerate(f['drift_cm'][:]):
            widths.append([label,f'{length:g}']+[f"{f[gas+'/sigma_'+a+'_mm'][i]:.4f}" for a in ('T','L')])
    error=max(float(f[g+'/max_Float64_error_'+a+'_mm'][()]) for g in ('pure_xe','xe_he') for a in ('T','L'))
intro=(f'Operating pressure: {pressure:g} bar. Xe-He helium fraction: {helium:.0f}% (user design inputs). '
       'Petit prescription: sigma_T,L [mm] = d_T,L [mm sqrt(bar/cm)] sqrt(L_cm/P_bar). '
       'The operating-pressure coefficients below multiply sqrt(L_cm). '
       'Pure-Xe reference coefficients are 3.5/0.9; Xe-He coefficients are 1.6/0.75. '
       f'Widths use Float32; maximum absolute difference from Float64 evaluation is {error:.3g} mm. '
       'The electron formula uses the supplied coefficients and pressure; it has no explicit electric-field argument.')
def htmltable(heads,rows):
    return '<table><tr>'+''.join('<th>'+v+'</th>' for v in heads)+'</tr>'+''.join('<tr>'+''.join('<td>'+v+'</td>' for v in r)+'</tr>' for r in rows)+'</table>'
body='<h1>Electron diffusion at operating pressure</h1><p>'+intro+'</p>'+htmltable(['Gas','dT at pressure [mm/√cm]','dL at pressure [mm/√cm]'],coefficients)+htmltable(['Gas','Drift [cm]','σT [mm]','σL [mm]'],widths)
(OUT/'electron_diffusion.html').write_text('<!doctype html><meta charset="utf-8"><title>Electron diffusion</title><style>body{font:16px system-ui;max-width:950px;margin:40px auto;line-height:1.5}table{border-collapse:collapse;margin:25px 0}td,th{border:1px solid #ccc;padding:10px}</style>'+body)
tex=r'\documentclass{article}\usepackage[margin=22mm]{geometry}\usepackage{booktabs}\setlength{\emergencystretch}{2em}\begin{document}\title{Electron diffusion at operating pressure}\author{}\date{\today}\maketitle'+'\n'+intro.replace('%',r'\%').replace('_',r'\_')+'\n'
for heads,rows in [(['Gas',r'$d_T(15\,\mathrm{bar})$ [mm/$\sqrt{\mathrm{cm}}$]',r'$d_L(15\,\mathrm{bar})$ [mm/$\sqrt{\mathrm{cm}}$]'],coefficients),(['Gas','Drift [cm]',r'$\sigma_T$ [mm]',r'$\sigma_L$ [mm]'],widths)]:
    tex+=r'\begin{center}\begin{tabular}{'+'l'*len(heads)+r'}\toprule '+' & '.join(heads)+r'\\\midrule'+'\n'
    tex+='\n'.join(' & '.join(v.replace('%',r'\%') for v in r)+r'\\' for r in rows)+r'\bottomrule\end{tabular}\end{center}'+'\n'
tex+=r'\end{document}';(TEX/'electron_diffusion.tex').write_text(tex)
subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error','electron_diffusion.tex'],cwd=TEX,check=True)
print(coefficients);print(widths)
