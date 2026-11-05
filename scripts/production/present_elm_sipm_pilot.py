#!/usr/bin/env python3
"""Generate ELM SiPM pilot figures and HTML/LaTeX/PDF from Julia reductions."""
from pathlib import Path
import html
import shutil
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/documentation'))
from build_report import build_report
OUT=ROOT/'results/elm_sipm_pilot'
LATEX=ROOT/'latex_generated/electron_readout/elm_sipm_pilot'
GASES=(('pure_xe','Pure Xe'),('xe_he15','Xe–15% He'))
TAGS=(('12p5','12.5 ns'),('25','25 ns'))
def fmt(x,n=3):return f'{float(x):,.{n}f}'
def table(heads,rows):
    return '<table><tr>'+''.join(f'<th>{html.escape(str(x))}</th>' for x in heads)+'</tr>'+''.join('<tr>'+''.join(f'<td>{html.escape(str(x))}</td>' for x in row)+'</tr>' for row in rows)+'</table>'
def save(fig,name):
    fig.tight_layout();fig.savefig(OUT/name,dpi=160);plt.close(fig);shutil.copy2(OUT/name,LATEX/name)
def esc(s):
    return str(s).replace('%',r'\%').replace('_',r'\_').replace('&',r'\&').replace('×',r'$\times$').replace('²',r'$^2$').replace('−',r'$-$').replace('–','--').replace('π',r'$\pi$').replace('τ',r'$\tau$').replace('M[1$-$(1$-$1/M)^N]',r'$M[1-(1-1/M)^N]$')
def tex_table(heads,rows):
    return r'\begin{center}\resizebox{\textwidth}{!}{\begin{tabular}{'+'l'+'r'*(len(heads)-1)+r'}\toprule'+'\n'+' & '.join(map(esc,heads))+r'\\\midrule'+'\n'+'\n'.join(' & '.join(map(esc,row))+r'\\' for row in rows)+r'\bottomrule\end{tabular}}\end{center}'+'\n'
def main():
    LATEX.mkdir(parents=True,exist_ok=True)
    with h5py.File(OUT/'elm_sipm_pilot.h5','r') as f:
        assert bool(f['metadata/complete'][()])
        meta={k:f['metadata'][k][()] for k in f['metadata']}
        data={tag:{gas:{k:f[f'tau_{tag}/{gas}'][k][()] for k in f[f'tau_{tag}/{gas}']} for gas,_ in GASES} for tag,_ in TAGS}
        pars={tag:{k:f[f'tau_{tag}/parameters'][k][()] for k in f[f'tau_{tag}/parameters']} for tag,_ in TAGS}
    p=pars['12p5'];rows=[];widths=[];pulses=[];worst=[];checks=[]
    for tag,tlabel in TAGS:
      for gas,glabel in GASES:
        d=data[tag][gas]
        rows.append([glabel,tlabel,fmt(d['detected_mean_total_pe'],0),fmt(d['sipm_charge_lost_percent'],4),fmt(d['asic_charge_lost_percent'],4),fmt(d['clipped_event_percent'],1),fmt(d['clipped_pulse_percent'],5)])
        widths.append([glabel,tlabel,*[fmt(d[k],4) for k in ('corrected_detected_rms_percent','corrected_sensor_rms_percent','corrected_unclipped_rms_percent','corrected_readout_rms_percent')]])
        pulses.append([glabel,tlabel,fmt(d['event_mean_average_pe'],2),fmt(d['event_max_average_pe'],1),fmt(d['largest_detected_pulse'],0),fmt(d['largest_sensor_pulse'],1),fmt(d['largest_expected_occupancy_percent'],2)])
        checks.append([glabel,tlabel,int(d['ledger_checks']),int(d['state_repairs']),int(d['variance_limits']),fmt(d['pde_mean_z'],2),fmt(d['quiet_noise_rms_pe'],4),fmt(1e6*d['maximum_state_repair_fraction'],3)])
        for i in np.argsort(d['sipm_loss_percent']+d['clip_loss_percent'])[-5:][::-1]:
            worst.append([glabel,tlabel,int(d['event_id'][i]),fmt(d['electron_drift_cm'][i],2),fmt(d['event_max_pe'][i],0),fmt(d['event_max_occupancy_percent'][i],2),fmt(d['sipm_loss_percent'][i],3),fmt(d['clip_loss_percent'][i],3)])
    heads=['Gas','Recovery','Mean detected PE/event','SiPM loss (%)','ASIC loss (%)','Events clipped (%)','Pulses clipped (%)']
    wh=['Gas','Recovery','Detected RMS (%)','SiPM RMS (%)','Noisy unclipped RMS (%)','ASIC RMS (%)']
    ph=['Gas','Recovery','Event-mean avg PE','Event-max avg PE','Largest count','Largest charge PE','Max occupancy (%)']
    worsth=['Gas','Recovery','Event ID','Drift (cm)','Max count','Occupancy (%)','SiPM loss (%)','ASIC loss (%)']
    ch=['Gas','Recovery','Ledger checks','State repairs','Variance limits','PDE z','Quiet RMS PE','Max repair (ppm)']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(3,2,figsize=(12,11))
    for col,(gas,label) in enumerate(GASES):
      d=data['12p5'][gas]
      bins=np.unique(np.r_[[-.5,.5],np.geomspace(1.5,max(d['inclusive_detected_pe'])+.5,80)])
      for key,lab in [('inclusive_detected_pe','Optical triggers'),('inclusive_optical_charge_pe','SiPM charge')]:
        axes[0,col].hist(d[key],bins=bins,histtype='step',label=lab)
      axes[0,col].set_xscale('symlog',linthresh=1);axes[0,col].set_yscale('log');axes[0,col].legend()
      axes[0,col].set(title=label+' · 12.5 ns recovery',xlabel='Counts or charge-equivalent PE / SiPM–250 ns',ylabel='Input cells')
      for row,key in [(1,'event_mean_pe'),(2,'event_max_pe')]:
        for tag,lab in TAGS:axes[row,col].hist(data[tag][gas][key],bins=25,histtype='step',label=lab)
        axes[row,col].set(xlabel=('Event mean' if row==1 else 'Event maximum')+' optical triggers per bucket',ylabel='Events');axes[row,col].legend()
        if row==2:axes[row,col].axvline(p['clip_threshold_pe'],color='crimson',ls='--')
    save(fig,'elm_sipm_pulses.png')
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for row,(tag,tlabel) in enumerate(TAGS):
      for col,(gas,label) in enumerate(GASES):
        d=data[tag][gas]
        keys=[('detected_corrected','Detected'),('sensor_corrected','After SiPM'),('readout_corrected','After ASIC')]
        edges=np.histogram_bin_edges(np.concatenate([d[k] for k,_ in keys]),bins=30)
        for key,lab in keys:axes[row,col].hist(d[key],bins=edges,histtype='step',label=lab)
        axes[row,col].set(title=label+' · '+tlabel,xlabel='Attachment-corrected PE / deposited MeV',ylabel='Events');axes[row,col].legend()
    save(fig,'elm_sipm_energy.png')
    fig,axes=plt.subplots(3,2,figsize=(12,11))
    for col,(gas,label) in enumerate(GASES):
      for tag,lab in TAGS:
        d=data[tag][gas]
        for row,key in [(0,'sipm_loss_percent'),(1,'clip_loss_percent'),(2,'event_max_occupancy_percent')]:
            axes[row,col].scatter(d['electron_drift_cm'],d[key],s=23,label=lab,alpha=.75)
      for row,y in [(0,'Event optical charge lost in SiPM (%)'),(1,'Event charge lost in ASIC (%)'),(2,'Maximum expected distinct-cell occupancy (%)')]:
        axes[row,col].set(xlabel='Electron reference drift (cm)',ylabel=y,title=label);axes[row,col].legend()
    save(fig,'elm_sipm_drift.png')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for ax,(gas,label) in zip(axes,GASES):
      for tag,lab in TAGS:
        ax.hist(data[tag][gas]['inclusive_occupancy_percent'],bins=np.geomspace(.001,100,90),histtype='step',label=lab)
      ax.set(xscale='log',yscale='log',xlabel='Expected distinct optical cells / illuminated cells (%)',ylabel='Input cells',title=label);ax.legend()
    save(fig,'elm_sipm_occupancy.png')
    description=('The same first 100 bb0nu tracks (source IDs 0–99) are replayed for both gases. Electron reference drift is sampled uniformly from 10 to 150 cm with a separate reproducible placement stream. Ion drift is its complement within 150 cm. This reference distance sets the common diffusion and attachment for the event hits; their source coordinates and primary-pair realization are retained. No ion histories are produced. Both recovery cases reuse the identical saved blue-photon histories; their SiPM stochastic draws are independent.')
    optics=('E/P = 2 kV/cm/bar at 15 bar; 3 mm diameter holes, 5 mm EL depth, 5 mm pitch and 250 ns buckets. Electron collection is unity. Uniform axial EL and absorbing walls give the analytical exit acceptance. TPB conversion 0.65, forward factor 0.5 and guide capture 1−1/1.5 apply once. The optical probability is '+fmt(data['12p5']['pure_xe']['optical_acceptance'],8)+'. Guide propagation and NOA 87 coupling are ideal; the xenon gap is 0.5 mm radially. The saved optical signal is blue photons at the SiPM entrance.')
    sensor=(f'S14160-3015PS: 3 × 3 mm², {int(p["physical_microcells"]):,} cells at 15 micrometre pitch, PDE 0.32 (peak at 460 nm used as a spectrum-average approximation), gain 360,000, +4 V and 25 C. Uniform illumination of the centered 3 mm circular guide footprint uses {int(p["illuminated_microcells"]):,} cells, round(39984 × π/4). The remaining 8,581 cells receive dark triggers only. The full-device 700 kcps dark rate is split by cell count; both populations have independent stationary recovery states. Provisional exponential recovery constants 12.5 and 25 ns are sensitivity parameters, not measured recharge constants for this device. Crosstalk, afterpulsing and gain spread remain deferred; PDE is independent of recovery.')
    model=('The validated statistical model retains first and second moments of cell recovery, uses exact conditional finite-cell occupancy moments, and a bounded moment-matched charge draw. Dense buckets use 128 substeps; sparse triggers use individual uniform times. Optical trigger counts are sampled PDE successes, and recurrent firing is allowed. Charge-equivalent PE is collected charge divided by e × gain. Distinct optical occupancy is the conditional expected number M[1−(1−1/M)^N] of cells reached by N optical triggers during a bucket. It counts each cell once; it is saved as an expectation rather than as sampled cell identities. See src/ion_signal/sipm_statistical.jl and the implementation tests.')
    electronics=(f'The ASIC adds {float(p["noise_rms_pe"]):.4f} PE RMS Gaussian noise to each acquired bucket and clips at 1 mA bucket-average current ({float(p["clip_threshold_pe"]):.2f} charge-equivalent PE). Signed values are retained. Pulse comparisons use all photon-positive input cells at every stage. Event totals sum the full acquisition interval on illuminated channels, including four quiet buckets at either end. SiPM loss is optical trigger count minus optical charge, divided by trigger count. ASIC loss is clipped charge divided by incoming sensor charge, including dark charge. Event widths are full-distribution RMS/mean after known attachment and deposited-energy correction. Optical and SiPM widths use optical quantities; ASIC widths subtract the expected full-area dark baseline before normalization. A correction for SiPM compression is not applied.')
    timing=(f'Photon replay/archive {float(meta["photon_wall_s"]):.2f} s; response 12.5 ns {float(pars["12p5"]["response_wall_s"]):.2f} s; response 25 ns {float(pars["25"]["response_wall_s"]):.2f} s; combined {float(meta["total_wall_s"]):.2f} s. Stage timings include compilation and I/O; process startup and presentation are separate.')
    validation=('Twenty dedicated tests cover placement reproducibility, uniform range, full-area dark rate, count/charge conservation, recurrent firing and 128-to-256 substep convergence. Saved-event checks verify primary checksums, PDE statistics, charge ledgers, response diagnostics and quiet-cell noise. Rare state feasibility repairs adjust only the second recovery moment, preserving sampled charge and the first moment. Their maximum size is reported relative to the smaller cell population; checks require less than 10 ppm and fewer than one repair per million batches. No charge-variance limits occurred.')
    figures=['elm_sipm_pulses.png','elm_sipm_energy.png','elm_sipm_drift.png','elm_sipm_occupancy.png']
    sections=[('Configuration',description+'\n\n'+optics+'\n\n'+sensor),('Statistical response and definitions',model+'\n\n'+electronics),('Timing and validation',timing+'\n\n'+validation)]
    body='<h1>ELM · finite-cell SiPM pilot</h1>'+''.join(f'<h2>{title}</h2>'+''.join('<p>'+html.escape(t)+'</p>' for t in text.split('\n\n')) for title,text in sections[:1])
    body+='<h2>Separate SiPM and ASIC losses</h2>'+table(heads,rows)+'<h2>Corrected event widths</h2>'+table(wh,widths)+'<h2>Pulse statistics</h2>'+table(ph,pulses)
    body+=''.join(f'<h2>{title}</h2>'+''.join('<p>'+html.escape(t)+'</p>' for t in text.split('\n\n')) for title,text in sections[1:])
    body+=''.join(f'<img src="{name}" alt="{name}">' for name in figures)+'<h2>Largest event losses</h2>'+table(worsth,worst)+'<h2>Validation checks</h2>'+table(ch,checks)
    cards=data['12p5']['pure_xe']['detector_card'].decode()+'\n'+pars['12p5']['sensor_card'].decode()+'\n'+pars['25']['sensor_card'].decode()+'\n'+p['electronics_card'].decode()
    body+='<p><a href="../../latex_generated/electron_readout/elm_sipm_pilot/elm_sipm_pilot.pdf">PDF</a> · <a href="files.toml">File index</a>. Complete workflow: <code>make elm-sipm-pilot</code>; presentation: <code>make elm-sipm-report</code>.</p><details><summary>Cards and provenance</summary><pre>'+html.escape(cards)+'</pre></details>'
    (OUT/'elm_sipm_pilot.html').write_text('<!doctype html><meta charset="utf-8"><title>ELM SiPM pilot</title><style>body{font:16px system-ui;line-height:1.6;max-width:1400px;margin:35px auto;padding:0 22px;color:#17334a}table{border-collapse:collapse;margin:18px 0;width:100%}td,th{padding:8px;text-align:right;border-bottom:1px solid #ccd6de}td:first-child,th:first-child{text-align:left}img{width:100%;margin:15px 0}pre{white-space:pre-wrap}</style>'+body)
    tex=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage[margin=22mm]{geometry}
\usepackage{graphicx,booktabs,amsmath,xurl,hyperref,fvextra}
\begin{document}\title{ELM: finite-cell SiPM pilot}\author{NAUSICAA}\date{}\maketitle
'''
    for title,text in sections[:1]:tex+=r'\section*{'+esc(title)+'}\n'+esc(text)+'\n'
    tex+=r'\section*{Charge losses and event widths}'+tex_table(heads,rows)+tex_table(wh,widths)+tex_table(ph,pulses)
    for title,text in sections[1:]:tex+=r'\section*{'+esc(title)+'}\n'+esc(text)+'\n'
    for name in figures:tex+=r'\clearpage\begin{figure}[p]\centering\includegraphics[width=\textwidth,height=.88\textheight,keepaspectratio]{'+name+r'}\end{figure}'+'\n'
    tex+=r'\clearpage\section*{Largest event losses}'+tex_table(worsth,worst)+r'\section*{Validation}'+tex_table(ch,checks)
    tex+=r'\section*{Reproduction}\path{results/elm_sipm_pilot/elm_sipm_pilot.h5} contains the Julia reductions; \path{results/elm_sipm_pilot/files.toml} resolves all histories and response cards. Run \texttt{make elm-sipm-pilot}; regenerate presentation with \texttt{make elm-sipm-report}.'+'\n'
    tex+=r'\section*{Input cards}\begin{Verbatim}[breaklines,breakanywhere,fontsize=\footnotesize]'+'\n'+cards+'\n'+r'\end{Verbatim}'+'\n'+r'\end{document}'+'\n'
    (LATEX/'elm_sipm_pilot.tex').write_text(tex);build_report(LATEX/'elm_sipm_pilot.tex')
    print(OUT/'elm_sipm_pilot.html');print(LATEX/'elm_sipm_pilot.pdf')
if __name__=='__main__':main()
