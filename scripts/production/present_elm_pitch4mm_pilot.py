#!/usr/bin/env python3
"""Present the 4 mm Xe–He ELM pilot and matched 5 mm reference reductions."""
from pathlib import Path
import html,json,sys,shutil
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/documentation'))
from build_report import build_report
from present_elm_sipm_pilot import table,tex_table,esc,fmt
OUT=ROOT/'results/elm_sipm_pitch4mm_pilot'
LATEX=ROOT/'latex_generated/electron_readout/elm_sipm_pitch4mm_pilot'
TAGS=(('12p5','12.5 ns'),('25','25 ns'))
def save(fig,name):
    fig.tight_layout();fig.savefig(OUT/name,dpi=160);plt.close(fig);shutil.copy2(OUT/name,LATEX/name)
def main():
    LATEX.mkdir(parents=True,exist_ok=True)
    with h5py.File(OUT/'elm_sipm_pitch4mm_pilot.h5','r') as f:
        assert f['metadata/complete'][()]
        meta={k:f['metadata'][k][()] for k in f['metadata']}
        data={tag:{k:f[f'tau_{tag}/xe_he15'][k][()] for k in f[f'tau_{tag}/xe_he15']} for tag,_ in TAGS}
        old={tag:{k:f[f'reference/tau_{tag}/xe_he15'][k][()] for k in f[f'reference/tau_{tag}/xe_he15']} for tag,_ in TAGS}
        pars={tag:{k:f[f'tau_{tag}/parameters'][k][()] for k in f[f'tau_{tag}/parameters']} for tag,_ in TAGS}
    p=pars['12p5'];lossrows=[];widthrows=[];pulserows=[];worst=[];checks=[]
    for tag,tlabel in TAGS:
      for pitch,d in [('5 mm / 3 mm hole',old[tag]),('4 mm / 2.5 mm hole',data[tag])]:
        lossrows.append([pitch,tlabel,fmt(d['sipm_charge_lost_percent'],4),fmt(d['asic_charge_lost_percent'],4),fmt(d['clipped_event_percent'],1),int(d['retained_events'])])
        widthrows.append([pitch,tlabel,fmt(d['corrected_detected_rms_percent'],4),fmt(d['corrected_sensor_rms_percent'],4),fmt(d['corrected_unclipped_rms_percent'],4),fmt(d['corrected_readout_rms_percent'],4),fmt(d['retained_readout_rms_percent'],4)])
        pulserows.append([pitch,tlabel,fmt(d['event_mean_average_pe'],2),fmt(d['event_max_average_pe'],1),fmt(d['largest_detected_pulse'],0),fmt(d['largest_sensor_pulse'],1)])
      d=data[tag]
      for i in np.argsort(d['clip_loss_percent']+d['sipm_loss_percent'])[-5:][::-1]:
        worst.append([tlabel,int(d['event_id'][i]),fmt(d['electron_drift_cm'][i],2),fmt(d['event_max_pe'][i],0),fmt(d['event_max_occupancy_percent'][i],2),fmt(d['sipm_loss_percent'][i],3),fmt(d['clip_loss_percent'][i],3)])
      checks.append([tlabel,int(d['ledger_checks']),int(d['state_repairs']),int(d['variance_limits']),fmt(1e6*d['maximum_state_repair_fraction'],3),fmt(d['pde_mean_z'],2),fmt(d['quiet_noise_rms_pe'],4)])
    lossheads=['Pitch / hole diameter','Recovery','SiPM loss (%)','ASIC loss (%)','Events clipped (%)','Retained / 100']
    widthheads=['Pitch / hole diameter','Recovery','Detected RMS (%)','SiPM RMS (%)','Unclipped RMS (%)','ASIC RMS (%)','Retained ASIC RMS (%)']
    pulseheads=['Pitch / hole diameter','Recovery','Event-mean avg PE','Event-max avg PE','Largest count','Largest charge PE']
    worstheads=['Recovery','Event ID','Drift (cm)','Maximum triggers','Occupancy (%)','SiPM loss (%)','ASIC loss (%)']
    checkheads=['Recovery','Ledger checks','State repairs','Variance limits','Max repair ppm','PDE z','Quiet RMS PE']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for col,(tag,label) in enumerate(TAGS):
      for row,key in [(0,'event_mean_pe'),(1,'event_max_pe')]:
        for pitch,d in [('5 mm',old[tag]),('4 mm',data[tag])]:
          axes[row,col].hist(d[key],bins=25,histtype='step',label=pitch)
        axes[row,col].set(title='Xe–15% He · '+label,xlabel=('Event mean' if row==0 else 'Event maximum')+' optical triggers / SiPM–250 ns',ylabel='Events');axes[row,col].legend()
        if row==1:axes[row,col].axvline(p['clip_threshold_pe'],color='crimson',ls='--')
    save(fig,'elm_pitch4_pulses.png')
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for col,(tag,label) in enumerate(TAGS):
      d=data[tag];keep=d['retained_event_mask'].astype(bool)
      for row,mask in [(0,np.ones(len(keep),dtype=bool)),(1,keep)]:
        edges=np.histogram_bin_edges(np.r_[d['detected_corrected'][mask],d['readout_corrected'][mask]],bins=25)
        for key,lab in [('detected_corrected','Detected'),('sensor_corrected','SiPM'),('readout_corrected','ASIC')]:axes[row,col].hist(d[key][mask],bins=edges,histtype='step',label=lab)
        axes[row,col].set(title=label+(' · all 100 events' if row==0 else f' · {keep.sum()} unclipped events'),xlabel='Attachment-corrected PE / deposited MeV',ylabel='Events');axes[row,col].legend()
    save(fig,'elm_pitch4_energy.png')
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for col,(tag,label) in enumerate(TAGS):
      for row,key in [(0,'sipm_loss_percent'),(1,'clip_loss_percent')]:
        for pitch,d in [('5 mm',old[tag]),('4 mm',data[tag])]:axes[row,col].scatter(d['electron_drift_cm'],d[key],s=25,alpha=.8,label=pitch)
        axes[row,col].set(title=label,xlabel='Electron reference drift (cm)',ylabel=('SiPM' if row==0 else 'ASIC')+' event charge loss (%)');axes[row,col].legend()
    save(fig,'elm_pitch4_drift.png')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for ax,(tag,label) in zip(axes,TAGS):
      d=data[tag];ax.hist(d['inclusive_occupancy_percent'],bins=np.geomspace(.001,100,90),histtype='step')
      ax.set(title=label,xscale='log',yscale='log',xlabel='Expected distinct optical cells / illuminated cells (%)',ylabel='Input SiPM–bucket cells')
    save(fig,'elm_pitch4_occupancy.png')
    description=('Xe–15% He only: the same 100 bb0nu source tracks and primary-pair realizations as the preceding pilot. Electron reference drift is sampled uniformly over 10–150 cm using the same deterministic placement stream; the saved reference and current events have identical IDs and drift distances. Ion distance is its complement within 150 cm. Original hit coordinates are retained, and each event uses its reference distance for diffusion and attachment. Only electron photon collection and response are replayed.')
    optics=(f'The square hole pitch is 4 mm; hole and guide diameter are 2.5 mm. EL depth is 5 mm at 15 bar and E/P = 2 kV/cm/bar, with unit electron collection. Guides have a 0.5 mm radial xenon gap. Uniform axial EL and absorbing walls give a recomputed exit acceptance. TPB conversion 0.65, forward factor 0.5 and guide trapping 1−1/1.5 apply once, giving emitted-VUV to saved-blue-photon probability {float(data["12p5"]["optical_acceptance"]):.8f}. Propagation and NOA 87 coupling are ideal.')
    sensor=(f'S14160-3015PS remains 3 × 3 mm² with 39,984 physical cells, PDE 0.32, gain 360,000, +4 V and 25 C. Uniform illumination of the 2.5 mm circular guide footprint uses {int(p["illuminated_microcells"]):,} cells, rounded from 39984 × π × 1.25² / 9. The remaining {int(p["physical_microcells"]-p["illuminated_microcells"]):,} cells receive dark triggers. The full-device 700 kcps dark rate is split by cell count. The validated statistical recovery model runs at 128 substeps per 250 ns bucket, with 12.5 and 25 ns provisional recovery constants. These are sensitivity inputs rather than measured recharge constants. PDE is the packaged peak-value approximation; crosstalk, afterpulsing and gain spread remain deferred.')
    definitions=(f'Optical trigger counts are sampled PDE successes; recurrent firing is allowed. Charge-equivalent PE is collected charge divided by e × gain. Distinct optical occupancy is the conditional expected number of cells reached by the triggers, divided by the illuminated count. Mean/max comparisons select the same photon-positive input cells at every stage. Electronics adds {float(p["noise_rms_pe"]):.4f} PE RMS Gaussian noise per bucket and clips at 1 mA ({float(p["clip_threshold_pe"]):.2f} charge-equivalent PE). Signed values are retained. SiPM loss compares optical charge with optical triggers; ASIC loss compares clipped charge with the incoming sensor charge. Event widths are full-distribution RMS/mean after known attachment and deposited-energy correction. Electronics totals subtract the expected dark baseline and include all acquired cells on illuminated channels plus four quiet buckets at either end.')
    selection=('The retained sample rejects an entire event when any acquired bucket clips at the ASIC. The table reports both its efficiency and the resulting conditional energy width. This pilot uses the simulated clipping flag; it is not an energy correction. No amplitude or signal-to-noise cut is applied. Both recovery cases reuse the same blue photons with independent response draws, and the 5 mm comparison uses saved results without new simulation.')
    timing=(f'Photon replay/archive: {float(meta["photon_wall_s"]):.2f} s; 12.5 ns response: {float(pars["12p5"]["response_wall_s"]):.2f} s; 25 ns response: {float(pars["25"]["response_wall_s"]):.2f} s; combined: {float(meta["total_wall_s"]):.2f} s. Timings include stage compilation and I/O; process startup and presentation are separate.')
    validation=('Thirty-one tests cover selected-gas cards, geometry, placement, illuminated cell count, dark-rate conservation, charge bounds, recurrent firing and substep convergence. Saved-event checks verify primary checksums, exact matching of IDs and placements against 5 mm, PDE sampling, charge ledgers and quiet noise. Rare second-moment feasibility repairs preserve sampled charge and the first moment; their maximum size and frequency are bounded and recorded. No charge-variance limits occurred.')
    sections=[('Setup',description+'\n\n'+optics+'\n\n'+sensor),('Definitions and event selection',definitions+'\n\n'+selection),('Timing and validation',timing+'\n\n'+validation)]
    figures=['elm_pitch4_pulses.png','elm_pitch4_energy.png','elm_pitch4_drift.png','elm_pitch4_occupancy.png']
    body='<h1>ELM · 4 mm pitch, 2.5 mm holes · Xe–He pilot</h1>'
    for title,text in sections[:1]:body+=f'<h2>{title}</h2>'+''.join('<p>'+html.escape(t)+'</p>' for t in text.split('\n\n'))
    body+='<h2>Charge losses and retained fraction</h2>'+table(lossheads,lossrows)+'<h2>Corrected event widths, before and after rejecting clipped events</h2>'+table(widthheads,widthrows)+'<h2>Pulse means and maxima</h2>'+table(pulseheads,pulserows)
    for title,text in sections[1:]:body+=f'<h2>{title}</h2>'+''.join('<p>'+html.escape(t)+'</p>' for t in text.split('\n\n'))
    body+=''.join(f'<img src="{name}" alt="{name}">' for name in figures)+'<h2>Largest event losses</h2>'+table(worstheads,worst)+'<h2>Validation</h2>'+table(checkheads,checks)
    cards=data['12p5']['detector_card'].decode()+'\n'+pars['12p5']['sensor_card'].decode()+'\n'+pars['25']['sensor_card'].decode()+'\n'+p['electronics_card'].decode()
    body+='<p><a href="../../latex_generated/electron_readout/elm_sipm_pitch4mm_pilot/elm_sipm_pitch4mm_pilot.pdf">PDF</a> · <a href="files.toml">Files and cards</a>. Run <code>make elm-pitch4-pilot</code>; regenerate presentation with <code>make elm-pitch4-report</code>.</p><details><summary>Cards</summary><pre>'+html.escape(cards)+'</pre></details>'
    (OUT/'elm_sipm_pitch4mm_pilot.html').write_text('<!doctype html><meta charset="utf-8"><title>ELM 4 mm pilot</title><style>body{font:16px system-ui;line-height:1.6;max-width:1400px;margin:35px auto;padding:0 22px;color:#17334a}table{border-collapse:collapse;margin:18px 0;width:100%}td,th{padding:8px;text-align:right;border-bottom:1px solid #ccd6de}td:first-child,th:first-child{text-align:left}img{width:100%;margin:15px 0}pre{white-space:pre-wrap}</style>'+body)
    tex=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage[margin=22mm]{geometry}
\usepackage{graphicx,booktabs,amsmath,xurl,hyperref,fvextra}
\begin{document}\title{ELM: 4 mm pitch, 2.5 mm holes\\Xe--He pilot}\author{NAUSICAA}\date{}\maketitle
'''
    for title,text in sections[:1]:tex+=r'\section*{'+esc(title)+'}\n'+esc(text)+'\n'
    tex+=r'\section*{Losses and retained fraction}'+tex_table(lossheads,lossrows)+r'\section*{Corrected event widths}'+tex_table(widthheads,widthrows)+tex_table(pulseheads,pulserows)
    for title,text in sections[1:]:tex+=r'\section*{'+esc(title)+'}\n'+esc(text)+'\n'
    for name in figures:tex+=r'\clearpage\begin{figure}[p]\centering\includegraphics[width=\textwidth,height=.88\textheight,keepaspectratio]{'+name+r'}\end{figure}'+'\n'
    tex+=r'\clearpage\section*{Largest event losses}'+tex_table(worstheads,worst)+r'\section*{Validation}'+tex_table(checkheads,checks)
    tex+=r'\section*{Reproduction}\path{results/elm_sipm_pitch4mm_pilot/elm_sipm_pitch4mm_pilot.h5}; file index \path{results/elm_sipm_pitch4mm_pilot/files.toml}. Run \texttt{make elm-pitch4-pilot}; presentation only \texttt{make elm-pitch4-report}.'+'\n'
    tex+=r'\section*{Input cards}\begin{Verbatim}[breaklines,breakanywhere,fontsize=\footnotesize]'+'\n'+cards+'\n'+r'\end{Verbatim}'+'\n'+r'\end{document}'+'\n'
    (LATEX/'elm_sipm_pitch4mm_pilot.tex').write_text(tex);build_report(LATEX/'elm_sipm_pitch4mm_pilot.tex')
    print(OUT/'elm_sipm_pitch4mm_pilot.html');print(LATEX/'elm_sipm_pitch4mm_pilot.pdf')
if __name__=='__main__':main()
