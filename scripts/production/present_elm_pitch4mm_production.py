#!/usr/bin/env python3
"""Present the full 4 mm Xe–He production from Julia event-level reductions."""
from pathlib import Path
import html,sys,shutil
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/documentation'))
from build_report import build_report
from present_elm_sipm_pilot import table,tex_table,esc,fmt
OUT=ROOT/'results/elm_sipm_pitch4mm_production'
LATEX=ROOT/'latex_generated/electron_readout/elm_sipm_pitch4mm_production'
def save(fig,name):
    fig.tight_layout();fig.savefig(OUT/name,dpi=160);plt.close(fig);shutil.copy2(OUT/name,LATEX/name)
def main():
    LATEX.mkdir(parents=True,exist_ok=True)
    with h5py.File(OUT/'elm_sipm_pitch4mm_production.h5','r') as f:
        assert f['metadata/complete'][()]
        meta={k:f['metadata'][k][()] for k in f['metadata']}
        d={k:f['tau_12p5/xe_he15'][k][()] for k in f['tau_12p5/xe_he15']}
        p={k:f['tau_12p5/parameters'][k][()] for k in f['tau_12p5/parameters']}
    n=int(d['event_count']);keep=d['retained_event_mask'].astype(bool)
    summary=[['Events',n],['Mean detected optical PE/event',fmt(d['detected_mean_total_pe'],1)],
        ['SiPM optical charge loss (%)',fmt(d['sipm_charge_lost_percent'],5)],
        ['ASIC charge loss (%)',fmt(d['asic_charge_lost_percent'],5)],
        ['Events with clipping (%)',fmt(d['clipped_event_percent'],3)],
        ['Input pulses with clipping (%)',fmt(d['clipped_pulse_percent'],5)],
        ['Retained events',int(d['retained_events'])],['Retained fraction (%)',fmt(100*d['retained_fraction'],3)],
        ['Largest optical trigger count/bucket',fmt(d['largest_detected_pulse'],0)],
        ['Largest optical SiPM charge/bucket (PE)',fmt(d['largest_sensor_pulse'],2)],
        ['Maximum expected distinct-cell occupancy (%)',fmt(d['largest_expected_occupancy_percent'],2)]]
    widths=[[name,fmt(d[f'corrected_{key}_rms_percent'],5),fmt(d[f'retained_{key}_rms_percent'],5)] for key,name in [('detected','Detected optical PE'),('sensor','After SiPM, optical charge'),('unclipped','Noisy, before clipping'),('readout','After ASIC')]]
    pulses=[]
    for stage,mean_key,max_key in [('Detected optical PE','event_mean_pe','event_max_pe'),('SiPM optical charge','event_mean_sensor_pe','event_max_sensor_pe'),('After ASIC','event_mean_after_pe','event_max_after_pe')]:
        pulses.append([stage,*[fmt(d[prefix+key],3) for key in [mean_key,max_key] for prefix in ['average_','median_']]])
    worst=[]
    for i in np.argsort(d['clip_loss_percent']+d['sipm_loss_percent'])[-15:][::-1]:
        worst.append([int(d['event_id'][i]),fmt(d['electron_drift_cm'][i],2),fmt(d['event_max_pe'][i],0),fmt(d['event_max_sensor_pe'][i],1),fmt(d['sipm_loss_percent'][i],4),fmt(d['clip_loss_percent'][i],4)])
    validation=[['Checked event ledgers',int(d['ledger_checks'])],['PDE estimate',fmt(d['pde_fraction'],7)],['PDE z-score',fmt(d['pde_mean_z'],3)],['Quiet noise RMS (PE)',fmt(d['quiet_noise_rms_pe'],5)],['Quiet noise mean (PE)',fmt(d['quiet_noise_mean_pe'],6)],['Second-moment state repairs',int(d['state_repairs'])],['Maximum repair / smaller population (ppm)',fmt(1e6*d['maximum_state_repair_fraction'],5)],['Largest next-trigger RMS impact bound (PE)',fmt(d['maximum_repair_single_trigger_rms_pe'],6)],['Charge-variance limits',int(d['variance_limits'])]]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for row,(mean_key,max_key) in enumerate([('event_mean_pe','event_max_pe'),('event_mean_after_pe','event_max_after_pe')]):
      for col,key in enumerate([mean_key,max_key]):
        axes[row,col].hist(d[key],bins=60,histtype='step',lw=1.5)
        axes[row,col].set(title=('Detected optical PE' if row==0 else 'After ASIC')+' · '+('event mean' if col==0 else 'event maximum'),xlabel='PE or charge-equivalent PE / SiPM–250 ns bucket',ylabel='Events',yscale='log')
        if col==1:axes[row,col].axvline(p['clip_threshold_pe'],color='crimson',ls='--',label='ASIC threshold');axes[row,col].legend()
    save(fig,'elm_full_pulses.png')
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for ax,mask,title in [(axes[0],np.ones(n,dtype=bool),f'All {n:,} events'),(axes[1],keep,f'{keep.sum():,} unclipped events')]:
      arrays=[d[k][mask] for k in ['detected_corrected','sensor_corrected','readout_corrected']]
      edges=np.histogram_bin_edges(np.concatenate(arrays),bins=70)
      for values,label in zip(arrays,['Detected','After SiPM','After ASIC']):ax.hist(values,bins=edges,histtype='step',label=label)
      ax.set(title=title,xlabel='Attachment-corrected PE / deposited MeV',ylabel='Events',yscale='log');ax.legend()
    save(fig,'elm_full_energy.png')
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for ax,key,label in zip(axes.flat,['sipm_loss_percent','clip_loss_percent','event_max_sensor_pe','event_max_occupancy_percent'],['SiPM optical charge loss (%)','ASIC charge loss (%)','Maximum optical SiPM charge (PE)','Maximum expected distinct-cell occupancy (%)']):
      ax.scatter(d['electron_drift_cm'],d[key],s=5,alpha=.35,rasterized=True)
      ax.set(xlabel='Electron reference drift (cm)',ylabel=label)
      if key=='event_max_sensor_pe':ax.axhline(p['clip_threshold_pe'],color='crimson',ls='--')
    save(fig,'elm_full_drift.png')
    setup=(f'Full available bb0nu sample: {n:,} events, Xe–15% He only. E/P = 2 kV/cm/bar at 15 bar; 4 mm square-grid pitch, 2.5 mm hole and light-guide diameter, 5 mm EL depth, and 0.5 mm radial xenon gap. Electron reference distances are sampled uniformly over 10–150 cm with deterministic event seeds. Ion distances are complementary within 150 cm. Source hits and primary-pair realizations are replayed; event reference distance determines diffusion and attachment. No ion or pure-Xe response was generated.')
    optics=(f'Electron collection is unity. Uniform axial EL with absorbing walls gives analytical hole-exit acceptance. TPB conversion 0.65, forward fraction 0.5 and meridional guide capture 1−1/1.5 apply once; emitted-VUV to saved-blue probability is {float(d["optical_acceptance"]):.8f}. Propagation and NOA 87 coupling are ideal. Saved blue photons are reused by the response stage.')
    sensor=(f'S14160-3015PS: 3 × 3 mm², 39,984 physical cells, 21,808 uniformly illuminated cells under the centered 2.5 mm guide. Packaged PDE is 0.32 (peak-value approximation), gain 360,000 at +4 V and 25 C. Total primary dark rate is 700 kcps over the full device, split by cell count between illuminated and dark-only populations. Statistical recovery uses a provisional 12.5 ns exponential constant and 128 substeps per 250 ns bucket. Recurrent firing is allowed; PDE remains constant with recovery. Crosstalk, afterpulsing and gain spread are deferred.')
    definitions=(f'Optical trigger counts are PDE successes. Charge-equivalent PE is charge divided by e × gain. Distinct optical occupancy is the conditional expected number of cells reached by triggers, divided by the illuminated count. Every event contributes one mean and maximum over its photon-positive input cells, with those same cells used at all stages. SiPM loss compares optical charge with trigger count; ASIC loss compares clipped charge with incoming sensor charge. Electronics adds {float(p["noise_rms_pe"]):.4f} PE RMS Gaussian noise and clips the bucket-average current at 1 mA ({float(p["clip_threshold_pe"]):.2f} PE). Signed outputs are retained.')
    energy=('Event sums use all acquired buckets on illuminated channels, including four quiet buckets at each end. Energy diagnostics correct known attachment survival and normalize by deposited MC energy. Optical and sensor widths use optical quantities; noisy totals subtract the expected dark baseline. Reported widths are full-distribution RMS/mean, including tails. The retained sample rejects an entire event if any acquired ASIC bucket clipped. Its width is conditional on that rejection and is reported with the retained fraction. No signal-amplitude cut or SiPM compression correction is applied.')
    timing=(f'Photon replay/archive: {float(meta["photon_wall_s"]):.2f} s; four-worker response and assembly: {float(meta["response_wall_s"]):.2f} s; total timed production: {float(meta["total_wall_s"]):.2f} s. Independent workers own disjoint 100-event shards and preserve event/channel seeds. Numerical reduction streams quiet-noise moments and retains event summaries; complete bucket histories remain in the independently reusable photon, sensor and ASIC files.')
    verification=('Validation checks source/primary ledgers, selected gas identity, count and charge bounds, PDE sampling, quiet-noise statistics and numerical recovery diagnostics. The streaming multi-shard reduction was checked against the pilot with identical event-level charge and resolution results. Rare feasibility repairs affect the second recovery moment while preserving charge and the first moment; their next-single-trigger variance effect is bounded by the second-moment correction divided by cell count. Its square root is checked against 1% of the electronics noise RMS, and the repair frequency against one per million batches.')
    sections=[('Setup',setup+'\n\n'+optics+'\n\n'+sensor),('Definitions',definitions+'\n\n'+energy),('Timing and validation',timing+'\n\n'+verification)]
    cards=d['detector_card'].decode()+'\n'+p['sensor_card'].decode()+'\n'+p['electronics_card'].decode()
    body='<h1>ELM · full Xe–He production · 4 mm pitch</h1>'+table(['Quantity','Value'],summary)+'<h2>Energy widths and clipping-event rejection</h2>'+table(['Stage','All-event RMS (%)','Retained RMS (%)'],widths)+'<h2>Pulse statistics</h2>'+table(['Stage','Event-mean average','Event-mean median','Event-max average','Event-max median'],pulses)
    for title,text in sections:body+=f'<h2>{title}</h2>'+''.join('<p>'+html.escape(t)+'</p>' for t in text.split('\n\n'))
    figures=['elm_full_pulses.png','elm_full_energy.png','elm_full_drift.png']
    body+=''.join(f'<img src="{name}" alt="{name}">' for name in figures)+'<h2>Largest event losses</h2>'+table(['Event ID','Drift (cm)','Max triggers','Max optical charge PE','SiPM loss (%)','ASIC loss (%)'],worst)+'<h2>Validation ledger</h2>'+table(['Check','Value'],validation)
    body+='<p><a href="../../latex_generated/electron_readout/elm_sipm_pitch4mm_production/elm_sipm_pitch4mm_production.pdf">PDF</a> · <a href="files.toml">Files and cards</a>. Workflow: <code>make elm-pitch4-production</code>; presentation: <code>make elm-pitch4-production-report</code>.</p><details><summary>Cards</summary><pre>'+html.escape(cards)+'</pre></details>'
    (OUT/'elm_sipm_pitch4mm_production.html').write_text('<!doctype html><meta charset="utf-8"><title>ELM full Xe–He production</title><style>body{font:16px system-ui;line-height:1.6;max-width:1400px;margin:35px auto;padding:0 22px;color:#17334a}table{border-collapse:collapse;margin:18px 0;width:100%}td,th{padding:8px;text-align:right;border-bottom:1px solid #ccd6de}td:first-child,th:first-child{text-align:left}img{width:100%;margin:15px 0}pre{white-space:pre-wrap}</style>'+body)
    tex=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage[margin=22mm]{geometry}
\usepackage{graphicx,booktabs,amsmath,xurl,hyperref,fvextra}
\begin{document}\title{ELM: full Xe--He production\\4 mm pitch, 2.5 mm holes}\author{NAUSICAA}\date{}\maketitle
'''
    tex+=tex_table(['Quantity','Value'],summary)+tex_table(['Stage','All-event RMS (%)','Retained RMS (%)'],widths)+tex_table(['Stage','Mean average','Mean median','Max average','Max median'],pulses)
    for title,text in sections:tex+=r'\section*{'+esc(title)+'}\n'+esc(text)+'\n'
    for name in figures:tex+=r'\clearpage\begin{figure}[p]\centering\includegraphics[width=\textwidth,height=.88\textheight,keepaspectratio]{'+name+r'}\end{figure}'+'\n'
    tex+=r'\clearpage\section*{Largest event losses}'+tex_table(['Event ID','Drift (cm)','Max triggers','Max charge PE','SiPM loss (%)','ASIC loss (%)'],worst)+r'\section*{Validation ledger}'+tex_table(['Check','Value'],validation)
    tex+=r'\section*{Reproduction}\path{results/elm_sipm_pitch4mm_production/elm_sipm_pitch4mm_production.h5}; file index \path{results/elm_sipm_pitch4mm_production/files.toml}. Run \texttt{make elm-pitch4-production}; regenerate presentation with \texttt{make elm-pitch4-production-report}.'+'\n'
    tex+=r'\section*{Input cards}\begin{Verbatim}[breaklines,breakanywhere,fontsize=\footnotesize]'+'\n'+cards+'\n'+r'\end{Verbatim}'+'\n'+r'\end{document}'+'\n'
    (LATEX/'elm_sipm_pitch4mm_production.tex').write_text(tex);build_report(LATEX/'elm_sipm_pitch4mm_production.tex')
    print(OUT/'elm_sipm_pitch4mm_production.html');print(LATEX/'elm_sipm_pitch4mm_production.pdf')
if __name__=='__main__':main()
