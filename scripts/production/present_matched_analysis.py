#!/usr/bin/env python3
"""Present matched ion/CSA and Xe–He/SiPM production from saved Julia reductions."""
from pathlib import Path
import html,sys,shutil,tomllib
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/documentation'))
from build_report import build_report
from present_elm_sipm_pilot import table,tex_table,esc,fmt
OUT=ROOT/'results/elm_matched_production'
LATEX=ROOT/'latex_generated/electron_readout/elm_matched_production'
def save(fig,name):
    fig.tight_layout();fig.savefig(OUT/name,dpi=160);plt.close(fig);shutil.copy2(OUT/name,LATEX/name)
def main():
    LATEX.mkdir(parents=True,exist_ok=True)
    with h5py.File(OUT/'matched_analysis.h5') as f:
        assert f['metadata/complete'][()]
        ion={k:f['ions'][k][()] for k in f['ions']}
        validation={k:f['validation'][k][()] for k in f['validation']}
        electron_branch=f["metadata/electron_branch"][()].decode()
        timings=tomllib.loads(f['metadata/run_timings_toml'][()].decode())
    with h5py.File(OUT/'electron_analysis.h5') as f:
        assert f['metadata/complete'][()]
        ele={k:f[electron_branch][k][()] for k in f[electron_branch]}
        pars={k:f[electron_branch.split('/')[0]+'/parameters'][k][()] for k in f[electron_branch.split('/')[0]+'/parameters']}
    graph=tomllib.loads((OUT/'manifest.toml').read_text());n=int(ion['event_count'])
    assert n==int(ele['event_count'])==graph['event_count']
    summary=[['Matched events',n],
      ['Corrected CSA maximum fast-model error (ENC)',fmt(validation['maximum_remaining_control_error_ENC'],5)],
      ['CSA maximum solver-refinement error (ENC)',fmt(validation['maximum_solver_refinement_error_ENC'],5)],['Mean initial pairs/event',fmt(ion['average_initial_ions'],2)],
      ['Mean detected optical PE/event',fmt(ele['detected_mean_total_pe'],2)],
      ['Ion CSA total noiseless bias (%)',fmt(ion['noiseless_total_bias_percent'],5)],
      ['Ion CSA total signed reconstructed bias (%)',fmt(ion['reconstructed_total_bias_percent'],5)],
      ['Events with a 50% compressed CSA history (%)',fmt(ion['compression50_event_percent'],3)],
      ['SiPM optical charge loss (%)',fmt(ele['sipm_charge_lost_percent'],5)],
      ['Electron ASIC charge loss (%)',fmt(ele['asic_charge_lost_percent'],5)],
      ['Electron events with clipping (%)',fmt(ele['clipped_event_percent'],3)],
      ['Input-positive electron pulses clipped (%)',fmt(ele['clipped_pulse_percent'],5)],
      ['Electron events retained',int(ele['retained_events'])],
      ['Electron retained fraction (%)',fmt(100*ele['retained_fraction'],3)]]
    pulses=[]
    for label,unit,mean,maxi,source,prefix in [
      ('Ion input','ions','event_mean_input','event_max_input',ion,'average_'),
      ('Noiseless CSA','ion equivalents','event_mean_noiseless','event_max_noiseless',ion,'average_'),
      ('Signed CSA output','ion equivalents','event_mean_reconstructed','event_max_reconstructed',ion,'average_'),
      ('Blue photons at SiPM','photons','event_mean_photons','event_max_photons',ele,'average_'),
      ('PDE optical triggers','PE','event_mean_pe','event_max_pe',ele,'average_'),
      ('SiPM optical charge','PE equivalents','event_mean_sensor_pe','event_max_sensor_pe',ele,'average_'),
      ('Signed electron ASIC','PE equivalents','event_mean_after_pe','event_max_after_pe',ele,'average_')]:
        pulses.append([label,unit,fmt(source[prefix+mean],3),fmt(source['median_'+mean],3),fmt(source[prefix+maxi],3),fmt(source['median_'+maxi],3)])
    worst=[[int(ele['event_id'][i]),fmt(ele['electron_drift_cm'][i],2),fmt(ele['event_max_pe'][i],0),fmt(ele['sipm_loss_percent'][i],4),fmt(ele['clip_loss_percent'][i],4)] for i in np.argsort(ele['clip_loss_percent'])[-10:][::-1]]
    widths=[[label,fmt(ele[f'corrected_{key}_rms_percent'],5),fmt(ele[f'retained_{key}_rms_percent'],5)] for key,label in [('detected','Detected optical PE'),('sensor','After SiPM'),('unclipped','Noisy, before clipping'),('readout','After electron ASIC')]]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    axes[0].hist(ion['ion_drift_cm'],bins=np.linspace(5,145,29),histtype='step',label='Ions')
    axes[0].hist(ion['electron_drift_cm'],bins=np.linspace(5,145,29),histtype='step',label='Electrons')
    axes[0].set(xlabel='Reference drift (cm)',ylabel='Events');axes[0].legend()
    axes[1].scatter(ion['ion_drift_cm'],ion['electron_drift_cm'],s=4,alpha=.25)
    axes[1].set(xlabel='Ion reference drift (cm)',ylabel='Electron reference drift (cm)',title='Shared placement: sum = 150 cm')
    save(fig,'matched_placement.png')
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for row,stat in enumerate(['mean','max']):
      for col,(source,key1,key2,label,unit) in enumerate([(ion,f'event_{stat}_input',f'event_{stat}_reconstructed','Ions / CSA','Ion-charge units'),(ele,f'event_{stat}_pe',f'event_{stat}_after_pe','Electrons / ASIC','PE / charge-equivalent PE')]):
        edges=np.histogram_bin_edges(np.r_[source[key1],source[key2]],bins=60)
        axes[row,col].hist(source[key1],bins=edges,histtype='step',label='Input')
        axes[row,col].hist(source[key2],bins=edges,histtype='step',label='Signed electronics output')
        axes[row,col].set(title=label+' · event '+('mean' if stat=='mean' else 'maximum'),xlabel=unit+' per channel/bucket',ylabel='Events',yscale='log');axes[row,col].legend()
    save(fig,'matched_pulses.png')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for ax,key,label in zip(axes,['event_mean_photons','event_max_photons'],['Event mean','Event maximum']):
      ax.hist(ele[key],bins=65,histtype='step');ax.set(xlabel=label+' blue photons per SiPM / 250 ns',ylabel='Events',yscale='log')
    save(fig,'matched_photons.png')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5));keep=ele['retained_event_mask'].astype(bool)
    for ax,mask,label in [(axes[0],np.ones(n,dtype=bool),'All events'),(axes[1],keep,'Events without electron ASIC clipping')]:
      values=[ele[k][mask] for k in ['detected_corrected','sensor_corrected','readout_corrected']]
      edges=np.histogram_bin_edges(np.concatenate(values),bins=65)
      for v,name in zip(values,['Detected optical PE','SiPM','Electron ASIC']):ax.hist(v,bins=edges,histtype='step',label=name)
      ax.set(title=label,xlabel='Attachment-corrected PE / deposited MeV',ylabel='Events',yscale='log');ax.legend()
    save(fig,'matched_electron_energy.png')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    axes[0].scatter(ion['ion_drift_cm'],100*(ion['noiseless_ratio']-1),s=4,alpha=.3)
    axes[0].set(xlabel='Ion reference drift (cm)',ylabel='Noiseless CSA total charge bias (%)')
    axes[1].scatter(ele['electron_drift_cm'],ele['clip_loss_percent'],s=4,alpha=.3)
    axes[1].set(xlabel='Electron reference drift (cm)',ylabel='Electron ASIC charge loss (%)')
    save(fig,'matched_losses.png')
    setup=(f'{n:,} bb0nu events share one saved realization. Ion reference drift is uniform over 5–145 cm in a 150 cm chamber; electron reference drift is its exact complement. Primary-pair draws are anchored to the existing reference card, so changing placement does not change the pair realization. Event IDs, primary checksums and both distances are validated across common realization, ion collection, photon collection, corrected CSA, SiPM charge and electron ASIC output. The 4 mm ion pitch and 40 ms frames correspond to 4 mm drift buckets. Electron readout uses a 4 mm hole pitch and 250 ns (0.25 mm) buckets.')
    optical=(f'Xe–15% He at 15 bar, E/P = 2 kV/cm/bar, 2.5 mm diameter holes/guides and 5 mm EL depth. A 0.5 mm radial xenon gap surrounds each guide. Unit electron collection and absorbing EL walls are assumed. Analytical hole-exit acceptance, TPB conversion 0.65, forward factor 0.5 and meridional capture 1−1/1.5 give total optical probability {float(ele["optical_acceptance"]):.8f}. Guide propagation and NOA 87 coupling are ideal. Saved photons are blue photons at the SiPM; PDE is applied once afterward.')
    sipm=(f'S14160-3015PS: 3 × 3 mm², 39,984 cells, 21,808 illuminated by the centered 2.5 mm disk, PDE 0.32, gain 360,000 and total dark rate 700 kcps. Illuminated and dark-only cell populations retain independent recovery states. The baseline uses a provisional 12.5 ns recovery constant with 128 statistical substeps per bucket. Electron electronics adds {float(pars["noise_rms_pe"]):.4f} PE RMS and clips at 1 mA ({float(pars["clip_threshold_pe"]):.2f} PE per bucket). Crosstalk, afterpulsing and gain spread remain deferred.')
    csa=('The ion arm uses the calibrated CSA v3 circuit, 5 ms averaging and the 40 ms response/noise covariance. The linear fast path is validated against transistor controls at 0.5 ENC tolerance; histories below the 0.63 V guard receive direct LTspice treatment, with adaptive controls and solver refinement. The simulator binary and circuit calibration are reused only with matching content hashes. The signed three-frame inverse reconstructs ion charge, including quiet/recovery frames and temporal noise correlations. A 50% compression diagnostic identifies events containing such a pad history; it is separate from electron ASIC clipping.')
    definitions=('Pulse tables use all input-positive channel/bucket cells, preserving those same cells through electronics. Each event contributes one mean and maximum. Signed response values remain intact. Ion total charge bias sums complete stored histories relative to collected ions. Electron energy widths are full-distribution RMS/mean after correcting known attachment survival and dividing by deposited MC energy; noisy totals subtract expected dark charge. The retained electron sample rejects an event when any electron ASIC bucket clips.')
    timing=(f'Detector production: {timings["detector_wall_s"]:.2f} s; concurrent SiPM/ASIC: {timings["sipm_wall_s"]:.2f} s; concurrent CSA pipeline: {timings["csa_wall_s"]:.2f} s; complete paired production and validation: {timings["total_wall_s"]:.2f} s. Reduction and presentation are separate.') if timings else 'Workflow timing was not recorded for this run.'
    sections=[('Shared production',setup),('Electron optical and response model',optical+'\n\n'+sipm),('Ion electronics',csa),('Definitions',definitions),('Timing',timing)]
    body='<h1>ELM · matched ion and Xe–He response</h1>'+table(['Quantity','Value'],summary)+'<h2>Per-channel/bucket pulse statistics</h2>'+table(['Stage','Unit','Mean average','Mean median','Maximum average','Maximum median'],pulses)+'<h2>Electron energy widths</h2>'+table(['Stage','All-event RMS (%)','Retained RMS (%)'],widths)
    body+='<h2>Events with the largest electron ASIC loss</h2>'+table(['Event ID','Electron drift (cm)','Maximum input PE','SiPM loss (%)','ASIC loss (%)'],worst)
    for title,text in sections:body+=f'<h2>{title}</h2>'+''.join('<p>'+html.escape(t)+'</p>' for t in text.split('\n\n'))
    figures=['matched_placement.png','matched_pulses.png','matched_photons.png','matched_electron_energy.png','matched_losses.png']
    body+=''.join(f'<img src="{name}" alt="{name}">' for name in figures)
    body+='<p><a href="../../latex_generated/electron_readout/elm_matched_production/elm_matched_production.pdf">PDF report</a> · <a href="manifest.toml">Combined manifest</a> · <a href="../../config/specs/runs/elm_matched_response_viewer.toml">Viewer card</a>.</p>'
    cards=ele['detector_card'].decode()+'\n'+pars['sensor_card'].decode()+'\n'+pars['electronics_card'].decode()+'\n'+ion['response_card'].decode()+'\n'+ion['reconstruction_card'].decode()+'\n'+ele['realisation_card'].decode()
    body+='<details><summary>Resolved cards</summary><pre>'+html.escape(cards)+'</pre></details>'
    (OUT/'elm_matched_production.html').write_text('<!doctype html><meta charset="utf-8"><title>Matched ELM response</title><style>body{font:16px system-ui;line-height:1.6;max-width:1400px;margin:35px auto;padding:0 22px;color:#17334a}table{border-collapse:collapse;margin:18px 0;width:100%}td,th{padding:8px;text-align:right;border-bottom:1px solid #ccd6de}td:first-child,th:first-child{text-align:left}img{width:100%;margin:15px 0}pre{white-space:pre-wrap}</style>'+body)
    tex=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage[margin=22mm]{geometry}
\usepackage{graphicx,booktabs,amsmath,xurl,hyperref,fvextra}
\begin{document}\title{ELM: matched ion and Xe--He response\\4 mm pitch, shared 5--145 cm placement}\author{NAUSICAA}\date{}\maketitle
'''
    tex+=tex_table(['Quantity','Value'],summary)+tex_table(['Stage','Unit','Mean average','Mean median','Max average','Max median'],pulses)+tex_table(['Stage','All-event RMS (%)','Retained RMS (%)'],widths)
    tex+=tex_table(['Event ID','Electron drift (cm)','Maximum input PE','SiPM loss (%)','ASIC loss (%)'],worst)
    for title,text in sections:tex+=r'\section*{'+esc(title)+'}\n'+esc(text)+'\n'
    for name in figures:tex+=r'\clearpage\begin{figure}[p]\centering\includegraphics[width=\textwidth,height=.88\textheight,keepaspectratio]{'+name+r'}\end{figure}'+'\n'
    tex+=r'\clearpage\section*{Reproduction}\path{results/elm_matched_production/manifest.toml} resolves all physical histories and their matching validation. Numerical reductions are \path{matched_analysis.h5} and \path{electron_analysis.h5} in that directory. Full workflow: \texttt{make elm-matched-production}; presentation: \texttt{make elm-matched-report}.'+'\n'
    tex+=r'\section*{Resolved cards}\begin{Verbatim}[breaklines,breakanywhere,fontsize=\footnotesize]'+'\n'+cards+'\n'+r'\end{Verbatim}'+'\n'+r'\end{document}'+'\n'
    (LATEX/'elm_matched_production.tex').write_text(tex);build_report(LATEX/'elm_matched_production.tex')
    print(OUT/'elm_matched_production.html');print(LATEX/'elm_matched_production.pdf')
if __name__=='__main__':main()
