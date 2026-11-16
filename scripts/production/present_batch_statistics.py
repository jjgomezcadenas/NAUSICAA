#!/usr/bin/env python3
"""Statistics page of a run kept as per-batch hit files.

  python scripts/production/present_batch_statistics.py results/readout/<run_label> [classification.toml] [--generated=N] [--comparison=DIR]

Reads <run_label>/<output_subdir>/statistics/events.h5, written by
prepare_batch_statistics.jl, and writes plots and index.html beside it. No
physical calculation is repeated: histograms are counts of the saved per-event
values, and the energy fits are read from the file. The optional TOML holds
the Monte Carlo classification of detached components of another sample.
--generated starts the selection table from the generated events; --comparison
names the folder of present_species_comparison.py, whose E2 thresholds close it.
"""
import os,html,sys,tomllib
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm

BLUE,ORANGE,RED='#287ba0','#b16d31','#a63f35'
LABEL={'bb0nu':'bb0nu','single_electrons':'single electrons','bi214':'²¹⁴Bi gammas'}

def text(v): return v.decode() if isinstance(v,bytes) else v
def stats_box(ax,values,unit='',loc='right'):
    v=np.asarray(values,float);v=v[np.isfinite(v)]
    ax.text(.97 if loc=='right' else .03,.96,f"N = {len(v)}\nmean = {np.mean(v):.4g} {unit}\nstd = {np.std(v):.4g} {unit}\nmedian = {np.median(v):.4g} {unit}",
            ha=loc,va='top',transform=ax.transAxes,fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='#ccd6dc'))
def histogram(ax,values,edges,xlabel,title,unit='',color=BLUE,loc='right'):
    v=np.asarray(values,float);v=v[np.isfinite(v)]
    ax.hist(v,bins=edges,color=color,alpha=.6);ax.set(xlabel=xlabel,ylabel='Events / bin',title=title,xlim=(edges[0],edges[-1]));ax.grid(alpha=.15)
    out=int(np.count_nonzero((v<edges[0])|(v>edges[-1])))
    stats_box(ax,v,unit,loc)
    if out: ax.text(.97 if loc=='right' else .03,.60,f"outside the axis: {out}",ha=loc,va='top',transform=ax.transAxes,fontsize=8,color='#526776')

def load(run_dir):
    run_dir=Path(run_dir);found=sorted(run_dir.glob('*/statistics/events.h5'))
    if len(found)!=1: raise SystemExit(f'Expected one statistics/events.h5 under {run_dir}, found {len(found)}')
    f=h5py.File(found[0]);e={k:f['events'][k][:] for k in f['events']}
    for arm in ('itrack','etrack'): e[arm+'_extent_r_mm']=np.maximum(e[arm+'_extent_x_mm'],e[arm+'_extent_y_mm'])
    meta=dict(path=found[0],out=found[0].parent,species=text(f['species'][()]),run_label=text(f['run_label'][()]),batches=[text(b) for b in f['batches'][:]],
              tc=tomllib.loads(text(f['track_analysis'][()])),dc=tomllib.loads(text(f['detector'][()])),
              ion_cut=float(f['ion_hit_cut_keV'][()]),electron_cut=float(f['electron_hit_cut_keV'][()]),
              total={k:int(v[()]) for k,v in f['selection/total'].items()},batch={k:v[:] for k,v in f['selection/batch'].items()},
              settings={k:text(v[()]) for k,v in f['selection/settings'].items()},
              fits={k:{n:float(v[()]) for n,v in g.items()} for k,g in f['energy_fits'].items()},
              by_drift={k:{n:v[:] for n,v in g.items()} for k,g in f['energy_fits_by_drift'].items()})
    return e,meta

def render(run_dir,classification=None,generated=None,comparison=None):
    e,m=load(run_dir);out=m['out'];plots=out/'plots';plots.mkdir(exist_ok=True)
    n=len(e['event_number']);L=m['dc']['chamber_length_mm'];R=m['dc']['detector_radius_mm'];tc=m['tc'];sample=LABEL.get(m['species'],m['species'])
    side=float(np.median(e['frame_side_mm']))

    # 1. Where the accepted events are
    fig,axes=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
    histogram(axes[0],e['vertex_z_mm'],np.linspace(0,L,53),'Vertex Z (mm)','Vertex along the drift · Z = 0 at the TopMetal plane','mm')
    histogram(axes[1],(np.hypot(e['vertex_x_mm'],e['vertex_y_mm'])/R)**2,np.linspace(0,1,41),'(vertex radius / detector radius)²','Vertex radius · uniform density is flat')
    histogram(axes[2],e['frame_shift_mm'],np.linspace(0,side/2,41),'Shift (mm)','Frame moved from the track centre to stay inside the detector','mm')
    axes[2].set_yscale('log')
    fig.suptitle(f'{sample} · accepted events · frame of {side:g} mm · detector length {L:g} mm, radius {R:g} mm')
    fig.savefig(plots/'event_position.png',dpi=140);plt.close(fig)

    # 2. Event energy before the hit cut
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for ax,key,title in ((axes[0],'ion_energy_keV','Ion energy estimator'),(axes[1],'electron_energy_keV','Electron energy estimator')):
        fit=m['fits'][key];mu,sigma=fit['mu'],fit['sigma'];edges=np.linspace(mu-5*sigma,mu+5*sigma,81);v=e[key][np.isfinite(e[key])]
        ax.hist(v,bins=edges,color=BLUE,alpha=.5);x=np.linspace(edges[0],edges[-1],400)
        window=norm.cdf((fit['upper']-mu)/sigma)-norm.cdf((fit['lower']-mu)/sigma)
        ax.plot(x,fit['entries']/window*(edges[1]-edges[0])*norm.pdf(x,mu,sigma),color=RED,lw=2,label='Gaussian fit')
        ax.set(xlim=(edges[0],edges[-1]),xlabel='Energy (keV)',ylabel='Events / bin',title=title);ax.grid(alpha=.15);ax.legend(loc='upper left')
        ax.text(.97,.96,f"Fit N = {fit['entries']:.0f} / {fit['total_entries']:.0f}\nμ = {mu:.2f} ± {fit['mu_error']:.2f} keV\nσ = {sigma:.2f} ± {fit['sigma_error']:.2f} keV\nσ/μ = {100*sigma/mu:.3f}%\nFWHM = {235.5*sigma/mu:.2f}%\nbelow the axis: {np.count_nonzero(v<edges[0])}",
                ha='right',va='top',transform=ax.transAxes,fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='#ccd6dc'))
    fig.suptitle(f'{sample} · full signed charge · before the hit cut');fig.savefig(plots/'event_energy.png',dpi=140);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for ax,key,title,label in ((axes[0],'ion_energy_keV','Ion energy resolution against ion drift','Ion drift (cm)'),(axes[1],'electron_energy_keV','Electron energy resolution against electron drift','Electron drift (cm)')):
        d=m['by_drift'][key];c=(d['edges'][:-1]+d['edges'][1:])/2;ok=np.isfinite(d['sigma'])
        ax.errorbar(c[ok],100*d['sigma'][ok]/d['mu'][ok],yerr=100*d['sigma_error'][ok]/d['mu'][ok],xerr=np.diff(d['edges'])[ok]/2,fmt='o',color=BLUE,capsize=3)
        for x,s,mu,k in zip(c[ok],d['sigma'][ok],d['mu'][ok],d['entries'][ok]): ax.annotate(f'{100*s/mu:.3f}%\nN={k}',(x,100*s/mu),textcoords='offset points',xytext=(0,9),ha='center',fontsize=8)
        ax.set(xlabel=label,ylabel='σ/μ (%)',title=title,xlim=(d['edges'][0],d['edges'][-1]));ax.grid(alpha=.15);ax.margins(y=.35)
    fig.suptitle(f'{sample} · Gaussian fit of the energy in six drift intervals');fig.savefig(plots/'energy_by_drift.png',dpi=140);plt.close(fig)

    # 3. Size of the main track
    fig,axes=plt.subplots(2,3,figsize=(15,8.5),layout='constrained')
    for row,arm,color in ((0,'itrack',BLUE),(1,'etrack',ORANGE)):
        histogram(axes[row,0],e[arm+'_length_mm'],np.linspace(0,500,51),'Track length (mm)',f'{arm} · path between the endpoints','mm',color)
        histogram(axes[row,1],e[arm+'_extent_r_mm'],np.linspace(0,side,41),'Extension in R (mm)',f'{arm} · largest of the X and Y ranges','mm',color)
        histogram(axes[row,2],e[arm+'_extent_z_mm'],np.linspace(0,side,41),'Extension in Z (mm)',f'{arm} · Z range','mm',color)
    fig.suptitle(f'{sample} · largest connected component of each track · all {n} accepted events');fig.savefig(plots/'track_size.png',dpi=140);plt.close(fig)

    # 4. Hits and components
    fig,axes=plt.subplots(2,3,figsize=(15,8.5),layout='constrained')
    for row,arm,color,cut in ((0,'itrack',BLUE,tc['ion_hit_cut_keV']),(1,'etrack',ORANGE,tc['electron_hit_cut_keV'])):
        histogram(axes[row,0],e[arm+'_hits'],np.linspace(0,max(400,np.percentile(e[arm+'_hits'],99.5)),51),'Hits per event',f'{arm} · hits above {cut:g} keV','',color)
        c=e[arm+'_components'];edges=np.arange(-.5,max(6,c.max())+1.5)
        axes[row,1].hist(c,bins=edges,color=color,alpha=.6);axes[row,1].set(xlabel='Connected components',ylabel='Events / bin',title=f'{arm} · components per event');axes[row,1].grid(alpha=.15);axes[row,1].set_yscale('log')
        axes[row,1].text(.97,.96,f"N = {n}\none component = {np.count_nonzero(c==1)} ({100*np.mean(c==1):.1f}%)\nmean = {c.mean():.2f}",ha='right',va='top',transform=axes[row,1].transAxes,fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='#ccd6dc'))
        multi=e[arm+'_second_keV'][c>1]
        histogram(axes[row,2],multi,np.linspace(0,300,61),'Energy of the second component (keV)',f'{arm} · events with more than one component','keV',color)
    fig.suptitle(f'{sample} · all {n} accepted events');fig.savefig(plots/'hits_components.png',dpi=140);plt.close(fig)

    # 5. Endpoint blobs
    fig,axes=plt.subplots(2,3,figsize=(15,8.5),layout='constrained')
    for row,arm,color in ((0,'itrack',BLUE),(1,'etrack',ORANGE)):
        e1,e2=e[arm+'_E1_keV'],e[arm+'_E2_keV'];ok=np.isfinite(e1)&np.isfinite(e2);edges=np.linspace(0,1000,101)
        h=axes[row,0].hist2d(e1[ok],e2[ok],bins=[edges,edges],cmap='Blues' if arm=='itrack' else 'Oranges',cmin=1)
        fig.colorbar(h[3],ax=axes[row,0],label='Events / bin');axes[row,0].plot([0,1000],[0,1000],color='#526776',lw=.8,ls='--')
        axes[row,0].set(xlabel='E₁ (keV)',ylabel='E₂ (keV)',title=f'{arm} · endpoint spheres, radius {tc["ion_blob_radius_mm" if arm=="itrack" else "electron_blob_radius_mm"]:g} mm')
        axes[row,0].text(.03,.96,f"N = {ok.sum()}\nmean E₁ = {e1[ok].mean():.1f} keV\nmean E₂ = {e2[ok].mean():.1f} keV",ha='left',va='top',transform=axes[row,0].transAxes,fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='#ccd6dc'))
        histogram(axes[row,1],e1,np.linspace(0,1000,81),'E₁ (keV)',f'{arm} · larger endpoint energy','keV',color)
        histogram(axes[row,2],e2,np.linspace(0,800,81),'E₂ (keV)',f'{arm} · smaller endpoint energy','keV',color)
    fig.suptitle(f'{sample} · endpoint energies of the largest component · all {n} accepted events');fig.savefig(plots/'blob_energies.png',dpi=140);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    edges=np.linspace(0,L/10,14);c=(edges[:-1]+edges[1:])/2
    for ax,arm,drift,color in ((axes[0],'itrack','ion_drift_cm',BLUE),(axes[1],'etrack','electron_drift_cm',ORANGE)):
        for key,marker,name in ((arm+'_E1_keV','o','E₁'),(arm+'_E2_keV','s','E₂')):
            med=[np.nanmedian(e[key][(e[drift]>=a)&(e[drift]<b)]) if np.any((e[drift]>=a)&(e[drift]<b)) else np.nan for a,b in zip(edges[:-1],edges[1:])]
            ax.plot(c,med,marker=marker,color=color,alpha=1 if name=='E₁' else .55,label=f'{name}: {np.nanmin(med):.0f}–{np.nanmax(med):.0f} keV')
        ax.set(xlabel=f'{"Ion" if arm=="itrack" else "Electron"} drift (cm)',ylabel='Median endpoint energy (keV)',title=f'{arm} · median endpoint energy against its own drift',ylim=(0,None));ax.grid(alpha=.15);ax.legend()
    fig.suptitle(f'{sample} · {n} events in 13 drift intervals');fig.savefig(plots/'blob_energy_by_drift.png',dpi=140);plt.close(fig)

    # Page
    t=m['total'];n0=t['input_events']
    rows=[('Events read (after the Geant4 energy window)',t['input_events']),('Rejected: energy deposit outside the detector',t['outside_detector_events']),
          (f"Rejected: main track within {float(m['settings'].get('el_exclusion_mm',0)):g} mm of the EL",t.get('main_track_near_el_events',0)),
          ('Rejected: track larger than the frame',t['track_larger_than_frame_events']),('Rejected: frame not inside the detector',t['frame_outside_detector_events']),
          ('Rejected: gamma deposit above the cut outside the frame',t['outside_frame_gamma_events']),('Accepted',t['accepted_events'])]
    selection=''.join(f'<tr><td>{a}</td><td>{b}</td><td>{100*b/n0:.1f}%</td></tr>' for a,b in rows)
    # Events remaining after each step; the E2 thresholds come from the comparison of the two samples
    fid=n0-t['outside_detector_events']-t.get('main_track_near_el_events',0)-t['track_larger_than_frame_events']-t['frame_outside_detector_events']
    common=([('Generated',generated)] if generated else [])+[('Energy deposited in the gas inside the Geant4 window',n0),
        (f"Fiducial: track inside the detector, {float(m['settings'].get('el_exclusion_mm',0)):g} mm from the EL, not larger than the frame; frame inside the detector",fid),
        ('Gamma veto: no gamma deposit above the cut outside the frame',t['accepted_events'])]
    threshold={};roi=None
    if comparison:
        with h5py.File(Path(comparison)/'comparison.h5') as c:
            if 'roi' in c: roi=(float(c['roi/lower_keV'][()]),float(c['roi/upper_keV'][()]))
            for arm in ('itrack','etrack'):
                w=c[f"E2_cut/{arm}/{'selected' if roi else 'one_component'}/working_points"];threshold[arm]=float(w['threshold_keV'][int(np.argmin(np.abs(w['target_signal_efficiency'][:]-0.8)))])
    def flow(arm):
        # With a comparison that defines the ROI: ROI cut, then the track selection of the track
        # (the lattice selection for the itrack when the table has it, one component otherwise)
        keep=np.ones(n,bool);steps=list(common)
        if roi: keep=(e['electron_energy_keV']>roi[0])&(e['electron_energy_keV']<roi[1]);steps.append((f'ROI cut: electron energy within {roi[0]:.1f}–{roi[1]:.1f} keV',int(keep.sum())))
        if arm=='itrack' and 'itrack_selected' in e: one=keep&(e['itrack_selected']==1);steps.append(('Track selection (itrack): no detached deposit; one segment, or up to three joined',int(one.sum())))
        else: one=keep&(e[arm+'_components']==1);steps.append((f'One component in the {arm}',int(one.sum())))
        if arm in threshold: steps.append((f'E₂ > {threshold[arm]:.0f} keV ({arm})',int((one&(e[arm+"_E2_keV"]>threshold[arm])).sum())))
        body=''.join(f'<tr><td>{a}</td><td>{v}</td><td>{"" if j==0 else f"{100*v/steps[j-1][1]:.1f}%"}</td><td>{100*v/steps[0][1]:.2f}%</td></tr>' for j,(a,v) in enumerate(steps))
        return f'<h3>{arm}</h3><table><tr><th>Step</th><th>Events</th><th>Of previous step</th><th>Of first step</th></tr>{body}</table>'
    flow_note='The first steps use the Monte Carlo truth and are common to both tracks.'+(' The ROI is one FWHM of the bb0nu electron-energy peak around Qbb, and the E₂ threshold keeps 80% of the bb0nu events that reach that step (<a href="'+html.escape(os.path.relpath(Path(comparison)/'index.html',out))+'">comparison of the two samples</a>).' if comparison else '')
    cutflow=f'<h2>Selection</h2><p>Events remaining after each step, applied in this order. {flow_note}</p>{flow("itrack")}{flow("etrack")}'
    b=m['batch'];per_batch=''.join(f"<tr><td>{tag}</td><td>{b['input_events'][i]}</td><td>{b['accepted_events'][i]}</td><td>{100*b['accepted_events'][i]/b['input_events'][i]:.1f}%</td></tr>" for i,tag in enumerate(m['batches']))
    fi,fe=m['fits']['ion_energy_keV'],m['fits']['electron_energy_keV']
    def line(arm):
        c=e[arm+'_components'];return (f"<tr><td>{arm}</td><td>{np.nanmedian(e[arm+'_length_mm']):.0f}</td><td>{np.nanmedian(e[arm+'_extent_r_mm']):.0f}</td><td>{np.nanmedian(e[arm+'_extent_z_mm']):.0f}</td>"
            f"<td>{np.median(e[arm+'_hits']):.0f}</td><td>{100*np.mean(c==1):.1f}%</td><td>{np.nanmean(e[arm+'_E1_keV']):.0f}</td><td>{np.nanmean(e[arm+'_E2_keV']):.0f}</td></tr>")
    truth=''
    if classification and Path(classification).exists():
        c=tomllib.loads(Path(classification).read_text());body=''.join(f"<tr><td>{r['track']}</td><td>{r['link_mm']:g}</td><td>{r['events']}</td><td>{r['one_component']} ({100*r['one_component']/r['events']:.1f}%)</td><td>{r['floating_gammas']} ({100*r['floating_gammas']/r['events']:.1f}%)</td><td>{r['broken']} ({100*r['broken']/r['events']:.1f}%)</td></tr>" for r in c['rows'] if r['sample']==m['species'])
        truth=(f"<h2>Detached components against the Monte Carlo truth</h2><p>From the {html.escape(c['sample_description'])}, with the same selection and settings. {html.escape(c['rule'])}</p>"
               f"<table><tr><th>Track</th><th>Connection radius (mm)</th><th>Events</th><th>One component</th><th>Only floating gammas</th><th>Broken track</th></tr>{body}</table>")
    viewer=Path(run_dir)/'batches'/m['batches'][0]/out.parent.name/'presentation'/'index.html'
    viewer_link=f'<p>Event viewer: <a href="../../batches/{m["batches"][0]}/{out.parent.name}/presentation/index.html#viewer">events of batch {m["batches"][0]}</a>.</p>' if viewer.exists() else ''
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ITACA · {sample} · statistics</title>
<style>body{{font:17px/1.5 system-ui;color:#203343;max-width:1500px;margin:auto;padding:24px}}img{{max-width:100%}}table{{border-collapse:collapse;margin:10px 0 18px}}th,td{{padding:7px 14px;border-bottom:1px solid #ddd;text-align:left}}th{{background:#eef4f7}}a{{color:#185b83}}h2{{margin-top:34px}}code{{background:#eef3f6;padding:2px 5px;border-radius:3px}}</style></head>
<body><a href="../../../../../notes/html/index.html">ITACA index</a><h1>{sample} · statistics of the accepted events</h1>
<p>Run <code>{html.escape(m['run_label'])}</code>: {len(m['batches'])} batches, {n} accepted events. Frame of {side:g} mm; itrack hits above {tc['ion_hit_cut_keV']:g} keV, voxels of {tc['ion_voxel_width_mm']:g} mm, connection radius {tc['ion_connection_mm']:g} mm; etrack hits above {tc['electron_hit_cut_keV']:g} keV, voxels of {tc['electron_voxel_width_mm']:g} mm, connection radius {tc['electron_connection_mm']:g} mm. Saved hit cuts: {m['ion_cut']:g} keV (ions), {m['electron_cut']:g} keV (electrons).</p>{viewer_link}
{cutflow}
<h3>Rejected events by cause</h3><table><tr><th>Step</th><th>Events</th><th>Fraction of events read</th></tr>{selection}</table>
<details><summary>Counts per batch</summary><table><tr><th>Batch</th><th>Events read</th><th>Accepted</th><th>Fraction</th></tr>{per_batch}</table></details>
<img src="plots/event_position.png" alt="Vertex position and frame shift">
<h2>Event energy before the hit cut</h2><p>Sum of all signed reconstructed charges. Ions: μ = {fi['mu']:.1f} keV, σ/μ = {100*fi['sigma']/fi['mu']:.2f}% ({235.5*fi['sigma']/fi['mu']:.2f}% FWHM). Electrons: μ = {fe['mu']:.1f} keV, σ/μ = {100*fe['sigma']/fe['mu']:.3f}% ({235.5*fe['sigma']/fe['mu']:.2f}% FWHM). Incomplete electron estimates: {int(np.count_nonzero(~np.isfinite(e['electron_energy_keV'])))}; incomplete ion estimates: {int(np.count_nonzero(~np.isfinite(e['ion_energy_keV'])))}.</p>
<img src="plots/event_energy.png" alt="Ion and electron energy"><img src="plots/energy_by_drift.png" alt="Energy resolution against drift">
<h2>Main track</h2><p>The main track is the connected component with most voxels. Length: sum of the distances between consecutive voxels of the tree path joining its two endpoints. Extension in R: the larger of its X and Y ranges; extension in Z: its Z range; each range includes one voxel width.</p>
<table><tr><th>Track</th><th>Length, median (mm)</th><th>Extension in R, median (mm)</th><th>Extension in Z, median (mm)</th><th>Hits, median</th><th>One component</th><th>Mean E₁ (keV)</th><th>Mean E₂ (keV)</th></tr>{line('itrack')}{line('etrack')}</table>
<img src="plots/track_size.png" alt="Track length and extensions"><img src="plots/hits_components.png" alt="Hits and components">{truth}
<h2>Endpoint energies</h2><p>E₁ and E₂ sum the voxel energies in spheres around the two endpoints of the main track, with E₁ ≥ E₂.</p>
<img src="plots/blob_energies.png" alt="Endpoint energies"><img src="plots/blob_energy_by_drift.png" alt="Endpoint energies against drift">
<p>Numerical results: <a href="events.h5">events.h5</a> (one row per accepted event).</p></body></html>'''
    (out/'index.html').write_text(page);print(out/'index.html')

if __name__=='__main__':
    options={k:v for k,v in (a[2:].split('=',1) for a in sys.argv[1:] if a.startswith('--'))};args=[a for a in sys.argv[1:] if not a.startswith('--')]
    if len(args) not in (1,2) or set(options)-{'generated','comparison'}: raise SystemExit(__doc__)
    render(args[0],args[1] if len(args)>1 else None,int(options['generated']) if 'generated' in options else None,options.get('comparison'))
