#!/usr/bin/env python3
"""Comparison page of two samples: distributions overlaid and the E2 cut.

  python scripts/production/present_species_comparison.py OUTPUT_DIR [GENERATED_SIGNAL GENERATED_BACKGROUND]

Reads OUTPUT_DIR/comparison.h5 (prepare_species_comparison.jl) and the two
statistics files it names; writes plots and index.html in OUTPUT_DIR.
Histograms are normalised to unit area so the shapes can be compared.
With the two numbers of generated events the selection tables start from them.
"""
import html,sys,tomllib
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SIG,BKG='#287ba0','#b16d31'
LABEL={'bb0nu':'bb0nu','single_electrons':'single electrons','bi214':'²¹⁴Bi gammas'}
def text(v): return v.decode() if isinstance(v,bytes) else v

def events(path):
    with h5py.File(path) as f:
        e={k:f['events'][k][:] for k in f['events']}
        e['selection']={k:int(v[()]) for k,v in f['selection/total'].items()}
    for arm in ('itrack','etrack'): e[arm+'_extent_r_mm']=np.maximum(e[arm+'_extent_x_mm'],e[arm+'_extent_y_mm'])
    return e

def overlay(ax,s,b,edges,xlabel,title,names,unit=''):
    for v,color,name in ((s,SIG,names[0]),(b,BKG,names[1])):
        v=np.asarray(v,float);v=v[np.isfinite(v)]
        ax.hist(v,bins=edges,density=True,histtype='stepfilled',alpha=.35,color=color)
        ax.hist(v,bins=edges,density=True,histtype='step',color=color,lw=1.4,label=f'{name}: N = {len(v)}, mean {np.mean(v):.4g}, std {np.std(v):.4g}, median {np.median(v):.4g} {unit}')
    ax.set(xlabel=xlabel,ylabel='Fraction of events / unit',title=title,xlim=(edges[0],edges[-1]));ax.grid(alpha=.15);ax.legend(fontsize=8,loc='upper right');ax.margins(y=.3)

def render(out,generated=None):
    out=Path(out);plots=out/'plots';plots.mkdir(exist_ok=True)
    with h5py.File(out/'comparison.h5') as f:
        meta={side:{k:text(f[side][k][()]) for k in f[side]} for side in ('signal','background')};tc=tomllib.loads(text(f['track_analysis'][()]))
        cuts={arm:{sel:{k:(v[()] if isinstance(v,h5py.Dataset) else {n:w[:] for n,w in v.items()}) for k,v in f[f'E2_cut/{arm}/{sel}'].items()} for sel in ('all','one_component')} for arm in ('itrack','etrack')}
    s=events(out/meta['signal']['file']);b=events(out/meta['background']['file'])
    names=(LABEL.get(meta['signal']['species'],meta['signal']['species']),LABEL.get(meta['background']['species'],meta['background']['species']))
    title=f'{names[0]} ({len(s["event_number"])} events) and {names[1]} ({len(b["event_number"])} events) · unit area'

    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    for row,arm in ((0,'itrack'),(1,'etrack')):
        overlay(axes[row,0],s[arm+'_length_mm'],b[arm+'_length_mm'],np.linspace(0,500,51),'Track length (mm)',f'{arm} · path between the endpoints',names,'mm')
        overlay(axes[row,1],s[arm+'_extent_r_mm'],b[arm+'_extent_r_mm'],np.linspace(0,160,41),'Extension in R (mm)',f'{arm} · largest of the X and Y ranges',names,'mm')
        overlay(axes[row,2],s[arm+'_extent_z_mm'],b[arm+'_extent_z_mm'],np.linspace(0,160,41),'Extension in Z (mm)',f'{arm} · Z range',names,'mm')
    fig.suptitle('Main track · '+title);fig.savefig(plots/'track_size.png',dpi=140);plt.close(fig)

    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    for row,arm in ((0,'itrack'),(1,'etrack')):
        overlay(axes[row,0],s[arm+'_E1_keV'],b[arm+'_E1_keV'],np.linspace(0,1000,81),'E₁ (keV)',f'{arm} · larger endpoint energy',names,'keV')
        overlay(axes[row,1],s[arm+'_E2_keV'],b[arm+'_E2_keV'],np.linspace(0,800,81),'E₂ (keV)',f'{arm} · smaller endpoint energy',names,'keV')
        overlay(axes[row,2],s[arm+'_hits'],b[arm+'_hits'],np.linspace(0,max(400,np.percentile(np.r_[s[arm+'_hits'],b[arm+'_hits']],99.5)),51),'Hits per event',f'{arm} · hits above the display cut',names)
    fig.suptitle('Endpoint energies and hits · '+title);fig.savefig(plots/'endpoint_energies.png',dpi=140);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    overlay(axes[0],s['ion_energy_keV'],b['ion_energy_keV'],np.linspace(2330,2580,81),'Energy (keV)','Ion energy estimator',names,'keV')
    overlay(axes[1],s['electron_energy_keV'],b['electron_energy_keV'],np.linspace(2430,2490,81),'Energy (keV)','Electron energy estimator',names,'keV')
    fig.suptitle('Event energy before the hit cut · '+title);fig.savefig(plots/'event_energy.png',dpi=140);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(13,5.6),layout='constrained')
    for ax,sel,label in ((axes[0],'all','all accepted events'),(axes[1],'one_component','events with one component in that track')):
        for arm,color in (('itrack',SIG),('etrack',BKG)):
            c=cuts[arm][sel];w=c['working_points'];i=int(np.argmin(np.abs(w['target_signal_efficiency']-0.8)))
            ax.plot(c['background_acceptance'],c['signal_efficiency'],color=color,lw=2,
                    label=f"{arm}: area {float(c['area']):.3f}; at {100*w['signal_efficiency'][i]:.0f}% signal, background {100*w['background_acceptance'][i]:.1f} ± {100*w['background_acceptance_error'][i]:.1f}% (E₂ > {w['threshold_keV'][i]:.0f} keV)")
            ax.plot(w['background_acceptance'],w['signal_efficiency'],'o',color=color,ms=4)
        ax.plot([0,1],[0,1],color='#526776',lw=.8,ls='--');ax.set(xlabel=f'Fraction of {names[1]} accepted',ylabel=f'Fraction of {names[0]} accepted',title=f'Cut E₂ > threshold · {label}',xlim=(0,1),ylim=(0,1));ax.grid(alpha=.15);ax.legend(fontsize=8,loc='lower right')
    fig.suptitle(f'Classical separation by the smaller endpoint energy · {names[0]} against {names[1]}');fig.savefig(plots/'e2_cut.png',dpi=140);plt.close(fig)

    def table(sel):
        rows=''
        for arm in ('itrack','etrack'):
            c=cuts[arm][sel];w=c['working_points']
            for k in range(len(w['target_signal_efficiency'])):
                rows+=(f"<tr><td>{arm}</td><td>{100*w['signal_efficiency'][k]:.1f}%</td><td>{w['threshold_keV'][k]:.0f}</td><td>{100*w['background_acceptance'][k]:.1f} ± {100*w['background_acceptance_error'][k]:.1f}%</td>"
                       f"<td>{w['signal_efficiency'][k]/max(w['background_acceptance'][k],1e-12):.2f}</td></tr>")
        return f"<table><tr><th>Track</th><th>{names[0]} accepted</th><th>E₂ threshold (keV)</th><th>{names[1]} accepted</th><th>Ratio</th></tr>{rows}</table>"
    def flow(arm,target=0.8):
        w=cuts[arm]['one_component']['working_points'];i=int(np.argmin(np.abs(w['target_signal_efficiency']-target)));thr=float(w['threshold_keV'][i])
        cols=[]
        for k,e in enumerate((s,b)):
            t=e['selection'];n=t['input_events'];fid=n-t['outside_detector_events']-t['main_track_near_el_events']-t['track_larger_than_frame_events']-t['frame_outside_detector_events']
            one=e[arm+'_components']==1
            cols.append(([generated[k]] if generated else [])+[n,fid,fid-t['outside_frame_gamma_events'],int(one.sum()),int((one&(e[arm+'_E2_keV']>thr)).sum())])
        steps=(['Generated'] if generated else [])+['Energy deposited in the gas inside the Geant4 window','Fiducial: track inside the detector, 50 mm from the EL, not larger than the frame; frame inside the detector',
            'Gamma veto: no gamma deposit above the cut outside the frame',f'One component in the {arm}',f'E₂ > {thr:.0f} keV ({arm})']
        def cell(c,j): return f'<td>{c[j]}</td><td>{"" if j==0 else f"{100*c[j]/c[j-1]:.1f}%"}</td><td>{100*c[j]/c[0]:.2f}%</td>'
        rows=''.join(f'<tr><td>{step}</td>{cell(cols[0],j)}{cell(cols[1],j)}</tr>' for j,step in enumerate(steps))
        return (f'<h3>{arm}</h3><table><tr><th rowspan="2">Step</th><th colspan="3">{names[0]}</th><th colspan="3">{names[1]}</th></tr>'
                f'<tr><th>Events</th><th>Of previous step</th><th>Of first step</th><th>Events</th><th>Of previous step</th><th>Of first step</th></tr>{rows}</table>')
    def frac(arm): c=cuts[arm]['one_component'];return f"{arm}: {100*float(c['signal_selected_fraction']):.1f}% of {names[0]} and {100*float(c['background_selected_fraction']):.1f}% of {names[1]}"
    areas='; '.join(f"{arm} {float(cuts[arm]['all']['area']):.3f}" for arm in ('itrack','etrack'))
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ITACA · {names[0]} and {names[1]}</title>
<style>body{{font:17px/1.5 system-ui;color:#203343;max-width:1500px;margin:auto;padding:24px}}img{{max-width:100%}}table{{border-collapse:collapse;margin:10px 0 18px}}th,td{{padding:7px 14px;border-bottom:1px solid #ddd;text-align:left}}th{{background:#eef4f7}}a{{color:#185b83}}h2{{margin-top:34px}}code{{background:#eef3f6;padding:2px 5px;border-radius:3px}}</style></head>
<body><a href="../../../notes/html/index.html">ITACA index</a><h1>{names[0]} and {names[1]} · comparison</h1>
<p>Runs <code>{html.escape(meta['signal']['run_label'])}</code> ({len(s['event_number'])} accepted events) and <code>{html.escape(meta['background']['run_label'])}</code> ({len(b['event_number'])} accepted events), with the same selection and track settings. Reports: <a href="{html.escape(str(Path(meta['signal']['file']).parent/'index.html'))}">{names[0]}</a> · <a href="{html.escape(str(Path(meta['background']['file']).parent/'index.html'))}">{names[1]}</a>.</p>
<h2>Selection</h2><p>Events remaining after each step, applied in this order. The first three steps use the Monte Carlo truth and are common to both tracks. The E₂ threshold keeps 80% of the {names[0]} events with one component.</p>{flow('itrack')}{flow('etrack')}
<h2>Classical separation by the smaller endpoint energy</h2><p>An event is kept when E₂ of the main track exceeds a threshold. The curves give the fraction of each sample kept as the threshold varies; the area under the curve is 0.5 for no separation and 1 for complete separation. Areas with all accepted events: {areas}. Errors are binomial.</p>
<img src="plots/e2_cut.png" alt="Signal efficiency against background acceptance for a cut on E2">
<h3>All accepted events</h3>{table('all')}
<h3>Events with one component in that track</h3><p>Fraction of the accepted events with one component: {frac('itrack')}; {frac('etrack')}. The efficiencies below are relative to those events.</p>{table('one_component')}
<h2>Endpoint energies and hits</h2><img src="plots/endpoint_energies.png" alt="Endpoint energies and hits">
<h2>Main track</h2><img src="plots/track_size.png" alt="Track length and extensions">
<h2>Event energy</h2><img src="plots/event_energy.png" alt="Event energy">
<p>Numerical results: <a href="comparison.h5">comparison.h5</a>.</p></body></html>'''
    (out/'index.html').write_text(page);print(out/'index.html')

if __name__=='__main__':
    if len(sys.argv) not in (2,4): raise SystemExit(__doc__)
    render(sys.argv[1],tuple(int(v) for v in sys.argv[2:4]) if len(sys.argv)==4 else None)
