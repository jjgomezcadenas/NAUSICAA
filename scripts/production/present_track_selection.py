#!/usr/bin/env python3
"""Page of the track selection on the readout lattice.

  python scripts/production/present_track_selection.py OUTPUT_DIR SIGNAL_RUN BACKGROUND_RUN... [--truth=RUN_LABEL:TABLE.txt]

Each RUN is results/readout/<run_label> with the run statistics (events.h5),
the selection (in the feature tables, or from prepare_track_selection.jl)
and, optionally, the components of prepare_second_pass.jl. The first run is the signal. --truth gives, for a
gamma sample, the table of xray_truth_distance.py. The ROI is one FWHM of the
signal electron-energy peak around Qbb. Writes plots and index.html.
"""
import html,sys,tomllib
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2];LABEL={'bb0nu':'bb0nu','single_electrons':'single electrons','bi214':'²¹⁴Bi gammas'};COLORS=['#287ba0','#b16d31','#3d8f5a','#7a5195']
def text(v): return v.decode() if isinstance(v,bytes) else v

def load(run):
    run=Path(run);sub=next(p.parent.parent for p in run.glob('*/statistics/events.h5'))
    with h5py.File(sub/'statistics/events.h5') as f:
        e={k:f['events'][k][:] for k in f['events']};species=text(f['species'][()]);label=text(f['run_label'][()]);sigma=float(f['energy_fits/electron_energy_keV/sigma'][()])
    defaults=dict(ion_segment_connection_mm=5.5,ion_segment_min_hits=6,ion_deposit_min_keV=10.0,ion_max_segments=3)
    if 'itrack_selected' in e:       # feature tables with the selection
        with h5py.File(sub/'statistics/events.h5') as f: tc=tomllib.loads(text(f['track_analysis'][()]))
        s={'segments':e['itrack_segments'],'deposits':e['itrack_deposits'],'deposit_keV':e['itrack_deposit_keV'],'segments_joined':e['itrack_segments_joined'],'selected':e['itrack_selected']}
        settings={k:tc.get(k,v) for k,v in defaults.items()};settings['ion_connection_mm']=tc['ion_connection_mm'];settings['ion_hit_cut_keV']=tc['ion_hit_cut_keV']
    else:
        with h5py.File(sub/'track_selection/events.h5') as f:
            assert np.array_equal(f['event_number'][:],e['event_number']);s={k:f[k][:] for k in ('segments','deposits','deposit_keV','segments_joined','selected')};settings={k:f['settings'][k][()] for k in f['settings']}
    comp=None
    if (sub/'second_pass/components_lattice.h5').exists():
        with h5py.File(sub/'second_pass/components_lattice.h5') as f:
            assert np.array_equal(f['event_number'][:],e['event_number']);radii=list(f['radii_mm'][:]);k=radii.index(float(settings['ion_segment_connection_mm']))+1
            comp={n:f[f'radius_{k}'][n][:] for n in ('offsets','voxels','energy_keV')}
    return dict(dir=run,events=e,sel=s,settings=settings,comp=comp,species=species,name=LABEL.get(species,species),run_label=run.name,sigma=sigma)

def render(out,runs,truth):
    out=Path(out);plots=out/'plots';plots.mkdir(parents=True,exist_ok=True);R=[load(r) for r in runs];sig=R[0]
    qbb=tomllib.loads((ROOT/'config/detector/energy.toml').read_text())['qbb_keV'];fwhm=2.3548*sig['sigma'];st=sig['settings']
    rows=[];flow={}
    for r in R:
        e=r['events'];s=r['sel'];n=len(e['event_number']);roi=np.abs(e['electron_energy_keV']-qbb)<fwhm/2
        old=(e['itrack_components']==1)&(e['etrack_components']==1);new=s['selected']==1;dep=s['deposits']>0
        apart=~dep&(s['segments']>1)&(s['segments_joined']==0);many=~dep&(s['segments']>int(st['ion_max_segments']))&(s['segments_joined']==1);none=~dep&(s['segments']==0)
        r.update(roi=roi,old=old,new=new);flow[r['name']]=dict(accepted=n,roi=int(roi.sum()),old=int((roi&old).sum()),new=int((roi&new).sum()),dep=int((roi&dep).sum()),apart=int((roi&apart).sum()),many=int((roi&(many|none)).sum()),
            one=int((roi&new&(s['segments']==1)).sum()),joined=int((roi&new&(s['segments']>1)).sum()))
    def cells(key,ref=None): return ''.join(f"<td>{flow[r['name']][key]}</td><td>{'' if ref is None else f'{100*flow[r['name']][key]/flow[r['name']][ref]:.1f}%'}</td>" for r in R)
    head=''.join(f'<th colspan="2">{r["name"]}</th>' for r in R);sub=''.join('<th>Events</th><th>Of the ROI</th>' for r in R)
    table=(f'<table><tr><th rowspan="2">Step</th>{head}</tr><tr>{sub}</tr>'
           f'<tr><td>Accepted: fiducial selection and gamma veto</td>{cells("accepted")}</tr>'
           f'<tr><td>ROI: electron energy within {qbb-fwhm/2:.1f}–{qbb+fwhm/2:.1f} keV</td>{cells("roi","roi")}</tr>'
           f'<tr><td>Rejected: a detached deposit</td>{cells("dep","roi")}</tr><tr><td>Rejected: segments that stay apart</td>{cells("apart","roi")}</tr><tr><td>Rejected: no segment, or more than {int(st["ion_max_segments"])}</td>{cells("many","roi")}</tr>'
           f'<tr><td><strong>Selected</strong></td>{cells("new","roi")}</tr><tr><td>&nbsp;&nbsp;one segment</td>{cells("one","roi")}</tr><tr><td>&nbsp;&nbsp;two or three segments joined at {float(st["ion_connection_mm"]):g} mm</td>{cells("joined","roi")}</tr>'
           f'<tr><td>For comparison: one component in each track, {float(st["ion_connection_mm"]):g} mm on 3 mm cells</td>{cells("old","roi")}</tr></table>')
    fs=flow[sig['name']];gain=''
    for r in R[1:]:
        fb=flow[r['name']];g_roi=(fs['roi']/fs['accepted'])/np.sqrt(fb['roi']/fb['accepted']);g_new=(fs['new']/fs['roi'])/np.sqrt(fb['new']/fb['roi']);g_old=(fs['old']/fs['roi'])/np.sqrt(fb['old']/fb['roi'])
        gain+=f"<tr><td>{r['name']}</td><td>{g_roi:.2f}</td><td>{g_old:.2f}</td><td>{g_new:.2f}</td><td>{g_new/g_old:.3f}</td></tr>"
    gain=f'<table><tr><th>Background</th><th>ROI cut</th><th>One component at 8 mm</th><th>Selection on the lattice</th><th>Ratio of the last two</th></tr>{gain}</table>'

    # Components beyond the largest, by hits and energy
    C=[r for r in R if r['comp'] is not None]
    fig,axes=plt.subplots(1,max(len(C),1),figsize=(6.2*max(len(C),1),5),layout='constrained',squeeze=False);hb=np.arange(.5,21.5,1);eb=np.logspace(np.log10(5),np.log10(1000),40)
    for ax,r in zip(axes[0],C):
        c=r['comp'];off=c['offsets'];first=np.zeros(len(c['voxels']),bool);first[off[:-1][np.diff(off)>0]]=True;ev=np.repeat(np.arange(len(off)-1),np.diff(off));m=~first&r['roi'][ev]
        h=ax.hist2d(np.clip(c['voxels'][m],1,20),c['energy_keV'][m],bins=[hb,eb],cmin=1,cmap='viridis');fig.colorbar(h[3],ax=ax,label='Components')
        ax.axvline(float(st['ion_segment_min_hits'])-.5,color='k',ls='--',lw=1);ax.plot([.5,float(st['ion_segment_min_hits'])-.5],[float(st['ion_deposit_min_keV'])]*2,color='k',ls=':',lw=1)
        ax.set_yscale('log');ax.set(xlabel='Hits in the component (20: 20 or more)',ylabel='Energy of the component (keV)',title=f"{r['name']} · {int(m.sum())} components in {int(r['roi'].sum())} ROI events")
    fig.suptitle(f"Components beyond the largest, connected at {float(st['ion_segment_connection_mm']):g} mm · dashed: segments to the right; dotted: detached deposits above");fig.savefig(plots/'components.png',dpi=140);plt.close(fig)

    fig,ax=plt.subplots(figsize=(9,5),layout='constrained');bins=np.arange(10,202,4)
    for r,c in zip(R,COLORS):
        v=r['sel']['deposit_keV'][r['roi']&(r['sel']['deposits']>0)];w=np.full(len(v),100/max(r['roi'].sum(),1))
        ax.hist(np.clip(v,0,199.9),bins=bins,weights=w,histtype='step',color=c,lw=1.6,label=f"{r['name']}: {len(v)} of {int(r['roi'].sum())} ROI events ({100*len(v)/max(r['roi'].sum(),1):.1f}%), median {np.median(v):.1f} keV")
    ax.set(xlabel='Energy of the largest detached deposit (keV); last bin: 198 keV or more',ylabel='ROI events with a detached deposit (% of the ROI events) / 4 keV',title='Largest detached deposit of the events that have one');ax.legend(fontsize=8);ax.grid(alpha=.15)
    fig.savefig(plots/'deposit_energy.png',dpi=140);plt.close(fig)

    truth_html=''
    for r in R:
        if r['run_label'] not in truth: continue
        t=np.loadtxt(truth[r['run_label']]);row={int(a):i for i,a in enumerate(t[:,0])};i=np.array([row.get(int(a),-1) for a in r['events']['event_number']])
        if (i<0).any(): raise SystemExit(f"Events of {r['run_label']} missing in the truth table")
        d=t[i,1];code=t[i,7].astype(int);x=(code==0)&np.isfinite(d);new=r['new'];old=r['old'];dep=r['sel']['deposits']>0
        names={0:'Photoelectric absorption with a K X-ray',1:'Photoelectric absorption without',2:'Compton scattering',3:'Pair production'}
        comp=''.join(f"<tr><td>{names[k]}</td><td>{int((code==k).sum())}</td><td>{int(((code==k)&old).sum())} ({100*((code==k)&old).sum()/max((code==k).sum(),1):.1f}%)</td><td>{int(((code==k)&new).sum())} ({100*((code==k)&new).sum()/max((code==k).sum(),1):.1f}%)</td></tr>" for k in names if (code==k).any())
        edges=np.array([0,2,4,6,8,10,12,16,24,60]);mid=(edges[:-1]+edges[1:])/2;fig,axes=plt.subplots(1,2,figsize=(15,5.2),layout='constrained');bins=np.arange(0,62,2)
        axes[0].hist(np.clip(d[x],0,59.9),bins=bins,histtype='step',color='#526776',lw=1.5,label=f'all: N = {int(x.sum())}, median {np.median(d[x]):.1f} mm, mean {d[x].mean():.1f} mm, std {d[x].std():.1f} mm')
        axes[0].hist(np.clip(d[x&~old],0,59.9),bins=bins,histtype='step',color=COLORS[1],lw=1.5,label=f'rejected by one component at 8 mm: N = {int((x&~old).sum())}')
        axes[0].hist(np.clip(d[x&~new],0,59.9),bins=bins,color=COLORS[0],alpha=.45,label=f'rejected by the selection on the lattice: N = {int((x&~new).sum())}')
        axes[0].set(xlabel='Distance from the X-ray deposit to the nearest main-track hit, Monte Carlo (mm); last bin: 58 mm or more',ylabel='Events / 2 mm',title='K X-ray of the photoelectric absorption');axes[0].legend(fontsize=8)
        for sel,c,name in ((~old,COLORS[1],'one component at 8 mm'),(~new,COLORS[0],'selection on the lattice')):
            f=[];er=[]
            for a,b in zip(edges[:-1],edges[1:]):
                m=x&(d>=a)&(d<b);p=sel[m].mean() if m.sum() else np.nan;f.append(100*p);er.append(100*np.sqrt(p*(1-p)/m.sum()) if m.sum() else 0)
            axes[1].errorbar(mid,f,er,marker='o',color=c,capsize=3,label=name)
        axes[1].set(xlabel='Distance to the nearest main-track hit (mm)',ylabel='Events rejected (%)',title='Fraction of the events with a K X-ray that is rejected',ylim=(0,105));axes[1].legend();axes[1].grid(alpha=.15)
        fig.suptitle(f"{r['name']} · {int(x.sum())} accepted photoelectric events with a K X-ray");fig.savefig(plots/f"xray_{r['species']}.png",dpi=140);plt.close(fig)
        truth_html+=(f"<h2>{r['name']}: against the Monte Carlo truth</h2><p>All {len(d)} accepted events, before the ROI cut. The main electron comes from the first interaction listed; the percentages are the fraction of each class that passes.</p>"
            f"<table><tr><th>Interaction that produced the main electron</th><th>Accepted</th><th>One component at 8 mm</th><th>Selection on the lattice</th></tr>{comp}</table>"
            f"<p>For the photoelectric events with a K X-ray, the distance from the X-ray deposit to the nearest hit of the main track decides whether it can be seen. Within about 4 mm it shares cells with the track.</p><img src=\"plots/xray_{r['species']}.png\" alt=\"K X-ray distance and rejection\">")

    back='../'*len(out.resolve().relative_to(ROOT).parts)
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ITACA · track selection</title>
<style>body{{font:17px/1.5 system-ui;color:#203343;max-width:1500px;margin:auto;padding:24px}}img{{max-width:100%}}table{{border-collapse:collapse;margin:10px 0 18px}}th,td{{padding:7px 14px;border-bottom:1px solid #ddd;text-align:left}}th{{background:#eef4f7}}a{{color:#185b83}}h2{{margin-top:34px}}code{{background:#eef3f6;padding:2px 5px;border-radius:3px}}</style></head>
<body><a href="{back}notes/html/index.html">ITACA index</a><h1>Track selection on the readout lattice</h1>
<p>Runs {', '.join(f"<code>{html.escape(r['run_label'])}</code>" for r in R)}. The selection and its reasons: <a href="{back}notes/html/event_analysis/track_selection.html">track selection</a>.</p>
<h2>Selection</h2><p>Itrack hits above {float(st['ion_hit_cut_keV']):g} keV are connected at {float(st['ion_segment_connection_mm']):g} mm. A component with at least {int(st['ion_segment_min_hits'])} hits is a segment; a smaller one with at least {float(st['ion_deposit_min_keV']):g} keV is a detached deposit. An event is selected with no detached deposit and one segment, or up to {int(st['ion_max_segments'])} segments that join at {float(st['ion_connection_mm']):g} mm. The ROI is one FWHM ({fwhm:.1f} keV) of the {sig['name']} electron-energy peak around Qbb = {qbb:g} keV.</p>{table}
<h3>Gain in signal over the square root of background</h3><p>Each number is the signal efficiency of the step divided by the square root of the background efficiency of the same step.</p>{gain}
<h2>What the components are</h2>{'<img src="plots/components.png" alt="Hits and energy of the components beyond the largest">' if C else ''}<img src="plots/deposit_energy.png" alt="Energy of the largest detached deposit">
{truth_html}</body></html>'''
    (out/'index.html').write_text(page);print(out/'index.html')

if __name__=='__main__':
    args=[a for a in sys.argv[1:] if not a.startswith('--')];truth=dict(a[8:].split(':',1) for a in sys.argv[1:] if a.startswith('--truth='))
    if len(args)<3: raise SystemExit(__doc__)
    render(args[0],args[1:],truth)
