#!/usr/bin/env python3
"""Page of a round of network trainings.

  python scripts/network/present_network.py ROUND_DIR

ROUND_DIR holds one folder per training with the training.h5 of
train_network.py. Writes plots and index.html in ROUND_DIR. The classical cut
on E2 is evaluated on the test events of the first training, which all
trainings of a round share. No training is repeated.
"""
import html,sys
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).parent))
from itaca_network import separation

COLOR={'itrack':'#287ba0','etrack':'#b16d31','both':'#3d8f5a'};TARGETS=(0.9,0.8,0.7,0.6)
def text(v): return v.decode() if isinstance(v,bytes) else v

def load(path):
    with h5py.File(path) as f:
        r=dict(name=path.parent.name,input=text(f['input'][()]),settings={k:(text(v[()]) if isinstance(v,h5py.Dataset) else None) for k,v in f['settings'].items()},
               history={k:v[:] for k,v in f['history'].items()},test={k:v[()] for k,v in f['test'].items() if isinstance(v,h5py.Dataset)},
               species=tuple(text(f[k][()]).replace('_',' ') for k in ('signal_species','background_species')),runs=(text(f['signal_run_label'][()]),text(f['background_run_label'][()])))
    r['label']=r['test']['label']>0;r['separation']=separation(r['test']['score'],r['label']);return r

def curve(score,label):
    """Signal efficiency and background acceptance for every threshold, from the highest score down."""
    order=np.argsort(-score,kind='stable');l=label[order];return np.cumsum(l)/l.sum(),np.cumsum(~l)/(~l).sum()

def at(sep,target): return int(np.argmin(np.abs(sep['target_signal_efficiency']-target)))

def by_drift(score,label,drift,edges,target=0.8):
    """Background acceptance in drift intervals, with one threshold set on all signal events."""
    thr=np.quantile(score[label],1-target);acc=[];err=[];eff=[]
    for a,b in zip(edges[:-1],edges[1:]):
        m=(drift>=a)&(drift<b);bk=score[m&~label]>thr;sg=score[m&label]>thr
        acc.append(bk.mean() if len(bk) else np.nan);err.append(np.sqrt(bk.mean()*(1-bk.mean())/len(bk)) if len(bk) else np.nan);eff.append(sg.mean() if len(sg) else np.nan)
    return np.array(acc),np.array(err),np.array(eff)

def render(root):
    root=Path(root);runs=[load(p) for p in sorted(root.glob('*/training.h5'))]
    if not runs: raise SystemExit(f'No training.h5 under {root}')
    plots=root/'plots';plots.mkdir(exist_ok=True);ref=runs[0];n_test=len(ref['label'])
    for r in runs[1:]:
        if not np.array_equal(r['test']['event_number'],ref['test']['event_number']) or not np.array_equal(r['label'],ref['label']): raise SystemExit(f"{r['name']} was tested on other events than {ref['name']}")
    classical={arm:dict(score=ref['test'][arm+'_E2_keV'].astype(float),sep=separation(ref['test'][arm+'_E2_keV'],ref['label'])) for arm in ('itrack','etrack')}
    top=max(float(r['settings']['fraction']) for r in runs);full=[r for r in runs if float(r['settings']['fraction'])==top];style=lambda r:'-' if r['settings']['stem']=='full' else '--'
    names=f"{ref['species'][0]} against {ref['species'][1]} · {int(ref['label'].sum())} + {int((~ref['label']).sum())} test events"

    fig,axes=plt.subplots(1,2,figsize=(15,6.2),layout='constrained')
    for ax,logx in ((axes[0],False),(axes[1],True)):
        for r in full:
            e,b=curve(r['test']['score'].astype(float),r['label']);k=at(r['separation'],0.8)
            ax.plot(b,e,color=COLOR[r['input']],ls=style(r),lw=1.8,label=f"{r['name']}: area {r['separation']['area']:.4f}; at 80%: {100*r['separation']['background_acceptance'][k]:.2f} ± {100*r['separation']['background_acceptance_error'][k]:.2f}%")
        for arm,c in classical.items():
            e,b=curve(c['score'],ref['label']);k=at(c['sep'],0.8)
            ax.plot(b,e,color=COLOR[arm],ls=':',lw=1.6,label=f"E₂ cut, {arm}: area {c['sep']['area']:.4f}; at 80%: {100*c['sep']['background_acceptance'][k]:.2f} ± {100*c['sep']['background_acceptance_error'][k]:.2f}%")
        ax.axhline(0.8,color='#526776',lw=.6);ax.set(xlabel=f"Fraction of {ref['species'][1]} accepted",ylabel=f"Fraction of {ref['species'][0]} accepted",ylim=(0,1));ax.grid(alpha=.15)
        if logx: ax.set_xscale('log');ax.set_xlim(1e-3,1);ax.set_title('Logarithmic background axis')
        else: ax.set_xlim(0,1);ax.set_title('Solid: full-resolution first stage; dashed: half resolution; dotted: classical cut');ax.legend(fontsize=7.5,loc='lower right')
    fig.suptitle(f'Separation on the test events · {names}');fig.savefig(plots/'separation.png',dpi=140);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(15,5.4),layout='constrained')
    for r in runs:
        ep=np.arange(1,len(r['history']['train_loss'])+1);c=COLOR[r['input']];a=1.0 if float(r['settings']['fraction'])==top else 0.45
        axes[0].plot(ep,r['history']['train_loss'],color=c,ls=style(r),lw=1,alpha=a*.6);axes[0].plot(ep,r['history']['validation_loss'],color=c,ls=style(r),lw=1.8,alpha=a,label=f"{r['name']}: selected epoch {int(r['settings']['best_epoch'])}, validation loss {r['history']['validation_loss'].min():.4f}")
        axes[1].plot(ep,100*r['history']['validation_background_at_80'],color=c,ls=style(r),lw=1.8,alpha=a,label=r['name'])
    axes[0].set(xlabel='Epoch',ylabel='Loss (binary cross-entropy)',title='Thick: validation events; thin: training events');axes[0].set_yscale('log');axes[0].grid(alpha=.15);axes[0].legend(fontsize=7.5)
    axes[1].set(xlabel='Epoch',ylabel='Background accepted at 80% signal (%)',title='Validation events');axes[1].set_yscale('log');axes[1].grid(alpha=.15,which='both')
    fig.suptitle(f"Training history · {int(ref['settings']['validation_events'])} validation events");fig.savefig(plots/'history.png',dpi=140);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(15,5.4),layout='constrained');edges=np.linspace(0,np.ceil(ref['test']['ion_drift_cm'].max()/10)*10,7);mid=(edges[:-1]+edges[1:])/2
    for r in full:
        acc,err,eff=by_drift(r['test']['score'].astype(float),r['label'],ref['test']['ion_drift_cm'],edges)
        axes[0].errorbar(mid,100*acc,100*err,color=COLOR[r['input']],ls=style(r),marker='o',ms=4,capsize=2,label=r['name']);axes[1].plot(mid,100*eff,color=COLOR[r['input']],ls=style(r),marker='o',ms=4,label=r['name'])
    for arm,c in classical.items():
        acc,err,eff=by_drift(c['score'],ref['label'],ref['test']['ion_drift_cm'],edges)
        axes[0].errorbar(mid,100*acc,100*err,color=COLOR[arm],ls=':',marker='s',ms=4,capsize=2,label=f'E₂ cut, {arm}');axes[1].plot(mid,100*eff,color=COLOR[arm],ls=':',marker='s',ms=4,label=f'E₂ cut, {arm}')
    axes[0].set(xlabel='Ion drift (cm); electron drift is the chamber length minus this',ylabel='Background accepted (%)',title='One threshold per method, set for 80% of all signal events',ylim=(0,None));axes[0].grid(alpha=.15);axes[0].legend(fontsize=7.5)
    axes[1].set(xlabel='Ion drift (cm)',ylabel='Signal accepted (%)',title='Signal efficiency with the same threshold');axes[1].grid(alpha=.15)
    fig.suptitle(f'Separation against drift · {names} in {len(mid)} drift intervals');fig.savefig(plots/'drift.png',dpi=140);plt.close(fig)

    groups={}
    for r in runs: groups.setdefault((r['input'],r['settings']['stem']),[]).append(r)
    curves={k:sorted(v,key=lambda r:int(r['settings']['train_events'])) for k,v in groups.items() if len(v)>1};learning=''
    if curves:
        fig,ax=plt.subplots(figsize=(8,5.2),layout='constrained')
        for (inp,stem),v in curves.items():
            x=[int(r['settings']['train_events']) for r in v];y=[100*r['separation']['background_acceptance'][at(r['separation'],0.8)] for r in v];e=[100*r['separation']['background_acceptance_error'][at(r['separation'],0.8)] for r in v]
            ax.errorbar(x,y,e,color=COLOR[inp],ls='-' if stem=='full' else '--',marker='o',capsize=3,label=f"{inp}, {stem}: "+', '.join(f'{a}: {b:.2f}%' for a,b in zip(x,y)))
        ax.set(xlabel='Training events (both species)',ylabel='Background accepted at 80% signal (%)',title=f'Learning curve · {n_test} test events');ax.set_xscale('log');ax.grid(alpha=.15,which='both');ax.legend(fontsize=8)
        fig.savefig(plots/'learning_curve.png',dpi=140);plt.close(fig);learning='<h2>Learning curve</h2><p>The same network trained on part of the training events and tested on the same test events.</p><img src="plots/learning_curve.png" alt="Background acceptance against the number of training events">'

    fig,axes=plt.subplots(1,len(full),figsize=(5.2*len(full),4.4),layout='constrained',squeeze=False)
    for ax,r in zip(axes[0],full):
        p=1/(1+np.exp(-r['test']['score'].astype(float)));bins=np.linspace(0,1,51)
        for sel,c,name in ((r['label'],COLOR['itrack'],ref['species'][0]),(~r['label'],COLOR['etrack'],ref['species'][1])): ax.hist(p[sel],bins=bins,histtype='step',color=c,lw=1.4,label=f'{name}: N = {int(sel.sum())}, median {np.median(p[sel]):.3f}')
        ax.set(xlabel='Network output (probability of signal)',ylabel='Events / 0.02',title=r['name']);ax.set_yscale('log');ax.legend(fontsize=7.5,loc='upper center')
    fig.suptitle('Network output on the test events');fig.savefig(plots/'scores.png',dpi=140);plt.close(fig)

    def row(name,inp,stem,train,epoch,sep,minutes):
        cells=''.join(f"<td>{100*sep['background_acceptance'][at(sep,t)]:.2f} ± {100*sep['background_acceptance_error'][at(sep,t)]:.2f}%</td>" for t in TARGETS)
        return f"<tr><td>{html.escape(name)}</td><td>{inp}</td><td>{stem}</td><td>{train}</td><td>{epoch}</td><td>{sep['area']:.4f}</td>{cells}<td>{minutes}</td></tr>"
    rows=''.join(row(r['name'],r['input'],r['settings']['stem'],int(r['settings']['train_events']),f"{int(r['settings']['best_epoch'])} of {len(r['history']['train_loss'])}",r['separation'],f"{float(r['settings']['time_s'])/60:.0f}") for r in runs)
    rows+=''.join(row(f'E₂ cut, {arm}',arm,'','','',c['sep'],'') for arm,c in classical.items())
    head=''.join(f'<th>Background at {int(100*t)}% signal</th>' for t in TARGETS);s=ref['settings']
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ITACA · classification network</title>
<style>body{{font:17px/1.5 system-ui;color:#203343;max-width:1500px;margin:auto;padding:24px}}img{{max-width:100%}}table{{border-collapse:collapse;margin:10px 0 18px}}th,td{{padding:7px 14px;border-bottom:1px solid #ddd;text-align:left}}th{{background:#eef4f7}}a{{color:#185b83}}h2{{margin-top:34px}}code{{background:#eef3f6;padding:2px 5px;border-radius:3px}}</style></head>
<body><a href="{'../'*len(root.resolve().relative_to(Path(__file__).resolve().parents[2]).parts)}notes/html/index.html">ITACA index</a><h1>Classification network · {ref['species'][0]} against {ref['species'][1]}</h1>
<p>Samples <code>{html.escape(ref['runs'][0])}</code> and <code>{html.escape(ref['runs'][1])}</code>, events with one component in each track. Test batches {list(map(int,np.atleast_1d(s['test_batches'])))} ({int(ref['label'].sum())} + {int((~ref['label']).sum())} events), validation batches {list(map(int,np.atleast_1d(s['validation_batches'])))}; the other batches train. Network and procedure: <a href="{'../'*len(root.resolve().relative_to(Path(__file__).resolve().parents[2]).parts)}notes/html/network/classification_network.html">classification network</a>.</p>
<h2>Separation on the test events</h2><p>An event is accepted when the network output exceeds a threshold; the threshold of each column keeps that fraction of the signal. The last rows apply the classical cut on E₂, the smaller endpoint energy of the main track, to the same test events. Errors are binomial.</p>
<table><tr><th>Training</th><th>Input</th><th>First stage</th><th>Training events</th><th>Selected epoch</th><th>Area</th>{head}<th>Time (min)</th></tr>{rows}</table>
<img src="plots/separation.png" alt="Signal efficiency against background acceptance">
<h2>Separation against drift</h2><img src="plots/drift.png" alt="Background acceptance and signal efficiency against drift">
{learning}
<h2>Training history</h2><img src="plots/history.png" alt="Loss and validation acceptance per epoch">
<h2>Network output</h2><img src="plots/scores.png" alt="Network output of signal and background test events">
<p>Numerical results: {' · '.join(f'<a href="{r["name"]}/training.h5">{r["name"]}</a>' for r in runs)}.</p></body></html>'''
    (root/'index.html').write_text(page);print(root/'index.html')

if __name__=='__main__':
    if len(sys.argv)!=2: raise SystemExit(__doc__)
    render(sys.argv[1])
