#!/usr/bin/env python3
"""Page of a round of network trainings.

  python scripts/network/network_presenter.py ROUND_DIR
  python scripts/network/network_presenter.py ROUND_DIR --only=NAME,NAME,... [--page=selected.html] [--comparison=comparison.h5] [--remark=NAME=text]

ROUND_DIR holds one folder per training with the training.h5 of
network_trainer.py. Writes plots and index.html in ROUND_DIR. The classical cut
on E2 is evaluated on the test events of the first training, which all
trainings of a round share. No training is repeated.

--only writes a page of the named trainings alone, in the order given, as
ROUND_DIR/PAGE with its plots in ROUND_DIR/plots/ under the name of the page;
index.html is not touched. That page also marks the classical cut on E2 of
all one-component events of the two samples, read from the comparison.h5 of
prepare_species_comparison.jl (default results/readout/<ROUND_DIR name>_comparison/comparison.h5).
--remark, repeatable, adds a note to the epochs of a training (for example
both_full=cancelled).
"""
import argparse,html,sys
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).parent))
from network_classifier import separation

EXTRA_COLORS=('#7a4fa3','#b23a48','#6b7480');LOSS_TICKS=(0.1,0.2,0.5,1,2,5,10);PERCENT_TICKS=(1,2,5,10,20,50,100)
COLOR={'itrack':'#287ba0','etrack':'#b16d31','both':'#3d8f5a'};TARGETS=(0.9,0.8,0.7,0.6)
def text(v): return v.decode() if isinstance(v,bytes) else v

def load(path):
    with h5py.File(path) as f:
        r=dict(name=path.parent.name,input=text(f['input'][()]),settings={k:(text(v[()]) if isinstance(v,h5py.Dataset) else None) for k,v in f['settings'].items()},
               history={k:v[:] for k,v in f['history'].items()},test={k:v[()] for k,v in f['test'].items() if isinstance(v,h5py.Dataset)},
               species=tuple(text(f[k][()]).replace('_',' ') for k in ('signal_species','background_species')),runs=(text(f['signal_run_label'][()]),text(f['background_run_label'][()])),
               selection=text(f['signal_selection'][()]) if 'signal_selection' in f else 'one component in each track')
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
    classical={arm:dict(score=ref['test'][arm+'_E2_keV'].astype(float),sep=separation(ref['test'][arm+'_E2_keV'],ref['label'])) for arm in ('itrack','etrack') if arm+'_E2_keV' in ref['test']}   # iItaca files have no etrack
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
<p>Samples <code>{html.escape(ref['runs'][0])}</code> and <code>{html.escape(ref['runs'][1])}</code>; {html.escape(ref['selection'])}. Test batches {list(map(int,np.atleast_1d(s['test_batches'])))} ({int(ref['label'].sum())} + {int((~ref['label']).sum())} events), validation batches {list(map(int,np.atleast_1d(s['validation_batches'])))}; the other batches train. Network and procedure: <a href="{'../'*len(root.resolve().relative_to(Path(__file__).resolve().parents[2]).parts)}notes/html/network/classification_network.html">classification network</a>.</p>
<h2>Separation on the test events</h2><p>An event is accepted when the network output exceeds a threshold; the threshold of each column keeps that fraction of the signal. The last rows apply the classical cut on E₂, the smaller endpoint energy of the main track, to the same test events. Errors are binomial.</p>
<table><tr><th>Training</th><th>Input</th><th>First stage</th><th>Training events</th><th>Selected epoch</th><th>Area</th>{head}<th>Time (min)</th></tr>{rows}</table>
<img src="plots/separation.png" alt="Signal efficiency against background acceptance">
<h2>Separation against drift</h2><img src="plots/drift.png" alt="Background acceptance and signal efficiency against drift">
{learning}
<h2>Training history</h2><img src="plots/history.png" alt="Loss and validation acceptance per epoch">
<h2>Network output</h2><img src="plots/scores.png" alt="Network output of signal and background test events">
<p>Numerical results: {' · '.join(f'<a href="{r["name"]}/training.h5">{r["name"]}</a>' for r in runs)}.</p></body></html>'''
    (root/'index.html').write_text(page);print(root/'index.html')

def learning_rate(r,points=2000):
    """Learning rate against epoch: the OneCycleLR of network_trainer.py (cosine, 30% rise, start max/25, end max/25e4), one step per batch."""
    s=r['settings'];n=int(s['train_events']);b=int(s['batch_size']);per_epoch=int(np.ceil(n/b));total=int(s['epochs'])*per_epoch;top=float(s['learning_rate'])
    done=len(r['history']['train_loss'])*(per_epoch-(n%b==1));k=np.unique(np.linspace(0,done-1,points).astype(int));rise=0.3*total-1;start=top/25;end=start/1e4
    cos=lambda a,z,p:z+(a-z)/2*(np.cos(np.pi*p)+1)
    return k/(per_epoch-(n%b==1)),np.where(k<=rise,cos(start,top,k/rise),cos(top,end,(k-rise)/(total-1-rise)))

def classical_cut(path,ref):
    """E2 cut on all one-component events of the two samples, from the comparison.h5 of prepare_species_comparison.jl."""
    with h5py.File(path) as f:
        labels=(text(f['signal/run_label'][()]),text(f['background/run_label'][()]))
        if labels!=ref['runs']: raise SystemExit(f'{path} compares {labels}, the trainings use {ref["runs"]}')
        out={}
        for arm in ('itrack','etrack'):
            g=f[f'E2_cut/{arm}/one_component'];w=g['working_points']
            out[arm]=dict(area=float(g['area'][()]),efficiency=g['signal_efficiency'][:],acceptance=g['background_acceptance'][:],signal_events=int(g['signal_events'][()]),background_events=int(g['background_events'][()]),
                          target_signal_efficiency=w['target_signal_efficiency'][:],background_acceptance=w['background_acceptance'][:],background_acceptance_error=w['background_acceptance_error'][:])
    return out

def render_only(root,names,page,comparison,remarks={}):
    root=Path(root);missing=[n for n in names if not (root/n/'training.h5').exists()]
    if missing: raise SystemExit(f'No training.h5 for {missing} under {root}')
    runs=[load(root/n/'training.h5') for n in names];ref=runs[0];plots=root/'plots';used=set()
    for r in runs:      # the colour of the input; another colour for a second training with the same input
        r['color']=COLOR[r['input']] if r['input'] not in used else next(c for c in EXTRA_COLORS if c not in used);used|={r['input'],r['color']};r['model']=r['settings'].get('model') or 'dense';plots.mkdir(exist_ok=True);tag=Path(page).stem;up='../'*len(root.resolve().relative_to(Path(__file__).resolve().parents[2]).parts)
    for r in runs[1:]:
        if not np.array_equal(r['test']['event_number'],ref['test']['event_number']) or not np.array_equal(r['label'],ref['label']) or r['runs']!=ref['runs']: raise SystemExit(f"{r['name']} was tested on other events than {ref['name']}")
    sig,bkg=ref['species'];n_sig=int(ref['label'].sum());n_bkg=int((~ref['label']).sum());targets=ref['separation']['target_signal_efficiency'];k80=at(ref['separation'],0.8)
    comparison=Path(comparison) if comparison else Path(__file__).resolve().parents[2]/'results/readout'/f'{root.resolve().name}_comparison/comparison.h5'
    if not comparison.exists(): raise SystemExit(f'No comparison file {comparison}: give --comparison')
    cut=classical_cut(comparison,ref);cut_test={arm:separation(ref['test'][arm+'_E2_keV'],ref['label']) for arm in ('itrack','etrack')}
    done=lambda r:len(r['history']['train_loss']);planned=lambda r:int(r['settings']['epochs'])
    epochs=lambda r:f"{done(r)} of {planned(r)}";remark=lambda r:f", {remarks[r['name']]}" if r['name'] in remarks else '';pc=lambda sep,k:f"{100*sep['background_acceptance'][k]:.2f} ± {100*sep['background_acceptance_error'][k]:.2f}%"

    # Loss per epoch, one panel per training
    cols=len(runs) if len(runs)<=3 else int(np.ceil(len(runs)/2));nrows=int(np.ceil(len(runs)/cols));fig,axes=plt.subplots(nrows,cols,figsize=(max(5.4,16.2/cols)*cols,5.2*nrows),layout='constrained',sharey=True,squeeze=False)
    for ax in axes.flat[len(runs):]: ax.set_visible(False)
    for ax,r in zip(axes.flat,runs):
        h=r['history'];ep=np.arange(1,done(r)+1);c=r['color'];best=int(np.argmin(h['validation_loss']));spikes=ep[h['validation_loss']>2*h['train_loss']]
        ax.plot(ep,h['train_loss'],color=c,lw=1.2,ls='--',marker='o',ms=3,label=f"training events: last {h['train_loss'][-1]:.4f}")
        ax.plot(ep,h['validation_loss'],color=c,lw=1.9,marker='o',ms=4,label=f"validation events: minimum {h['validation_loss'][best]:.4f} at epoch {best+1}, maximum {h['validation_loss'].max():.2f}")
        ax.plot(best+1,h['validation_loss'][best],marker='o',ms=11,mfc='none',mec='#203343',ls='none',label=f"selected weights (epoch {int(r['settings']['best_epoch'])})")
        ax.set_yscale('log');ax.set_xlim(0.5,planned(r)+0.5);ax.set_xticks(range(2,planned(r)+1,2));ax.grid(alpha=.2,which='both');ax.set_xlabel('Epoch')
        ax.set_yticks(LOSS_TICKS);ax.set_yticklabels([f'{v:g}' for v in LOSS_TICKS]);ax.minorticks_off()
        ax.set_title(f"{r['name']} · {r['model']} · {epochs(r)} epochs{remark(r)}\nvalidation loss above twice the training loss at {'epochs '+', '.join(map(str,spikes)) if len(spikes) else 'no epoch'}",fontsize=10);ax.legend(fontsize=8,loc='upper right')
    for ax in axes[:,0]: ax.set_ylabel('Loss (binary cross-entropy per event)')
    axes[0][0].set_ylim(0.1,12*max(r['history']['validation_loss'].max() for r in runs))
    fig.suptitle(f"Loss per epoch · {int(ref['settings']['train_events'])} training and {int(ref['settings']['validation_events'])} validation events · {sig} against {bkg}");fig.savefig(plots/f'{tag}_loss.png',dpi=140);plt.close(fig)

    # Validation acceptance and learning rate per epoch
    fig,axes=plt.subplots(1,2,figsize=(15,5.4),layout='constrained',width_ratios=(1.5,1))
    for r in runs:
        h=r['history'];y=100*h['validation_background_at_80'];axes[0].plot(np.arange(1,done(r)+1),y,color=r['color'],lw=1.8,marker='o',ms=4,label=f"{r['name']} ({epochs(r)} epochs{remark(r)}): minimum {y.min():.2f}% at epoch {int(np.argmin(y))+1}, last {y[-1]:.2f}%, maximum {y.max():.1f}%")
    axes[0].set(xlabel='Epoch',ylabel=f'{bkg.capitalize()} accepted at 80% {sig} efficiency (%)',title=f"{int(ref['settings']['validation_events'])} validation events\nthreshold set each epoch to keep 80% of the {sig} validation events",xlim=(0.5,max(map(planned,runs))+0.5))
    axes[0].set_yscale('log');axes[0].set_xticks(range(2,max(map(planned,runs))+1,2));axes[0].grid(alpha=.2,which='both');axes[0].legend(fontsize=8,loc='upper right');axes[0].set_yticks(PERCENT_TICKS);axes[0].set_yticklabels([f'{v:g}' for v in PERCENT_TICKS]);axes[0].minorticks_off();axes[0].set_ylim(1,600)
    seen={}
    for r in runs:
        x,lr=learning_rate(r);key=(int(r['settings']['train_events']),int(r['settings']['batch_size']),planned(r),float(r['settings']['learning_rate']))
        if key in seen and done(r)<=seen[key]: continue
        seen[key]=done(r);axes[1].plot(x,lr,color='#203343' if len(seen)==1 else r['color'],lw=1.8,label=f"{int(np.ceil(key[0]/key[1]))} steps per epoch\nstarts at {lr[0]:.1e}\nmaximum {lr.max():.1e} at epoch {x[np.argmax(lr)]:.1f}\n{lr[-1]:.1e} at the last step of epoch {done(r)}")
    same=len({(int(r['settings']['train_events']),int(r['settings']['batch_size']),planned(r),float(r['settings']['learning_rate'])) for r in runs})==1
    axes[1].set(xlabel='Epoch',ylabel='Learning rate',title=('The same schedule in the %d trainings'%len(runs) if same else 'Schedule of each training')+'\nOneCycleLR, one step per batch',xlim=(0,max(map(planned,runs))));axes[1].set_yscale('log');axes[1].set_xticks(range(0,max(map(planned,runs))+1,2));axes[1].grid(alpha=.2,which='both');axes[1].legend(fontsize=8.5,loc='lower center')
    fig.suptitle(f'Validation acceptance and learning rate per epoch · {sig} against {bkg}');fig.savefig(plots/f'{tag}_validation.png',dpi=140);plt.close(fig)

    # Background acceptance against signal efficiency on the test events
    fig,ax=plt.subplots(figsize=(11,7.2),layout='constrained');floor=1/n_bkg;low=1.0
    for r in runs:
        e,b=curve(r['test']['score'].astype(float),r['label']);low=min(low,max(b[e>=0.3].min(),floor));ax.plot(e,np.maximum(b,floor/2),color=r['color'],lw=1.9,label=f"{r['name']} ({epochs(r)} epochs{remark(r)}): area {r['separation']['area']:.4f}; at 80%: {pc(r['separation'],k80)}")
    for arm,c in cut.items():
        k=int(np.argmin(np.abs(c['target_signal_efficiency']-0.8)));ax.plot(c['efficiency'],c['acceptance'],color=COLOR[arm],ls=':',lw=1.5)
        ax.plot(c['target_signal_efficiency'][k],c['background_acceptance'][k],color=COLOR[arm],marker='s',ms=9,mec='#203343',ls=':',lw=1.5,label=f"E₂ cut, {arm}, all {c['signal_events']} + {c['background_events']} events with one {arm} component: area {c['area']:.4f}; at 80%: {pc(c,k)}")
    ax.axvline(0.8,color='#526776',lw=.7);ax.set_yscale('log');ax.set(xlim=(0.3,1),ylim=(10**np.floor(np.log10(low)),1),xlabel=f'Fraction of {sig} accepted',ylabel=f'Fraction of {bkg} accepted');ax.grid(alpha=.2,which='both');ax.legend(fontsize=8.5,loc='upper left')
    ax.set_title(f'Separation on the test events · {n_sig} {sig} + {n_bkg} {bkg}\nSolid: network; dotted with square at 80%: classical cut on E₂; one test event of {bkg} is {floor:.1e}',fontsize=11);fig.savefig(plots/f'{tag}_roc.png',dpi=140);plt.close(fig)

    cells=lambda sep:''.join(f"<td>{pc(sep,int(np.argmin(np.abs(sep['target_signal_efficiency']-t))))}</td>" for t in targets)
    hours=lambda t:f'{int(t//3600)} h {int(round(t%3600/60)):02d} min'
    rows=''.join(f"<tr><td>{html.escape(r['name'])}</td><td>{r['input']}</td><td>{r['model']}</td><td>{int(r['settings']['parameters']):,}</td><td>{int(r['settings']['train_events']):,}</td><td>{epochs(r)}{remark(r)}; selected {int(r['settings']['best_epoch'])}</td><td>{hours(float(r['settings']['time_s']))}</td><td>{r['history']['epoch_s'].mean():.0f}</td><td>{r['separation']['area']:.4f}</td>{cells(r['separation'])}</tr>" for r in runs)
    rows+=''.join(f"<tr><td>E₂ cut, {arm}, test events</td><td>{arm}</td><td></td><td></td><td></td><td></td><td></td><td></td><td>{c['area']:.4f}</td>{cells(c)}</tr>" for arm,c in cut_test.items())
    rows+=''.join(f"<tr><td>E₂ cut, {arm}, all {c['signal_events']:,} + {c['background_events']:,} events with one {arm} component</td><td>{arm}</td><td></td><td></td><td></td><td></td><td></td><td></td><td>{c['area']:.4f}</td>{cells(c)}</tr>" for arm,c in cut.items())
    head=''.join(f'<th>{int(round(100*t))}%</th>' for t in targets);s=ref['settings'];unfinished=[r for r in runs if done(r)<planned(r)]
    note=''.join(f"<p><code>{html.escape(r['name'])}</code> has completed {epochs(r)} epochs{remark(r)}; its numbers are those of the weights with the lowest validation loss among them.</p>" for r in unfinished)
    stems=sorted({r['settings']['stem'] for r in runs})
    page_text=f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ITACA · classification network · {', '.join(html.escape(n) for n in names)}</title>
<style>body{{font:17px/1.5 system-ui;color:#203343;max-width:1500px;margin:auto;padding:24px}}img{{max-width:100%}}table{{border-collapse:collapse;margin:10px 0 18px}}th,td{{padding:7px 14px;border-bottom:1px solid #ddd;text-align:left}}th{{background:#eef4f7}}a{{color:#185b83}}h2{{margin-top:34px}}code{{background:#eef3f6;padding:2px 5px;border-radius:3px}}</style></head>
<body><a href="{up}notes/html/index.html">ITACA index</a><h1>Classification network · {sig} against {bkg} · {', '.join(stems)}-resolution first stage</h1>
<p>Samples <code>{html.escape(ref['runs'][0])}</code> and <code>{html.escape(ref['runs'][1])}</code>; {html.escape(ref['selection'])}. Test batches {list(map(int,np.atleast_1d(s['test_batches'])))} ({n_sig:,} {sig} + {n_bkg:,} {bkg}), validation batches {list(map(int,np.atleast_1d(s['validation_batches'])))} ({int(s['validation_events']):,} events); the other batches train. Batch size {int(s['batch_size'])}, maximum learning rate {float(s['learning_rate']):g}, weight decay {float(s['weight_decay']):g}, seed {int(s['seed'])}. Network and procedure: <a href="{up}notes/html/network/classification_network.html">classification network</a>.</p>{note}
<h2>Separation on the test events</h2><p>An event is accepted when the network output exceeds a threshold; the threshold of each column keeps that fraction of the {sig} test events. The weights are those of the epoch with the lowest validation loss. The last rows apply the classical cut on E₂, the smaller endpoint energy of the main track: to the same test events, and to all events of the two samples with one component in that track (<a href="{up}results/readout/{comparison.resolve().parent.name}/index.html">comparison and classical separation</a>). Errors are binomial.</p>
<table><tr><th rowspan="2">Training</th><th rowspan="2">Input</th><th rowspan="2">Model</th><th rowspan="2">Parameters</th><th rowspan="2">Training events</th><th rowspan="2">Epochs</th><th rowspan="2">Time</th><th rowspan="2">Time per epoch (s)</th><th rowspan="2">Area under the curve</th><th colspan="{len(targets)}">{bkg.capitalize()} accepted at a {sig} efficiency of</th></tr><tr>{head}</tr>{rows}</table>
<img src="plots/{tag}_roc.png" alt="Background acceptance against signal efficiency on the test events">
<h2>Loss per epoch</h2><p>Training loss: mean over the epoch, with the random symmetry operations and shifts. Validation loss: after the epoch, without them.</p><img src="plots/{tag}_loss.png" alt="Training and validation loss per epoch">
<h2>Validation acceptance and learning rate</h2><img src="plots/{tag}_validation.png" alt="Validation background acceptance at 80% signal efficiency and learning rate per epoch">
<p>Numerical results: {' · '.join(f'<a href="{r["name"]}/training.h5">{r["name"]}</a>' for r in runs)} · <a href="{up}results/readout/{comparison.resolve().parent.name}/comparison.h5">comparison.h5</a>.</p></body></html>"""
    (root/page).write_text(page_text);print(root/page)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter);p.add_argument('round_dir');p.add_argument('--only',default='');p.add_argument('--page',default='selected.html');p.add_argument('--comparison',default='');p.add_argument('--remark',action='append',default=[])
    a=p.parse_args();only=[n for n in a.only.split(',') if n]
    if only:
        if Path(a.page).name!=a.page or a.page=='index.html': raise SystemExit('--page is a file name in ROUND_DIR other than index.html')
        render_only(a.round_dir,only,a.page,a.comparison,dict(v.split('=',1) for v in a.remark))
    else: render(a.round_dir)
