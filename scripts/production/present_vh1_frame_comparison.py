#!/usr/bin/env python3
"""Present the VH1 frame comparison saved by reduce_vh1_frame_comparison.jl:
one table and the distributions of the run against its reference, same events.

    python scripts/production/present_vh1_frame_comparison.py results/readout/<run>/vh1_frame_comparison/comparison.h5
"""
from pathlib import Path
import sys,html,os
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
COLORS={'run':'#c37b34','reference':'#287ba0'}

def load(g):
    e=g['events'];v={k:e[k][:] for k in e}
    h=g['histograms'];hist={k:h[k][:] for k in h}
    hits={k:(g['hits'][k][()] if g['hits'][k].shape==() else g['hits'][k][:]) for k in g['hits']}
    cond={}
    for k in g['conditioning']:
        c=g['conditioning'][k];cond[int(c['buckets'][()])]={n:(c[n][()] if c[n].shape==() else c[n][:]) for n in c}
    meta={k:g[k][()] for k in ('run_label','frame_ms','conversion_window_ms','ion_drift_velocity_cm_s','ion_bucket_ms')}
    meta={k:(x.decode() if isinstance(x,bytes) else float(x)) for k,x in meta.items()}
    return dict(events=v,hist=hist,hits=hits,cond=cond,meta=meta)

def fwhm_gauss(x):return 2.355*np.std(x)

def summary(d):
    v=d['events'];h=d['hits']
    n=len(v['event_number']);qn=np.sum(v['empty_cells']);qm=np.sum(v['empty_sum'])/qn
    noise=np.sqrt(np.sum(v['empty_sumsq'])/qn-qm*qm)
    cells=np.sum(v['signal_pad_buckets'])
    e=h['ion_energy_keV']
    return dict(events=n,pads=int(np.sum(v['ion_pads'])),cells=int(cells),
        noise_observed=noise,noise_predicted=float(np.mean(v['predicted_bucket_noise'])),
        snr=float(np.sum(v['signal_snr_sum'])/cells),input_per_cell=float(np.sum(v['signal_input_ions'])/cells),
        saturated=int(np.sum(v['saturated_pads'])),adc_limit=int(np.sum(v['adc_limit_pads'])),invalid=int(np.sum(v['invalid_pads'])),
        rms=float(np.mean(v['bucket_residual_rms'])),signal_rms=float(np.mean(v['signal_residual_rms'])),
        energy_mean=float(np.mean(e)),energy_fwhm=float(fwhm_gauss(e)),energy_fwhm_pct=float(100*fwhm_gauss(e)/np.mean(e)),
        hits_above=float(np.mean(h['itrack_hits_above_cut'])),hits_total=float(np.mean(h['itrack_hits'])),
        ions=float(np.mean(h['ions'])),time=float(np.mean(v['vh1_s'])),
        cond=d['cond'])

def centres(edges):return .5*(edges[1:]+edges[:-1])

def plots(run,ref,out,cut):
    plots=out/'plots';plots.mkdir(exist_ok=True)
    labels={'run':f"{run['meta']['frame_ms']:.0f} ms frames",'reference':f"{ref['meta']['frame_ms']:.0f} ms frames"}
    pair=[('run',run),('reference',ref)]
    def stats(x):return f"n {len(x)}, mean {np.mean(x):.2f}, std {np.std(x):.2f}, median {np.median(x):.2f}"
    def hstats(counts,edges):
        c=centres(edges);n=counts.sum();m=(c*counts).sum()/n;s=np.sqrt(((c-m)**2*counts).sum()/n)
        return f"n {n}, mean {m:.2f}, std {s:.2f}"
    fig,ax=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    for key,d in pair:
        h=d['hist'];e=h['residual_edges']
        ax[0,0].stairs(h['empty_residual_counts'],e,label=f"{labels[key]}: {hstats(h['empty_residual_counts'],e)}",color=COLORS[key],lw=1.6)
        ax[0,1].stairs(h['signal_residual_counts'],e,label=f"{labels[key]}: {hstats(h['signal_residual_counts'],e)}",color=COLORS[key],lw=1.6)
        ax[1,0].stairs(h['snr_counts'],h['snr_edges'],label=f"{labels[key]}: {hstats(h['snr_counts'],h['snr_edges'])}",color=COLORS[key],lw=1.6)
        x=d['hits']['ion_energy_keV']
        ax[1,1].hist(x,bins=np.linspace(2350,2550,41),histtype='step',lw=1.6,color=COLORS[key],label=f"{labels[key]}: {stats(x)} keV, FWHM(2.355 std) {fwhm_gauss(x):.1f} keV")
    ax[0,0].set(title='Noise after reconstruction: valid pad-buckets without input ions',xlabel='reconstructed ions',ylabel='pad-buckets')
    ax[0,1].set(title='Residual on valid pad-buckets with input ions',xlabel='reconstructed minus input ions',ylabel='pad-buckets')
    ax[1,0].set(title='S/N per pad-bucket with input ions (last bin: overflow)',xlabel='reconstructed ions / predicted noise',ylabel='pad-buckets')
    ax[1,1].set(title='Ion energy per event, all accepted events',xlabel='keV',ylabel='events');ax[1,1].set_ylim(top=1.45*ax[1,1].get_ylim()[1])
    for a in ax.flat:a.grid(alpha=.15);a.legend(fontsize=8,loc='upper left');a.title.set_fontsize(11)
    fig.savefig(plots/'distributions.png',dpi=140);plt.close(fig)
    fig,ax=plt.subplots(1,3,figsize=(15,4.6),layout='constrained')
    for key,d in pair:
        x=d['hits']['itrack_hits_above_cut'];y=d['hits']['itrack_hits']
        ax[0].hist(x,bins=np.arange(-.5,max(np.max(run['hits']['itrack_hits_above_cut']),np.max(ref['hits']['itrack_hits_above_cut']))+1.5,2),histtype='step',lw=1.6,color=COLORS[key],label=f"{labels[key]}: {stats(x)}")
        ax[1].hist(y,bins=40,histtype='step',lw=1.6,color=COLORS[key],label=f"{labels[key]}: {stats(y)}")
        nb=max(d['cond']);c=d['cond'][nb]
        ax[2].plot(np.arange(1,len(c['reading_drop_noise_ions'])+1),c['reading_drop_noise_ions'],'.-',color=COLORS[key],ms=4,label=f"{labels[key]}: {nb} buckets, {len(c['reading_drop_noise_ions'])} readings; bucket noise after fit {np.median(c['predicted_bucket_noise_ions']):.1f} ions")
    ax[0].set(title=f'itrack hits at or above {cut:g} keV per event',xlabel='hits',ylabel='events')
    ax[1].set(title='itrack hits written per event',xlabel='hits',ylabel='events')
    ax[2].set(title='Noise of each reading drop before the fit, longest event',xlabel='reading',ylabel='ions')
    for a in ax.flat:a.grid(alpha=.15);a.legend(fontsize=8);a.title.set_fontsize(11)
    fig.savefig(plots/'hits_and_readings.png',dpi=140);plt.close(fig)

def render(path):
    out=path.parent
    with h5py.File(path) as f:
        run=load(f['run']);ref=load(f['reference']);batches=f['batches'][:];n=int(f['events'][()])
    cut=float(run['hits']['hit_cut_keV']);output_cut=float(run['hits']['output_cut_keV'])
    plots(run,ref,out,cut)
    a=summary(run);b=summary(ref);mr=run['meta'];mf=ref['meta']
    def row(name,fa,fb):return f'<tr><th>{html.escape(name)}</th><td>{fa}</td><td>{fb}</td></tr>'
    def both(name,fmt,key):return row(name,fmt.format(a[key]),fmt.format(b[key]))
    nb=sorted(set(a['cond'])&set(b['cond']));nbmax=max(nb)
    rows=''.join([
        row('Run card',html.escape(mr['run_label']),html.escape(mf['run_label'])),
        row('Ion drift velocity / bucket / frame',f"{mr['ion_drift_velocity_cm_s']:.1f} cm/s · {mr['ion_bucket_ms']:.0f} ms · {mr['frame_ms']:.0f} ms, 2 × {mr['conversion_window_ms']:.1f} ms",
            f"{mf['ion_drift_velocity_cm_s']:.1f} cm/s · {mf['ion_bucket_ms']:.0f} ms · {mf['frame_ms']:.0f} ms, 2 × {mf['conversion_window_ms']:.1f} ms"),
        both('Events',"{}",'events'),both('Pads read',"{:,}",'pads'),both('Pad-buckets with input ions',"{:,}",'cells'),
        both('Noise after reconstruction, observed on pad-buckets without input (ions, RMS)',"{:.2f}",'noise_observed'),
        both('Noise after reconstruction, predicted from the covariance (ions)',"{:.2f}",'noise_predicted'),
        both('Mean input ions per pad-bucket with input',"{:.1f}",'input_per_cell'),
        both('Signal-to-noise per pad-bucket with input, mean',"{:.2f}",'snr'),
        both('Saturated pads (noiseless output left the window)',"{}",'saturated'),both('Pads with a converter code at a range end',"{}",'adc_limit'),both('Invalid pads',"{}",'invalid'),
        both('Bucket residual RMS, valid pad-buckets (ions)',"{:.2f}",'rms'),both('Bucket residual RMS, pad-buckets with input (ions)',"{:.2f}",'signal_rms'),
        both('Ion energy per event: mean (keV)',"{:.1f}",'energy_mean'),both('Ion energy per event: FWHM as 2.355 × std (keV)',"{:.1f}",'energy_fwhm'),both('Ion energy per event: FWHM / mean',"{:.2f} %",'energy_fwhm_pct'),
        both('Reconstructed ions per event, mean',"{:,.0f}",'ions'),
        both(f'itrack hits at or above {cut:g} keV per event, mean',"{:.1f}",'hits_above'),both(f'itrack hits written per event, mean (output cut {output_cut:g} keV)',"{:.1f}",'hits_total'),
        row(f'Response matrix, {nbmax} buckets: readings',a['cond'][nbmax]['readings'],b['cond'][nbmax]['readings']),
        row(f'Condition number of the response matrix, {nbmax} buckets',f"{a['cond'][nbmax]['condition_response']:.1f}",f"{b['cond'][nbmax]['condition_response']:.1f}"),
        row(f'Condition number after noise whitening, {nbmax} buckets',f"{a['cond'][nbmax]['condition_whitened']:.1f}",f"{b['cond'][nbmax]['condition_whitened']:.1f}"),
        row('Noise of the first reading drop, before the fit (ions)',f"{a['cond'][nbmax]['reading_drop_noise_ions'][0]:.1f}",f"{b['cond'][nbmax]['reading_drop_noise_ions'][0]:.1f}"),
        both('VH1 electronics and reconstruction (s/event)',"{:.3f}",'time')])
    cond_rows=''.join(f"<tr><td>{k}</td><td>{a['cond'][k]['readings']}</td><td>{a['cond'][k]['condition_whitened']:.1f}</td><td>{np.median(a['cond'][k]['predicted_bucket_noise_ions']):.2f}</td><td>{b['cond'][k]['condition_whitened']:.1f}</td><td>{np.median(b['cond'][k]['predicted_bucket_noise_ions']):.2f}</td></tr>" for k in nb)
    index=os.path.relpath(ROOT/'notes/html/index.html',out);pilots=os.path.relpath(ROOT/'results/readout/vh1_pilots_100/index.html',out)
    page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>VH1 frames: {mr['frame_ms']:.0f} ms against {mf['frame_ms']:.0f} ms</title><style>body{{font:17px/1.5 system-ui;max-width:1150px;margin:2em auto;color:#203343}}table{{border-collapse:collapse;margin:1em 0}}th,td{{text-align:left;padding:7px 12px;border-bottom:1px solid #ddd;vertical-align:top}}img{{max-width:100%}}</style>
<a href="{index}">ITACA index</a> · <a href="{pilots}">VH1 pilots</a>
<h1>VH1 with {mr['frame_ms']:.0f} ms frames against {mf['frame_ms']:.0f} ms: the same {n} bb0nu events</h1>
<p>Batches {', '.join(str(x) for x in batches)} of the pilot NEXUS files, read with the same seeds, 3 mm pads and the same VH1 circuit (20 fF, device changes). The run changes only the ion drift velocity (<code>drift.ion_drift_velocity_cm_s</code>), the ion bucket (<code>ions.ion_bucket_ms</code>) and the VH1 card, so one frame spans 3.0 mm in both and the ions per pad-bucket are statistically the same. The electron arm ran unchanged and is not shown.</p>
<p>The circuit's SPICE results carry over: the gain, the charge fraction and the output window are AC results at 1 Hz, the noise is the stationary output spectrum, which the readout integrates over each conversion window, and the capacitor swap settles in microseconds, far below a {mr['conversion_window_ms']:g} ms window. No SPICE rerun was needed for the shorter frame.</p>
<table><tr><th></th><th>Run</th><th>Reference</th></tr>{rows}</table>
<p>Noise after reconstruction is the RMS of the reconstructed ions on valid pad-buckets without input ions; the prediction is the mean of the covariance diagonal. The bucket residual compares reconstructed with input ions. The ion energy and the hits are read from the chained hit files of both runs; hits include the noise hits of the fixed frame.</p>
<img src="plots/distributions.png" alt="Noise, residual, S/N and energy distributions">
<img src="plots/hits_and_readings.png" alt="Hits per event and reading noise">
<h2>Response matrix by event length</h2>
<table><tr><th>Buckets</th><th>Readings</th><th>Run: condition (whitened)</th><th>Run: bucket noise median (ions)</th><th>Reference: condition (whitened)</th><th>Reference: bucket noise median (ions)</th></tr>{cond_rows}</table>
<p>Numerical results: <a href="comparison.h5">comparison.h5</a> (schema itaca.vh1_frame_comparison.v1), made by <code>scripts/production/reduce_vh1_frame_comparison.jl</code>; this page by <code>present_vh1_frame_comparison.py</code>.</p></html>'''
    (out/'index.html').write_text(page)
    print(out/'index.html')
    for k in ('noise_observed','noise_predicted','snr','saturated','adc_limit','rms','energy_mean','energy_fwhm','hits_above'):
        print(f"{k}: run {a[k]:.3f} reference {b[k]:.3f}" if isinstance(a[k],float) else f"{k}: run {a[k]} reference {b[k]}")

if __name__=='__main__':render(Path(sys.argv[1]))
