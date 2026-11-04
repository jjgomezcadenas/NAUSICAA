#!/usr/bin/env python3
"""Present saved Julia ion-CSA response diagnostics as HTML, figures and LaTeX."""
from pathlib import Path
import argparse,csv,hashlib,json,os,tomllib
import h5py,numpy as np
ROOT=Path(__file__).resolve().parents[2]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/ion_csa_matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plotly.offline import get_plotlyjs

def plain(x):
 if isinstance(x,bytes):return x.decode()
 if isinstance(x,np.ndarray):return [plain(v) for v in x]
 return x.item() if isinstance(x,np.generic) else x

def tree(g):
 return {'meta':{k:plain(v) for k,v in g.attrs.items()},**{k:tree(v) if isinstance(v,h5py.Group) else plain(v[()]) for k,v in g.items()}}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',required=True,type=Path);a=p.parse_args();directory=a.results.resolve()
 config=tomllib.loads((directory/'manifest.toml').read_text());run=tomllib.loads((ROOT/config['run_card']).read_text())
 out=ROOT/run['report_directory']/config['study_id'];out.mkdir(parents=True,exist_ok=True)
 with h5py.File(directory/'summary.h5') as h:
  names=[plain(n) for n in h['groups'][:]];data={n:tree(h[n]) for n in names}
  calibration={k:tree(g) for k,g in h['calibration'].items() if isinstance(g,h5py.Group)}
  validation=tree(h['validation'])
 first_meta=data[names[0]]['meta'];pitch_mm=float(first_meta['pitch_mm']);frame_ms=round(1000*float(first_meta['frame_s']));frame_key=f'{frame_ms}ms'
 geometry=f'{pitch_mm:g} mm / {frame_ms} ms'
 colors=['#176589','#bd6333','#5f7e43'];labels={'bb0nu':'Double beta','single_electrons':'Single electron'}
 fig,axes=plt.subplots(2,2,figsize=(11,8))
 for i,name in enumerate(names):
  d=data[name];label=labels[d['meta']['sample']];color=colors[i]
  edges=np.array(d['snr']['edges']);count=np.array(d['snr']['count']);axes[0,0].stairs(count/sum(count),edges,label=label,color=color)
  frac=np.array(d['fraction_events_with_compression']);severities=np.array(d['compression_fractions'])
  axes[0,1].plot(100*severities,100*frac,'o-',color=color,label=label)
  for key,style in (('linear_bias','--'),('mean_bias','-')):
   e=np.array(d[key]['edges']);n=np.array(d[key]['count']);axes[1,1].stairs(n/sum(n),100*e,color=color,linestyle=style,label=label+(' linear' if style=='--' else ' CSA'))
 axes[0,0].set_xscale('log');axes[0,0].set_xlabel('Expected S/N per positive pixel–frame');axes[0,0].set_ylabel('Fraction of cells / bin')
 axes[0,1].set_xlabel('Local response loss criterion (%)');axes[0,1].set_ylabel('Events with an affected pad (%)')
 for i,(key,c) in enumerate(sorted(calibration.items())):
  cv=np.array(c['covariance_e2']);axes[1,0].plot(np.arange(9),cv[:9]/cv[0],'o-',color=colors[i],label=key.replace('ms',' ms'))
 axes[1,0].set_xlabel('Frame lag');axes[1,0].set_ylabel('Noise correlation')
 axes[1,1].set_xlabel('Summed noiseless charge bias (%)');axes[1,1].set_ylabel('Fraction of events / bin')
 bias_values=np.concatenate([100*np.array(data[n][key]) for n in names for key in ('linear_charge_bias_fraction','noiseless_charge_bias_fraction')])
 axes[1,1].set_xlim(float(bias_values.min())-.1,float(bias_values.max())+.1)
 for ax in axes.flat:ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=8)
 fig.suptitle(('CSA v3 corrected pilot' if config.get('reconstruction_method') else 'CSA v3 pilot')+f' · {geometry} · 100 events per sample',fontsize=14);fig.tight_layout()
 for ext in ('png','pdf'):fig.savefig(out/f'diagnostics.{ext}',dpi=170)
 plt.close(fig)
 corrected=config.get('reconstruction_method')=='three_frame_inverse'
 comparison={}
 if corrected:
  with h5py.File(ROOT/config['parent_study_directory']/'summary.h5') as original:
   comparison={name:tree(original[name])['meta'] for name in names}
 rows=[]
 for name in names:
  d=data[name];m=d['meta'];j=2
  rows.append([m['sample'],len(d['event_id']),m['positive_cells'],m['pads'],m['direct_spice_pads'],m['ENC_e_rms'],m['inclusive_mean_snr'],m['event_weighted_mean_snr'],d['mean_compressed_histories_per_event'][j],100*d['fraction_events_with_compression'][j],100*m['mean_compression_deficit_fraction'],100*m['mean_noiseless_charge_bias_fraction'],100*m['mean_noisy_charge_bias_fraction']])
 with (directory/'summary.csv').open('w',newline='') as f:
  writer=csv.writer(f,lineterminator='\n');writer.writerow(['sample','events','positive_cells','pad_histories','direct_spice_histories','ENC_ions_rms','cell_weighted_mean_snr','event_weighted_mean_snr','mean_histories_loss50_per_event','events_loss50_percent','mean_positive_cell_compression_deficit_percent','mean_noiseless_total_charge_bias_percent','mean_noisy_total_charge_bias_percent']);writer.writerows(rows)
 performance=json.loads((directory/'acquisition_benchmark.json').read_text()) if (directory/'acquisition_benchmark.json').exists() else {}
 payload={'groups':data,'calibration':calibration,'validation':validation,'study_id':config['study_id'],'performance':performance,'reconstruction':corrected,'comparison':comparison}
 body=HTML.replace('@@PLOTLY@@',get_plotlyjs()).replace('@@DATA@@',json.dumps(payload,allow_nan=False,separators=(',',':')))
 body=body.replace('3 mm pitch · 30 ms frames',f'{pitch_mm:g} mm pitch · {frame_ms} ms frames')
 if corrected:
  body=body.replace('Ion CSA v3 · response pilot','Ion CSA v3 · corrected response pilot').replace('<span class="badge">S/N measures quality', '<span class="badge">Three-frame inverse reconstruction</span><span class="badge">S/N measures quality')
  body=body.replace("['Frame ENC',", "['Interior frame ENC',")
  body=body.replace('<th>Gain (µV/ion)</th>', '<th>Filtered gain (µV/ion)</th>')
  body=body.replace('The central-frame packet calibration leaves neighbouring-frame contributions in the total; its predicted offset is listed below.', 'The inverse removes the original packet-normalization bias; the comparison table shows the change.')
 if corrected:
  with (directory/'comparison.csv').open('w',newline='') as stream:
   writer=csv.writer(stream,lineterminator='\n');writer.writerow(['sample','ENC_before_ions_rms','ENC_after_ions_rms','mean_snr_before','mean_snr_after','mean_total_bias_before_percent','mean_total_bias_after_percent'])
   for name in names:
    old=comparison[name];new=data[name]['meta'];writer.writerow([new['sample'],old['ENC_e_rms'],new['ENC_e_rms'],old['inclusive_mean_snr'],new['inclusive_mean_snr'],100*old['mean_noiseless_charge_bias_fraction'],100*new['mean_noiseless_charge_bias_fraction']])
 (directory/'ion_csa.html').write_text(body)
 outputs=[{'sample':d['meta']['sample'],'pitch_mm':d['meta']['pitch_mm'],'frame_s':d['meta']['frame_s'],'selected_event_ids':d['event_id'],'path':d['meta']['response_path'],'sha256':d['meta']['response_sha256'],'bytes':Path(d['meta']['response_path']).stat().st_size} for d in data.values()]
 (directory/'outputs.json').write_text(json.dumps({'study_id':config['study_id'],'schema':'nausicaa.ion_csa_response.pilot.v2' if corrected else 'nausicaa.ion_csa_response.pilot.v1','matrix_axes_hdf5':['pad','frame'],'matrix_axes_julia':['frame','pad'],'responses':outputs},indent=2)+'\n')
 table='\n'.join(f"{labels[r[0]]} & {r[1]} & {r[5]:.2f} & {r[6]:.3f} & {r[7]:.3f} & {r[9]:.2f} & {r[10]:.3f} & {r[11]:+.3f} "+r'\\' for r in rows)
 caltable='\n'.join(f"{key[:-2]} & {c['meta']['tau_s']:.4f} & {1e6*c['meta']['gain_V_per_ion']:.3f} & {c['meta']['ENC_e_rms']:.2f} & {c['covariance_e2'][1]/c['covariance_e2'][0]:.4f} "+r'\\' for key,c in sorted(calibration.items()))
 v=validation['meta'];maxref=max(validation['refinement_error_ENC'])
 predicted_bias=100*calibration[frame_key]['meta']['predicted_linear_total_bias_fraction']
 tex=r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm]{geometry}
\usepackage{amsmath,booktabs,graphicx,hyperref}
\title{Ion CSA v3 response pilot}
\date{\today}
\begin{document}\maketitle
The pilot processes 100 deterministically spaced event IDs from each 8000-event
sample at 3 mm pitch and 30 ms frame width. Every pixel with charge retains its
complete event history, including empty frames, two quiet leading frames and a
three-second recovery window. There is no S/N cut. Signed reconstructed counts,
including negative noise fluctuations, are retained.

The circuit is the provisional v3 design, nominal 20 s discharge, with 5 ms
centred voltage averages, 1.25 pF external capacitance and 3 fF feedback. The input
uses uniform mean collection within each frame and the established central
0.5 mm square-pad weighting benchmark in an infinite grounded plane. The pad
weighting potential is one, with zero on the surrounding plane and at infinity,
in homogeneous dielectric. The ion speed is 3 m/s through the finite 5 mm gap.
The charge integral is $e[1-\phi_w(g)]$ per collected ion.

The same transistor circuit was recalibrated at 20, 30 and 40 ms. Noise uses its
one-sided stationary operating-point spectrum from 0.0001 Hz to 1 MHz, filtered
by the 5 ms average and the decay-corrected frame difference. Gaussian samples
retain temporal covariance; pad channels are independent. This stationary noise
approximation is also applied to compressed mean signals.
\begin{center}\begin{tabular}{rrrrr}\toprule
Frame (ms) & Fitted $\tau$ (s) & Gain ($\mu$V/ion) & ENC (ions RMS) & Lag-1 correlation\\\midrule
'''+caltable+r'''
\bottomrule\end{tabular}\end{center}

S/N is the noiseless recovered signal divided by frame ENC, evaluated over all
positive input pixel--frame cells. The cell-weighted and event-weighted means
are shown separately. Compression compares the noiseless transistor response
with the calibrated linear prediction; the flags require a deficit of at least
one frame ENC. A pad history is counted once per event if any positive-input
bucket passes the indicated local-loss criterion.

\begin{center}\scriptsize\begin{tabular}{lrrrrrrr}\toprule
Sample & Events & ENC & Mean S/N & Event mean S/N & Loss50 events (\%) & Deficit (\%) & Total bias (\%)\\\midrule
'''+table+r'''
\bottomrule\end{tabular}\end{center}

The positive-cell compression deficit is
$\sum_{N_n>0}\max(\widehat N_{\rm linear,n}-\widehat N_{\rm CSA,n},0)/\sum N_n$.
The signed total bias is $(\sum\widehat N_{\rm CSA,n}-\sum N_n)/\sum N_n$,
including quiet frames and recovery. The linear prediction's own calibration
bias is saved separately, as is the noisy summed estimate. The isolated-packet
calibration sets the central-frame response to one ion; its earlier and later
frame contributions remain in the complete history. Their saved three-tap
kernel predicts a total linear charge bias of @BIAS@\%. This identifies the
common positive offset in the present scalar estimator. These diagnostics
retain that existing estimator and apply no detection threshold.

\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{diagnostics.pdf}
\caption{Saved pilot quality, compression, covariance and charge-bias diagnostics.
Dashed charge-bias curves use the linear estimator; solid curves include
transistor response where directly simulated.}
\end{figure}

The fast response is checked against control histories spanning predicted
output voltage. Direct SPICE is expanded until remaining control errors fall
below the configured 0.1 ENC. The final voltage guard is @GUARD@ V and the
largest remaining control error is @CONTROL@ ENC. Independent solver refinement
changes the sampled response by at most @REFINE@ ENC. Independent checks cover
finite-path charge, spectral integration, covariance convergence and positive
definiteness, sampled zero-input mean/variance/covariance, input ion ledgers and
Float32 storage errors relative to Float64 calculations.

Response files use \texttt{nausicaa.ion\_csa\_response.pilot.v1}. Each saves the
selected IDs, parent file/hash, exact cards, code identity, stable per-pad seed,
frame covariance, input counts, mean voltages, linear/noiseless/noisy recovered
ions and compression flags. Files are written atomically and marked complete
after the requested pilot ledgers pass.

The companion \texttt{ion\_csa.html} provides sample and compression-severity
selectors, diagnostic distributions and typical/weak/compressed pad histories.
S/N is a quality observable; final analysis determines acceptable performance.
Spatial noise occupancy requires a finite readout area; the present controls
use zero-input histories on the selected channels.

From the repository root:
\begin{verbatim}
work/event-viewer-venv/bin/python scripts/production/run_ion_csa_pilot.py
# Presentation only:
work/event-viewer-venv/bin/python scripts/production/present_ion_csa_pilot.py \
  --results RESULTS_DIRECTORY
\end{verbatim}
\end{document}
'''
 tex=tex.replace('3 mm pitch and 30 ms frame width',f'{pitch_mm:g} mm pitch and {frame_ms} ms frame width')
 tex=tex.replace('@BIAS@',f'{predicted_bias:.3f}').replace('@GUARD@',f"{v['guard_V']:.2f}").replace('@CONTROL@',f"{v['maximum_remaining_control_error_ENC']:.4g}").replace('@REFINE@',f'{maxref:.4g}')
 if corrected:
  start=tex.index("The isolated-packet")
  stop=tex.index("\n\n\\begin{figure}",start)
  tex=tex[:start]+r"The corrected reconstruction solves $M\widehat{\mathbf N}=\mathbf y$ using the saved three-frame response. The same banded solve is applied to the noiseless signal and the existing noise realisation. Its covariance is $M^{-1}C_yM^{-T}$. Quiet boundary frames are retained, and exact history-dependent frame variances enter S/N. This removes the linear normalization bias and temporal mixing while preserving the circuit compression."+tex[stop:]
  tex=tex.replace('nausicaa.ion\\_csa\\_response.pilot.v1','nausicaa.ion\\_csa\\_response.pilot.v2')
  tex=tex.replace('Ion CSA v3 response pilot','Ion CSA v3 corrected response pilot')
  tex=tex.replace('The fast response is checked','The inherited fast response is checked')
  tex=tex.replace('Direct SPICE is expanded until remaining control errors fall', 'The inherited voltage range keeps remaining transformed control errors')
  tex=tex.replace('python scripts/production/run_ion_csa_pilot.py','python scripts/production/run_ion_csa_correction.py')
  tex=tex.replace('Frame (ms) & Fitted $\\tau$ (s) & Gain ($\\mu$V/ion)', 'Frame (ms) & Fitted $\\tau$ (s) & Filtered gain ($\\mu$V/ion)')
  comparison_rows='\n'.join(f"{labels[data[n]['meta']['sample']]} & {comparison[n]['ENC_e_rms']:.2f} & {data[n]['meta']['ENC_e_rms']:.2f} & {100*comparison[n]['mean_noiseless_charge_bias_fraction']:+.3f} & {100*data[n]['meta']['mean_noiseless_charge_bias_fraction']:+.3f} "+r'\\' for n in names)
  comparison_tex=r"\begin{center}\small\begin{tabular}{lrrrr}\toprule Sample & ENC before & ENC after & Bias before (\%) & Bias after (\%)\\\midrule"+'\n'+comparison_rows+r"\bottomrule\end{tabular}\end{center}"
  tex=tex.replace('The same transistor circuit was recalibrated',comparison_tex+'\n\nThe same transistor circuit was recalibrated')
 if performance:
  benchmark=(f"A fresh acquisition benchmark of {performance['fresh_cases']} histories on "
    f"{performance['workers']} workers took {performance['elapsed_s']:.2f} s. "
    f"At that measured rate, the {sum(r[4] for r in rows)} selected transistor histories "
    f"require about {performance['estimated_acquisition_s_for_pilot_selected_histories']:.0f} s of acquisition. "
    "Current generation, reduction, noise sampling and presentation have separate timings.\n\n")
  tex=tex.replace('From the repository root:',benchmark+'From the repository root:')
 (out/'report.tex').write_text(tex)
 (out/'sources.json').write_text(json.dumps({'summary':str((directory/'summary.h5').relative_to(ROOT)),'summary_sha256':hashlib.sha256((directory/'summary.h5').read_bytes()).hexdigest(),'presenter_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
 print(directory/'ion_csa.html');print(out/'report.tex')

HTML=r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Ion CSA v3 · pilot</title><script>@@PLOTLY@@</script><style>
:root{--ink:#203440;--muted:#60727d;--line:#dce4e8;--blue:#176589;--paper:#f3f6f7}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.5 system-ui,sans-serif}main{max-width:1440px;margin:auto;padding:28px}header{margin-bottom:22px}h1{font-size:32px;margin:0 0 8px}h2{font-size:20px;margin:0 0 12px}.sub{color:var(--muted)}.badge{display:inline-block;background:#e1edf2;color:#145b7c;border-radius:16px;padding:4px 11px;margin:9px 6px 0 0;font-size:13px}.panel{background:white;border:1px solid var(--line);border-radius:12px;padding:20px;margin:16px 0}.controls{display:flex;gap:20px;align-items:end;flex-wrap:wrap}label{font-size:13px;color:var(--muted);display:block}select{display:block;padding:9px 12px;min-width:170px;border:1px solid #bdcdd6;border-radius:6px;background:white;color:var(--ink);font:inherit}.cards{display:grid;grid-template-columns:repeat(5,1fr);gap:12px;margin-top:18px}.card{background:#f6f9fa;border-radius:8px;padding:14px}.value{font-size:27px;font-weight:650}.caption{font-size:12px;color:var(--muted)}.grid{display:grid;grid-template-columns:1fr 1fr;gap:16px}.plot{height:330px;min-width:0}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:right;border-bottom:1px solid var(--line);padding:10px}th:first-child,td:first-child{text-align:left}.scroll{overflow:auto}.note{background:#eaf2f6;border-left:3px solid var(--blue);padding:12px 15px;margin-top:14px}#story{height:620px}details{color:var(--muted)}footer{font-size:12px;color:var(--muted);margin:20px 0}@media(max-width:900px){.grid{grid-template-columns:1fr}.cards{grid-template-columns:repeat(2,1fr)}main{padding:15px}}
</style></head><body><main><header><h1>Ion CSA v3 · response pilot</h1><div class="sub">Complete pad histories, transistor compression and correlated frame noise</div><span class="badge">3 mm pitch · 30 ms frames</span><span class="badge">100 events per sample</span><span class="badge">20 s nominal retention · 5 ms averages</span><span class="badge">S/N measures quality · no S/N cuts</span></header>
<section class="panel"><div class="controls"><div><label for="sample">Event sample</label><select id="sample"></select></div><div><label for="severity">Compression diagnostic</label><select id="severity"><option value="0">≥10% local response loss</option><option value="1">≥25% local response loss</option><option value="2" selected>≥50% local response loss</option><option value="3">≥75% local response loss</option></select></div><div><label for="snrscale">S/N view</label><select id="snrscale"><option value="log">Full positive range · log axis</option><option value="low">Low S/N · signed linear axis</option></select></div><div><label for="normalization">Histogram display</label><select id="normalization"><option value="probability">Fraction per bin</option><option value="count">Entries per bin</option></select></div></div><div class="cards" id="cards"></div><div class="note">S/N = noiseless recovered ions / frame ENC. All positive input pixel–frame cells enter the quality distributions. Negative noisy estimates remain in the response files. Compression is measured from the noiseless response, with a deficit of at least one ENC.</div></section>
<section class="panel" id="comparisonPanel" style="display:none"><h2>Before and after reconstruction</h2><div id="comparisonTable" class="scroll"></div><div class="caption">Same circuit, voltage filter, selected events and noise realisation. The inverse unmixes frame contributions and transforms their covariance.</div></section><section class="panel"><h2>Samples and accounting</h2><div class="scroll" id="table"></div></section>
<div class="grid"><section class="panel"><h2>Pixel–frame S/N</h2><div class="plot" id="snr"></div><div class="caption" id="snrcaption"></div></section><section class="panel"><h2>Mean S/N per event</h2><div class="plot" id="eventsnr"></div><div class="caption">Each event contributes one mean over its positive input cells.</div></section><section class="panel"><h2>Compressed pad histories per event</h2><div class="plot" id="stories"></div><div class="caption">One history is one event/pixel over its complete frame sequence.</div></section><section class="panel"><h2>Summed charge bias</h2><div class="plot" id="bias"></div><div class="caption">Signed bias relative to input ions over the saved event and recovery window. Linear calibration bias, transistor response and one noisy realisation are shown separately. The central-frame packet calibration leaves neighbouring-frame contributions in the total; its predicted offset is listed below.</div></section><section class="panel"><h2>Temporal noise covariance</h2><div class="plot" id="cov"></div><div class="caption">Independent pad channels; temporally correlated Gaussian noise from the v3 operating-point spectrum.</div></section><section class="panel"><h2>Zero-input noise control</h2><div class="plot" id="noise"></div><div class="caption">Sampled control histories generated with the same covariance.</div></section></div>
<section class="panel"><h2>Inspect a pad history</h2><div class="controls"><div><label for="example">History</label><select id="example"></select></div><div><label for="window">Time window</label><select id="window"><option value="event">Whole event interval</option><option value="all">Event plus recovery</option></select></div></div><div id="story"></div><div class="caption" id="storycaption"></div></section>
<section class="panel"><h2>Calibration and numerical controls</h2><div id="calibration" class="scroll"></div><p id="performance" class="sub"></p><details><summary>Model, interpretation and provenance</summary><p>The uniform mean ion current uses the finite 5 mm gap above a central 0.5 mm square pad in an infinite grounded plane, with homogeneous dielectric and 3 m/s approach. Each ion contributes e[1−φw(g)]. The transistor circuit uses the selected v3 geometry and biases. Noise is the stationary quiescent-point model, applied consistently to every mean response.</p><p>Direct SPICE replaces the calibrated fast response for histories selected by voltage and for validation controls. The range expands until remaining control differences are below 0.1 ENC. This pilot diagnoses electronics quality; track and blob analysis will determine acceptable performance.</p><pre id="provenance"></pre></details></section><footer>Generated from versioned HDF5 response files and Julia reductions. No detector charge cut or S/N acceptance cut is applied.</footer></main>
<script>const DATA=@@DATA@@;
const labels={bb0nu:'Double beta',single_electrons:'Single electron'};const colors=['#176589','#bd6333','#5f7e43'];const $=id=>document.getElementById(id);let groupKeys=Object.keys(DATA.groups);let exampleKeys=[];
for(const k of groupKeys){let o=document.createElement('option');o.value=k;o.textContent=labels[DATA.groups[k].meta.sample];$('sample').append(o)}
const fmt=(x,n=2)=>Number(x).toLocaleString(undefined,{maximumFractionDigits:n,minimumFractionDigits:n});
function layout(x,y,extra={}){return {margin:{l:63,r:16,t:15,b:57},font:{family:'system-ui',size:12,color:'#203440'},paper_bgcolor:'white',plot_bgcolor:'white',xaxis:{title:{text:x},gridcolor:'#eef1f2'},yaxis:{title:{text:y},gridcolor:'#eef1f2'},legend:{orientation:'h',y:1.16,font:{size:11}},...extra}}
function hist(h,name,color,scale=1){let n=h.count;let den=$('normalization').value==='probability'?n.reduce((a,b)=>a+b,0)+h.meta.underflow+h.meta.overflow:1;let edges=h.edges;return {x:n.map((_,i)=>(edges[i]+edges[i+1])*.5*scale),y:n.map(v=>v/den),name,mode:'lines',line:{color,width:2,shape:'hvh'},type:'scatter'}}
function renderExamples(){let d=DATA.groups[$('sample').value];$('example').innerHTML='';exampleKeys=Object.keys(d.examples).filter(k=>k!=='meta');for(let k of exampleKeys){let m=d.examples[k].meta;let o=document.createElement('option');o.value=k;o.textContent=`${m.label} · event ${m.event_id} · pad (${m.pixel_ix}, ${m.pixel_iy})`;$('example').append(o)}renderStory()}
function renderStory(){let d=DATA.groups[$('sample').value];let e=d.examples[$('example').value];if(!e)return;let t=e.edge_time_s;let x=e.counts.map((_,i)=>(t[i]+t[i+1])/2);let range=$('window').value==='event'?e.event_window_s:[t[0],t[t.length-1]];let traces=[{x,y:e.counts,name:'Input ions',type:'bar',marker:{color:'#176589'},yaxis:'y'},{x:t,y:e.output_at_edges_V.map(v=>1000*v),name:'CSA output',mode:'lines',line:{color:'#744b83'},yaxis:'y2'},{x,y:e.recovered_ions,name:'Recovered + noise',mode:'markers',marker:{color:'#8796a2',size:4},yaxis:'y3'},{x,y:e.noiseless_estimated_ions,name:'Noiseless CSA',mode:'lines',line:{color:'#bd6333',width:2},yaxis:'y3'},{x,y:e.linear_estimated_ions,name:'Linear reference',mode:'lines',line:{color:'#176589',dash:'dot'},yaxis:'y3'}];let l=layout('Event-relative time (s)','');l.xaxis.range=range;l.yaxis={domain:[.73,1],title:{text:'Input ions'}};l.yaxis2={domain:[.39,.65],title:{text:'Vout (mV)'},anchor:'x'};l.yaxis3={domain:[0,.3],title:{text:'Recovered ions'},anchor:'x',zeroline:true,zerolinecolor:'#aabbc5'};l.legend.y=1.09;l.margin.t=35;Plotly.react('story',traces,l,{responsive:true,displaylogo:false});$('storycaption').textContent=`${e.model_used===1?'Direct transistor simulation':'Validated calibrated response'}; original event and signed frame alignment retained. The recovery selector shows the complete saved post-event window.`}
function render(){let d=DATA.groups[$('sample').value],m=d.meta,j=Number($('severity').value),n=d.event_id.length;let cards=[['Frame ENC',fmt(m.ENC_e_rms,2),'ion charges RMS'],['Mean S/N',fmt(m.inclusive_mean_snr,3),'positive cells, cell weighted'],['Affected events',fmt(100*d.fraction_events_with_compression[j],1)+'%','selected local-loss criterion'],['Histories / event',fmt(d.mean_compressed_histories_per_event[j],2),'compressed pad histories'],['Compression deficit',fmt(100*m.mean_compression_deficit_fraction,3)+'%','positive-input-cell diagnostic']];$('cards').innerHTML=cards.map(c=>`<div class="card"><label>${c[0]}</label><div class="value">${c[1]}</div><div class="caption">${c[2]}</div></div>`).join('');let probability=$('normalization').value==='probability',y=probability?'Fraction per bin':'Entries per bin';let low=$('snrscale').value==='low';let h=low?d.snr_low:d.snr;let snrtr=[hist(h,'CSA expected S/N',colors[0])];if(!low)snrtr.unshift(hist(d.input_snr,'Input ions / ENC',colors[1]));let l=layout('S/N',y);if(!low)l.xaxis.type='log';Plotly.react('snr',snrtr,l,{responsive:true,displaylogo:false});$('snrcaption').textContent=`${m.positive_cells.toLocaleString()} positive input cells; ${m.nonpositive_snr_cells} nonpositive expected S/N values. The response files retain all signed estimates.`;let el=layout('Mean S/N',y);el.xaxis.type='log';Plotly.react('eventsnr',[hist(d.event_snr,'Events',colors[0])],el,{responsive:true,displaylogo:false});Plotly.react('stories',[{...hist(d['stories_'+(j+1)],'Events',colors[1]),type:'bar',marker:{color:colors[1]}}],layout('Compressed pad histories',y),{responsive:true,displaylogo:false});Plotly.react('bias',[hist(d.linear_bias,'Linear calibration',colors[0],100),hist(d.mean_bias,'Noiseless CSA',colors[1],100),hist(d.noisy_bias,'With noise',colors[2],100)],layout('Summed charge bias (%)',y),{responsive:true,displaylogo:false});let c=DATA.calibration[Math.round(m.frame_s*1000)+'ms'],cv=c.covariance_e2;Plotly.react('cov',[{x:cv.slice(0,13).map((_,i)=>i),y:cv.slice(0,13).map(v=>v/cv[0]),type:'bar',marker:{color:colors[0]},name:'Noise correlation'}],layout('Frame lag','Correlation'),{responsive:true,displaylogo:false});let nh={edges:c.noise_control_hist_edges_e,count:c.noise_control_hist_count,meta:{underflow:c.meta.noise_control_underflow,overflow:c.meta.noise_control_overflow}};Plotly.react('noise',[hist(nh,'Zero-input controls',colors[2])],layout('Recovered ions (zero input)',y),{responsive:true,displaylogo:false});renderStory()}
let rows=groupKeys.map(k=>{let d=DATA.groups[k],m=d.meta;return `<tr><td>${labels[m.sample]}</td><td>${d.event_id.length}</td><td>${m.positive_cells.toLocaleString()}</td><td>${fmt(m.mean_ions_per_positive_cell,1)}</td><td>${fmt(m.inclusive_mean_snr,3)}</td><td>${fmt(m.event_weighted_mean_snr,3)}</td><td>${m.direct_spice_pads}/${m.pads}</td><td>${fmt(100*m.mean_noiseless_charge_bias_fraction,3)}</td></tr>`}).join('');$('table').innerHTML='<table><tr><th>Sample</th><th>Events</th><th>Positive cells</th><th>Mean ions / cell</th><th>Cell mean S/N</th><th>Event mean S/N</th><th>SPICE / pad histories</th><th>Mean charge bias (%)</th></tr>'+rows+'</table>';
$('calibration').innerHTML='<table><tr><th>Frame</th><th>Fitted tau (s)</th><th>Gain (µV/ion)</th><th>ENC (ions RMS)</th><th>Lag-1 correlation</th><th>Linear total-charge bias (%)</th></tr>'+Object.entries(DATA.calibration).sort().map(([k,c])=>`<tr><td>${k}</td><td>${fmt(c.meta.tau_s,4)}</td><td>${fmt(c.meta.gain_V_per_ion*1e6,3)}</td><td>${fmt(c.meta.ENC_e_rms,2)}</td><td>${fmt(c.covariance_e2[1]/c.covariance_e2[0],4)}</td><td>${fmt(100*c.meta.predicted_linear_total_bias_fraction,3)}</td></tr>`).join('')+'</table>';
if(DATA.reconstruction){$('comparisonPanel').style.display='block';$('comparisonTable').innerHTML='<table><tr><th>Sample</th><th>ENC before → after</th><th>Mean S/N before → after</th><th>Mean charge bias before → after (%)</th></tr>'+groupKeys.map(k=>{let m=DATA.groups[k].meta,p=DATA.comparison[k];return `<tr><td>${labels[m.sample]}</td><td>${fmt(p.ENC_e_rms)} → ${fmt(m.ENC_e_rms)}</td><td>${fmt(p.inclusive_mean_snr,3)} → ${fmt(m.inclusive_mean_snr,3)}</td><td>${fmt(100*p.mean_noiseless_charge_bias_fraction,3)} → ${fmt(100*m.mean_noiseless_charge_bias_fraction,3)}</td></tr>`}).join('')+'</table>';}
if(DATA.performance.fresh_cases)$('performance').textContent=`Fresh LTspice benchmark: ${DATA.performance.fresh_cases} histories on ${DATA.performance.workers} workers in ${fmt(DATA.performance.elapsed_s,2)} s. Estimated acquisition for this pilot’s selected histories: ${fmt(DATA.performance.estimated_acquisition_s_for_pilot_selected_histories,0)} s; other pipeline stages are separate.`;
$('provenance').textContent=JSON.stringify({study_id:DATA.study_id,validation:DATA.validation.meta,solver_refinement_max_error_ENC:Math.max(...DATA.validation.refinement_error_ENC)},null,2);$('sample').onchange=()=>{renderExamples();render()};for(let id of ['severity','snrscale','normalization'])$(id).onchange=render;$('example').onchange=renderStory;$('window').onchange=renderStory;renderExamples();render();
</script></body></html>'''
if __name__=='__main__':main()
