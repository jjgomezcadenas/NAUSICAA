#!/usr/bin/env python3
"""Present Julia event-charge estimates, positive-hit energies and projections."""
from pathlib import Path
import sys,json,tomllib
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def hist(ax,g,title,unit,log=False):
    ax.stairs(g['counts'][:],g['edges'][:],fill=True,alpha=.65)
    ax.set(title=title,xlabel=unit,ylabel='Entries');ax.grid(alpha=.15)
    if log:ax.set_yscale('log')
    text=f"Entries = {g['entries'][()]}\nμ = {g['mean'][()]:.4g}\nσ = {g['std'][()]:.4g}\nσ/μ = {100*g['std'][()]/g['mean'][()]:.3g}%"
    if g['overflow'][()]:text+=f"\nAbove range = {g['overflow'][()]}"
    ax.text(.97,.96,text,ha='right',va='top',transform=ax.transAxes,fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))

def render(path):
 from plotly.offline import get_plotlyjs
 out=path.parent;plots=out/'plots';plots.mkdir(exist_ok=True)
 assets=out/'viewer';assets.mkdir(exist_ok=True)
 (assets/'plotly.min.js').write_text(get_plotlyjs())
 with h5py.File(path) as f:
    dc=tomllib.loads(f['cards/detector_card'][()].decode());rc=tomllib.loads(f['card'][()].decode());tc=rc['track_analysis']
    e=f['events'];ids=[str(v) for v in e['event_number'][:]];summary=f['track_summary']
    if any(id not in f['track_examples'] or 'backbones' not in f['track_examples'][id]['itrack'] for id in ids):
        raise ValueError('Prepare projections and backbones from the saved hit file first.')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for ax,key,title in [(axes[0],'ion_energy_keV','Ion energy estimator'),(axes[1],'electron_energy_keV','Electron energy estimator')]:
        hist(ax,summary[key],title,'Energy (keV)')
    fig.suptitle('Full signed charge · before the hit cut')
    fig.savefig(plots/'event_energy.png',dpi=150);plt.close(fig)
    def clean(value):
        if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
        if isinstance(value,list):return [clean(v) for v in value]
        if isinstance(value,float) and not np.isfinite(value):return None
        return value
    for j,id in enumerate(ids):
        event=f['track_examples'][id]
        payload={'event':id,'ordinal':j+1,'energy':{a:float(e[key][j]) for a,key in [('itrack','ion_energy_keV'),('etrack','electron_energy_keV')]},'tracks':{},'histograms':{}}
        for arm,histarm in [('itrack','ions'),('etrack','pe')]:
            g=event[arm]
            payload['tracks'][arm]={k:g[k][()].tolist() for k in ['x_edges_mm','y_edges_mm','z_edges_mm','xy_keV','xz_keV','yz_keV']}
            payload['tracks'][arm]['backbones']=[{k:g['backbones'][c][k][()].tolist() for k in ['x_mm','y_mm','z_mm','hit_count']} for c in sorted(g['backbones'],key=int)]
            payload['tracks'][arm]['blobs']={scale:{k:(v[()].T.tolist() if k=='centres_mm' else v[()].tolist()) for k,v in group.items()} for scale,group in g.get('blobs',{}).items()}
            h=f['hit_histograms'][id][histarm]
            payload['histograms'][arm]={k:h[k][()].tolist() for k in ['edges','counts','entries','mean','std']}
        (assets/f'event_{id}.js').write_text('receiveEvent('+json.dumps(clean(payload),separators=(',',':'),allow_nan=False)+');')
    options=''.join(f'<option value="{id}">{id}</option>' for id in ids)
    panels=''.join(f'<section class="track-panel" id="panel-{arm}"><div id="{arm}-energy"></div><div class="views">'+''.join(f'<div class="plot" id="{arm}-{view}"></div>' for view in ['xy','xz','yz','3d'])+f'</div><div class="hist" id="{arm}-hist"></div></section>' for arm in ['itrack','etrack'])
    invalid_i=int(np.count_nonzero(e['invalid_ion_buckets'][:]));invalid_e=int(np.count_nonzero(e['invalid_pe_buckets'][:]))
    page='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ITACA · Event energy and tracks</title>
<style>body{font:17px/1.5 system-ui;color:#203343;max-width:1800px;margin:auto;padding:24px}img{max-width:100%}.views{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}.plot{height:440px;min-width:0}.hist{height:330px}.viewer-frame{height:1250px;min-height:1250px;border:1px solid #ccd6dc;border-radius:8px;padding:12px;box-sizing:border-box;overflow:hidden}.track-panel{height:100%;display:none;grid-template-rows:30px minmax(0,1fr) 200px}.track-panel.active{display:grid}.track-panel .views{min-height:0;grid-template-rows:repeat(2,minmax(0,1fr))}.track-panel .plot,.track-panel .hist{height:100%;min-height:0}.tabs button[aria-selected="true"]{background:#203343;color:white}select,button{font:inherit;padding:7px}nav{background:#fff;padding:12px 0;border-bottom:1px solid #ddd}a{color:#185b83}@media(max-width:750px){body{padding:10px}.plot{font-size:12px}}</style>
<script src="viewer/plotly.min.js"></script></head><body><a href="../../../../../notes/html/index.html">ITACA index</a><h1>bb0nu · Event energy and tracks</h1>
'''+f'''<p>{len(ids)} events. TopMetal pitch {dc['ion_pitch_mm']:g} mm; time buckets {dc['ion_bucket_ms']:g} ms. SiPM pitch {dc['EL_pitch_mm']:g} mm; reconstructed time buckets {dc['electron_bucket_us']*rc['electrons']['rebin_factor']:g} μs.</p>
<h2>Event energy before the hit cut</h2><p>Sum all signed reconstructed charges. Ions include CSA deconvolution; electrons include SiPM and mean attachment corrections.</p><img src="plots/event_energy.png" alt="Ion and electron energy estimates before the hit cut">
<p>Incomplete estimates are excluded: {invalid_i} ion events and {invalid_e} electron events.</p>
<nav><button id="prev">Previous</button> <label>MC event <select id="event">{options}</select></label> <button id="next">Next</button> <span id="status"></span></nav>
<p>Displayed hits: itrack &gt; {tc.get('ion_hit_cut_keV',0):g} keV; etrack &gt; {tc.get('electron_hit_cut_keV',0):g} keV. Positions in mm. Colour shows summed hit energy, with one colour scale shared by the three projections. Circles and spheres mark endpoint blobs. Lines follow Petit’s tree paths. Spatial limits are fixed around the tracks with 25% extra space. Drag the 3D view to rotate.</p>
<p>kNN settings: itrack k = {tc.get('ion_graph_k',tc['graph_k'])}, {tc.get('ion_edge_sigma_scale',tc['edge_sigma_scale']):g}σₜ; etrack k = {tc.get('electron_graph_k',tc['graph_k'])}, {tc['edge_sigma_scale']:g}σₜ.</p>
<p><label>Sphere radius <select id="blob-radius"><option value="2">2σₜ</option><option value="3" selected>3σₜ</option></select></label> <label><input id="show-path" type="checkbox" checked> Show connecting paths</label></p><div class="tabs" role="tablist"><button id="tab-itrack" role="tab" aria-selected="true">itrack</button> <button id="tab-etrack" role="tab" aria-selected="false">etrack</button></div><div class="viewer-frame">{panels}</div><p><a href="controls.h5">Numerical results</a></p>
'''
    page+='''<script>
const selector=document.getElementById('event');let eventScript;
let activeArm='itrack',currentData;
for(const id of ['blob-radius','show-path'])document.getElementById(id).onchange=()=>{if(currentData)receiveEvent(currentData);};
function showArm(arm){activeArm=arm;for(const a of ['itrack','etrack']){document.getElementById('panel-'+a).classList.toggle('active',a===arm);document.getElementById('tab-'+a).setAttribute('aria-selected',String(a===arm));}requestAnimationFrame(()=>{for(const v of ['xy','xz','yz','3d','hist']){const p=document.getElementById(arm+'-'+v);if(p.data)Plotly.Plots.resize(p);}});}
for(const arm of ['itrack','etrack'])document.getElementById('tab-'+arm).onclick=()=>showArm(arm);showArm(activeArm);
const config={responsive:true,scrollZoom:false,doubleClick:false,modeBarButtonsToRemove:['zoom2d','pan2d','select2d','lasso2d','zoomIn2d','zoomOut2d','autoScale2d','resetScale2d','zoom3d','pan3d','resetCameraDefault3d','resetCameraLastSave3d'],displaylogo:false};
const centres=e=>e.slice(1).map((v,i)=>(v+e[i])/2);
function loadEvent(){document.getElementById('status').textContent='Loading…';if(eventScript)eventScript.remove();eventScript=document.createElement('script');eventScript.src='viewer/event_'+selector.value+'.js';eventScript.onerror=()=>document.getElementById('status').textContent='Could not read event';document.body.appendChild(eventScript);}
function receiveEvent(data){if(data.event!==selector.value)return;currentData=data;
 document.getElementById('status').textContent='Event '+data.ordinal+' / '+selector.options.length;
 const ranges={};for(const axis of ['x','y','z']){const values=['itrack','etrack'].flatMap(a=>data.tracks[a][axis+'_edges_mm']);const lo=Math.min(...values),hi=Math.max(...values),margin=Math.max(1,(hi-lo)*.125);ranges[axis]=[lo-margin,hi+margin];}
 for(const arm of ['itrack','etrack']){const g=data.tracks[arm];
 const blob=g.blobs?.[document.getElementById('blob-radius').value];const showPath=document.getElementById('show-path').checked;
 const energyMax=Math.max(1,...['xy','xz','yz'].flatMap(v=>g[v+'_keV'].flat()));
 document.getElementById(arm+'-energy').textContent='Event energy before the hit cut: '+(data.energy[arm]===null?'incomplete':data.energy[arm].toFixed(2)+' keV');
 if(blob)document.getElementById(arm+'-energy').textContent+=' · r = '+blob.radius_mm.toFixed(2)+' mm · E₁ = '+blob.energy_keV[0].toFixed(1)+' keV · E₂ = '+blob.energy_keV[1].toFixed(1)+' keV · overlap = '+blob.overlap_keV.toFixed(1)+' keV';
 for(const [view,a,b] of [['xy','x','y'],['xz','x','z'],['yz','y','z']]){
 const matrix=g[view+'_keV'];
 const lines=g.backbones.filter(p=>showPath&&p.x_mm.length>1).map((p,i)=>({type:'scatter',mode:'lines',x:p[a+'_mm'],y:p[b+'_mm'],line:{color:i===0?'#a63f35':'#697985',width:i===0?1.1:.65},showlegend:false,hoverinfo:'skip'}));
 if(blob)for(let q=0;q<2;q++){const ia='xyz'.indexOf(a),ib='xyz'.indexOf(b),angles=Array.from({length:65},(_,i)=>2*Math.PI*i/64);lines.push({type:'scatter',mode:'lines+markers',x:angles.map(t=>blob.centres_mm[ia][q]+blob.radius_mm*Math.cos(t)),y:angles.map(t=>blob.centres_mm[ib][q]+blob.radius_mm*Math.sin(t)),line:{color:q?'#c07822':'#217b52',width:2},marker:{size:0},showlegend:false,hoverinfo:'skip'});}
 Plotly.react(arm+'-'+view,[{type:'heatmap',x:centres(g[a+'_edges_mm']),y:centres(g[b+'_edges_mm']),z:matrix,zmin:0,zmax:energyMax,colorscale:[[0,'#ffffff'],[.01,'#e1edf3'],[.06,'#b4d2e2'],[.25,'#64a3c7'],[.6,'#327da8'],[1,'#145178']],colorbar:{title:{text:'keV'},thickness:12,len:.8,outlinewidth:0},hovertemplate:a.toUpperCase()+': %{x:.2f} mm<br>'+b.toUpperCase()+': %{y:.2f} mm<br>%{z:.3f} keV<extra></extra>'},...lines],{dragmode:false,plot_bgcolor:'#ffffff',title:arm+' · '+view.toUpperCase()+' · MC '+data.event,margin:{l:65,r:80,t:50,b:60},xaxis:{title:{text:a.toUpperCase()+' (mm)'},range:ranges[a],fixedrange:true,constrain:'domain',showgrid:false,zeroline:false,showline:true,linecolor:'#a8b4bd',ticks:'outside'},yaxis:{title:{text:b.toUpperCase()+' (mm)'},range:ranges[b],fixedrange:true,constrain:'domain',showgrid:false,zeroline:false,showline:true,linecolor:'#a8b4bd',ticks:'outside',scaleanchor:'x',scaleratio:1},uirevision:data.event},config);
 }
 const paths=g.backbones.filter(p=>showPath&&p.x_mm.length>1).map((p,i)=>({type:'scatter3d',mode:'lines',x:p.x_mm,y:p.y_mm,z:p.z_mm,line:{color:i===0?'#a63f35':'#697985',width:i===0?2.3:1.2},showlegend:false,hoverinfo:'skip'}));
 if(blob)for(let q=0;q<2;q++){const x=[],y=[],z=[];for(let i=0;i<=16;i++){const xx=[],yy=[],zz=[],t=Math.PI*i/16;for(let j=0;j<=32;j++){const phi=2*Math.PI*j/32;xx.push(blob.centres_mm[0][q]+blob.radius_mm*Math.sin(t)*Math.cos(phi));yy.push(blob.centres_mm[1][q]+blob.radius_mm*Math.sin(t)*Math.sin(phi));zz.push(blob.centres_mm[2][q]+blob.radius_mm*Math.cos(t));}x.push(xx);y.push(yy);z.push(zz);}const colour=q?'#c07822':'#217b52';paths.push({type:'surface',x,y,z,opacity:.25,colorscale:[[0,colour],[1,colour]],showscale:false,hoverinfo:'skip'});paths.push({type:'scatter3d',mode:'markers',x:[blob.centres_mm[0][q]],y:[blob.centres_mm[1][q]],z:[blob.centres_mm[2][q]],marker:{size:4,color:colour},showlegend:false,hoverinfo:'skip'});}
 Plotly.react(arm+'-3d',paths,{title:arm+' · 3D · MC '+data.event,margin:{l:0,r:0,t:50,b:0},scene:{bgcolor:'#ffffff',dragmode:'orbit',xaxis:{title:'X (mm)',range:ranges.x,showbackground:false,showgrid:false,zeroline:false,showline:true,linecolor:'#a8b4bd',ticks:'outside'},yaxis:{title:'Y (mm)',range:ranges.y,showbackground:false,showgrid:false,zeroline:false,showline:true,linecolor:'#a8b4bd',ticks:'outside'},zaxis:{title:'Z (mm)',range:ranges.z,showbackground:false,showgrid:false,zeroline:false,showline:true,linecolor:'#a8b4bd',ticks:'outside'},aspectmode:'data'},uirevision:data.event},config).then(()=>{
 const panel=document.getElementById(arm+'-3d');
 if(!panel.cameraLock){panel.cameraLock=true;panel.on('plotly_relayout',change=>{
  if(panel.restoringCamera||!Object.keys(change).some(k=>k.startsWith('scene.camera')))return;
  const camera=panel.layout.scene.camera,eye=camera?.eye;if(!eye)return;
  const distance=Math.hypot(eye.x,eye.y,eye.z),fixed=Math.sqrt(3)*1.25;
  const center=camera.center||{x:0,y:0,z:0};
  if(Math.abs(distance-fixed)<1e-6&&Math.hypot(center.x,center.y,center.z)<1e-6)return;
  panel.restoringCamera=true;
  Plotly.relayout(panel,{'scene.camera':{...camera,eye:{x:eye.x*fixed/distance,y:eye.y*fixed/distance,z:eye.z*fixed/distance},center:{x:0,y:0,z:0}}}).finally(()=>panel.restoringCamera=false);
 });}
 });
 const hist=data.histograms[arm];const ratio=hist.mean?100*hist.std/hist.mean:0;
 Plotly.react(arm+'-hist',[{type:'bar',x:centres(hist.edges),y:hist.counts,width:hist.edges.slice(1).map((v,i)=>v-hist.edges[i]),marker:{color:arm==='itrack'?'#286c91':'#c26c20'}}],{title:arm+' · hit energies · MC '+data.event,dragmode:false,xaxis:{fixedrange:true,title:{text:'Hit energy (keV)'}},yaxis:{fixedrange:true,title:{text:'Entries'},type:'log'},margin:{l:65,r:30,t:55,b:60},annotations:[{xref:'paper',yref:'paper',x:1,y:1,showarrow:false,align:'right',bgcolor:'white',text:'N = '+hist.entries+'<br>μ = '+Number(hist.mean).toPrecision(4)+' keV<br>σ = '+Number(hist.std).toPrecision(4)+' keV<br>σ/μ = '+ratio.toPrecision(3)+'%'}],uirevision:data.event},config);
 }
}
for(const arm of ['itrack','etrack']){
 const panel=document.getElementById(arm+'-3d');
 panel.addEventListener('wheel',event=>{event.preventDefault();event.stopImmediatePropagation();},{capture:true,passive:false});
 panel.addEventListener('touchmove',event=>{if(event.touches.length>1){event.preventDefault();event.stopImmediatePropagation();}},{capture:true,passive:false});
}
selector.onchange=loadEvent;document.getElementById('prev').onclick=()=>{if(selector.selectedIndex>0){selector.selectedIndex--;loadEvent();}};document.getElementById('next').onclick=()=>{if(selector.selectedIndex<selector.options.length-1){selector.selectedIndex++;loadEvent();}};loadEvent();
</script></body></html>'''
    (out/'index.html').write_text(page)
 print(out/'index.html')
if __name__=='__main__':render(Path(sys.argv[1]))
