#!/usr/bin/env python3
"""Present event energy, graph components, endpoint blobs and track projections."""
from pathlib import Path
import sys,json,tomllib
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm

def render(path):
 from plotly.offline import get_plotlyjs
 out=path.parent;plots=out/'plots';plots.mkdir(exist_ok=True)
 assets=out/'viewer';assets.mkdir(exist_ok=True)
 (assets/'plotly.min.js').write_text(get_plotlyjs())
 with h5py.File(path) as f:
    dc=tomllib.loads(f['cards/detector_card'][()].decode());rc=tomllib.loads(f['card'][()].decode());tc=rc['track_analysis']
    species=f['species'][()].decode() if 'species' in f else 'bb0nu'
    sample_label={'bb0nu':'bb0nu','single_electrons':'single electrons'}.get(species,species)
    e=f['events'];all_ids=[str(v) for v in e['event_number'][:]];positions={id:j for j,id in enumerate(all_ids)}
    ids=sorted(f['track_examples'],key=lambda id:positions[id]);summary=f['track_summary']
    ion_setting=(f"Fixed itrack voxel width {tc['ion_voxel_width_mm']:g} mm and connection radius {tc['ion_connection_mm']:g} mm."
                 if tc.get('ion_graph_mode')=='fixed' else
                 'Itrack uses Petit’s event-dependent voxel width and connection radius.')
    electron_setting=(f"Etrack voxel width {tc['electron_voxel_width_mm']:g} mm and connection radius {tc['electron_connection_mm']:g} mm."
                      if tc.get('electron_graph_mode')=='fixed' else
                      'Etrack uses Petit’s KDT settings: below 50 cm, voxel width 3σₜ and connection radius 6σₜ; above 50 cm, 1.5σₜ and 2.25σₜ.')
    short_ion_radius=int(np.count_nonzero(f['kdt_connection_radius_mm/itrack'][:]<dc['ion_pitch_mm']))
    short_radius_note=(f'The ion connection radius is below the {dc["ion_pitch_mm"]:g} mm pad pitch in {short_ion_radius} events. ' if short_ion_radius else '')
    blob_settings=(f'Itrack sphere radius: {tc["ion_blob_radius_mm"]:g} mm; etrack sphere radius: {tc["electron_blob_radius_mm"]:g} mm.'
                   if 'ion_blob_radius_mm' in tc and 'electron_blob_radius_mm' in tc else
                   'Itrack sphere radius: 5σₜ; etrack sphere radius: 3σₜ.')
    single_ion=int(summary['itrack_components/counts'][1])
    comparison=''
    earlier=out.parent/'presentation_kdt_v2'
    if tc.get('ion_graph_mode')=='fixed' and (earlier/'controls.h5').exists():
        with h5py.File(earlier/'controls.h5') as old:
            previous=int(old['track_summary/itrack_components/counts'][1])
        comparison=f'<p>Single itrack component: {single_ion}/{len(all_ids)} with the fixed 3 mm / 6 mm settings; <a href="../presentation_kdt_v2/index.html">{previous}/{len(all_ids)} with Petit’s diffusion-scaled settings</a>.</p>'
    prior_cut=out.parent/'presentation_ionflat_components'/'controls.h5'
    if tc['ion_hit_cut_keV']==5 and prior_cut.exists():
        with h5py.File(prior_cut) as old:
            previous=int(np.count_nonzero(old['events/itrack_component_count'][:]>1))
        comparison=f'<p>Multi-component itracks: {previous}/{len(all_ids)} at 3 keV; {len(all_ids)-single_ion}/{len(all_ids)} at 5 keV. Voxel width and connection radius are 3 and 6 mm in both runs.</p>'
    prior_geometry=out.parent/'presentation_ionflat_5keV'/'controls.h5'
    if tc.get('ion_connection_mm')==7 and tc.get('electron_connection_mm')==8 and prior_geometry.exists():
        with h5py.File(prior_geometry) as old:
            pi=int(np.count_nonzero(old['events/itrack_component_count'][:]>1))
            pe=int(np.count_nonzero(old['events/etrack_component_count'][:]>1))
        comparison=f'<p>Before the geometry change: {pi} multi-component itracks and {pe} multi-component etracks. With 3 mm / 7 mm itrack and 4 mm / 8 mm etrack settings: {len(all_ids)-single_ion} and {np.count_nonzero(e["etrack_component_count"][:]>1)}, respectively.</p>'
    if any(id not in f['track_examples'] or 'backbones' not in f['track_examples'][id]['itrack'] for id in ids):
        raise ValueError('Prepare projections and backbones from the saved hit file first.')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for ax,key,title in [(axes[0],'ion_energy_keV','Ion energy estimator'),(axes[1],'electron_energy_keV','Electron energy estimator')]:
        h=summary[key];fit=f['energy_fits'][key];mu=float(fit['mu'][()]);sigma=float(fit['sigma'][()]);edges=h['edges'][:]
        ax.stairs(h['counts'][:],edges,fill=True,alpha=.45,color='#287ba0')
        x=np.linspace(edges[0],edges[-1],400)
        norm_in_window=norm.cdf((fit['upper'][()]-mu)/sigma)-norm.cdf((fit['lower'][()]-mu)/sigma)
        expected=fit['entries'][()]/norm_in_window*np.diff(edges).mean()*norm.pdf(x,mu,sigma)
        ax.plot(x,expected,color='#a63f35',lw=2,label='Gaussian fit')
        ax.set(xlim=(edges[0],edges[-1]),xlabel='Energy (keV)',ylabel='Events / bin',title=title)
        ax.text(.97,.96,f"Fit N = {fit['entries'][()]} / {fit['total_entries'][()]}\nμ = {mu:.2f} keV\nσ = {sigma:.2f} keV\nσ/μ = {100*sigma/mu:.3f}%",ha='right',va='top',transform=ax.transAxes,bbox=dict(facecolor='white',edgecolor='#a8b4bd',alpha=.92))
        ax.legend(loc='upper left');ax.grid(alpha=.15)
    fig.suptitle('Full signed charge · before the hit cut')
    fig.savefig(plots/'event_energy.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for ax,arm in zip(axes,['itrack','etrack']):
        h=summary[arm+'_components'];counts=h['counts'][:];centres=(h['edges'][:-1]+h['edges'][1:])/2
        ax.bar(centres,counts,width=.8,color='#287ba0' if arm=='itrack' else '#b16d31')
        ax.set(xlabel='Connected components',ylabel='Events',title=arm+' · KDT graph')
        ax.text(.98,.97,f"N = {h['entries'][()]}\nμ = {h['mean'][()]:.2f}\nSingle track = {counts[1] if len(counts)>1 else 0}",ha='right',va='top',transform=ax.transAxes)
        ax.grid(alpha=.15)
    fig.savefig(plots/'track_components.png',dpi=150);plt.close(fig)
    failure_note=''
    if 'itrack_component_count' in e:
        components=e['itrack_component_count'][:]
        failure_note=f'<p>{np.count_nonzero(components>1)} events have more than one itrack component.</p>'
        truth=out/'component_mc_classification.h5'
        if truth.exists():
            with h5py.File(truth) as mc:
                clear=int(np.count_nonzero(mc['all_secondary_gamma'][:]))
                broken=int(len(mc['all_secondary_gamma'])-clear)
                broken_ids=', '.join(map(str,mc['event_number'][:][~mc['all_secondary_gamma'][:].astype(bool)]))
            failure_note=f'<p><strong>MC classification at {tc["ion_hit_cut_keV"]:g} keV:</strong> {clear} events have secondary components matching gamma deposits; {broken} have a disconnected primary track'+(f' (MC events {broken_ids})' if broken else '')+f'. A secondary component matches a gamma deposit when at least 80% of its energy lies within {tc["ion_connection_mm"]:g} mm of MC gamma hits and nearer to those hits than to other MC hits. <a href="component_mc_classification.h5">Numerical classification</a>.</p>'+failure_note
        electron_truth=out/'etrack_component_mc_classification.h5'
        if electron_truth.exists():
            with h5py.File(electron_truth) as mc:
                clear_e=int(np.count_nonzero(mc['all_secondary_gamma'][:]))
                broken_e=int(len(mc['all_secondary_gamma'])-clear_e)
            failure_note+=f'<p><strong>Etrack MC classification:</strong> {clear_e} events have secondary components matching gamma deposits; {broken_e} have a disconnected primary track. The matching distance is {tc["electron_connection_mm"]:g} mm. <a href="etrack_component_mc_classification.h5">Numerical classification</a>.</p>'
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for ax,arm in zip(axes,['itrack','etrack']):
        scale=next(iter(f['blob_energies'][arm]));b=f['blob_energies'][arm][scale];e1=b['E1_keV'][:];e2=b['E2_keV'][:];good=np.isfinite(e1)&np.isfinite(e2)
        ax.scatter(e1[good],e2[good],s=9,alpha=.35,color='#287ba0' if arm=='itrack' else '#b16d31')
        radius_label=(f'{tc["ion_blob_radius_mm" if arm=="itrack" else "electron_blob_radius_mm"]:g} mm' if scale=='fixed' else scale.replace('sigma','')+'σₜ')
        ax.set(xlabel='E₁ (keV)',ylabel='E₂ (keV)',title=arm+' · endpoint spheres\nradius '+radius_label)
        ax.text(.98,.97,f"N = {sum(good)}\nμ E₁ = {np.mean(e1[good]):.1f} keV\nμ E₂ = {np.mean(e2[good]):.1f} keV",ha='right',va='top',transform=ax.transAxes,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
        ax.grid(alpha=.15)
    fig.savefig(plots/'blob_energies.png',dpi=150);plt.close(fig)
    def clean(value):
        if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
        if isinstance(value,list):return [clean(v) for v in value]
        if isinstance(value,float) and not np.isfinite(value):return None
        return value
    source=out/f['hit_file'][()].decode()
    with h5py.File(source) as hits_file:
        ion_drift_cm=hits_file['events/ion_drift_cm'][:]
        electron_drift_cm=hits_file['events/electron_drift_cm'][:]
    for id in ids:
        j=positions[id]
        event=f['track_examples'][id]
        payload={'event':id,'ordinal':j+1,'drift_cm':{'itrack':float(ion_drift_cm[j]),'etrack':float(electron_drift_cm[j])},'energy':{a:float(e[key][j]) for a,key in [('itrack','ion_energy_keV'),('etrack','electron_energy_keV')]},'tracks':{}}
        for arm in ['itrack','etrack']:
            g=event[arm]
            payload['tracks'][arm]={k:g[k][()].tolist() for k in ['x_edges_mm','y_edges_mm','z_edges_mm','xy_keV','xz_keV','yz_keV']}
            payload['tracks'][arm]['backbones']=[{k:g['backbones'][c][k][()].tolist() for k in ['x_mm','y_mm','z_mm','hit_count']} for c in sorted(g['backbones'],key=int)]
            payload['tracks'][arm]['blobs']={scale:{k:(v[()].T.tolist() if k=='centres_mm' else v[()].tolist()) for k,v in group.items()} for scale,group in g.get('blobs',{}).items()}
        (assets/f'event_{id}.js').write_text('receiveEvent('+json.dumps(clean(payload),separators=(',',':'),allow_nan=False)+');')
    options=''.join(f'<option value="{id}">{id}</option>' for id in ids)
    panels=''.join(f'<section class="track-panel" id="panel-{arm}"><div id="{arm}-energy"></div><div class="views">'+''.join(f'<div class="plot" id="{arm}-{view}"></div>' for view in ['xy','xz','yz','3d'])+'</div></section>' for arm in ['itrack','etrack'])
    invalid_i=int(np.count_nonzero(e['invalid_ion_buckets'][:]));invalid_e=int(np.count_nonzero(e['invalid_pe_buckets'][:]))
    page='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ITACA · Event energy and tracks</title>
<style>body{font:17px/1.5 system-ui;color:#203343;max-width:1800px;margin:auto;padding:24px}img{max-width:100%}.views{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}.plot{height:440px;min-width:0}.viewer-frame{height:980px;min-height:980px;border:1px solid #ccd6dc;border-radius:8px;padding:12px;box-sizing:border-box;overflow:hidden}.track-panel{height:100%;display:none;grid-template-rows:36px minmax(0,1fr)}.track-panel.active{display:grid}.track-panel .views{min-height:0;grid-template-rows:repeat(2,minmax(0,1fr))}.track-panel .plot{height:100%;min-height:0}.tabs button[aria-selected="true"]{background:#203343;color:white}select,button{font:inherit;padding:7px}nav{background:#fff;padding:12px 0;border-bottom:1px solid #ddd}a{color:#185b83}@media(max-width:750px){body{padding:10px}.plot{font-size:12px}}</style>
<script src="viewer/plotly.min.js"></script></head><body><a href="../../../../../notes/html/index.html">ITACA index</a><h1>bb0nu · Event energy and tracks</h1>
'''+f'''<p>{len(all_ids)} events; {len(ids)} event displays (one per {tc.get('example_every',1)}). TopMetal pitch {dc['ion_pitch_mm']:g} mm; time buckets {dc['ion_bucket_ms']:g} ms. SiPM pitch {dc['EL_pitch_mm']:g} mm; reconstructed time buckets {dc['electron_bucket_us']*rc['electrons']['rebin_factor']:g} μs.</p>
<h2 id="results">Event energy before the hit cut</h2><p>Sum all signed reconstructed charges. Ions include CSA deconvolution; electrons include SiPM and mean attachment corrections.</p><img src="plots/event_energy.png" alt="Ion and electron energy estimates before the hit cut">
<h2>Track components and endpoint blobs</h2><p>Both distributions use all {len(all_ids)} events. {ion_setting} {electron_setting} {short_radius_note}{single_ion} events have one itrack component. E₁ and E₂ sum voxel energies around the endpoints of the largest component, with E₁ ≥ E₂. {blob_settings}</p>{comparison}<img src="plots/track_components.png" alt="Connected component counts for itracks and etracks">{failure_note}<img src="plots/blob_energies.png" alt="Endpoint blob energy correlations for itracks and etracks">
<p>Incomplete estimates are excluded: {invalid_i} ion events and {invalid_e} electron events.</p>
<nav id="viewer"><button id="prev">Previous</button> <label>MC event <select id="event">{options}</select></label> <button id="next">Next</button> <span id="status"></span></nav>
<p>Displayed hits: itrack &gt; {tc.get('ion_hit_cut_keV',0):g} keV; etrack &gt; {tc.get('electron_hit_cut_keV',0):g} keV. Positions in mm. Colour shows summed hit energy, with one colour scale shared by the three projections. Circles and spheres mark endpoint blobs. Lines follow Petit’s tree paths. Spatial limits are fixed around the tracks with 25% extra space. Drag the 3D view to rotate.</p>
<p><label><input id="show-path" type="checkbox" checked> Show connecting paths</label></p><div class="tabs" role="tablist"><button id="tab-itrack" role="tab" aria-selected="true">itrack</button> <button id="tab-etrack" role="tab" aria-selected="false">etrack</button></div><div class="viewer-frame">{panels}</div><p>Data: <a href="controls.h5">graph and plot results</a> · <a href="../hits.h5">reconstructed event hits and energies</a></p>
'''
    page+='''<script>
const selector=document.getElementById('event');let eventScript;
let activeArm='itrack',currentData;
document.getElementById('show-path').onchange=()=>{if(currentData)receiveEvent(currentData);};
function showArm(arm){activeArm=arm;for(const a of ['itrack','etrack']){document.getElementById('panel-'+a).classList.toggle('active',a===arm);document.getElementById('tab-'+a).setAttribute('aria-selected',String(a===arm));}requestAnimationFrame(()=>{for(const v of ['xy','xz','yz','3d']){const p=document.getElementById(arm+'-'+v);if(p.data)Plotly.Plots.resize(p);}});}
for(const arm of ['itrack','etrack'])document.getElementById('tab-'+arm).onclick=()=>showArm(arm);showArm(activeArm);
const config={responsive:true,scrollZoom:false,doubleClick:false,modeBarButtonsToRemove:['zoom2d','pan2d','select2d','lasso2d','zoomIn2d','zoomOut2d','autoScale2d','resetScale2d','zoom3d','pan3d','resetCameraDefault3d','resetCameraLastSave3d'],displaylogo:false};
const centres=e=>e.slice(1).map((v,i)=>(v+e[i])/2);
function loadEvent(){document.getElementById('status').textContent='Loading…';if(eventScript)eventScript.remove();eventScript=document.createElement('script');eventScript.src='viewer/event_'+selector.value+'.js';eventScript.onerror=()=>document.getElementById('status').textContent='Could not read event';document.body.appendChild(eventScript);}
function receiveEvent(data){if(data.event!==selector.value)return;currentData=data;
 document.getElementById('status').textContent='Event '+data.ordinal+' / '+selector.options.length+' · L(ions) = '+data.drift_cm.itrack.toFixed(1)+' cm · L(electrons) = '+data.drift_cm.etrack.toFixed(1)+' cm';
 const ranges={};for(const axis of ['x','y','z']){const values=['itrack','etrack'].flatMap(a=>data.tracks[a][axis+'_edges_mm']);const lo=Math.min(...values),hi=Math.max(...values),margin=Math.max(1,(hi-lo)*.125);ranges[axis]=[lo-margin,hi+margin];}
 for(const arm of ['itrack','etrack']){const g=data.tracks[arm];
 const blob=Object.values(g.blobs||{})[0];const showPath=document.getElementById('show-path').checked;
 const energyMax=Math.max(1,...['xy','xz','yz'].flatMap(v=>g[v+'_keV'].flat()));
 document.getElementById(arm+'-energy').textContent='L(drift) = '+data.drift_cm[arm].toFixed(1)+' cm · Event energy before the hit cut: '+(data.energy[arm]===null?'incomplete':data.energy[arm].toFixed(2)+' keV');
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
 }
}
for(const arm of ['itrack','etrack']){
 const panel=document.getElementById(arm+'-3d');
 panel.addEventListener('wheel',event=>{event.preventDefault();event.stopImmediatePropagation();},{capture:true,passive:false});
 panel.addEventListener('touchmove',event=>{if(event.touches.length>1){event.preventDefault();event.stopImmediatePropagation();}},{capture:true,passive:false});
}
selector.onchange=loadEvent;document.getElementById('prev').onclick=()=>{if(selector.selectedIndex>0){selector.selectedIndex--;loadEvent();}};document.getElementById('next').onclick=()=>{if(selector.selectedIndex<selector.options.length-1){selector.selectedIndex++;loadEvent();}};loadEvent();
</script></body></html>'''
    (out/'index.html').write_text(page.replace('<h1>bb0nu ·',f'<h1>{sample_label} ·'))
 print(out/'index.html')
if __name__=='__main__':render(Path(sys.argv[1]))
