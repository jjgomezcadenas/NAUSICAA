/* Event navigation and view state. Numerical maps come from the Julia snapshot. */
'use strict';
const $=id=>document.getElementById(id), planes=['xy','xz','yz','xyz'];
let current=null,busy=false,rendering=false,rotating=false,rotationTimer=null,saveTimer=null;
let selectedComponent=1;
let browseDirectory=null,parentDirectory=null;
const layerNames=['energy','nodes','edges','backbone','endpoints'];
const defaultCamera={eye:{x:1.7,y:1.7,z:1.4},up:{x:0,y:0,z:1},center:{x:0,y:0,z:0}};
const clone=x=>JSON.parse(JSON.stringify(x));
const config={responsive:true,displaylogo:false,scrollZoom:true,modeBarButtonsToRemove:['select2d','lasso2d','toImage'],toImageButtonOptions:{format:'png'}};

function message(text,error=false){$('status').textContent=text;$('status').className=error?'error':busy?'busy':'';}
function setBusy(value){
  busy=value;document.body.classList.toggle('loading',value);
  document.querySelectorAll('button,input,select').forEach(e=>e.disabled=value||!current);
  document.querySelectorAll('[data-file-control]').forEach(e=>e.disabled=value);
  if(current&&!value){$('drift-distance').disabled=$('drift-mode').value!=='fixed';$('nbins').disabled=['transport','pads'].includes(current.active_tab);$('apply-bins').disabled=['transport','pads'].includes(current.active_tab);$('previous').disabled=current.position===1;$('next').disabled=current.position===current.event_count;}
}
async function api(path,data){
  const response=await fetch(path,data===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)});
  const result=await response.json();if(!response.ok)throw Error(result.error||response.statusText);return result;
}
function captureView(){
  if(!current)return {};
  const view={camera:clone($('xyz').layout.scene.camera||defaultCamera),ranges:{},
    stage:current.active_tab,selected_component:selectedComponent,layers:Object.fromEntries(layerNames.map(l=>[l,$('layer-'+l).checked]))};
  for(const p of ['xy','xz','yz'])view.ranges[p]={xaxis:clone($(p).layout.xaxis.range),yaxis:clone($(p).layout.yaxis.range)};
  return view;
}
async function saveView(){
  clearTimeout(saveTimer);
  if(!current||rendering)return;
  try{await api('/api/view',{sample:current.sample,event_id:current.event_id,view:captureView(),opacity:Number($('opacity').value)});}
  catch(e){console.warn('View state:',e.message);}
}
function queueSave(){if(!rendering&&!rotating){clearTimeout(saveTimer);saveTimer=setTimeout(saveView,350);}}

function graphPanel(data){
  const active=data.active_tab==='tracks',g=data.graph;
  document.body.classList.toggle('tracks',active);$('graph-panel').hidden=!active;
  $('stage-title').textContent=data.active_tab==='pads'?'Pad charge':data.active_tab==='transport'?'Transported hits':active?'Monte Carlo tracks':'Monte Carlo energy';
  transportPanel(data);
  document.querySelectorAll('[data-integrated]').forEach(e=>e.textContent=(data.active_tab==='pads'?'Charge':data.active_tab==='transport'?'Ions':'Energy')+' integrated over '+(data.active_tab==='pads'&&e.dataset.integrated==='z'?'time':e.dataset.integrated));
  for(const tab of ['energy','tracks','transport','pads']){
    $('tab-'+tab).classList.toggle('active',data.active_tab===tab);
    $('tab-'+tab).setAttribute('aria-selected',String(data.active_tab===tab));
  }
  selectedComponent=data.view?.selected_component??1;
  if(g.status==='complete'&&!g.components.some(c=>c.component_id===selectedComponent))selectedComponent=g.components[0]?.component_id||0;
  if(data.view?.selected_component===0)selectedComponent=0;
  for(const l of layerNames)$('layer-'+l).checked=data.view?.layers?.[l]??(l!=='nodes');
  $('component-rows').replaceChildren();
  if(!active||g.status!=='complete')return;
  const s=g.settings;
  $('graph-settings').textContent=`Applied: ${s.method} · voxels ${s.voxel_width_mm} mm · distance ${s.resolved_distance_mm} mm · threshold ≥ ${s.energy_threshold_keV} keV`+(s.method==='KDT'?'':` · k ${s.k}`);
  $('graph-summary').textContent=`${g.component_count} components · retained ${g.retained_energy_MeV.toFixed(5)} MeV · below threshold ${g.below_threshold_energy_MeV.toFixed(5)} MeV`+(g.component_count?' · select a component to highlight it in all views':' · no voxels pass the graph threshold');
  for(const c of g.components){
    const row=document.createElement('tr');row.dataset.component=c.component_id;
    const first=document.createElement('td'),button=document.createElement('button'),dot=document.createElement('span');
    dot.className='component-dot';dot.style.background=c.colour;button.append(dot,document.createTextNode(String(c.component_id)));
    button.setAttribute('aria-label',`Highlight component ${c.component_id}`);
    button.onclick=()=>selectComponent(c.component_id);first.append(button);row.append(first);
    const endpoints=c.endpoints_mm.map(p=>'('+p.map(v=>v.toFixed(2)).join(', ')+')').join(' → ')||'—';
    for(const value of [c.energy_MeV.toFixed(5),(100*c.energy_fraction).toFixed(2),c.node_count,c.edge_count,c.length_mm.toFixed(2),endpoints,c.confidence.toFixed(3)]){
      const cell=document.createElement('td');cell.textContent=value;row.append(cell);
    }
    $('component-rows').append(row);
  }
}

function transportPanel(data){
  const active=['transport','pads'].includes(data.active_tab),t=data.transport,c=data.transport_controls;
  $('transport-panel').hidden=!active;document.body.classList.toggle('transported',active);
  $('xz-title').textContent=data.active_tab==='pads'?'X–time':'XZ';$('yz-title').textContent=data.active_tab==='pads'?'Y–time':'YZ';
  $('pad-ledger').hidden=data.active_tab!=='pads';
  if(data.active_tab==='pads'){const p=data.pads;$('pad-ledger').textContent=`Collected: ${p.total_ions.toLocaleString()} ions = ${p.total_charge_fC.toFixed(3)} fC · pitch ${p.pitch_mm.toFixed(2)} mm · slices ${(p.z_bin_width_mm/p.drift_speed_mm_ms).toFixed(1)} ms · grid phase (${p.grid_offset_mm.map(v=>v.toFixed(3)).join(', ')}) mm · full transported cloud`; }
  if(!active)return;
  $('drift-mode').value=c.mode;$('drift-distance').value=c.fixed_cm;
  $('drift-distance').min=c.min_cm;$('drift-distance').max=c.max_cm;
  $('transport-label').textContent=`Transported distance: ${t.drift_cm.toFixed(2)} cm · σT: ${t.sigma_t_mm.toFixed(3)} mm · σL: ${t.sigma_l_mm.toFixed(3)} mm`;
  $('transport-ledger').textContent=`Produced: ${t.sampled_ions.toLocaleString()} ions · retained: ${t.retained_ions.toLocaleString()} · outside histogram: ${t.outside_ions.toLocaleString()} · realisation ${t.realisation} · diffusion grid ${c.nbins}³ · padding ${c.nsigma}σ · thermal Gaussian model`;
}

function styledTraces(data,plane){
  const traces=clone(data.figures[plane].data);
  if(data.active_tab==='tracks'&&data.graph.status==='complete')traces.push(...clone(data.graph.overlays[plane]));
  return traces.map((t,i)=>({...t,...traceStyle(t,i,plane)}));
}
function traceStyle(t,i,plane){
  const active=current.active_tab==='tracks';
  if(i===0)return {visible:!active||$('layer-energy').checked,opacity:plane==='xyz'?current.opacity:1};
  return {visible:active&&$('layer-'+t.meta.layer).checked,
    opacity:selectedComponent===0||t.meta.component_id===selectedComponent?1:.22};
}
async function updateLayers(){
  if(!current)return;
  await Promise.all(planes.map(p=>{
    const styles=$(p).data.map((t,i)=>traceStyle(t,i,p));
    return Plotly.restyle($(p),{visible:styles.map(s=>s.visible),opacity:styles.map(s=>s.opacity)});
  }));
  document.querySelectorAll('#component-rows tr').forEach(r=>r.classList.toggle('selected',Number(r.dataset.component)===selectedComponent));
}
async function selectComponent(id){if(busy)return;selectedComponent=id;await updateLayers();queueSave();}

async function render(data){
  rendering=true;
  const old=current;current=data;
  try{
    graphPanel(data);
    for(const p of planes){
      const f=data.figures[p],layout=clone(f.layout);
      layout.uirevision=data.sample+':'+data.event_id+':'+data.nbins;
      if(p==='xyz')layout.scene.camera=(!data.view?.stage||data.view.stage===data.active_tab)&&data.view?.camera||clone(defaultCamera);
      else if((!data.view?.stage||data.view.stage===data.active_tab)&&data.view?.ranges?.[p]){
        layout.xaxis.range=data.view.ranges[p].xaxis;layout.yaxis.range=data.view.ranges[p].yaxis;
      }
      // Reuse the 3D WebGL context across events and files. Rebuild the 2D
      // heatmaps to refresh their edge caches when the histogram grid changes.
      if(p==='xyz'&&$(p).data){
        await Plotly.react($(p),styledTraces(data,p),layout,config);
      }else{
        Plotly.purge($(p));
        await Plotly.newPlot($(p),styledTraces(data,p),layout,config);
      }
      // Each panel owns its zoom/pan state; persist edits without relaying them.
      $(p).removeAllListeners('plotly_relayout');
      $(p).on('plotly_relayout',queueSave);
    }
    $('sample').replaceChildren(...data.samples.map(key=>{
      const option=document.createElement('option');option.value=key;
      option.textContent=(data.sample_files[key]||key).split('/').pop();option.title=data.sample_files[key]||key;return option;
    }));
    $('sample').value=data.sample;$('event-id').value=data.event_id;$('nbins').value=data.nbins;
    $('file-path').value=data.file_path.startsWith(data.file_directory+'/')?data.file_path.slice(data.file_directory.length+1):data.file_path;
    $('file-path').title=data.file_path;
    $('file-directory').textContent='Default directory: '+data.file_directory;
    browseDirectory=data.file_path.slice(0,data.file_path.lastIndexOf('/'));
    $('energy').textContent=data.energy_MeV.toFixed(5)+' MeV';
    $('hit-count').textContent=data.hit_count.toLocaleString()+' hits · '+data.occupied_cells.toLocaleString()+' occupied cells';
    $('position').textContent=`Event ${data.position.toLocaleString()} of ${data.event_count.toLocaleString()}`;
    $('opacity').value=data.opacity;
    await updateLayers();
    const fmt=x=>Number(x.toPrecision(3)).toString();
    $('grid-info').textContent=`Display cells ${data.bin_counts.x} × ${data.bin_counts.y} × ${data.bin_counts.z} · widths x/y/z: `+
      ['x','y','z'].map(a=>fmt(data.bin_widths[a][0])).join(' / ')+' mm';
    if(data.active_tab==='pads')$('grid-info').textContent=`Pad/time grid ${data.bin_counts.x} × ${data.bin_counts.y} × ${data.bin_counts.z} · x/y pitch ${fmt(data.pads.pitch_mm)} mm · time slices ${fmt(data.bin_widths.z[0])} ms`;
    $('footnote').textContent=data.active_tab==='pads'?'Ideal hole-to-pad assignment · all retained ions · charge in fC · relative arrival time · independent panel zoom':data.active_tab==='transport'?'Source coordinates · transported ions · linear ion-count colours · independent panel zoom':`${data.frame==='source_frame'?'Source coordinates':'Vertex coordinates'} · all MC hits · linear energy colours · independently scaled axes · independent panel zoom`;
    const sampleLabel=data.file_path.split('/').pop();
    document.title=`Event ${data.event_id} · ${sampleLabel} · NAUSICAA`;
  }catch(e){current=old;throw e;}finally{rendering=false;}
}

async function navigate(request){
  if(busy)return;
  stopRotation();await saveView();setBusy(true);message(['transport','pads'].includes(request.tab)||['transport','pads'].includes(current?.active_tab)?'Sampling Gaussian ion transport…':request.action==='file'?'Opening file · preparing the hit table if needed…':request.tab==='tracks'||current?.active_tab==='tracks'?'Preparing MC graph and energy maps…':'Loading event and preparing energy maps…');
  try{
    const data=await api('/api/navigate',request);await render(data);
    message(`Event ${data.event_id} · ${data.elapsed_seconds.toFixed(2)} s reader + presentation · hover for ${data.active_tab==='pads'?'pad charge':data.active_tab==='transport'?'ion counts':'bin energies'}`);
  }catch(e){message(e.message,true);if(current){$('sample').value=current.sample;$('event-id').value=current.event_id;$('nbins').value=current.nbins;}}
  finally{setBusy(false);}
}

async function listFiles(directory){
  $('picker-status').textContent='Reading directory…';
  try{
    const data=await api('/api/files',{directory});browseDirectory=data.directory;parentDirectory=data.parent;
    $('picker-directory').value=data.directory;$('file-list').replaceChildren();
    for(const entry of data.entries){
      const button=document.createElement('button');button.type='button';button.className='file-entry';
      button.textContent=(entry.directory?'Folder · ':'')+entry.name;button.title=entry.path;
      button.onclick=()=>{
        if(entry.directory)listFiles(entry.path);
        else{$('file-path').value=entry.path;$('file-picker').close();navigate({action:'file',path:entry.path});}
      };
      $('file-list').append(button);
    }
    $('picker-status').textContent=data.entries.length?'Choose an HDF5 file. NEXUS files are converted on first use.':'No subdirectories or HDF5 files in this directory.';
  }catch(e){$('picker-status').textContent=e.message;}
}

function stopRotation(){rotating=false;clearTimeout(rotationTimer);$('rotate').textContent='Rotate';}
async function orbit(){
  if(!rotating)return;
  if(!document.hidden){
    const camera=clone($('xyz').layout.scene.camera||defaultCamera),angle=.035;
    const {x,y}=camera.eye;camera.eye.x=x*Math.cos(angle)-y*Math.sin(angle);camera.eye.y=x*Math.sin(angle)+y*Math.cos(angle);
    await Plotly.relayout($('xyz'),{'scene.camera':camera});
  }
  if(rotating)rotationTimer=setTimeout(orbit,65);
}
async function resetViews(){
  stopRotation();rendering=true;
  try{
    await Promise.all(planes.map(p=>p==='xyz'?Plotly.relayout($(p),{'scene.camera':clone(defaultCamera)}):
      Plotly.relayout($(p),{'xaxis.range':current.figures[p].layout.xaxis.range,'yaxis.range':current.figures[p].layout.yaxis.range,'xaxis.autorange':false,'yaxis.autorange':false})));
  }finally{rendering=false;queueSave();}
}

$('next').onclick=()=>navigate({action:'next'});$('previous').onclick=()=>navigate({action:'previous'});
for(const tab of ['energy','tracks','transport','pads'])$('tab-'+tab).onclick=()=>{if(current?.active_tab!==tab)navigate({action:'tab',tab});};
for(const layer of layerNames)$('layer-'+layer).onchange=async()=>{await updateLayers();queueSave();};
$('all-components').onclick=()=>selectComponent(0);
$('transport-form').onsubmit=e=>{e.preventDefault();navigate({action:'transport',drift_mode:$('drift-mode').value,fixed_cm:Number($('drift-distance').value)});};
$('drift-mode').onchange=()=>{$('drift-distance').disabled=$('drift-mode').value!=='fixed';};
$('new-realisation').onclick=()=>navigate({action:'realisation'});
$('sample').onchange=()=>navigate({action:'sample',sample:$('sample').value});
$('file-form').onsubmit=e=>{e.preventDefault();const path=$('file-path').value.trim();if(path)navigate({action:'file',path});};
$('browse-files').onclick=()=>{$('file-picker').showModal();listFiles(browseDirectory);};
$('close-picker').onclick=()=>$('file-picker').close();
$('directory-form').onsubmit=e=>{e.preventDefault();listFiles($('picker-directory').value.trim());};
$('parent-directory').onclick=()=>listFiles(parentDirectory);
$('event-form').onsubmit=e=>{e.preventDefault();navigate({action:'id',event_id:$('event-id').value});};
$('bins-form').onsubmit=e=>{e.preventDefault();navigate({action:'bins',nbins:Number($('nbins').value)});};
$('reset').onclick=resetViews;
$('opacity').oninput=async()=>{if(current){current.opacity=Number($('opacity').value);await Plotly.restyle($('xyz'),{opacity:current.opacity},[0]);queueSave();}};
$('rotate').onclick=()=>{if(rotating){stopRotation();queueSave();}else{rotating=true;$('rotate').textContent='Stop';orbit();}};
document.querySelectorAll('[data-camera]').forEach(button=>button.onclick=async()=>{
  stopRotation();const mode=button.dataset.camera;
  const camera=mode==='xy'?{eye:{x:0,y:0,z:2.5},up:{x:0,y:1,z:0}}:
    mode==='xz'?{eye:{x:0,y:-2.5,z:0},up:{x:0,y:0,z:1}}:{eye:{x:2.5,y:0,z:0},up:{x:0,y:0,z:1}};
  await Plotly.relayout($('xyz'),{'scene.camera':{...camera,center:{x:0,y:0,z:0}}});
});
document.addEventListener('keydown',e=>{
  if(e.key==='Enter'&&!$('file-picker').open&&!['INPUT','SELECT','TEXTAREA','BUTTON'].includes(document.activeElement.tagName)){
    e.preventDefault();if(!busy&&current&&current.position<current.event_count)navigate({action:'next'});
  }
});
window.addEventListener('pagehide',()=>{
  if(current&&!rendering)navigator.sendBeacon('/api/view',new Blob([JSON.stringify({sample:current.sample,event_id:current.event_id,view:captureView(),opacity:Number($('opacity').value)})],{type:'application/json'}));
});
window.addEventListener('error',e=>message('Browser rendering error: '+e.message,true));
window.addEventListener('unhandledrejection',e=>message('Browser rendering error: '+(e.reason?.message||e.reason),true));
(async()=>{
  setBusy(true);
  try{const data=await api('/api/event');await render(data);message(`Event ${data.event_id} · ready · Return advances to the next event`);}
  catch(e){message(e.message,true);}
  finally{setBusy(false);if(current?.auto_rotate){rotating=true;$('rotate').textContent='Stop';orbit();}}
})();
