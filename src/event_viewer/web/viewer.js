/* Event navigation and view state. Numerical maps come from the Julia snapshot. */
'use strict';
const $=id=>document.getElementById(id), planes=['xy','xz','yz','xyz'];
let current=null,busy=false,rendering=false,rotating=false,rotationTimer=null,saveTimer=null;
const defaultCamera={eye:{x:1.7,y:1.7,z:1.4},up:{x:0,y:0,z:1},center:{x:0,y:0,z:0}};
const clone=x=>JSON.parse(JSON.stringify(x));
const config={responsive:true,displaylogo:false,scrollZoom:true,modeBarButtonsToRemove:['select2d','lasso2d','toImage'],toImageButtonOptions:{format:'png'}};

function message(text,error=false){$('status').textContent=text;$('status').className=error?'error':busy?'busy':'';}
function setBusy(value){
  busy=value;document.body.classList.toggle('loading',value);
  document.querySelectorAll('button,input,select').forEach(e=>e.disabled=value||!current);
  if(current&&!value){$('previous').disabled=current.position===1;$('next').disabled=current.position===current.event_count;}
}
async function api(path,data){
  const response=await fetch(path,data===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)});
  const result=await response.json();if(!response.ok)throw Error(result.error||response.statusText);return result;
}
function captureView(){
  if(!current)return {};
  const view={camera:clone($('xyz').layout.scene.camera||defaultCamera),ranges:{}};
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

async function render(data){
  rendering=true;
  const old=current;current=data;
  try{
    for(const p of planes){
      const f=data.figures[p],layout=clone(f.layout);
      layout.uirevision=data.sample+':'+data.event_id+':'+data.nbins;
      if(p==='xyz')layout.scene.camera=data.view?.camera||clone(defaultCamera);
      else if(data.view?.ranges?.[p]){
        layout.xaxis.range=data.view.ranges[p].xaxis;layout.yaxis.range=data.view.ranges[p].yaxis;
      }
      // Recreate traces when their event/grid changes; release the old WebGL
      // scene and heatmap edge caches. Camera/layer edits use restyle/relayout.
      Plotly.purge($(p));
      await Plotly.newPlot($(p),clone(f.data),layout,config);
      // Each panel owns its zoom/pan state; persist edits without relaying them.
      $(p).on('plotly_relayout',queueSave);
    }
    $('sample').value=data.sample;$('event-id').value=data.event_id;$('nbins').value=data.nbins;
    $('energy').textContent=data.energy_MeV.toFixed(5)+' MeV';
    $('hit-count').textContent=data.hit_count.toLocaleString()+' hits · '+data.occupied_cells.toLocaleString()+' occupied cells';
    $('position').textContent=`Event ${data.position.toLocaleString()} of ${data.event_count.toLocaleString()}`;
    $('opacity').value=data.opacity;
    await Plotly.restyle($('xyz'),{opacity:data.opacity});
    const fmt=x=>Number(x.toPrecision(3)).toString();
    $('grid-info').textContent=`Display cells ${data.bin_counts.x} × ${data.bin_counts.y} × ${data.bin_counts.z} · widths x/y/z: `+
      ['x','y','z'].map(a=>fmt(data.bin_widths_mm[a][0])).join(' / ')+' mm';
    $('footnote').textContent=`${data.frame==='source_frame'?'Source coordinates':'Vertex coordinates'} · all MC hits · linear energy colours · independently scaled axes · independent panel zoom`;
    document.title=`Event ${data.event_id} · ${data.sample==='bb0nu'?'0νββ':'electrons'} · NAUSICAA`;
  }catch(e){current=old;throw e;}finally{rendering=false;}
}

async function navigate(request){
  if(busy)return;
  stopRotation();await saveView();setBusy(true);message('Loading event and preparing energy maps…');
  try{
    const data=await api('/api/navigate',request);await render(data);
    message(`Event ${data.event_id} · ${data.elapsed_seconds.toFixed(2)} s reader + presentation · hover for bin energies`);
  }catch(e){message(e.message,true);if(current){$('sample').value=current.sample;$('event-id').value=current.event_id;$('nbins').value=current.nbins;}}
  finally{setBusy(false);}
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
$('sample').onchange=()=>navigate({action:'sample',sample:$('sample').value});
$('event-form').onsubmit=e=>{e.preventDefault();navigate({action:'id',event_id:$('event-id').value});};
$('bins-form').onsubmit=e=>{e.preventDefault();navigate({action:'bins',nbins:Number($('nbins').value)});};
$('reset').onclick=resetViews;
$('opacity').oninput=async()=>{if(current){await Plotly.restyle($('xyz'),{opacity:Number($('opacity').value)});queueSave();}};
$('rotate').onclick=()=>{if(rotating){stopRotation();queueSave();}else{rotating=true;$('rotate').textContent='Stop';orbit();}};
document.querySelectorAll('[data-camera]').forEach(button=>button.onclick=async()=>{
  stopRotation();const mode=button.dataset.camera;
  const camera=mode==='xy'?{eye:{x:0,y:0,z:2.5},up:{x:0,y:1,z:0}}:
    mode==='xz'?{eye:{x:0,y:-2.5,z:0},up:{x:0,y:0,z:1}}:{eye:{x:2.5,y:0,z:0},up:{x:0,y:0,z:1}};
  await Plotly.relayout($('xyz'),{'scene.camera':{...camera,center:{x:0,y:0,z:0}}});
});
document.addEventListener('keydown',e=>{
  if(e.key==='Enter'&&!['INPUT','SELECT','TEXTAREA'].includes(document.activeElement.tagName)){
    e.preventDefault();if(!busy&&current&&current.position<current.event_count)navigate({action:'next'});
  }
});
window.addEventListener('pagehide',()=>{
  if(current&&!rendering)navigator.sendBeacon('/api/view',new Blob([JSON.stringify({sample:current.sample,event_id:current.event_id,view:captureView(),opacity:Number($('opacity').value)})],{type:'application/json'}));
});
(async()=>{
  setBusy(true);
  try{const data=await api('/api/event');await render(data);message(`Event ${data.event_id} · ready · Return advances to the next event`);}
  catch(e){message(e.message,true);}
  finally{setBusy(false);if(current?.auto_rotate){rotating=true;$('rotate').textContent='Stop';orbit();}}
})();
