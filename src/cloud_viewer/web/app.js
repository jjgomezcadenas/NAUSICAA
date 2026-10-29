const board=document.getElementById('board');
const input=document.getElementById('event');
const sampleInput=document.getElementById('sample');
const gasInput=document.getElementById('gas');
const status=document.getElementById('status');
const counter=document.getElementById('counter');
const names={ion:'Ions',pure_xe:'PE · pure Xe',xe_he:'PE · Xe + 15% He'};
const sampleNames={bb0nu:'Double beta (0νββ)',electrons:'Single electron'};
let current=null,loading=false,pending=null,eventCounts={};
const config={responsive:true,displaylogo:false,scrollZoom:true,modeBarButtonsToRemove:['select2d','lasso2d']};

function centres(edges){return edges.slice(0,-1).map((v,i)=>(v+edges[i+1])/2)}
function bounds(edges){const lo=edges[0],hi=edges[edges.length-1],mid=(lo+hi)/2,half=Math.max((hi-lo)*.7,1);return[mid-half,mid+half]}
function makeCard(panel,index){
  const card=document.createElement('section');card.className='card';card.id=`card-${index}`;
  const head=document.createElement('div');head.className='card-head';
  const title=document.createElement('div');title.className='card-title';title.textContent=names[panel.carrier];
  const meta=document.createElement('div');meta.className='card-meta';
  meta.textContent=`event ${panel.event_number} · L(ions) = ${panel.ion_drift_cm.toFixed(1)} cm · L(elec) = ${panel.electron_drift_cm.toFixed(1)} cm · ${panel.signal_count.toLocaleString()} ${panel.signal_unit}`;
  head.append(title,meta);card.append(head);
  const grid=document.createElement('div');grid.className='plots';
  for(const plane of ['xy','xz','yz','xyz']){
    const wrap=document.createElement('div');wrap.className='plot-wrap';
    const label=document.createElement('div');label.className='plot-label';
    label.textContent=plane==='xyz'?'3D · drag to rotate':
      plane==='xy'?'XY':`${plane[0].toUpperCase()}–time`;
    const plot=document.createElement('div');plot.className='plot';plot.id=`plot-${index}-${plane}`;
    wrap.append(label,plot);grid.append(wrap);
  }
  card.append(grid);return card;
}
function figure2d(panel,plane){
  const a=plane[0],b=plane[1],x=centres(panel.edges[a]),y=centres(panel.edges[b]);
  const z=panel.projections[plane].map(row=>row.map(v=>v||null));
  const max=Math.max(1,...panel.projections[plane].map(row=>Math.max(...row)));
  return [{type:'heatmap',x,y,z,colorscale:'Viridis',zmin:0,zmax:max,hoverongaps:false,
    colorbar:{thickness:7,len:.67,title:{text:panel.signal_unit,side:'right'}},
    hovertemplate:`${panel.axis_names[a]}: %{x:.2f} ${panel.axis_units[a]}<br>${panel.axis_names[b]}: %{y:.2f} ${panel.axis_units[b]}<br>${panel.signal_unit}: %{z:,}<extra></extra>`}];
}
function figure3d(panel){
  const c=panel.cloud3d,values=c.ions,max=Math.max(1,...values);
  return [{type:'scatter3d',mode:'markers',x:c.x,y:c.y,z:c.z,
    marker:{size:2.6,color:values,colorscale:'Viridis',cmin:0,cmax:max,opacity:.78,
      colorbar:{thickness:7,len:.65,title:{text:panel.signal_unit,side:'right'}}},
    hovertemplate:`x: %{x:.2f} mm<br>y: %{y:.2f} mm<br>${panel.axis_names.z}: %{z:.2f} ${panel.axis_units.z}<br>${panel.signal_unit}: %{marker.color:,}<extra></extra>`}];
}
function layout(panel,plane){
  const base={margin:{l:43,r:39,t:4,b:38},paper_bgcolor:'#fff',plot_bgcolor:'#fbfdff',font:{size:10,color:'#315066'},showlegend:false};
  if(plane==='xyz')return {...base,margin:{l:0,r:22,t:0,b:0},scene:{
    xaxis:{title:'x [mm]',range:bounds(panel.edges.x)},
    yaxis:{title:'y [mm]',range:bounds(panel.edges.y)},
    zaxis:{title:`${panel.axis_names.z} [${panel.axis_units.z}]`,range:bounds(panel.edges.z)},
    aspectmode:'cube',dragmode:'orbit',camera:{eye:{x:1.5,y:1.5,z:1.2}}}};
  const a=plane[0],b=plane[1];
  return {...base,xaxis:{title:`${panel.axis_names[a]} [${panel.axis_units[a]}]`,range:bounds(panel.edges[a]),zeroline:false},
    yaxis:{title:`${panel.axis_names[b]} [${panel.axis_units[b]}]`,range:bounds(panel.edges[b]),zeroline:false}};
}
async function render(payload){
  board.replaceChildren();board.className='board';
  const h=document.createElement('h2');h.className='sample-heading';
  h.textContent=`${sampleNames[payload.sample]} · event ${payload.event} of ${payload.event_count}`;
  board.append(h);
  for(let i=0;i<2;i++){
    const panel=payload.panels[i];
    board.append(makeCard(panel,i));
  }
  for(let i=0;i<2;i++){
    for(const plane of ['xy','xz','yz','xyz']){
      const panel=payload.panels[i];
      await Plotly.newPlot(`plot-${i}-${plane}`,plane==='xyz'?figure3d(panel):figure2d(panel,plane),layout(panel,plane),config);
    }
    status.textContent=`Rendered ${i+1} of 2 panels`;
  }
  status.textContent=`Event ${payload.event} of ${payload.event_count} · ${sampleNames[payload.sample]} · ${names[payload.gas]}`;
}
async function load(event,sample=sampleInput.value,gas=gasInput.value){
  if(loading){pending={event,sample,gas};return}
  loading=true;status.textContent=`Loading event ${event}…`;
  try{
    const query=new URLSearchParams({event,sample,gas});
    const response=await fetch(`/api/event?${query}`);
    const payload=await response.json();
    if(!response.ok)throw Error(payload.error||`HTTP ${response.status}`);
    await render(payload);current=payload.event;input.value=current;
    input.max=payload.event_count;counter.textContent=`of ${payload.event_count}`;
    document.getElementById('previous').disabled=current===1;
    document.getElementById('next').disabled=current===payload.event_count;
    history.replaceState(null,'',`?${new URLSearchParams({sample,gas,event:current})}`);
  }catch(error){status.textContent=`Error: ${error.message}`;console.error(error)}
  finally{loading=false;if(pending){const next=pending;pending=null;load(next.event,next.sample,next.gas)}}
}
document.getElementById('previous').onclick=()=>load(current-1);
document.getElementById('next').onclick=()=>load(current+1);
document.getElementById('go').onclick=()=>load(Number(input.value));
function choose(){
  const event=Math.min(current||1,eventCounts[sampleInput.value]||Infinity);
  load(event);
}
sampleInput.addEventListener('change',choose);
gasInput.addEventListener('change',choose);
input.addEventListener('keydown',e=>{if(e.key==='Enter')load(Number(input.value))});
document.addEventListener('keydown',e=>{if(['INPUT','SELECT','TEXTAREA'].includes(document.activeElement.tagName))return;
  if(e.key==='ArrowLeft'&&current>1)load(current-1);
  if(e.key==='ArrowRight'&&current!==null&&current<(eventCounts[sampleInput.value]||Infinity))load(current+1)});
const query=new URLSearchParams(location.search);
if(['bb0nu','electrons'].includes(query.get('sample')))sampleInput.value=query.get('sample');
if(['pure_xe','xe_he'].includes(query.get('gas')))gasInput.value=query.get('gas');
fetch('/api/info').then(response=>response.json()).then(info=>{
  eventCounts=info.event_counts;
  load(Number(query.get('event')||1));
}).catch(error=>{status.textContent=`Error: ${error.message}`});
