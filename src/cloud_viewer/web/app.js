const board=document.getElementById('board');
const input=document.getElementById('event');
const status=document.getElementById('status');
const counter=document.getElementById('counter');
const names={ion:'Ions',pure_xe:'Electrons · pure Xe',xe_he:'Electrons · Xe + 15% He'};
const sampleNames={bb0nu:'Double beta (0νββ)',electrons:'Single electron'};
const carrierName=panel=>panel.carrier==='ion'?'ions':'electrons';
let current=null,loading=false;
const config={responsive:true,displaylogo:false,scrollZoom:true,modeBarButtonsToRemove:['select2d','lasso2d']};

function centres(edges){return edges.slice(0,-1).map((v,i)=>(v+edges[i+1])/2)}
function bounds(edges){const lo=edges[0],hi=edges[edges.length-1],mid=(lo+hi)/2,half=Math.max((hi-lo)*.7,1);return[mid-half,mid+half]}
function makeCard(panel,index){
  const card=document.createElement('section');card.className='card';card.id=`card-${index}`;
  const head=document.createElement('div');head.className='card-head';
  const title=document.createElement('div');title.className='card-title';title.textContent=names[panel.carrier];
  const meta=document.createElement('div');meta.className='card-meta';
  meta.textContent=`event ${panel.event_id} · L(ions) = ${panel.ion_drift_cm.toFixed(1)} cm · L(elec) = ${panel.electron_drift_cm.toFixed(1)} cm · ${panel.retained_ions.toLocaleString()} ${carrierName(panel)}`;
  head.append(title,meta);card.append(head);
  const grid=document.createElement('div');grid.className='plots';
  for(const plane of ['xy','xz','yz','xyz']){
    const wrap=document.createElement('div');wrap.className='plot-wrap';
    const label=document.createElement('div');label.className='plot-label';label.textContent=plane==='xyz'?'3D · drag to rotate':plane.toUpperCase();
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
    colorbar:{thickness:7,len:.67,title:{text:carrierName(panel),side:'right'}},
    hovertemplate:`${a}: %{x:.2f} mm<br>${b}: %{y:.2f} mm<br>${carrierName(panel)}: %{z:,}<extra></extra>`}];
}
function figure3d(panel){
  const c=panel.cloud3d,values=c.ions,max=Math.max(1,...values);
  return [{type:'scatter3d',mode:'markers',x:c.x,y:c.y,z:c.z,
    marker:{size:2.6,color:values,colorscale:'Viridis',cmin:0,cmax:max,opacity:.78,
      colorbar:{thickness:7,len:.65,title:{text:carrierName(panel),side:'right'}}},
    hovertemplate:`x: %{x:.2f} mm<br>y: %{y:.2f} mm<br>z: %{z:.2f} mm<br>${carrierName(panel)}: %{marker.color:,}<extra></extra>`}];
}
function layout(panel,plane){
  const base={margin:{l:43,r:39,t:4,b:38},paper_bgcolor:'#fff',plot_bgcolor:'#fbfdff',font:{size:10,color:'#315066'},showlegend:false};
  if(plane==='xyz')return {...base,margin:{l:0,r:22,t:0,b:0},scene:{
    xaxis:{title:'x [mm]',range:bounds(panel.edges.x)},
    yaxis:{title:'y [mm]',range:bounds(panel.edges.y)},
    zaxis:{title:'z [mm]',range:bounds(panel.edges.z)},
    aspectmode:'cube',dragmode:'orbit',camera:{eye:{x:1.5,y:1.5,z:1.2}}}};
  const a=plane[0],b=plane[1];
  return {...base,xaxis:{title:`${a} [mm]`,range:bounds(panel.edges[a]),zeroline:false},
    yaxis:{title:`${b} [mm]`,range:bounds(panel.edges[b]),zeroline:false}};
}
async function render(payload){
  board.replaceChildren();board.className='board';
  for(let i=0;i<6;i++){
    const panel=payload.panels[i];
    if(i%3===0){const h=document.createElement('h2');h.className='sample-heading';
      h.textContent=sampleNames[panel.sample]+` · event ${panel.event_id}`;board.append(h)}
    board.append(makeCard(panel,i));
  }
  for(let i=0;i<6;i++){
    for(const plane of ['xy','xz','yz','xyz']){
      const panel=payload.panels[i];
      await Plotly.newPlot(`plot-${i}-${plane}`,plane==='xyz'?figure3d(panel):figure2d(panel,plane),layout(panel,plane),config);
    }
    status.textContent=`Rendered ${i+1} of 6 cloud panels`;
  }
  status.textContent=`Event ${payload.event} · six clouds · 3 projections and 3D per cloud`;
}
async function load(event){
  if(loading)return;
  loading=true;status.textContent=`Loading event ${event}…`;
  try{
    const response=await fetch(`/api/event?event=${encodeURIComponent(event)}`);
    const payload=await response.json();
    if(!response.ok)throw Error(payload.error||`HTTP ${response.status}`);
    await render(payload);current=payload.event;input.value=current;
    input.max=payload.event_count-1;counter.textContent=`of ${payload.event_count-1}`;
    document.getElementById('previous').disabled=current===0;
    document.getElementById('next').disabled=current===payload.event_count-1;
    history.replaceState(null,'',`?event=${current}`);
  }catch(error){status.textContent=`Error: ${error.message}`;console.error(error)}
  finally{loading=false}
}
document.getElementById('previous').onclick=()=>load(current-1);
document.getElementById('next').onclick=()=>load(current+1);
document.getElementById('go').onclick=()=>load(Number(input.value));
input.addEventListener('keydown',e=>{if(e.key==='Enter')load(Number(input.value))});
document.addEventListener('keydown',e=>{if(['INPUT','TEXTAREA'].includes(document.activeElement.tagName))return;
  if(e.key==='ArrowLeft'&&current>0)load(current-1);
  if(e.key==='ArrowRight'&&current!==null)load(current+1)});
load(Number(new URLSearchParams(location.search).get('event')||0));
