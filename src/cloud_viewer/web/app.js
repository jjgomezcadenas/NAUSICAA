const board=document.getElementById('board');
const input=document.getElementById('event');
const sampleInput=document.getElementById('sample');
const gasInput=document.getElementById('gas');
const status=document.getElementById('status');
const counter=document.getElementById('counter');
const names={ion:'Ions',pure_xe:'PE · pure Xe',xe_he:'PE · Xe + 15% He'};
const sampleNames={bb0nu:'Double beta (0νββ)',electrons:'Single electron'};
let current=null,loading=false,exporting=false,pending=null,eventCounts={},productionInfo=null;
const config={responsive:true,displaylogo:false,scrollZoom:true,modeBarButtonsToRemove:['select2d','lasso2d']};

function centres(edges){return edges.slice(0,-1).map((v,i)=>(v+edges[i+1])/2)}
function bounds(edges){const lo=edges[0],hi=edges[edges.length-1],mid=(lo+hi)/2,half=(hi-lo)*.625;return[mid-half,mid+half]}
function planeLabel(panel,plane){return plane==='xyz'?'3D':plane==='xy'?'XY':panel.axis_units.z==='mm'?plane.toUpperCase():`${plane[0].toUpperCase()}–time`}
function makeCard(panel,index){
  const card=document.createElement('section');card.className='card';card.id=`card-${index}`;
  const head=document.createElement('div');head.className='card-head';
  const title=document.createElement('div');title.className='card-title';title.textContent=panel.title||names[panel.carrier];
  const meta=document.createElement('div');meta.className='card-meta';
  meta.textContent=`event ${panel.event_number} · L(ions) = ${panel.ion_drift_cm.toFixed(1)} cm · L(elec) = ${panel.electron_drift_cm.toFixed(1)} cm · ${panel.signal_count.toLocaleString()} ${panel.signal_unit}`;
  if(panel.source_event_id!==undefined)meta.textContent+=` · MC ID ${panel.source_event_id}`;
  if(panel.clipped_buckets!==undefined)meta.textContent+=panel.clipped_buckets?
    ` · electronics clipping: ${panel.clipped_buckets} buckets, ${panel.clip_loss_pe.toFixed(1)} PE lost`:
    ' · no electronics clipping';
  head.append(title,meta);card.append(head);
  const grid=document.createElement('div');grid.className='plots';
  for(const plane of ['xy','xz','yz','xyz']){
    const wrap=document.createElement('div');wrap.className='plot-wrap';
    const label=document.createElement('div');label.className='plot-label';
    label.textContent=plane==='xyz'?'3D · drag to rotate':planeLabel(panel,plane);
    const plot=document.createElement('div');plot.className='plot';plot.id=`plot-${index}-${plane}`;
    wrap.append(label,plot);grid.append(wrap);
  }
  card.append(grid);return card;
}
function figure2d(panel,plane){
  const a=plane[0],b=plane[1],x=centres(panel.edges[a]),y=centres(panel.edges[b]);
  const z=panel.projections[plane].map(row=>row.map(v=>v||null));
  const max=Math.max(1,...panel.projections[plane].map(row=>Math.max(...row.map(Math.abs))));
  return [{type:'heatmap',x,y,z,colorscale:panel.signed?'RdBu':'Viridis',zmin:panel.signed?-max:0,zmax:max,hoverongaps:false,
    colorbar:{thickness:7,len:.67,title:{text:panel.signal_unit,side:'right'}},
    hovertemplate:`${panel.axis_names[a]}: %{x:.2f} ${panel.axis_units[a]}<br>${panel.axis_names[b]}: %{y:.2f} ${panel.axis_units[b]}<br>${panel.signal_unit}: %{z:,}<extra></extra>`}];
}
function figure3d(panel){
  const c=panel.cloud3d,values=c.ions,max=values.reduce((a,v)=>Math.max(a,Math.abs(v)),1);
  return [{type:'scatter3d',mode:'markers',x:c.x,y:c.y,z:c.z,
    customdata:c.arrival_time,
    marker:{size:2.6,color:values,colorscale:panel.signed?'RdBu':'Viridis',cmin:panel.signed?-max:0,cmax:max,opacity:.78,
      colorbar:{thickness:7,len:.65,title:{text:panel.signal_unit,side:'right'}}},
    hovertemplate:`x: %{x:.2f} mm<br>y: %{y:.2f} mm<br>${panel.axis_names.z}: %{z:.2f} ${panel.axis_units.z}<br>`+
      (c.arrival_time?`relative arrival time: %{customdata:.2f} ${panel.arrival_time_unit}<br>`:'')+
      `${panel.signal_unit}: %{marker.color:,}<extra></extra>`}];
}
function layout(panel,plane){
  const ranges=panel.view_ranges||Object.fromEntries(Object.entries(panel.edges).map(([axis,edges])=>[axis,bounds(edges)]));
  const spatial=panel.axis_units.z==='mm';
  const base={margin:{l:43,r:39,t:4,b:38},paper_bgcolor:'#fff',plot_bgcolor:'#fbfdff',font:{size:10,color:'#315066'},showlegend:false};
  if(plane==='xyz')return {...base,margin:{l:0,r:22,t:0,b:0},scene:{
    xaxis:{title:{text:'x [mm]'},range:ranges.x},
    yaxis:{title:{text:'y [mm]'},range:ranges.y},
    zaxis:{title:{text:`${panel.axis_names.z} [${panel.axis_units.z}]`},range:ranges.z},
    aspectmode:spatial?'data':'cube',dragmode:'orbit',camera:{eye:{x:1.5,y:1.5,z:1.2}}}};
  const a=plane[0],b=plane[1];
  const equal=spatial||plane==='xy';
  return {...base,xaxis:{title:{text:`${panel.axis_names[a]} [${panel.axis_units[a]}]`},range:ranges[a],zeroline:false,constrain:'domain'},
    yaxis:{title:{text:`${panel.axis_names[b]} [${panel.axis_units[b]}]`},range:ranges[b],zeroline:false,constrain:'domain',...(equal?{scaleanchor:'x',scaleratio:1}:{})}};
}
async function render(payload){
  board.replaceChildren();board.className='board';
  const h=document.createElement('h2');h.className='sample-heading';
  h.textContent=`${payload.label?payload.label+' · ':''}${sampleNames[payload.sample]} · event ${payload.event} of ${payload.event_count}`;
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
  status.textContent=`Event ${payload.event} of ${payload.event_count} · ${sampleNames[payload.sample]} · ${payload.panels[1].signal_unit} · ${payload.gas==='pure_xe'?'pure Xe':'Xe + 15% He'}`;
}
async function load(event,sample=sampleInput.value,gas=gasInput.value){
  if(exporting)return;
  if(loading){pending={event,sample,gas};return}
  loading=true;status.textContent=`Loading event ${event}…`;
  try{
    const query=new URLSearchParams({event,sample,gas});
    for(const key of ["pitch","stage","electron_stage","operating_point"])query.set(key,document.getElementById(key).value);
    const response=await fetch(`/api/event?${query}`);
    const payload=await response.json();
    if(!response.ok)throw Error(payload.error||`HTTP ${response.status}`);
    await render(payload);current=payload.event;input.value=current;
    input.max=payload.event_count;counter.textContent=`of ${payload.event_count}`;
    document.getElementById('previous').disabled=current===1;
    document.getElementById('next').disabled=current===payload.event_count;
    history.replaceState(null,'',`?${query}`);
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
for(const key of ['pitch','stage','electron_stage'])document.getElementById(key).addEventListener('change',choose);
function fillChoices(element, values, preferred, labels={}){
  element.replaceChildren(...values.map(v=>new Option(labels[v]||v,v)));
  element.value=values.includes(preferred)?preferred:values[0];
  element.disabled=values.length===1;
}
function configureBranches(info, preferred={}){
  eventCounts=info.event_counts;
  fillChoices(sampleInput,info.samples,preferred.sample||info.samples[0],sampleNames);
  fillChoices(gasInput,info.gases,preferred.gas||info.default_gas,{pure_xe:'Pure Xe',xe_he:'Xe + 15% He'});
  if(info.pitches)fillChoices(document.getElementById('pitch'),info.pitches,preferred.pitch||info.default_pitch);
}
document.getElementById('operating_point').addEventListener('change',()=>{
  const info=productionInfo.operating_points.find(p=>p.id===document.getElementById('operating_point').value);
  configureBranches({...productionInfo,...info},{sample:sampleInput.value,gas:gasInput.value,pitch:document.getElementById('pitch').value});
  choose();
});
input.addEventListener('keydown',e=>{if(e.key==='Enter')load(Number(input.value))});
document.addEventListener('keydown',e=>{if(['INPUT','SELECT','TEXTAREA'].includes(document.activeElement.tagName))return;
  if(e.key==='ArrowLeft'&&current>1)load(current-1);
  if(e.key==='ArrowRight'&&current!==null&&current<(eventCounts[sampleInput.value]||Infinity))load(current+1)});
const query=new URLSearchParams(location.search);
if(['bb0nu','electrons'].includes(query.get('sample')))sampleInput.value=query.get('sample');
if(['pure_xe','xe_he'].includes(query.get('gas')))gasInput.value=query.get('gas');
fetch('/api/info').then(response=>response.json()).then(info=>{
  productionInfo=info;
  const points=document.getElementById('operating_point');
  if(info.operating_points?.length){
    fillChoices(points,info.operating_points.map(p=>p.id),query.get('operating_point')||info.default_operating_point,
      Object.fromEntries(info.operating_points.map(p=>[p.id,p.label])));
    configureBranches({...info,...info.operating_points.find(p=>p.id===points.value)},Object.fromEntries(query));
  }else configureBranches(info,Object.fromEntries(query));
  document.getElementById("production-controls").hidden=!info.production;
  if(info.production)document.getElementById("footnote").textContent="Source coordinates in mm; both drift directions mapped to the same z axis. Shared initial frame with 25% padding and equal spatial scale. Each plot zooms independently; signed signals retained.";
  for(const key of ['stage','electron_stage']){
    const select=document.getElementById(key);
    if([...select.options].some(o=>o.value===query.get(key)))select.value=query.get(key);
  }
  load(Number(query.get('event')||1));
}).catch(error=>{status.textContent=`Error: ${error.message}`});

// Export current Plotly figures, including their zoom ranges and 3D camera.
function pngFilename(){
  const point=document.getElementById('operating_point');
  const label=point.selectedOptions[0]?.text||point.value;
  const configuration=label.replace(/E\/P\s*/gi,'EP').replace(/(\d)\s+(mm|cm)/g,'$1$2')
    .toLowerCase().replace(/[^a-z0-9]+/g,'_').replace(/^_|_$/g,'');
  const event=String(current).padStart(5,'0');
  return `ITACA_${sampleInput.value}_${configuration}_ion${document.getElementById('pitch').value}mm_${gasInput.value}_event${event}_${document.getElementById('stage').value}_${document.getElementById('electron_stage').value}.png`;
}
document.getElementById('produce-png').onclick=async()=>{
  const button=document.getElementById('produce-png');
  if(loading||current===null)return;
  exporting=true;
  const controls=[...document.querySelectorAll('.controls button,.controls select,.controls input')];
  const disabled=controls.map(c=>c.disabled);controls.forEach(c=>c.disabled=true);
  try{
    const canvas=document.createElement('canvas');canvas.width=2400;canvas.height=1530;
    const ctx=canvas.getContext('2d');ctx.fillStyle='#ffffff';ctx.fillRect(0,0,2400,1530);
    ctx.fillStyle='#17334a';ctx.font='bold 28px sans-serif';
    ctx.fillText(`ITACA · ${sampleNames[sampleInput.value]} · event ${current}`,30,42);
    ctx.font='20px sans-serif';
    ctx.fillText(document.getElementById('operating_point').selectedOptions[0].text,30,78);
    for(let i=0;i<2;i++){
      const x=30+i*1200,card=document.getElementById(`card-${i}`);
      ctx.font='bold 21px sans-serif';ctx.fillText(card.querySelector('.card-title').textContent,x,118);
      ctx.font='15px sans-serif';ctx.fillText(card.querySelector('.card-meta').textContent,x,148,1140);
      for(const [j,plane] of ['xy','xz','yz','xyz'].entries()){
        const px=x+(j%2)*580,py=185+Math.floor(j/2)*600;
        const div=document.getElementById(`plot-${i}-${plane}`);
        const uri=await Plotly.toImage(div,{format:'png',width:565,height:550,scale:1});
        const img=new Image();await new Promise((ok,fail)=>{img.onload=ok;img.onerror=fail;img.src=uri;});
        ctx.font='bold 18px sans-serif';ctx.fillText(plane==='xyz'?'3D':card.querySelectorAll('.plot-label')[j].textContent,px,py);
        ctx.drawImage(img,px,py+10);
      }
    }
    ctx.font='17px sans-serif';ctx.fillText(document.getElementById('footnote').textContent,30,1490);
    const link=document.createElement('a');link.download=pngFilename();
    link.href=canvas.toDataURL('image/png');link.click();
    status.textContent=`PNG produced: ${link.download}`;
  }catch(error){status.textContent=`PNG export failed: ${error.message}`;console.error(error)}
  finally{controls.forEach((c,i)=>c.disabled=disabled[i]);exporting=false}
};
