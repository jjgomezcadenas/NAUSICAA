// Exercise event switching and Plotly inputs without a browser rendering engine.
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const root=process.argv[2],html=fs.readFileSync(root+'/index.html','utf8');
const nav=html.match(/<nav\b[^>]*>([\s\S]*?)<\/nav>/);assert(nav);
const ids=[...nav[1].matchAll(/<option value="(\d+)">/g)].map(m=>m[1]);assert(ids.length>0);
const elements={};function el(id){return elements[id]??={value:ids[0],selectedIndex:0,options:ids,classList:{toggle(){}},setAttribute(){},addEventListener(){},on(){}};}
const plots={};const ctx={requestAnimationFrame:f=>f(),document:{getElementById:el,createElement:()=>({remove(){}}),body:{appendChild(){}}},Plotly:{Plots:{resize(){}},react:(id,data,layout,config)=>{plots[id]={data,layout,config};el(id).layout=layout;el(id).data=data;return Promise.resolve();},relayout:()=>Promise.resolve()}};vm.createContext(ctx);
el('show-path').checked=true;
vm.runInContext(html.match(/<script>([\s\S]*)<\/script>/)[1],ctx);
for(const id of [ids[0],ids[Math.floor(ids.length/2)],ids.at(-1)]){
 el('event').value=id;vm.runInContext(fs.readFileSync(root+'/viewer/event_'+id+'.js','utf8'),ctx);
 assert.equal(Object.keys(plots).length,8);
 assert(el('status').textContent.includes('L(ions) ='));
 assert(el('status').textContent.includes('L(electrons) ='));
 for(const arm of ['itrack','etrack']){
  for(const view of ['xy','xz','yz']){
   const p=plots[arm+'-'+view],h=p.data[0];assert(p.layout.title.includes(id));assert.equal(p.config.scrollZoom,false);
   for(const axis of ['xaxis','yaxis']){assert.equal(p.layout[axis].fixedrange,true);assert.equal(p.layout[axis].constrain,'domain');assert(p.layout[axis].title.text.includes('(mm)'));}
   assert.equal(h.z.length,h.y.length);assert(h.z.every(row=>row.length===h.x.length));
   assert(p.data.slice(1).every(t=>['lines','lines+markers'].includes(t.mode)));assert(h.z.every(row=>row.every(Number.isFinite)));
  }
  const p=plots[arm+'-3d'];assert.equal(p.data.filter(t=>t.type==='surface').length,2);for(const sphere of p.data.filter(t=>t.type==='surface'))for(const a of ['x','y','z'])assert(sphere[a].flat().every(Number.isFinite));assert.equal(p.config.scrollZoom,false);
  for(const axis of ['xaxis','yaxis','zaxis']){assert.equal(p.layout.scene[axis].showbackground,false);assert.equal(p.layout.scene[axis].showgrid,false);}
  assert(!plots[arm+'-hist']);
 }
}
el('show-path').checked=false;el('show-path').onchange();assert(!plots['itrack-3d'].data.some(t=>t.mode==='lines'));
el('tab-etrack').onclick();el('tab-itrack').onclick();
// Zoom to the main track: every view takes the saved zoom cube, which contains the main-track boxes; a second press restores the full frame.
const axisRanges=()=>Object.fromEntries(['itrack','etrack'].flatMap(arm=>[
 ...[['xy','x','y'],['xz','x','z'],['yz','y','z']].map(([v,a,b])=>[arm+'-'+v,{[a]:plots[arm+'-'+v].layout.xaxis.range,[b]:plots[arm+'-'+v].layout.yaxis.range}]),
 [arm+'-3d',{x:plots[arm+'-3d'].layout.scene.xaxis.range,y:plots[arm+'-3d'].layout.scene.yaxis.range,z:plots[arm+'-3d'].layout.scene.zaxis.range}]]));
const close=(u,v)=>Math.abs(u-v)<1e-2;   // Float32 positions: 0.01 mm
let zoomChecked=0;
for(const id of ids){
 const text=fs.readFileSync(root+'/viewer/event_'+id+'.js','utf8'),data=JSON.parse(text.slice('receiveEvent('.length,-2));
 if(!data.zoom_mm)continue;
 el('event').value=id;vm.runInContext(text,ctx);
 const full=axisRanges();assert.equal(el('zoom').textContent||'Zoom to track','Zoom to track');
 el('zoom').onclick();assert.equal(el('zoom').textContent,'Full frame');
 const zoom=axisRanges(),side=data.zoom_mm[1][0]-data.zoom_mm[0][0];
 for(const [i,axis] of ['x','y','z'].entries()){
  assert(close(data.zoom_mm[1][i]-data.zoom_mm[0][i],side),'zoom region is a cube');
  for(const r of Object.values(zoom))if(r[axis]){assert(close(r[axis][0],data.zoom_mm[0][i])&&close(r[axis][1],data.zoom_mm[1][i]),'all views take the zoom range on '+axis);}
  for(const arm of ['itrack','etrack']){const box=data.tracks[arm].main_box_mm;if(!box)continue;
   assert(box[0][i]>=data.zoom_mm[0][i]-1e-6&&box[1][i]<=data.zoom_mm[1][i]+1e-6,'zoom region contains the '+arm+' main-track box on '+axis);}
 }
 el('zoom').onclick();assert.equal(el('zoom').textContent,'Zoom to track');
 assert.deepEqual(axisRanges(),full,'second press restores the full frame');
 zoomChecked++;
}
assert(zoomChecked>0,'no event carries a zoom region');
console.log(`Passed: ${ids.length} selectable events; three event changes; axes, heatmaps, overlays, 3D lines and tabs; zoom to track and back on ${zoomChecked} events.`);
