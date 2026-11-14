// Exercise event switching and Plotly inputs without a browser rendering engine.
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const root=process.argv[2],html=fs.readFileSync(root+'/index.html','utf8');
const ids=[...html.split('<nav>')[1].split('</nav>')[0].matchAll(/<option value="(\d+)">/g)].map(m=>m[1]);assert(ids.length>0);
const elements={};function el(id){return elements[id]??={value:ids[0],selectedIndex:0,options:ids,classList:{toggle(){}},setAttribute(){},addEventListener(){},on(){}};}
const plots={};const ctx={requestAnimationFrame:f=>f(),document:{getElementById:el,createElement:()=>({remove(){}}),body:{appendChild(){}}},Plotly:{Plots:{resize(){}},react:(id,data,layout,config)=>{plots[id]={data,layout,config};el(id).layout=layout;el(id).data=data;return Promise.resolve();},relayout:()=>Promise.resolve()}};vm.createContext(ctx);
el('blob-radius').value='3';el('show-path').checked=true;
vm.runInContext(html.match(/<script>([\s\S]*)<\/script>/)[1],ctx);
for(const id of [ids[0],ids[Math.floor(ids.length/2)],ids.at(-1)]){
 el('event').value=id;vm.runInContext(fs.readFileSync(root+'/viewer/event_'+id+'.js','utf8'),ctx);
 assert.equal(Object.keys(plots).length,10);
 for(const arm of ['itrack','etrack']){
  for(const view of ['xy','xz','yz']){
   const p=plots[arm+'-'+view],h=p.data[0];assert(p.layout.title.includes(id));assert.equal(p.config.scrollZoom,false);
   for(const axis of ['xaxis','yaxis']){assert.equal(p.layout[axis].fixedrange,true);assert.equal(p.layout[axis].constrain,'domain');assert(p.layout[axis].title.text.includes('(mm)'));}
   assert.equal(h.z.length,h.y.length);assert(h.z.every(row=>row.length===h.x.length));
   assert(p.data.slice(1).every(t=>['lines','lines+markers'].includes(t.mode)));assert(h.z.every(row=>row.every(Number.isFinite)));
  }
  const p=plots[arm+'-3d'];assert.equal(p.data.filter(t=>t.type==='surface').length,2);for(const sphere of p.data.filter(t=>t.type==='surface'))for(const a of ['x','y','z'])assert(sphere[a].flat().every(Number.isFinite));assert.equal(p.config.scrollZoom,false);
  for(const axis of ['xaxis','yaxis','zaxis']){assert.equal(p.layout.scene[axis].showbackground,false);assert.equal(p.layout.scene[axis].showgrid,false);}
  assert(plots[arm+'-hist'].layout.title.includes(id));
 }
}
el('blob-radius').value='2';el('blob-radius').onchange();el('show-path').checked=false;el('show-path').onchange();assert(!plots['itrack-3d'].data.some(t=>t.mode==='lines'));
el('tab-etrack').onclick();el('tab-itrack').onclick();
console.log(`Passed: ${ids.length} selectable events; three event changes; axes, heatmaps, overlays, 3D lines and tabs.`);
