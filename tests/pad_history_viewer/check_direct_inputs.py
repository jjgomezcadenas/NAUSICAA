"""Validate direct source histories, CSA refinement and the updated viewer."""
from pathlib import Path
import argparse,json,sys,urllib.request,urllib.error
import h5py,numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from pad_history_viewer.engine import Engine

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--url');args=parser.parse_args()
 engine=Engine(card=ROOT/'config/pad_history_viewer_direct.toml',work=ROOT/'notes/physics/work/direct_history_test')
 checks=0
 try:
  for sample,event in (('p2_electrons',3790),('p4_bb0nu',3340)):
   pads=engine.pads(sample,event);pad=min(pads,key=lambda p:p['minimum_output_V']);px,py=pad['pixel_x'],pad['pixel_y']
   payload=engine.history(sample,event,px,py);n=payload['meta']['original_buckets']
   assert len(payload['edge_time_s'])==len(payload['counts'])+1;checks+=1
   with h5py.File(ROOT/engine.samples[sample]['path']) as source:
    assert source.attrs['schema'].decode()=='nausicaa.pixel_charge.direct.v2';checks+=1
    ids=source['events/event_id'][:];i=np.flatnonzero(ids==event)[0];sizes=source['events/cells_with_charge'][:];start=int(sizes[:i].sum());stop=start+int(sizes[i])
    x=source['cells/pixel_x'][start:stop];y=source['cells/pixel_y'][start:stop];sl=source['cells/time_slice'][start:stop];ions=source['cells/n_ions'][start:stop]
    expected=np.zeros(n,dtype=np.int64);use=(x==px)&(y==py);np.add.at(expected,sl[use],ions[use])
    np.testing.assert_array_equal(expected,payload['counts'][:n]);checks+=1
    np.testing.assert_allclose(source[f'event_axes/{event}/time_edges_ms'][:]*.001,payload['edge_time_s'][:n+1],atol=1e-6);checks+=1
   coarse=engine.simulate(sample,event,px,py);fine=engine.simulate(sample,event,px,py,refine=True)
   with h5py.File(coarse) as a,h5py.File(fine) as b:
    np.testing.assert_allclose(a['frame_signal_V'][:],b['frame_signal_V'][:],rtol=1e-4,atol=1e-7);checks+=1
    assert a.attrs['charge_error']<1e-7;checks+=1
    assert a.attrs['baseline_validation_max_error_V']<1e-9;checks+=1
  if args.url:
   def get(route):
    with urllib.request.urlopen(args.url+route) as r:return json.load(r)
   assert get('/api/config')['input_description']=='Direct carrier collection';checks+=1
   for sample in engine.samples:
    entries=get('/api/standard?sample='+sample)
    assert len(entries)==60 and len({v['event_id'] for v in entries})==60;checks+=1
    assert entries==get('/api/standard?sample='+sample);checks+=1
    with h5py.File(engine.directory/'statistics.h5') as stats:
     ids=stats[sample]['event_id'][:];flags=stats[sample]['criterion_2/affected_pads'][:]
     flagged=dict(zip(ids,flags));assert all(flagged[e['event_id']]==0 for e in entries);checks+=1

    e=entries[0]
    payload=get(f'/api/history?sample={sample}&event_id={e["event_id"]}&pixel_x={e["pixel_x"]}&pixel_y={e["pixel_y"]}')
    assert payload['spice'] is not None;checks+=1
   try:get('/api/history?sample=wrong&event_id=1&pixel_x=0&pixel_y=0')
   except urllib.error.HTTPError as error:assert error.code==400;checks+=1
   else:raise AssertionError('Invalid sample accepted')
   request=urllib.request.Request(args.url+'/api/simulate',data=b'{}',headers={'Origin':'https://example.invalid','Content-Type':'application/json'},method='POST')
   try:urllib.request.urlopen(request)
   except urllib.error.HTTPError as error:assert error.code==403;checks+=1
   else:raise AssertionError('Cross-origin mutation accepted')

 finally:engine.close()
 print(f'{checks} direct-source, refinement and viewer checks passed')
if __name__=='__main__':main()
