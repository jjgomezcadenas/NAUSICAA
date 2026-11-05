"""Verify real-history provenance, nonlinear agreement and the viewer HTTP interface."""
from pathlib import Path
import json,sys,urllib.request,urllib.error
import h5py,numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from pad_history_viewer.engine import Engine

def main():
 engine=Engine(work=ROOT/'notes/physics/work/pad_history_test')
 checks=0
 try:
  payload=engine.history('p4_bb0nu',3340,58,-55)
  assert len(payload['edge_time_s'])==len(payload['counts'])+1;checks+=1
  assert max(payload['counts'])==30633;checks+=1
  assert 0 in payload['counts'];checks+=1
  assert sum(payload['counts'])==sum(engine.history('p4_bb0nu',3340,58,-55)['counts']);checks+=1
  with h5py.File(ROOT/'results/pixel_charge_1mm_4mm/bb0nu.h5') as source:
   ids=source['events/event_id'][:];i=np.flatnonzero(ids==3340)[0];sizes=source['events/cells_with_charge'][:]
   start=int(sizes[:i].sum());stop=start+int(sizes[i]);x=source['cells/pixel_x'][start:stop];y=source['cells/pixel_y'][start:stop]
   expected=np.zeros(len(source['event_axes/3340/time_edges_ms'])-1,dtype=np.int64)
   sl=source['cells/time_slice'][start:stop];n=source['cells/n_ions'][start:stop]
   np.add.at(expected,sl[(x==58)&(y==-55)],n[(x==58)&(y==-55)])
   np.testing.assert_array_equal(expected,payload['counts'][:len(expected)]);checks+=1
  engine.simulate('p4_bb0nu',3340,58,-55)
  actual=engine.history('p4_bb0nu',3340,58,-55)
  assert actual['spice']['meta']['charge_error']<1e-7;checks+=1
  assert actual['spice']['meta']['baseline_validation_max_error_V']<1e-9;checks+=1
  with h5py.File(ROOT/'results/ltspice/pad_history_readout/spice_analysis.h5') as old:
   ref=old['p4_bb0nu_largest_bucket'];times=ref['edge_time_s'][:]
   ends=np.array(actual['edge_time_s'][1:]);use=(ends>=times[1])&(ends<=times[-1])
   expected=np.interp(ends[use],times[1:],ref['frame_signal_V'][:])
   new=np.array(actual['spice']['estimated_ions'])[use]*actual['meta']['gain_V_per_ion']
   np.testing.assert_allclose(new,expected,atol=2e-6,rtol=2e-4);checks+=1
  try:engine.history('p4_bb0nu',3340,999999,999999)
  except ValueError:checks+=1
  else:raise AssertionError('Absent pad accepted')
 finally:engine.close()
 if len(sys.argv)>1:
  base=sys.argv[1]
  def get(route):
   with urllib.request.urlopen(base+route) as r:return json.load(r)
  assert len(get('/api/config')['samples'])==4;checks+=1
  assert len(get('/api/events?sample=p4_bb0nu'))==8000;checks+=1
  for sample in ('p2_bb0nu','p2_electrons','p4_bb0nu','p4_electrons'):
   standard=get('/api/standard?sample='+sample)
   assert len(standard)==60 and len({v['event_id'] for v in standard})==60;checks+=1
   assert standard==get('/api/standard?sample='+sample);checks+=1
   with h5py.File(engine.directory/'statistics.h5') as stat:
    event_ids=stat[sample]['event_id'][:];flags=stat[sample]['criterion_2/affected_pads'][:]
    flag_lookup=dict(zip(event_ids,flags))
    assert all(flag_lookup[v['event_id']]==0 for v in standard);checks+=1

  live=get('/api/history?sample=p4_bb0nu&event_id=3340&pixel_x=58&pixel_y=-55')
  assert live['spice'] is not None;checks+=1
  try:get('/api/history?sample=wrong&event_id=1&pixel_x=0&pixel_y=0')
  except urllib.error.HTTPError as error:assert error.code==400;checks+=1
  else:raise AssertionError('Invalid sample accepted')
  request=urllib.request.Request(base+'/api/simulate',data=b'{}',headers={'Origin':'https://example.invalid','Content-Type':'application/json'},method='POST')
  try:urllib.request.urlopen(request)
  except urllib.error.HTTPError as error:assert error.code==403;checks+=1
  else:raise AssertionError('Cross-origin mutation accepted')
 print(f'{checks} real-history, SPICE and HTTP integration checks passed')
if __name__=='__main__':main()
