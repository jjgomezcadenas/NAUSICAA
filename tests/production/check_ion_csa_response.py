"""Independently check pilot HDF5 keys, ion ledgers and sampled-response accounting."""
from pathlib import Path
import argparse,json
import h5py,numpy as np

def text(x):return x.decode() if isinstance(x,bytes) else str(x)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);a=p.parse_args();checks=0;negatives=0;events=0;seeds=set()
 with h5py.File(a.results/'summary.h5') as summary:
  for raw in summary['groups'][:]:
   group=summary[text(raw)];path=Path(text(group.attrs['response_path']))
   with h5py.File(path) as h:
    assert text(h.attrs['schema'])=='nausicaa.ion_csa_response.pilot.v1';checks+=1
    assert bool(h['metadata/complete'][()]);checks+=1
    ids=h['metadata/selected_event_ids'][:];assert len(ids)==100 and len(set(ids))==100;checks+=1
    with h5py.File(text(h['metadata/source_path'][()])) as source:
     total=0
     for eid in ids:
      e=h[f'events/{eid}'];s=source[f'events/{eid}'];n=e['counts'][:].T;mu=e['noiseless_estimated_ions'][:].T;noise=e['noise_ions'][:].T;y=e['recovered_ions'][:].T
      enc=float(h['calibration'].attrs['ENC_e_rms']);frame=float(h['metadata/frame_s'][()]);edges=e['edge_time_s'][:];first=int(e['first_bucket'][()]);px=e['pixel_ix'][:];py=e['pixel_iy'][:]
      assert np.all(np.diff(edges)>0) and np.allclose(np.diff(edges),frame,rtol=0,atol=1e-12);checks+=1
      assert abs(edges[0]-first*frame)<1e-12;checks+=1
      assert n.shape==mu.shape==noise.shape==y.shape and len(edges)==n.shape[0]+1;checks+=1
      assert int(n.sum())==int(s['initial_pairs'][()])==int(e['initial_pairs'][()]);checks+=1;total+=int(n.sum())
      pads={tuple(v):i for i,v in enumerate(zip(px,py))};expected=np.zeros_like(n)
      for x,z,bucket,count in zip(s['pixel_ix'][:],s['pixel_iy'][:],s['bucket'][:],s['count'][:]):expected[int(bucket)-first,pads[(x,z)]]=count
      np.testing.assert_array_equal(n,expected);checks+=1
      np.testing.assert_allclose(y.astype(float),mu.astype(float)+noise.astype(float),rtol=1e-6,atol=.001*enc);checks+=1
      assert np.all(np.isfinite(y));checks+=1
      negative=int((y<0).sum());negatives+=negative
      assert text(e['pair_checksum'][()])==text(s['pair_checksum'][()]);checks+=1
      for seed in e['noise_seed'][:]:assert int(seed) not in seeds;seeds.add(int(seed))
      events+=1
     assert total==int(h['metadata/input_ion_total'][()]);checks+=1
 assert negatives>0;checks+=1
 result={'passed':True,'checks':checks,'events':events,'negative_estimates_retained':negatives,'distinct_pad_seeds':len(seeds)}
 (a.results/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
