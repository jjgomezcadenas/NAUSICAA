"""Independently verify saved corrected counts/covariance against a banded solve."""
from pathlib import Path
import argparse,json
import h5py,numpy as np
from scipy.linalg import solve_banded

def text(x):return x.decode() if isinstance(x,bytes) else str(x)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);a=p.parse_args();checks=0;events=0;negative=0
 with h5py.File(a.results/'summary.h5') as summary:
  for name in summary['groups'][:]:
   path=Path(text(summary[text(name)].attrs['response_path']))
   with h5py.File(path) as source:parent_path=text(source['metadata/parent_response_path'][()])
   with h5py.File(path) as h,h5py.File(parent_path) as parent:
    assert text(h.attrs['schema'])=='nausicaa.ion_csa_response.pilot.v2';checks+=1
    assert bool(h['metadata/complete'][()]);checks+=1
    np.testing.assert_array_equal(h['metadata/selected_event_ids'][:],parent['metadata/selected_event_ids'][:]);checks+=1
    taps=h['calibration/precorrection_linear_reconstruction_kernel'][:]
    for eid in h['metadata/selected_event_ids'][:]:
     e=h[f'events/{eid}'];old=parent[f'events/{eid}'];counts=e['counts'][:].T;n,pads=counts.shape
     bands=np.zeros((3,n));bands[0,1:]=taps[0];bands[1,:]=taps[1];bands[2,:-1]=taps[2]
     sigma=e['noise_sigma_ions'][:]
     np.testing.assert_array_equal(counts,old['counts'][:].T);checks+=1
     for field in ('linear_estimated_ions','noiseless_estimated_ions','noise_ions','recovered_ions'):
      expected=solve_banded((1,1),bands,old[field][:].T.astype(float));actual=e[field][:].T.astype(float)
      assert np.max(np.abs(actual-expected)/sigma[:,None])<.001;checks+=1
     for field in ('noise_seed','output_at_edges_V','average_at_edges_V','pixel_ix','pixel_iy','edge_time_s'):
      np.testing.assert_array_equal(e[field][:],old[field][:]);checks+=1
     matrix=np.diag(np.full(n,taps[1]))+np.diag(np.full(n-1,taps[0]),1)+np.diag(np.full(n-1,taps[2]),-1)
     c=h[f'covariance_by_frames/{n}/covariance_e2'][:];cv=h['calibration/precorrection_covariance_e2'][:]
     original=cv[np.abs(np.arange(n)[:,None]-np.arange(n)[None,:])]
     np.testing.assert_allclose(matrix@c@matrix.T,original,rtol=1e-10,atol=1e-8);checks+=1
     np.testing.assert_allclose(sigma**2,np.diag(c),rtol=1e-12);checks+=1
     assert abs(e['linear_estimated_ions'][:].sum()/counts.sum()-1)<1e-6;checks+=1
     negative+=int((e['recovered_ions'][:]<0).sum());events+=1
 assert events==200 and negative>0;checks+=1
 status={'passed':True,'events':events,'checks':checks,'negative_estimates_retained':negative}
 (a.results/'independent_validation.json').write_text(json.dumps(status,indent=2)+'\n');print(status)
if __name__=='__main__':main()
