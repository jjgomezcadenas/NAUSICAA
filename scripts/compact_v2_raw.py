#!/usr/bin/env python3
"""Losslessly shuffle/compress v2 HDF5 acquisitions for version control."""
from pathlib import Path
import h5py,numpy as np
ROOT=Path(__file__).resolve().parents[1]
def compact(path):
 with h5py.File(path) as src:
  need=[]
  src.visititems(lambda n,o:need.append(n) if isinstance(o,h5py.Dataset) and o.dtype.kind in 'fci' and not o.shuffle else None)
  if not need:return
  tmp=path.with_suffix('.compact.h5')
  with h5py.File(tmp,'w') as dst:
   def copy(a,b):
    b.attrs.update(dict(a.attrs))
    for k,o in a.items():
     if isinstance(o,h5py.Group):copy(o,b.create_group(k))
     else:
      data=o[:];d=b.create_dataset(k,data=data,compression='gzip',compression_opts=9,shuffle=True) if o.dtype.kind in 'fci' else b.create_dataset(k,data=data)
      d.attrs.update(dict(o.attrs));assert np.array_equal(d[:],data)
   copy(src,dst)
 tmp.replace(path)
 print('Compacted',path.relative_to(ROOT),flush=True)
def main():
 for name in ('rogelio_v2','rogelio_v2_retuned'):
  compact(ROOT/'results/ltspice'/name/'raw.h5')
if __name__=='__main__':main()
