#!/usr/bin/env python3
"""Cache transistor responses for the initial pads in the standard-event menus."""
from pathlib import Path
import concurrent.futures,sys
import h5py
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from pad_history_viewer.engine import Engine

def main():
 engine=Engine(work=ROOT/'notes/physics/work/standard_pad_cache')
 try:
  tasks=[]
  for name in engine.samples:
   path=engine.worker.request('standard',sample=name)
   with h5py.File(path) as h:
    tasks.extend((name,int(event),int(px),int(py)) for event,px,py in zip(h['event_id'][:],h['pixel_x'][:],h['pixel_y'][:]))
  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
   futures={pool.submit(engine.simulate,*task):task for task in tasks}
   for i,future in enumerate(concurrent.futures.as_completed(futures),1):
    future.result()
    if i%40==0 or i==len(tasks):print(f'{i}/{len(tasks)} standard histories cached',flush=True)
 finally:engine.close()
if __name__=='__main__':main()
