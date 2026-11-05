"""Serialize Julia calculations and acquire cached LTspice waveforms."""
from pathlib import Path
import concurrent.futures
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import threading
import tomllib
import uuid
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[2]

def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def scalar(value):
    if isinstance(value,bytes):return value.decode()
    if isinstance(value,np.generic):return value.item()
    return value

def write_toml(path,values):
    path.write_text(''.join(f'{key} = '+('true' if val is True else 'false' if val is False else json.dumps(val))+'\n' for key,val in values.items()))

class Worker:
    """One resident Julia calculator; requests are serialized by a thread lock."""
    def __init__(self,card,work,julia=None):
        self.work=Path(work);self.work.mkdir(parents=True,exist_ok=True);self.lock=threading.Lock()
        self.lockfile=(self.work/'worker.lock').open('w')
        try:fcntl.flock(self.lockfile,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('A Julia worker already uses this work directory')
        self.log=(self.work/'worker.log').open('a')
        env=os.environ.copy()
        env.setdefault('JULIA_DEPOT_PATH',str(ROOT/'notes/physics/work/julia-depot')+':'+str(Path.home()/'.julia'))
        executable=julia or env.get('JULIA') or shutil.which('julia') or str(Path.home()/'.juliaup/bin/julia')
        worker=ROOT/tomllib.loads(Path(card).read_text()).get('worker_script','src/pad_history_viewer/worker.jl')
        self.process=subprocess.Popen([executable,'--startup-file=no','--project='+str(ROOT),str(worker),str(self.work),str(card)],cwd=ROOT,env=env,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,text=True,bufsize=1)
        ready=self.process.stdout.readline().strip().split('\t')
        if len(ready)!=2 or ready[0]!='READY':
            self.close();raise RuntimeError(f'Julia worker failed; see {self.work}/worker.log')
        self.directory=Path(ready[1])
    def request(self,op,**fields):
        with self.lock:
            rid=uuid.uuid4().hex;req=self.work/(rid+'.request.toml');res=self.work/(rid+'.response.toml')
            write_toml(req,dict(op=op,**fields))
            try:
                self.process.stdin.write(rid+'\n');self.process.stdin.flush()
                if self.process.stdout.readline().strip()!=rid:raise RuntimeError('Julia worker stopped; see worker.log')
                answer=tomllib.loads(res.read_text())
                if answer['status']!='ok':raise ValueError(answer['message'])
                return Path(answer['path'])
            finally:
                req.unlink(missing_ok=True);res.unlink(missing_ok=True)
    def close(self):
        if getattr(self,'process',None) and self.process.poll() is None:
            self.process.stdin.close()
            try:self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:self.process.terminate()
        if getattr(self,'log',None):self.log.close()
        if getattr(self,'lockfile',None):self.lockfile.close()

class Engine:
    """Own the calculation worker, immutable case files and simulator selection."""
    def __init__(self,card=None,work=None,simulator=None,julia=None):
        self.card=Path(card or ROOT/'config/pad_history_viewer_direct.toml').resolve()
        self.config=tomllib.loads(self.card.read_text())
        self.readout=tomllib.loads((ROOT/self.config['readout_card']).read_text())
        self.samples={s['name']:s for s in self.readout['samples']}
        self.worker=Worker(self.card,work or ROOT/'work/pad_history_viewer',julia)
        self.directory=self.worker.directory
        source=ROOT/self.readout['source_decks']
        self.model=source/'corrected.inc';self.source=source/'40ms_single_N400_r2.cir'
        candidates=[ROOT/'notes/physics/work/ltspice-review/legacy/LTspice.pkg/Payload/Applications/LTspice.app/Contents/MacOS/LTspice',Path('/Applications/LTspice.app/Contents/MacOS/LTspice')]
        self.simulator=Path(simulator).resolve() if simulator else next((p for p in candidates if p.exists()),None)
        self.events_cache={};self.closed=False
    def validate(self,sample,event=None,px=None,py=None):
        if sample not in self.samples:raise ValueError('Unknown sample')
        if event is not None:
            event=int(event)
            if event not in self.events(sample):raise ValueError(f'Event {event} is absent')
        return {'sample':sample,**({} if event is None else {'event_id':event}),**({} if px is None else {'pixel_x':int(px),'pixel_y':int(py)})}
    def events(self,sample):
        if sample not in self.samples:raise ValueError('Unknown sample')
        if sample not in self.events_cache:
            path=self.worker.request('events',sample=sample)
            with h5py.File(path) as h:self.events_cache[sample]=h['event_id'][:].tolist()
            path.unlink()
        return self.events_cache[sample]
    def pads(self,sample,event):
        fields=self.validate(sample,event);path=self.worker.request('pads',**fields)
        with h5py.File(path) as h:
            columns={key:h[key][:].tolist() for key in h}
        path.unlink()
        return [dict(zip(columns,values)) for values in zip(*columns.values())]
    def linear(self,sample,event,px,py):
        return self.worker.request('history',**self.validate(sample,event,px,py))
    def history(self,sample,event,px,py):
        path=self.linear(sample,event,px,py);spice=path.with_name('spice.h5')
        with h5py.File(path) as h:
            result={'meta':{k:scalar(v) for k,v in h.attrs.items()},**{key:h[key][:].tolist() for key in h}}
        result['spice']=None
        if spice.exists():
            with h5py.File(spice) as h:
                result['spice']={'meta':{k:scalar(v) for k,v in h.attrs.items()},**{key:h[key][:].tolist() for key in ('time_s','output_V','output_at_edges_V','estimated_ions','relative_loss','deficit_e')}}
        return result
    def simulate(self,sample,event,px,py,refine=False):
        if self.simulator is None or not self.simulator.exists():raise ValueError('Set --simulator to the LTspice executable')
        fields=self.validate(sample,event,px,py)
        path=self.linear(sample,event,px,py);folder=path.parent;destination=folder/('refined.h5' if refine else 'spice.h5')
        with (folder/'simulation.lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            if destination.exists():return destination
            self.worker.request('prepare',**fields,refine=refine)
            suffix='refined' if refine else 'history'
            model=folder/'corrected.inc'
            if not model.exists():shutil.copyfile(self.model,model)
            if sha(model)!=sha(self.model):raise ValueError('Cached model differs from the selected circuit')
            body='\n'.join(x for x in self.source.read_text().splitlines() if not x.lower().startswith(('is1 ','.tran','.save','.end')))+'\n'
            inp=tomllib.loads((folder/(suffix+'_input.toml')).read_text())
            step=self.config['solver_step_s']/(2 if refine else 1)
            deck=folder/(suffix+'.cir')
            deck.write_text(body+f'Is1 0 in PWL FILE="{suffix}.pwl"\n.tran 0 {inp["stoptime_s"]:.17e} 0 {step:.17e}\n.save V(in) V(out) I(Is1)\n.end\n')
            from ion_signal.spice_acquisition import run
            raw=folder/(suffix+'_raw.h5')
            for attempt in range(3):
                if attempt:
                    step=min(step/4,self.config['current_step_s']/2**attempt)
                    deck.write_text(body+f'Is1 0 in PWL FILE="{suffix}.pwl"\n.tran 0 {inp["stoptime_s"]:.17e} 0 {step:.17e}\n.save V(in) V(out) I(Is1)\n.end\n')
                wave,_=run(self.simulator,deck)
                with h5py.File(raw,'w') as h:
                    for key,val in wave.items():h.create_dataset(key,data=np.asarray(val).real,compression='gzip',shuffle=True)
                    h.attrs.update(deck_sha256=sha(deck),model_sha256=sha(model),simulator_sha256=sha(self.simulator),acquisition_sha256=sha(Path(__file__)),solver_step_s=step,charge_resolution_retries=attempt)
                try:
                    answer=self.worker.request('reduce',**fields,refine=refine)
                    break
                except ValueError as error:
                    if 'charge_error' not in str(error) or attempt==2:raise
            # Reduced HDF5 retains the exact charge/extrema and all sampled values.
            for ext in ('.raw','.log','.op.raw'):
                deck.with_suffix(ext).unlink(missing_ok=True)
            raw.unlink()
            return answer
    def plan(self):return self.worker.request('plan')
    def aggregate(self):return self.worker.request('aggregate')
    def statistics(self):
        path=self.directory/'statistics.h5'
        if not path.exists():return None
        with h5py.File(path) as h:
            out={'meta':{k:scalar(v) for k,v in h.attrs.items()},'criteria':[dict(loss=float(r),noise_multiple=float(s)) for r,s in zip(h['criterion_loss_fraction'][:],h['criterion_noise_multiple'][:])],'samples':{}}
            for name in self.samples:
                if name not in h:continue
                g=h[name]
                out['samples'][name]={'meta':{k:scalar(v) for k,v in g.attrs.items()},'criteria':[{k:scalar(v) for k,v in g[f'criterion_{i+1}'].attrs.items()} for i in range(len(out['criteria']))],'control_pads_with_flag':g['control_pads_with_flag'][:].tolist()}
            return out
    def close(self):
        if not self.closed:self.worker.close();self.closed=True
