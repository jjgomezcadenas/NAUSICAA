"""Local MC event-viewer application and persistent Julia-worker client."""
from pathlib import Path
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse
import atexit
import fcntl
import json
import os
import subprocess
import threading
import time
import tomllib
import uuid
import webbrowser
import h5py
from plotly.offline import get_plotlyjs
from .presentation import snapshot_payload

ROOT=Path(__file__).resolve().parents[2]
WEB=Path(__file__).parent/'web'


def write_toml(path, values):
    """Atomically write flat application state; JSON quoting is TOML-compatible."""
    def scalar(v):
        if isinstance(v,bool):return 'true' if v else 'false'
        if isinstance(v,str):return json.dumps(v,ensure_ascii=False)
        return str(v)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(''.join(f'{k} = {scalar(v)}\n' for k,v in values.items()))
    temp.replace(path)


class Worker:
    """Serialize file-backed TOML requests through one resident Julia process."""
    def __init__(self, work, card):
        self.work=work
        self.log=(work/'worker.log').open('a')
        env=os.environ.copy()
        env.setdefault('JULIA_DEPOT_PATH',str(ROOT/'work/julia-depot')+':'+str(Path.home()/'.julia'))
        self.process=subprocess.Popen([env.get('JULIA','julia'),'--startup-file=no','--project='+str(ROOT),
            str(Path(__file__).with_name('worker.jl')),str(work),str(card)],cwd=ROOT,env=env,
            stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,text=True,bufsize=1)
        if self.process.stdout.readline().strip()!='READY':
            self.close();raise RuntimeError(f'Julia worker failed to start; see {work / "worker.log"}')

    def request(self, op, **fields):
        rid=uuid.uuid4().hex
        req=self.work/(rid+'.request.toml');res=self.work/(rid+'.response.toml')
        write_toml(req,dict(op=op,**fields))
        try:
            self.process.stdin.write(rid+'\n');self.process.stdin.flush()
            if self.process.stdout.readline().strip()!=rid:
                raise RuntimeError('Julia worker stopped; see worker.log')
            result=tomllib.loads(res.read_text())
            if result['status']!='ok':raise ValueError(result['message'])
            return result.get('path')
        finally:
            req.unlink(missing_ok=True);res.unlink(missing_ok=True)

    def close(self):
        if getattr(self,'process',None) and self.process.poll() is None:
            self.process.stdin.close()
            try:self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:self.process.terminate()
        if getattr(self,'log',None):self.log.close()


class Session:
    """Own successful navigation state, cached event IDs and the current figure payload."""
    def __init__(self, card, work_override=None):
        self.card=card;self.config=tomllib.loads(card.read_text())
        self.work=Path(work_override or ROOT/self.config['session']['work_dir']).resolve()
        self.work.mkdir(parents=True,exist_ok=True)
        self.lockfile=(self.work/'server.lock').open('w')
        try:fcntl.flock(self.lockfile,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('A viewer already uses this work directory; use --work-dir for another session')
        self.lock=threading.Lock();self.ids={};self.payload=None
        self.state_path=self.work/'session.toml';self.state={};self.view={}
        if self.config['session'].get('restore',True) and self.state_path.exists():
            self.state=tomllib.loads(self.state_path.read_text())
            self.view=json.loads(self.state.get('view_json','{}'))
        self.sources={key:str((ROOT/path).resolve()) for key,path in self.config['samples'].items()}
        self.files=json.loads(self.state.get('files_json','{}'))
        self.sources.update(self.files)
        self.default_directory=(ROOT/self.config['session'].get('file_directory','events')).resolve()
        tc=self.config.get('transport',{})
        self.transport_settings=json.loads(self.state.get('transport_settings_json',json.dumps(
            dict(mode=tc.get('mode','random'),fixed_cm=tc.get('fixed_cm',80.0)))))
        self.realisations=json.loads(self.state.get('realisations_json','{}'))
        self.display=dict(self.config['display'])
        self.display['nbins']=int(self.state.get('nbins',self.display['nbins']))
        self.display['voxel_opacity']=float(self.state.get('voxel_opacity',self.display['voxel_opacity']))
        self.worker=Worker(self.work,card)

    def catalogue(self,sample):
        if sample not in self.sources:raise ValueError('Unknown sample')
        if sample not in self.ids:
            path=self.worker.request('open_sample',sample=sample,source_path=self.sources[sample])
            with h5py.File(path) as f:self.ids[sample]=[str(i) for i in f['event_id'][:]]
        return self.ids[sample]

    def save(self):
        write_toml(self.state_path,{**self.state,'view_json':json.dumps(self.view,separators=(',',':')),
                                   'transport_settings_json':json.dumps(self.transport_settings),
                                   'realisations_json':json.dumps(self.realisations),
                                   'files_json':json.dumps(self.files,separators=(',',':'))})

    def resolve_file(self, path):
        """Resolve filenames relative to the configured event directory; accept absolute paths."""
        p=Path(path).expanduser()
        return (p if p.is_absolute() else self.default_directory/p).resolve()

    def list_files(self, directory=None):
        """List local subdirectories and HDF5 files for the browser's file picker."""
        path=self.resolve_file(directory) if directory else self.default_directory
        if not path.is_dir():raise ValueError(f'Directory does not exist: {path}')
        entries=[dict(name=p.name,path=str(p),directory=p.is_dir()) for p in path.iterdir()
                 if not p.name.startswith('.') and (p.is_dir() or p.suffix.lower()=='.h5')]
        entries.sort(key=lambda e:(not e['directory'],e['name'].lower()))
        return dict(directory=str(path),parent=str(path.parent),entries=entries)

    def navigate(self, request):
        with self.lock:
            current=self.payload
            action=request.get('action','current')
            if action=='current' and current is not None:return current
            tab=request.get('tab',self.state.get('active_tab','energy'))
            if tab not in ('energy','tracks','transport','pads'):raise ValueError('Unknown viewer tab')
            sample=request.get('sample',current['sample'] if current else self.state.get('sample',self.config['session']['initial_sample']))
            if action=='file':
                selected=self.resolve_file(request.get('path',''))
                if not selected.is_file():raise ValueError(f'File does not exist: {selected}')
                prepared=self.worker.request('prepare_file',path=str(selected))
                sample=next((key for key,path in self.sources.items() if path==prepared),'file:'+prepared)
                self.sources[sample]=prepared
            if sample not in self.sources:sample=self.config['session']['initial_sample']
            ids=self.catalogue(sample)
            if not ids:raise ValueError('The selected sample is empty')
            event_id=current['event_id'] if current and current['sample']==sample else self.state.get('event_id',ids[0])
            if event_id not in ids:event_id=ids[0]
            if action in ('sample','file'):event_id=ids[0]
            elif action in ('next','previous'):
                index=ids.index(event_id)+(1 if action=='next' else -1)
                if not 0<=index<len(ids):raise ValueError('Already at the sample boundary')
                event_id=ids[index]
            elif action=='id':
                event_id=str(int(request['event_id']))
                if event_id not in ids:raise ValueError(f'Event {event_id} is absent from this sample')
            elif action not in ('current','bins','tab','transport','realisation'):raise ValueError('Unknown navigation action')
            nbins=int(request.get('nbins',self.display['nbins']))
            if not 1<=nbins<=500:raise ValueError('Display bin target must be between 1 and 500')
            display={**self.display,'nbins':nbins}
            transport_settings=dict(self.transport_settings)
            if action=='transport':
                transport_settings=dict(mode=request['drift_mode'],fixed_cm=float(request['fixed_cm']))
            event_key=self.sources[sample]+'|'+event_id
            realisation=self.realisations.get(event_key,0)+(1 if action=='realisation' else 0)
            start=time.perf_counter()
            path=self.worker.request('pads' if tab=='pads' else 'transport' if tab=='transport' else 'reconstruct' if tab=='tracks' else 'load_event',sample=sample,
                                     source_path=self.sources[sample],event_id=event_id,nbins=nbins,
                                     drift_mode=transport_settings['mode'],fixed_cm=transport_settings['fixed_cm'],realisation=realisation)
            payload=snapshot_payload(path,display)
            same=current is not None and (sample,event_id)==(current['sample'],current['event_id'])
            restoring=current is None and (sample,event_id)==(self.state.get('sample'),self.state.get('event_id'))
            if not same and not restoring:
                self.view={'layers':self.view['layers']} if 'layers' in self.view else {}
            if sample.startswith('file:'):self.files[sample]=self.sources[sample]
            self.transport_settings=transport_settings;self.realisations[event_key]=realisation
            payload.update(samples=list(self.sources),sample_files=self.sources.copy(),file_path=self.sources[sample],
                           transport_controls={**self.config.get('transport',{}),**transport_settings,'realisation':realisation},
                           file_directory=str(self.default_directory),view=self.view,opacity=display['voxel_opacity'],active_tab=tab,
                           auto_rotate=bool(display.get('auto_rotate',False)),elapsed_seconds=time.perf_counter()-start)
            self.state.update(sample=sample,event_id=event_id,nbins=nbins,voxel_opacity=display['voxel_opacity'],active_tab=tab)
            self.display=display;self.payload=payload;self.save()
            return payload

    def save_view(self,data):
        with self.lock:
            if self.payload and (data.get('sample'),data.get('event_id'))==(self.payload['sample'],self.payload['event_id']):
                self.view=data.get('view',{})
                opacity=float(data.get('opacity',self.display['voxel_opacity']))
                if not 0<=opacity<=1:raise ValueError('Opacity must be between 0 and 1')
                self.display['voxel_opacity']=opacity;self.state['voxel_opacity']=opacity
                self.payload.update(view=self.view,opacity=opacity)
                self.save()

    def close(self):
        self.worker.close();self.lockfile.close()


def run(card, *, port=None, open_browser=True, work_dir=None):
    """Start the localhost viewer and keep the Julia reader alive until shutdown."""
    session=Session(card,work_dir);atexit.register(session.close)
    bundle=get_plotlyjs().encode()
    class Handler(BaseHTTPRequestHandler):
        def log_message(self,fmt,*args):pass
        def send_bytes(self,body,kind,status=200):
            self.send_response(status);self.send_header('Content-Type',kind)
            self.send_header('Content-Length',str(len(body)))
            self.send_header('Cache-Control','public, max-age=86400' if self.path=='/plotly.js' else 'no-store')
            self.end_headers()
            try:self.wfile.write(body)
            except (BrokenPipeError,ConnectionResetError):pass
        def json(self,data,status=200):
            self.send_bytes(json.dumps(data,allow_nan=False,separators=(',',':')).encode(),'application/json',status)
        def do_GET(self):
            path=urlparse(self.path).path
            try:
                if path=='/api/event':self.json(session.navigate({}));return
                if path=='/api/health':self.json(dict(status='ready',worker_pid=session.worker.process.pid));return
                if path=='/plotly.js':self.send_bytes(bundle,'application/javascript');return
                files={'/':('index.html','text/html; charset=utf-8'),'/viewer.js':('viewer.js','application/javascript'),'/viewer.css':('viewer.css','text/css')}
                if path not in files:self.send_error(404);return
                name,kind=files[path];self.send_bytes((WEB/name).read_bytes(),kind)
            except Exception as exc:self.json(dict(error=str(exc)),400)
        def do_POST(self):
            try:
                origin=self.headers.get('Origin')
                if origin and urlparse(origin).netloc!=self.headers.get('Host'):raise ValueError('Request origin differs from viewer')
                size=int(self.headers.get('Content-Length',0))
                if size>65536:raise ValueError('Request too large')
                data=json.loads(self.rfile.read(size))
                if self.path=='/api/navigate':self.json(session.navigate(data))
                elif self.path=='/api/files':self.json(session.list_files(data.get('directory')))
                elif self.path=='/api/view':session.save_view(data);self.json(dict(status='ok'))
                else:self.send_error(404)
            except Exception as exc:self.json(dict(error=str(exc)),400)
    host=session.config['server']['host']
    if host not in ('127.0.0.1','localhost'):raise ValueError('Viewer host must be localhost')
    server=ThreadingHTTPServer((host,session.config['server']['port'] if port is None else port),Handler)
    url=f'http://127.0.0.1:{server.server_port}/'
    (session.work/'url.txt').write_text(url+'\n')
    print(f'MC event viewer: {url}\nSession: {session.work}\nPress Ctrl-C to stop.',flush=True)
    if open_browser and session.config['server'].get('open_browser',True):webbrowser.open(url)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close();session.close();atexit.unregister(session.close)
