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
        self.display=dict(self.config['display'])
        self.display['nbins']=int(self.state.get('nbins',self.display['nbins']))
        self.display['voxel_opacity']=float(self.state.get('voxel_opacity',self.display['voxel_opacity']))
        self.worker=Worker(self.work,card)

    def catalogue(self,sample):
        if sample not in self.config['samples']:raise ValueError('Unknown sample')
        if sample not in self.ids:
            path=self.worker.request('open_sample',sample=sample)
            with h5py.File(path) as f:self.ids[sample]=[str(i) for i in f['event_id'][:]]
        return self.ids[sample]

    def save(self):
        write_toml(self.state_path,{**self.state,'view_json':json.dumps(self.view,separators=(',',':'))})

    def navigate(self, request):
        with self.lock:
            current=self.payload
            action=request.get('action','current')
            if action=='current' and current is not None:return current
            sample=request.get('sample',current['sample'] if current else self.state.get('sample',self.config['session']['initial_sample']))
            if sample not in self.config['samples']:sample=self.config['session']['initial_sample']
            ids=self.catalogue(sample)
            if not ids:raise ValueError('The selected sample is empty')
            event_id=current['event_id'] if current and current['sample']==sample else self.state.get('event_id',ids[0])
            if event_id not in ids:event_id=ids[0]
            if action=='sample':event_id=ids[0]
            elif action in ('next','previous'):
                index=ids.index(event_id)+(1 if action=='next' else -1)
                if not 0<=index<len(ids):raise ValueError('Already at the sample boundary')
                event_id=ids[index]
            elif action=='id':
                event_id=str(int(request['event_id']))
                if event_id not in ids:raise ValueError(f'Event {event_id} is absent from this sample')
            elif action not in ('current','bins'):raise ValueError('Unknown navigation action')
            nbins=int(request.get('nbins',self.display['nbins']))
            if not 1<=nbins<=500:raise ValueError('Display bin target must be between 1 and 500')
            display={**self.display,'nbins':nbins}
            start=time.perf_counter()
            path=self.worker.request('load_event',sample=sample,event_id=event_id,nbins=nbins)
            payload=snapshot_payload(path,display)
            same=current is not None and (sample,event_id)==(current['sample'],current['event_id'])
            restoring=current is None and (sample,event_id)==(self.state.get('sample'),self.state.get('event_id'))
            if not same and not restoring:self.view={}
            payload.update(samples=list(self.config['samples']),view=self.view,opacity=display['voxel_opacity'],
                           auto_rotate=bool(display.get('auto_rotate',False)),elapsed_seconds=time.perf_counter()-start)
            self.state.update(sample=sample,event_id=event_id,nbins=nbins,voxel_opacity=display['voxel_opacity'])
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
