"""Local HTTP interface for cached pad histories and asynchronous SPICE jobs."""
from pathlib import Path
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from urllib.parse import urlparse,parse_qs
import argparse,concurrent.futures,json,threading,uuid,webbrowser
import h5py
import numpy as np
from plotly.offline import get_plotlyjs
from .engine import Engine,ROOT
WEB=Path(__file__).with_name('web')
LABELS={'p2_bb0nu':'2 mm · 20 ms · double beta','p2_electrons':'2 mm · 20 ms · electrons','p4_bb0nu':'4 mm · 40 ms · double beta','p4_electrons':'4 mm · 40 ms · electrons'}

class Application:
    """Own data access and independently pollable SPICE calculation jobs."""
    def __init__(self,engine):
        self.engine=engine;self.pool=concurrent.futures.ThreadPoolExecutor(max_workers=2)
        self.jobs={};self.pending={};self.lock=threading.Lock()
    def config(self):
        return {'samples':[dict(name=k,label=LABELS.get(k,k)) for k in self.engine.samples],
                'initial_sample':self.engine.config['initial_sample'],
                'simulator_available':self.engine.simulator is not None,
                'directory':str(self.engine.directory),'input_description':self.engine.config.get('input_description','Historical binned collection')}
    def ranked(self,sample):
        path=self.engine.plan()
        with h5py.File(path) as h:
            g=h[sample];order=np.argsort(g['minimum_output_V'][:])[:60]
            cols={k:g[k][:] for k in ('event_id','pixel_x','pixel_y','minimum_output_V','max_bucket_ions')}
            return [{k:float(v[i]) if k=='minimum_output_V' else int(v[i]) for k,v in cols.items()} for i in order]
    def standard(self,sample):
        self.engine.validate(sample)
        path=self.engine.worker.request('standard',sample=sample)
        with h5py.File(path) as h:
            columns={key:h[key][:] for key in ('event_id','pixel_x','pixel_y','minimum_output_V','max_bucket_ions')}
            return [{k:float(v[i]) if k=='minimum_output_V' else int(v[i]) for k,v in columns.items()} for i in range(len(columns['event_id']))]
    def start(self,fields):
        sample=fields['sample'];event=int(fields['event_id']);px=int(fields['pixel_x']);py=int(fields['pixel_y'])
        self.engine.validate(sample,event,px,py)
        key=(sample,event,px,py)
        with self.lock:
            if key in self.pending:return {'job_id':self.pending[key]}
            rid=uuid.uuid4().hex;self.jobs[rid]={'state':'running'};self.pending[key]=rid
        def task():
            try:
                self.engine.simulate(*key)
                result={'state':'complete'}
            except Exception as error:result={'state':'error','message':str(error)}
            with self.lock:self.jobs[rid]=result;self.pending.pop(key,None)
        self.pool.submit(task)
        return {'job_id':rid}
    def close(self):
        self.pool.shutdown(wait=True);self.engine.close()

def handler(app):
    plotly=get_plotlyjs().encode()
    class Handler(BaseHTTPRequestHandler):
        def log_message(self,*args):pass
        def respond(self,status,payload,mime='application/json'):
            body=json.dumps(payload,allow_nan=False).encode() if mime=='application/json' else payload
            self.send_response(status);self.send_header('Content-Type',mime);self.send_header('Content-Length',str(len(body)))
            self.send_header('Cache-Control','no-store');self.send_header('X-Content-Type-Options','nosniff');self.end_headers();self.wfile.write(body)
        def do_GET(self):
            try:
                url=urlparse(self.path);q={k:v[0] for k,v in parse_qs(url.query).items()}
                if url.path in ('/','/app.js','/style.css'):
                    file={'/':'index.html','/app.js':'app.js','/style.css':'style.css'}[url.path]
                    mime={'/':'text/html; charset=utf-8','/app.js':'application/javascript','/style.css':'text/css'}[url.path]
                    return self.respond(200,(WEB/file).read_bytes(),mime)
                if url.path=='/plotly.js':return self.respond(200,plotly,'application/javascript')
                if url.path=='/api/config':result=app.config()
                elif url.path=='/api/events':result=app.engine.events(q['sample'])
                elif url.path=='/api/pads':result=app.engine.pads(q['sample'],q['event_id'])
                elif url.path=='/api/ranked':result=app.ranked(q['sample'])
                elif url.path=='/api/standard':result=app.standard(q['sample'])
                elif url.path=='/api/history':result=app.engine.history(q['sample'],q['event_id'],q['pixel_x'],q['pixel_y'])
                elif url.path=='/api/statistics':result=app.engine.statistics()
                elif url.path=='/api/job':
                    with app.lock:result=app.jobs.get(q['id'],{'state':'error','message':'Unknown job'})
                else:return self.respond(404,{'error':'Unknown route'})
                return self.respond(200,result)
            except (ValueError,KeyError) as error:self.respond(400,{'error':str(error)})
            except Exception as error:self.respond(500,{'error':str(error)})
        def do_POST(self):
            try:
                if self.headers.get('Origin') not in (None,'http://'+self.headers.get('Host','')):
                    return self.respond(403,{'error':'Use this application from its local address'})
                size=int(self.headers.get('Content-Length','0'))
                if not 0<size<4096:return self.respond(400,{'error':'Invalid request size'})
                fields=json.loads(self.rfile.read(size))
                if urlparse(self.path).path!='/api/simulate':return self.respond(404,{'error':'Unknown route'})
                self.respond(202,app.start(fields))
            except (ValueError,KeyError) as error:self.respond(400,{'error':str(error)})
            except Exception as error:self.respond(500,{'error':str(error)})
    return Handler

def run(card=None,port=None,open_browser=True,work=None,simulator=None,julia=None):
    engine=Engine(card,work,simulator,julia);app=Application(engine)
    try:server=ThreadingHTTPServer(('127.0.0.1',engine.config['port'] if port is None else port),handler(app))
    except Exception:app.close();raise
    server.daemon_threads=True
    url=f'http://127.0.0.1:{server.server_port}'
    print('Pad history viewer:',url,flush=True)
    (engine.worker.work/'server.json').write_text(json.dumps({'url':url,'directory':str(engine.directory)})+'\n')
    if open_browser:webbrowser.open(url)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close();app.close()
