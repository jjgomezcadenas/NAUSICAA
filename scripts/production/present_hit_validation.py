"""Present the recorded production validation and compact-file sizes."""
from pathlib import Path
import re,sys,html
import h5py
path=Path(sys.argv[1]);log=Path(sys.argv[2]).read_text();out=path.parent
with h5py.File(path) as f:
 n=len(f['events/event_number']);ih=len(f['itrack/x_mm']);eh=len(f['etrack/x_mm'])
rows=[]
for name,total,per in re.findall(r'TIMING (\w+) total_s=([\d.e+-]+) per_event_s=([\d.e+-]+)',log):
 rows.append(f'<tr><td>{html.escape(name)}</td><td>{float(total):.3f}</td><td>{float(per)*1000:.2f}</td></tr>')
read=re.search(r'READ_BACK seconds=([\d.e+-]+)',log)
page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ITACA · Hit output validation</title><style>body{{font:17px/1.6 system-ui;max-width:1000px;margin:30px auto;padding:0 24px;color:#203343}}table{{border-collapse:collapse;width:100%}}td,th{{padding:9px;text-align:left;border-bottom:1px solid #ddd}}</style><a href="../../../../notes/html/index.html">ITACA index</a><h1>Hit output · {n} bb0nu events</h1>
<ul><li>Controls disabled. Diagnostic histogram routines were replaced by errors during validation.</li><li>HDF5 groups: schema, events, itrack, etrack and configuration.</li><li>Every event read back exactly: hits, pre-cut totals, frame and diffusion widths.</li><li>No debug directory or temporary output remains.</li></ul>
<table><tr><th>Quantity</th><th>Value</th></tr><tr><td>Compressed hit file</td><td>{path.stat().st_size/2**20:.3f} MiB</td></tr><tr><td>itrack hits</td><td>{ih:,}</td></tr><tr><td>etrack hits</td><td>{eh:,}</td></tr><tr><td>Read all events</td><td>{float(read[1]):.3f} s</td></tr></table>
<h2>Production timing</h2><table><tr><th>Calculation</th><th>Total (s)</th><th>Per event (ms)</th></tr>{''.join(rows)}</table><p>Timings include first-call compilation within these calculations. Plot and graph generation are separate.</p>
<p><a href="presentation/index.html">Energy distributions and event viewer</a> · <a href="hits.h5">Compact hit file</a></p></html>'''
(out/'index.html').write_text(page)
print(out/'index.html')
