"""Present the saved kNN connectivity scan."""
from pathlib import Path
import sys,tomllib
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(sys.argv[1]);out=p.parent
with h5py.File(p) as f:
    tc=tomllib.loads(f['card'][()].decode())['track_analysis']
    ion_scale=tc.get('ion_edge_sigma_scale',tc['edge_sigma_scale']);electron_scale=tc['edge_sigma_scale']
    distance_label=f'itrack {ion_scale:g}σₜ; etrack {electron_scale:g}σₜ'
    n=int(f['events'][()]);ks=sorted(int(k[2:]) for k in f['ions']);rows=[];noise_rows=[];fixed_rows=[]
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    for arm,label in [('ions','itrack'),('pe','etrack')]:
        fractions=[100*f[f'{arm}/k_{k}/reference_broken_fraction'][()] for k in ks]
        ax.plot(ks,fractions,'o-',label=label)
        for k,v in zip(ks,fractions):ax.annotate(f'{v:.2f}%',(k,v),xytext=(0,8 if arm=='ions' else -17),textcoords='offset points',ha='center')
    ax.set(xlabel='k: nearest neighbours considered',ylabel='Previously connected signal tracks broken (%)',xticks=ks,ylim=(0,None));ax.legend();ax.grid(alpha=.2)
    ax.set_title(f'{n} events · 1 keV hit cut · {distance_label}')
    fig.savefig(out/'knn_scan.png',dpi=150);plt.close(fig)
    for k in ks:
        cells=[];noise=[];fixed=[]
        for arm in ['ions','pe']:
            g=f[f'{arm}/k_{k}'];cells.append(f'<td>{g["broken"][()]} / {g["baseline_connected"][()]} ({100*g["broken_fraction"][()]:.2f}%)</td>')
            fixed.append(f'<td>{g["reference_broken"][()]} / {g["reference_connected"][()]} ({100*g["reference_broken_fraction"][()]:.2f}%)</td>')
            noise.append(f'<td>{g["all_disconnected_after"][()]} / {n} ({100*g["all_disconnected_after"][()]/n:.1f}%)</td>')
        fixed_rows.append(f'<tr><td>{k}</td>'+''.join(fixed)+'</tr>')
        rows.append(f'<tr><td>{k}</td>'+''.join(cells)+'</tr>');noise_rows.append(f'<tr><td>{k}</td>'+''.join(noise)+'</tr>')
    comparison_text=''
    if 'comparison_6sigma' in f:
        c=f['comparison_6sigma'];den=int(c['reference_connected'][()]);cr=[]
        for k in ks:
            a=int(c[f'broken_6sigma_k_{k}'][()]);b=int(c[f'broken_7sigma_k_{k}'][()])
            cr.append(f'<tr><td>{k}</td><td>{a}/{den} ({100*a/den:.2f}%)</td><td>{b}/{den} ({100*b/den:.2f}%)</td></tr>')
        comparison_text=f'<h2>Ions: 6σₜ versus 7σₜ</h2><p>The same {den} signal tracks, connected before the cut at 6σₜ and k = 10.</p><table><tr><th>k</th><th>6σₜ</th><th>7σₜ</th></tr>{"".join(cr)}</table>'
    gap_text=''
    if (out/'ion_gap_check.h5').exists():
        with h5py.File(out/'ion_gap_check.h5') as gaps:
            e=gaps['events'];gaprows=''.join(f'<tr><td>{int(id)}</td><td>{limit:.3f}</td><td>{gap:.3f}</td><td>{int(count)}</td></tr>' for id,limit,gap,count in zip(e['event_number'][:],e['max_distance_mm'][:],e['minimum_component_gap_mm'][:],e['radius_components'][:]))
            gap_text=f'<h2>Why the 12 itracks remain broken</h2><p>A graph containing every possible connection within 6σₜ gives the same components as k = 30 in all 12 cases. The closest gap between components exceeds the distance limit.</p><table><tr><th>Event</th><th>6σₜ limit (mm)</th><th>Smallest component gap (mm)</th><th>Components</th></tr>{gaprows}</table>'
    page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>ITACA · kNN scan</title><style>body{{font:17px/1.6 system-ui;max-width:1050px;margin:35px auto;padding:0 20px;color:#203343}}table{{border-collapse:collapse;width:100%}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #ddd}}img{{max-width:100%}}</style><a href="../../../../notes/html/index.html">ITACA index</a><h1>kNN scan · {n} bb0nu events</h1><ul><li>Hit cut: 1 keV for itrack and etrack. Maximum connection distance: {distance_label}, using each event's ion or electron diffusion width.</li><li>Rebuild the existing kNN graph at every k, before and after the hit cut.</li><li>Broken fraction: signal-only events connected before the cut that have more than one component or no surviving component afterwards, divided by events connected before the cut at that k.</li><li>Signal-only hits have positive input signal; these labels are used only for this diagnostic.</li></ul>{comparison_text}<h2>Same initially connected tracks at every k</h2><p>Reference sample: signal-only tracks connected before the cut at k = 10. The denominator stays fixed.</p><table><tr><th>k</th><th>itrack</th><th>etrack</th></tr>{''.join(fixed_rows)}</table><h2>Using each k’s own initially connected sample</h2><table><tr><th>k</th><th>itrack: broken / initially connected</th><th>etrack: broken / initially connected</th></tr>{''.join(rows)}</table><img src="knn_scan.png"><h2>All retained hits, including noise</h2><p>Events with other than one connected component after the cut, divided by all {n} events.</p><table><tr><th>k</th><th>itrack</th><th>etrack</th></tr>{''.join(noise_rows)}</table>{gap_text}<p><a href="{p.name}">Numerical scan</a> · <a href="../tracks_500_cut1keV/index.html">1 keV track views</a></p></html>'''
    (out/'index.html').write_text(page)
print(out/'index.html')
