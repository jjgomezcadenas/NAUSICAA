"""Make the figures and the HTML page from results/fast_trigger_gamma_outside_frame.h5; run after fast_trigger_gamma_outside_frame.jl."""
from pathlib import Path
import html
import time
START=time.perf_counter()
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path(__file__).resolve().parent
OUT=BASE/'results'
PLOTS=OUT/'fast_trigger_gamma_outside_frame';PLOTS.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'figure.dpi':120})
sections=[]
def annotation(ax,s,unit):
    ax.text(.97,.97,f"N={int(s['n'][()])}\nmean={s['mean'][()]:.3g} {unit}\nrms={s['rms'][()]:.3g} {unit}",transform=ax.transAxes,ha='right',va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
def save(fig,name):
    fig.tight_layout();fig.savefig(PLOTS/name);plt.close(fig)
    return f'<img src="fast_trigger_gamma_outside_frame/{name}" loading="lazy">'
with h5py.File(OUT/'fast_trigger_gamma_outside_frame.h5') as f:
 for sample,label in [('bb0nu','bb0nu (two electrons)'),('single_electrons','Single electrons')]:
    g=f[sample];ev=g['events'];ga=g['gammas'];es=g['event_statistics'];gs=g['gamma_statistics']
    block=[f'<h2>{label}</h2>']
    sel=g['selection'];input_count=int(sel['input'][()])
    block.append(f'<ul><li><strong>Input:</strong> {input_count:,} events.</li><li><strong>Inside detector:</strong> {int(sel["contained"][()]):,} events.</li><li><strong>Vertex within 50–1450 mm:</strong> {len(ev["event_number"]):,} events, containing {len(ga["event_number"]):,} gamma deposits.</li></ul>')
    block.append(f'<p>Normalization: only the {sel["contained"][()]} events with all deposits inside the detector box enter efficiency numerators and denominators. The {sel["rejected_containment"][()]} events outside it are excluded from the sample.</p>')
    block.append('<h3>Selection</h3><table><tr><th>Step</th><th>Accepted</th><th>Rejected at this step</th><th>Efficiency</th></tr>')
    block.append(f'<tr><td>All deposits inside Z = 0–1500 mm</td><td>{sel["contained"][()]}</td><td>{sel["rejected_containment"][()]}</td><td>100% (normalization sample)</td></tr>')
    block.append(f'<tr><td>Primary vertex in Z = 50–1450 mm</td><td>{sel["accepted"][()]}</td><td>{sel["rejected_vertex"][()]}</td><td>{sel["fiducial_percent_of_contained"][()]:.2f}% of contained</td></tr>')
    block.append(f'<tr><td>Reject gamma veto (gamma deposits outside the frame above {sel["gamma_energy_cut_keV"][()]:g} keV)</td><td>{sel["accepted_gamma_cut"][()]}</td><td>{sel["rejected_gamma_cut"][()]}</td><td>{sel["final_percent"][()]:.2f}% of contained</td></tr></table>')
    inside=g['all_gammas_inside_frame'];stats=g['all_gammas_inside_frame_statistics']
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ax,key,symbol in zip(axes,['nearest_pair_r_mm','nearest_pair_z_mm'],['R','Z']):
        ax.scatter(inside[key][:],inside['energy_keV'][:],s=12,alpha=.55)
        ax.set_xscale('symlog',linthresh=1);ax.set_yscale('log')
        distance_label = {'D': 'Distance to primary track', 'Z': 'Z separation to primary track', 'R': 'R separation to primary track'}[symbol]
        ax.set_title(f'All gamma deposits inside frame\nEnergy vs {symbol} separation of closest hit pair');ax.set_xlabel(f'{distance_label} (mm)');ax.set_ylabel('Gamma deposited energy (keV)');ax.grid(alpha=.2)
        ax.text(.97,.97,f"mean E = {stats['energy_keV/mean'][()]:.3g} keV\nRMS = {stats['energy_keV/rms'][()]:.3g} keV\nmean {symbol} = {stats[key+'/mean'][()]:.3g} mm\nRMS = {stats[key+'/rms'][()]:.3g} mm",transform=ax.transAxes,ha='right',va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
    block.append(save(fig,sample+'_gammas.png'))
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ax,axis in zip(axes,'rz'):
        ax.scatter(ga['outside_frame_'+axis+'_mm'][:],ga['energy_keV'][:],s=12,alpha=.55);ax.set_xscale('symlog',linthresh=1);ax.set_yscale('log');ax.set_title(f'Before gamma veto: all fiducial events\nEnergy vs maximum distance outside frame, {axis.upper()}');ax.set_xlabel(f'Distance outside frame, {axis.upper()} (mm)');ax.set_ylabel('E (keV)');ax.text(.03,.97,f'N={len(ga["energy_keV"])}',transform=ax.transAxes,va='top');ax.grid(alpha=.2)
    block.append(save(fig,sample+'_scatter.png'))
    fig,axes=plt.subplots(1,3,figsize=(15,4))
    for ax,axis in zip(axes,'xyz'):
        for prefix,lab in [('primary','Primary electrons'),('gamma_free','Main track')]:
            s=es[f'{prefix}_d{axis}_mm'];ax.stairs(s['counts'][:],s['edges'][:],label=lab)
        ax.set_title(f'Event extent in {axis.upper()}');ax.set_xlabel(f'Δ{axis.upper()} (mm)');ax.set_ylabel('Events');annotation(ax,es[f'gamma_free_d{axis}_mm'],'mm');ax.legend(loc='lower right',fontsize=8)
    block.append(save(fig,sample+'_extents.png'))
    block.append('<p>Extent of the main track per event. The statistics box in each extent plot refers to the main track.</p><table><tr><th>Axis</th><th>Mean (mm)</th><th>RMS (mm)</th><th>95th percentile (mm)</th><th>Maximum (mm)</th></tr>')
    for axis in 'xyz':
        s=es[f'gamma_free_d{axis}_mm'];block.append('<tr><td>'+axis.upper()+'</td>'+''.join(f'<td>{s[k][()]:.2f}</td>' for k in ['mean','rms','p95','max'])+'</tr>')
    block.append('</table>')
    block.append('<h3>Outside-frame gamma deposits surviving the energy cut</h3>')
    if sel['surviving_outside_gamma_deposits'][()]:
        s=g['surviving_outside_gamma_statistics/energy_keV'];fig,ax=plt.subplots(figsize=(9,4))
        ax.stairs(s['counts'][:],s['edges'][:],fill=True,alpha=.6);ax.set_title('Outside-frame gamma deposits in surviving events');ax.set_xlabel('E (keV)');ax.set_ylabel('Gamma deposits');ax.set_xlim(0,sel['gamma_energy_cut_keV'][()])
        ax.text(.97,.97,f"N={s['n'][()]}\nmean E={s['mean'][()]:.3g} keV\nRMS={s['rms'][()]:.3g} keV",transform=ax.transAxes,ha='right',va='top')
        block.append(save(fig,sample+'_surviving_outside_energy.png'))
    else:
        block.append('<p>No outside-frame gamma deposits survive in this sample.</p>')
    for dataset,heading in [('example_event_numbers','Five examples before the gamma cut'),('surviving_outside_example_numbers','Surviving events with gamma deposits outside the frame')]:
      if len(g[dataset])==0:continue
      block.append(f'<h3>{heading}</h3>')
      for eid in g[dataset][:]:
        d=g['event_details'][str(eid)];xyz=np.array([d[k+'_mm'][:] for k in 'xyz']);primary=d['primary'][:].astype(bool);gamma=d['gamma_particle_id'][:]!=0
        fig,axes=plt.subplots(1,2,figsize=(12,5));lo=d['main_track_low_mm'][:]-20;hi=d['main_track_high_mm'][:]+20
        masks=[(primary,'Primary electron hits','navy'),(~primary&~gamma,'Delta rays','teal'),(gamma,'Gamma deposits','darkorange')]
        for ax,(a,b) in zip(axes,[(0,1),(0,2)]):
            for mask,lab,col in masks:ax.scatter(xyz[a,mask],xyz[b,mask],s=9,label=lab,c=col,alpha=.65)
            from matplotlib.patches import Rectangle
            ax.add_patch(Rectangle((lo[a],lo[b]),hi[a]-lo[a],hi[b]-lo[b],fill=False,linestyle='--',edgecolor='red'))
            ax.set_title('XYZ'[a]+'–'+'XYZ'[b]+' projection');ax.set_xlabel('XYZ'[a]+' (mm)');ax.set_ylabel('XYZ'[b]+' (mm)');ax.update_datalim([[lo[a],lo[b]],[hi[a],hi[b]]]);ax.autoscale_view();ax.margins(.12);ax.set_aspect('equal',adjustable='box')
        axes[0].legend(fontsize=7)
        if dataset=='surviving_outside_example_numbers':
            deposits=g['surviving_outside_gamma_deposits'];energies=deposits['energy_keV'][:][deposits['event_number'][:]==eid]
            block.append('<p>MC event '+str(eid)+': outside-frame gamma deposited energies = '+', '.join(f'{e:.2f}' for e in energies)+' keV.</p>')
        status='accepted' if eid in g['after_gamma_cut_events/event_number'][:] else 'rejected'
        fig.suptitle(f'{label}, MC event {eid} ({status}), vertex Z={d["vertex_z_mm"][()]:.1f} mm. Frame: main-track box + 20 mm each side')
        block.append(save(fig,f'{sample}_{dataset}_event_{eid}.png'))
    sections.append('\n'.join(block))
page='''<!doctype html><html lang="en"><meta charset="utf-8"><title>ITACA — MC fast-trigger prototype</title><style>ul{max-width:1100px;line-height:1.6}body{font:16px system-ui;max-width:1500px;margin:35px auto;padding:0 25px;color:#203047}h1,h2{color:#16365c}img{width:100%;margin:18px 0;border:1px solid #ddd}table{border-collapse:collapse}td,th{padding:9px 18px;border:1px solid #ddd}p{max-width:1100px;line-height:1.6}</style><h1>Fast trigger: full-sample pre-selection</h1><ul><li>Full bb0nu and single-electron NEXUS samples with 1 mm hits.</li><li>Sample the primary vertex uniformly in Z from 0 to 1500 mm (seed 42). Apply the same translation to every hit.</li><li>Require all deposited hits inside the detector Z interval.</li><li>Require the primary vertex in 50–1450 mm.</li><li>The first distributions use events passing these two cuts.</li><li>Apply the gamma veto: reject an event if any gamma deposit outside the 20 mm frame has total deposited energy above the threshold set by TOML.</li><li>The surviving-deposit distribution and examples apply all three cuts.</li></ul><p>Definitions.</p><ul><li>Primary track: hits of the primary electron or electrons.</li><li>Main track: the primary track plus its delta rays.</li><li>Gamma deposit: all the energy deposited by a gamma and its descendants.</li><li>Distance to track: minimum 3D distance between a gamma hit and a primary hit.</li><li>Distance outside frame: maximum separation beyond the frame in R or Z, with R = √(dx² + dy²).</li></ul><p>In the event plots the dashed frame is the main-track bounding box enlarged by 20 mm on each side.</p>'''+''.join(sections)+'</html>'
(OUT/'fast_trigger_gamma_outside_frame.html').write_text(page)
print(OUT/'fast_trigger_gamma_outside_frame.html')

print(f"Presentation seconds: {time.perf_counter()-START:.3f}")
