#!/usr/bin/env python3
"""Plot reconstructed buckets, grouped by input signal in each bucket."""
from pathlib import Path
import json, sys, tomllib
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'results/readout/bb0nu_full_ion3mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42/reconstruction_bucket_histograms_100/controls.h5'

def render(path):
    output=path.parent;plots=output/'plots';plots.mkdir(exist_ok=True)
    with h5py.File(path) as f:
        if f['schema'][()].decode()!='itaca.reconstruction_controls.v3':
            raise ValueError('Regenerate channel histograms with the current Julia calculation.')
        dc=tomllib.loads(f['cards/detector_card'][()].decode())
        rc=tomllib.loads(f['card'][()].decode())
        ion_bucket=dc['ion_bucket_ms'];pe_bucket=dc['electron_bucket_us']*rc['electrons']['rebin_factor']
        ids=[int(i) for i in f['events/event_number'][:]];images={}
        for event in ids:
            fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
            for col,(arm,device,unit,bucket) in enumerate([('ions','Pads','ions',f'{ion_bucket:g} ms'),('pe','SiPMs','PE',f'{pe_bucket:g} μs')]):
                for row,(group,label,color) in enumerate([('buckets_with_signal','with signal','#c26c20'),('buckets_without_signal','without signal','#176baf')]):
                    ax=axes[row,col];g=f[f'zs/{event}/{arm}/histograms/{group}']
                    entries=int(g['entries'][()]);channels=int(g['channels'][()]);invalid=int(g['invalid_buckets'][()])
                    if entries:
                        ax.stairs(g['counts'][:],g['edges'][:],color=color,fill=True,alpha=.65)
                        ax.set_yscale('log');ax.set_ylim(bottom=.8)
                        text=f"Buckets = {entries}\nMean = {g['mean'][()]:.3g} {unit}\nσ = {g['std'][()]:.3g} {unit}"
                    else:
                        text='No valid buckets'
                    if invalid:text+=f'\nInvalid buckets = {invalid}'
                    ax.text(.97,.96,text,ha='right',va='top',transform=ax.transAxes,fontsize=10,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
                    if group=='buckets_with_signal':
                        ax.set_xlim(right=1200 if arm=='ions' else 3000)
                    ax.set_title(f'{device} · buckets {label}');ax.set_xlabel(f'Reconstructed {unit} per {bucket} bucket');ax.set_ylabel('Buckets');ax.grid(alpha=.15)
            fig.suptitle(f'bb0nu · Event {event}',fontsize=16)
            name=f'event_{event}_channels.png';fig.savefig(plots/name,dpi=120);plt.close(fig);images[str(event)]=name
        options=''.join(f'<option value="{i}">{i}</option>' for i in ids);first=ids[0]
        page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>ITACA · Bucket signal and noise</title><style>body{{font:17px/1.6 system-ui;color:#203343;max-width:1300px;margin:25px auto;padding:0 20px}}img{{width:100%}}select,button{{font:inherit;padding:6px;margin:5px}}.controls{{position:sticky;top:0;background:white;padding:10px}}a{{color:#185b83}}</style>
<a href="../../../../notes/html/index.html">ITACA index</a><h1>bb0nu · Bucket signal and noise · {len(ids)} events</h1>
<ul><li>With signal: that pad–bucket or SiPM–bucket received input ions or photons.</li><li>Without signal: that pad–bucket or SiPM–bucket received no input ions or photons.</li><li>Histograms contain signed reconstructed values. Invalid buckets are reported separately.</li><li>TopMetal: {dc["ion_pitch_mm"]:g} mm pitch, {ion_bucket:g} ms buckets. SiPMs: {dc["EL_pitch_mm"]:g} mm pitch, reconstructed PE rebinned to {pe_bucket:g} μs. No ZS cut.</li></ul>
<div class="controls"><button id="previous">Previous</button><label>Event number <select id="event">{options}</select></label><button id="next">Next</button><a id="png" href="plots/{images[str(first)]}" download>Download PNG</a></div>
<img id="event-image" src="plots/{images[str(first)]}" alt="Event {first}: four bucket histograms">
<script>const images={json.dumps(images)};const select=document.getElementById('event');function show(){{const v=select.value;document.getElementById('event-image').src='plots/'+images[v];document.getElementById('event-image').alt='Event '+v+': four bucket histograms';document.getElementById('png').href='plots/'+images[v];}}select.addEventListener('change',show);document.getElementById('previous').onclick=()=>{{select.selectedIndex=Math.max(0,select.selectedIndex-1);show();}};document.getElementById('next').onclick=()=>{{select.selectedIndex=Math.min(select.options.length-1,select.selectedIndex+1);show();}};</script></html>'''
        (output/'index.html').write_text(page)
    print(output/'index.html')
if __name__=='__main__':render(Path(sys.argv[1]) if len(sys.argv)>1 else DEFAULT)
