"""Present saved Julia timing summaries; no detector calculations."""
from pathlib import Path
import sys, html
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(sys.argv[1]); out=p.parent
with h5py.File(p) as f:
    stages=[('Input and frame selection','input_s'),('Detector simulation','detector_total_s'),('CSA + ADC, including preparation','ion_electronics_total_s'),('SiPM + ASIC','sipm_asic_s'),('Ion reconstruction, including preparation','ion_reconstruction_total_s'),('SiPM reconstruction and rebinning','electron_reconstruction_total_s'),('Other event operations','other_s'),('Complete event','event_s')]
    rows=[]
    for label,key in stages:
        g=f['summary/'+key]
        rows.append(f'<tr><td>{label}</td>'+''.join(f'<td>{1000*g[n][()]:.3f}</td>' for n in ['mean','median','p90'])+'</tr>')
    g=f['summary/event_s']; mean=g['mean'][()]; med=g['median'][()]; p90=g['p90'][()]
    fig,ax=plt.subplots(figsize=(8,4.5));ax.hist(1000*f['events/event_s'][:],bins=20,color='#286c91')
    ax.set(xlabel='Complete event processing time (ms)',ylabel='Events')
    ax.text(.97,.96,f'Entries = 100\nMean = {1000*mean:.1f} ms\nMedian = {1000*med:.1f} ms\n90th percentile = {1000*p90:.1f} ms',ha='right',va='top',transform=ax.transAxes)
    fig.tight_layout();fig.savefig(out/'event_time.png',dpi=150);plt.close(fig)
    first=g['first_half_mean'][()];second=g['second_half_mean'][()]
    checks={n:int(f['summary/'+n+'_total'][()]) for n in ['failed_pads','nonlinear_pads','invalid_pe']}
    extras={n:float(f[n][()]) for n in ['warmup_s','input_setup_s','electronics_setup_s','elapsed_before_write_s']}
    gc=float(f['summary/gc_s/mean'][()])
    detail=''.join(f'<li>{label}: {1000*f["summary/"+key+"/mean"][()]:.3f} ms/event.</li>' for label,key in [('Ion detector simulation','ion_detector_s'),('Electron detector simulation','electron_detector_s'),('Pair counts and detector settings','pairs_settings_s'),('CSA settings preparation','electronics_preparation_s'),('CSA + ADC calculation','csa_adc_s'),('Ion reconstruction preparation','ion_reconstruction_preparation_s'),('Ion reconstruction calculation','ion_reconstruction_s'),('SiPM reconstruction','electron_reconstruction_s'),('Electron rebinning','electron_rebin_s')])
page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>ITACA · Readout timing</title><style>body{{font:17px/1.6 system-ui;max-width:1050px;margin:40px auto;padding:0 20px;color:#203343}}table{{border-collapse:collapse;width:100%}}td,th{{padding:9px;text-align:left;border-bottom:1px solid #ddd}}img{{max-width:100%}}</style><a href="../../../../notes/html/event_analysis/index.html">Event analysis</a><h1>Readout timing: 100 bb0nu events</h1>
<ul><li>TopMetal pitch: 3 mm; time buckets: 20 ms. SiPM pitch: 4 mm; time buckets: 250 ns.</li><li>Xe–15% He; frame margin: 20 mm; electron lifetime: 30 ms.</li><li>SiPM calculation: 64 intervals, fully recovered initial state. The limit of 16 avalanches is a computational choice.</li><li>One warm-up event, followed by 100 accepted events. One CPU thread for BLAS. All channels in the frame include noise.</li><li>Event numbers, positions and frames come from the saved fast-trigger selection. This benchmark starts by reading accepted-event hits; it does not repeat sample selection.</li></ul>
<h2>Time per event</h2><table><tr><th>Stage</th><th>Mean (ms)</th><th>Median (ms)</th><th>90th percentile (ms)</th></tr>{''.join(rows)}</table>
<p>Total event processing: {100*mean:.2f} s. First 50 events: {1000*first:.1f} ms/event; last 50: {1000*second:.1f} ms/event. Memory cleanup is included in these times ({1000*gc:.1f} ms/event on average).</p><img src="event_time.png" alt="Event processing times">
<h2>Calculation details</h2><ul>{detail}</ul>
<h2>Run preparation</h2><ul><li>Source checksum and saved event selection: {extras['input_setup_s']:.3f} s.</li><li>Common electronics preparation: {extras['electronics_setup_s']:.3f} s.</li><li>Warm-up event, including first-call compilation: {extras['warmup_s']:.3f} s.</li><li>Elapsed from script start to summary-file writing: {extras['elapsed_before_write_s']:.2f} s, including module loading and compilation. Julia process startup and report generation are additional.</li></ul>
<h2>Checks and output</h2><ul><li>Nonlinear reconstructed pad vectors: {checks['nonlinear_pads']}.</li><li>Pad vectors flagged as not converged: {checks['failed_pads']}.</li><li>Invalid rebinned PE cells: {checks['invalid_pe']}.</li><li>Detector vectors, ADC vectors and reconstructed vectors remain in memory. No zero suppression or event-file writing is timed.</li><li>Noise accuracy validation remains a separate study.</li></ul><p>Numerical results: <a href="controls.h5">controls.h5</a>.</p></html>'''
(out/'index.html').write_text(page)
print(out/'index.html');print(checks);print(extras);print('halves',first,second)
