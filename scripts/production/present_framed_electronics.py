#!/usr/bin/env python3
"""Plot saved Julia electronics controls; no detector or electronics calculation."""
from pathlib import Path
import sys, html
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'results/readout/bb0nu_full_ion3mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42/electronics_100events_efficient/controls.h5'
def main(path):
    output=path.parent; plots=output/'plots';plots.mkdir(exist_ok=True)
    blocks=[]
    with h5py.File(path) as f, h5py.File(output/'plot_controls.h5') as reduced:
        ev=f['events']; rows=[]
        for i,n in enumerate(ev['event_number']):
            rows.append('<tr>'+''.join(f'<td>{x}</td>' for x in [int(n),int(ev['ion_channels'][i]),int(ev['ion_samples'][i]),int(ev['sipm_channels'][i]),int(ev['sipm_buckets'][i]),f"{ev['ion_total_s' if 'ion_total_s' in ev else 'ion_s'][i]:.3f}",f"{ev['sipm_s'][i]:.3f}",int(ev['nonlinear_channels'][i]),int(ev['adc_below'][i]+ev['adc_above'][i]),int(ev['asic_clipped_buckets'][i])])+'</tr>')
        blocks.append('<h2>Events and timing</h2><p>Times cover electronics processing; the first event includes compilation. Every channel in the frame is simulated.</p><details><summary>Per-event timing and counts</summary><table><tr>'+''.join(f'<th>{x}</th>' for x in ['Event','Pads','ADC samples/pad','SiPMs','Buckets/SiPM','Complete ion calculation (s)','SiPM + ASIC (s)','Nonlinear pads','ADC out of range','ASIC clipped buckets'])+'</tr>'+''.join(rows)+'</table></details>')
        ion_timing='ion_total_s' if 'ion_total_s' in ev else 'ion_s'
        blocks.append(f'<p>Mean time after the first event: complete ion calculation <strong>{reduced["timing/"+ion_timing+"/mean_excluding_first"][()]:.3f} s/event</strong>; SiPM + ASIC {reduced["timing/sipm_s/mean_excluding_first"][()]:.3f} s/event.</p>')
        if 'efficiency' in reduced:
            q=reduced['efficiency']
            blocks.append('<h2>CSA calculation time</h2><ul>'
                +f'<li>One-time ion setup: {f["ion_setup_s"][()]:.3f} s; SiPM setup: {f["sipm_setup_s"][()]:.3f} s.</li>'
                +f'<li>First ion event, including compilation and initial noise preparation: {ev["ion_total_s"][0]:.3f} s.</li>'
                +f'<li>Same 100 events, excluding the first: {q["before_s"][()]:.3f} → {q["after_s"][()]:.3f} s/event ({q["speedup"][()]:.2f}× faster).</li>'
                +'<li>Ion time includes event preparation, current calculation, current-limit checks, CSA response, correlated noise and ADC conversion. The earlier preparation timer also included loading the SiPM cards.</li>'
                +f'<li>Estimated ion calculation for 8,000 events: {q["predicted_8000_s"][()]/60:.1f} min, plus one-time setup and compilation. Detector transport, SiPM processing and plots are separate.</li>'
                +f'<li>Mean allocations: {reduced["timing/ion_allocated_bytes/mean_excluding_first"][()]/2**20:.2f} MiB/event. Maximum retained CSA preparation and working arrays: {q["max_prepared_MiB"][()]:.2f} MiB.</li>'
                +f'<li>Maximum noiseless-voltage change in saved examples: {q["maximum_voltage_change_V"][()]:.3g} V. {q["changed_adc_samples"][()]} ADC samples changed, by at most {q["maximum_adc_change"][()]} count.</li>'
                +'<li>All 100 events retain the same linear/nonlinear choices and current-limit flags. SiPM examples are unchanged.</li></ul>')
        summary=reduced['summary']
        blocks.append(f'<p>{summary["event_count"][()]} events; {summary["nonlinear_channels"][()]} nonlinear pad vectors; {summary["adc_out_of_range"][()]} ADC samples outside range; {summary["asic_clipped_buckets"][()]} clipped ASIC buckets.</p>')
        if 'events_with_current_exceedances' in summary:
            blocks.append(f'<h2>Peak-current checks</h2><p>{summary["pads_with_current_exceedances"][()]} pads in {summary["events_with_current_exceedances"][()]} events exceed the calibration sample maximum. Largest ratio: {summary["maximum_current_ratio"][()]:.3f}. Events use the existing linear/nonlinear models. No SPICE runs are launched.</p>')
            if 'event_number' in f['current_exceedances']:
                g=f['current_exceedances'];cases=[]
                for j,event in enumerate(g['event_number']):
                    cases.append(f'<tr><td>{int(event)}</td><td>{int(g["pixel_ix"][j])}, {int(g["pixel_iy"][j])}</td><td>{g["peak_current_A"][j]:.5g}</td><td>{g["ratio"][j]:.3f}</td></tr>')
                blocks.append('<table><tr><th>Event</th><th>Pad indices</th><th>Peak current (A)</th><th>Ratio to calibration maximum</th></tr>'+''.join(cases)+'</table>')
        if 'comparison' in reduced:
            q=reduced['comparison/ion_s']
            blocks.append(f'<h2>Noise optimisation</h2><p>Mean CSA + ADC time after the first event: <strong>{q["before_s"][()]:.3f} → {q["after_s"][()]:.3f} s/event</strong> ({q["speedup"][()]:.2f}× faster). Same five events and channels.</p>')
            blocks.append('<p>The noiseless ion examples and SiPM examples are unchanged. Individual ion noise vectors differ because the sampler changed; their covariance is preserved.</p>')
            benchmark=output/'noise_benchmark.h5'
            if benchmark.exists():
                with h5py.File(benchmark) as b:
                    blocks.append(f'<p>Noise alone, {b["pads"][()]} pads with {b["samples"][()]} samples each: {b["reference_median_s"][()]:.4f} → {b["prepared_median_s"][()]:.4f} s ({b["speedup"][()]:.1f}× faster). FFT length: {b["reference_period"][()]} → {b["prepared_period"][()]}.</p>')
            blocks.append('<p><a href="../electronics_noise_reference/index.html">Previous electronics report</a></p>')
        blocks.append('<h2>Validation</h2><table><tr><th>Quantity</th><th>Value</th></tr>'+''.join(f'<tr><td>{label}</td><td>{f["validation/"+key][()]:.6g} {unit}</td></tr>' for key,label,unit in [
            ('linear_maximum_voltage_error_V','Maximum linear voltage residual against SPICE','V'),
            ('nonlinear_maximum_voltage_error_V','Maximum nonlinear voltage residual against SPICE','V'),
            ('nonlinear_cases','Nonlinear SPICE cases',''),
            ('noise_covariance_relative_error','Maximum short-lag covariance error / variance',''),
            ('voltage_noise_rms_V','CSA voltage-noise RMS before filtering','V')])+'</table>')
        if 'convolution_25us_error_V' in f['validation']:
            blocks.append(f'<p>Maximum linear-voltage difference from the exact recurrence: {f["validation/convolution_25us_error_V"][()]:.3g} V with 25 μs convolution steps; {f["validation/convolution_6p25us_error_V"][()]:.3g} V with 6.25 μs steps.</p>')
        blocks.append('<p>The raw voltage noise includes slowly varying components. Its RMS is not the reconstructed ion uncertainty.</p>')
        blocks.append('<h2>Example pad and SiPM vectors</h2><p>For each event, show the channel with the largest input and a channel with zero input, when present. Ion times are relative to the first ADC sample; SiPM times are relative to the first photon bucket. The vectors retain their event-t0 offsets in memory.</p>')
        for event in f['examples']:
            blocks.append(f'<h3>Event {event}</h3>')
            for name,g in f['examples'][event].items():
                fig,axes=plt.subplots(2,1,figsize=(10,6),layout='constrained')
                if name.startswith('ions'):
                    t=g['time_ms'][:];axes[0].step(g['input_time_ms'][:],g['input'][:],where='mid',label='Ions / 20 ms bucket')
                    axes[0].set_ylabel('Ions / bucket');axes[0].set_xlabel('Time (ms)')
                    axes[1].plot(t,g['adc_midpoint_V'][:],lw=.6,label='ADC code-centre voltage')
                    axes[1].plot(t,g['voltage_V'][:],lw=1.2,label='CSA without noise')
                    axes[1].set_ylabel('Voltage (V)');axes[1].set_xlabel('Time (ms)')
                    title=f'Event {event} · pad {int(g["channel"][()])}'
                else:
                    t=g['time_us'][:];axes[0].step(t,g['input'][:],where='mid',label='Photons / 250 ns bucket')
                    axes[0].set_ylabel('Photons / bucket');axes[0].set_xlabel('Time (μs)')
                    axes[1].step(t,g['response_pe'][:],where='mid',label='SiPM + ASIC')
                    axes[1].set_ylabel('Response (PE / bucket)');axes[1].set_xlabel('Time (μs)')
                    title=f'Event {event} · SiPM {int(g["channel"][()])}'
                for ax in axes:ax.grid(alpha=.2);ax.legend(loc='upper right')
                for ax,field in zip(axes,['input','adc_midpoint_V' if name.startswith('ions') else 'response_pe']):
                    stats=reduced[f'statistics/{event}/{name}/{field}']
                    unit='V' if field=='adc_midpoint_V' else ('PE' if field=='response_pe' else ('ions' if name.startswith('ions') else 'photons'))
                    ax.text(.02,.95,f"N = {stats['entries'][()]}\nMean = {stats['mean'][()]:.4g} {unit}\nσ = {stats['std'][()]:.4g} {unit}",transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
                fig.suptitle(title);filename=f'event_{event}_{name}.png';fig.savefig(plots/filename,dpi=130);plt.close(fig)
                blocks.append(f'<img src="plots/{filename}" alt="{html.escape(title)}" loading="lazy">')
        blocks.append('<h2>Noise examples</h2><p>ADC voltage minus simulated noiseless voltage includes CSA noise and quantisation. SiPM examples have zero input photons and include dark activity and ASIC noise.</p>')
        fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
        event=next(iter(reduced['examples']))
        for ax,prefix in zip(axes,['ions','sipm']):
            name=next(n for n in reduced[f'examples/{event}'] if n.startswith(prefix))
            g=reduced[f'examples/{event}/{name}/value'];ax.stairs(g['counts'][:],g['edges'][:],fill=True,alpha=.7)
            unit='V' if prefix=='ions' else 'PE'
            ax.set_xlabel(f'Noise response ({unit})');ax.set_ylabel('Samples' if prefix=='ions' else 'Buckets')
            ax.set_title(f'Event {event} · {name}')
            ax.text(.97,.95,f"N = {g['entries'][()]}\nMean = {g['mean'][()]:.4g} {unit}\nσ = {g['std'][()]:.4g} {unit}",transform=ax.transAxes,va='top',ha='right',fontsize=9)
        fig.savefig(plots/'noise.png',dpi=130);plt.close(fig)
        blocks.append('<img src="plots/noise.png" alt="Noise examples">')
        blocks.append('<h2>Output vector sizes</h2>')
        fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
        for ax,name,title in zip(axes,['adc','asic'],['Ion ADC vectors','SiPM ASIC vectors and flags']):
            g=reduced[f'vector_sizes/{name}/value']
            ax.stairs(g['counts'][:],g['edges'][:],fill=True,alpha=.65)
            ax.set_xlabel('Size per event (MiB)');ax.set_ylabel('Events');ax.set_title(title)
            ax.grid(alpha=.2)
            ax.text(.97,.95,f"N = {g['entries'][()]}\nMean = {g['mean'][()]:.2f} MiB\nMedian = {g['median'][()]:.2f} MiB\nMaximum = {g['maximum'][()]:.2f} MiB",transform=ax.transAxes,ha='right',va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
        fig.savefig(plots/'vector_sizes.png',dpi=130);plt.close(fig)
        blocks.append('<img src="plots/vector_sizes.png" alt="Distributions of ADC and ASIC vector size per event"><p>Full event vectors remain in memory. Saved controls contain summaries and selected channel examples.</p>')
    style='body{font:17px/1.6 system-ui;color:#203343;background:#f3f6f8;margin:0}main{max-width:1150px;margin:auto;padding:30px;background:white}a{color:#185b83}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:8px;text-align:left;border-bottom:1px solid #d8e3e8}img{max-width:100%}h1,h2{color:#185b83}'
    page='<!doctype html><html lang="en"><head><meta charset="utf-8"><title>ITACA · bb0nu electronics</title><style>'+style+'</style></head><body><main><a href="../../../../notes/html/index.html">ITACA index</a><h1>bb0nu · Electronics response</h1><ul><li>TopMetal pitch 3 mm; ion buckets 20 ms.</li><li>SiPMs pitch 4 mm; ASIC buckets 250 ns.</li><li>ADC: 12 bits, 0–1.8 V, 100 μs sampling; additional input noise zero.</li><li>Baseline 10 ms before induced current; final sampling 5 ms after collection.</li></ul>'+''.join(blocks)+'</main></body></html>'
    (output/'index.html').write_text(page)
    print(output/'index.html')
if __name__=='__main__':main(Path(sys.argv[1]) if len(sys.argv)>1 else DEFAULT)
