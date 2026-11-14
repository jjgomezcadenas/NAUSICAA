#!/usr/bin/env python3
"""Present saved Julia reconstruction/ZS controls. No signal calculation."""
from pathlib import Path
import html,json,sys
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'results/readout/bb0nu_full_ion3mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42/reconstruction/controls.h5'

def draw_hist(ax,g,label,color):
    if g['entries'][()]:
        ax.stairs(g['counts'][:],g['edges'][:],label=label,color=color,lw=1.3)

def annotation(g,unit):
    return f"N = {int(g['entries'][()])}\nMean = {g['mean'][()]:.3g} {unit}\nσ = {g['std'][()]:.3g} {unit}"

def style(ax,x,y):
    ax.set_xlabel(x);ax.set_ylabel(y);ax.grid(alpha=.2)

def render(path):
    output=path.parent;plots=output/'plots';plots.mkdir(exist_ok=True)
    with h5py.File(path) as f,h5py.File(output/'plot_controls.h5') as reduced:
        ids=[int(i) for i in f['events/event_number'][:]];summary=reduced['summary'];event_labels={}
        for event in ids:
            fig,axes=plt.subplots(3,2,figsize=(13,12),layout='constrained')
            for col,(arm,unit,title) in enumerate([('ions','ions','Pads · 20 ms buckets'),('pe','PE','SiPMs · 5 μs buckets')]):
                g=f[f'zs/{event}/{arm}'];hs=g['histograms'];ax=axes[0,col]
                draw_hist(ax,hs['all'],'All valid cells','#697c8c')
                draw_hist(ax,hs['signal'],'Input signal > 0','#c26c20')
                draw_hist(ax,hs['zero_input'],'Zero input','#176baf')
                ax.set_yscale('log');ax.set_title(title);style(ax,f'Reconstructed {unit} / channel–bucket','Cells')
                ax.legend(fontsize=8,loc='upper right')
                ax.text(.03,.95,annotation(hs['all'],unit),transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
                ax=axes[1,col]
                if 'zero_input' in g['noise_histograms']:
                    noise=g['noise_histograms/zero_input'];quiet=g['noise_histograms/quiet_channels']
                    draw_hist(ax,noise,'Zero-input cells','#176baf');draw_hist(ax,quiet,'Channels with no input in event','#a03085')
                    ax.set_yscale('log');ax.legend(fontsize=8,loc='upper right')
                    ax.text(.03,.95,annotation(noise,unit),transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
                ax.set_title('Noise only');style(ax,f'Reconstructed {unit} / channel–bucket','Cells')
                ax=axes[2,col];cuts=g['thresholds'][:];noise_n=g['noise_cells'][()];total=g['input_total'][()]
                if noise_n:ax.plot(cuts,100*g['noise_above'][:]/noise_n,label='Zero-input cells retained',color='#176baf')
                ax.plot(cuts,100*g['input_below'][:]/total,label='Input charge below cut',color='#c26c20')
                ax.legend(fontsize=8);style(ax,f'Threshold ({unit})','Fraction (%)');ax.set_ylim(bottom=0)
                invalid=100*g['input_invalid'][()]/total
                ax.text(.97,.95,f"Invalid cells = {int(g['invalid_cells'][()])}\nInput charge in invalid cells = {invalid:.3g}%",transform=ax.transAxes,ha='right',va='top',fontsize=9)
            fig.suptitle(f'Event {event} · Reconstruction and ZS controls · No cut applied',fontsize=15)
            filename=f'event_{event}_zs.png';fig.savefig(plots/filename,dpi=115);plt.close(fig)
            event_labels[str(event)]=filename
        fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
        for col,(arm,unit,title) in enumerate([('ions','ions','Reconstructed ions'),('pe','PE','Reconstructed PE')]):
            g=reduced[arm];ax=axes[0,col];h=g['event_relative_residual_percent']
            ax.stairs(h['counts'][:],h['edges'][:],fill=True,alpha=.6)
            ax.set_title(title);style(ax,'(Reconstructed − reference) / reference (%)','Events')
            ax.text(.03,.95,annotation(h,'%'),transform=ax.transAxes,va='top',fontsize=9)
            ax=axes[1,col];x=g['thresholds'][:]
            ax.plot(x,100*g['noise_fraction'][:],label='Zero-input cells retained',color='#176baf')
            ax.plot(x,100*g['quiet_fraction'][:],label='Empty-channel cells retained',color='#a03085',ls='--')
            ax.plot(x,100*g['input_lost_fraction'][:],label='Input charge below cut',color='#c26c20')
            ax.fill_between(x,100*g['loss_p16'][:],100*g['loss_p84'][:],alpha=.15,color='#c26c20',label='16–84% of event charge-loss fractions')
            style(ax,f'Threshold ({unit})','Fraction (%)');ax.legend(fontsize=8);ax.set_ylim(bottom=0)
            ax.text(.97,.95,f"Events = {len(ids)}\nInput charge in invalid cells = {100*g['invalid_input_fraction'][()]:.3g}%",transform=ax.transAxes,ha='right',va='top',fontsize=9)
        fig.savefig(plots/'summary.png',dpi=130);plt.close(fig)
        examples=[]
        for event in f['examples']:
            for pad,g in f['examples'][event].items():
                fig,axes=plt.subplots(2,1,figsize=(10,6),layout='constrained')
                axes[0].step(g['time_ms'][:],g['input'][:],where='mid',label='Input ions')
                axes[0].step(g['time_ms'][:],g['reconstructed'][:],where='mid',label='Reconstructed ions')
                style(axes[0],'Time from collection-frame start (ms)','Ions / bucket')
                axes[1].plot(g['average_time_ms'][:],g['measured_drop_V'][:],label='Measured voltage drop')
                axes[1].plot(g['average_time_ms'][:],g['fitted_drop_V'][:],label='Fitted voltage drop')
                style(axes[1],'Time from first ADC sample (ms)','Voltage drop (V)')
                stat=reduced[f'examples/{event}/{pad}']
                axes[0].text(.03,.95,f"Input sum = {stat['input_sum'][()]:.4g} ions\nReconstructed sum = {stat['reconstructed_sum'][()]:.4g} ions",transform=axes[0].transAxes,va='top',fontsize=9)
                axes[1].text(.03,.95,f"RMS residual = {stat['voltage_residual_rms_V'][()]:.3g} V",transform=axes[1].transAxes,va='top',fontsize=9)
                for ax in axes:ax.legend(loc='upper right')
                fig.suptitle(f'Event {event} · {pad}')
                name=f'event_{event}_{pad}.png';fig.savefig(plots/name,dpi=115);plt.close(fig)
                examples.append(f'<img src="plots/{name}" alt="Event {event} {pad}" loading="lazy">')
        rows=''.join(f'<tr><td>{label}</td><td>{summary[key][()]:.4g}</td></tr>' for key,label in [
            ('events','Events'),('nonlinear_pads','Pads reconstructed with nonlinear model'),('failed_pads','Failed ion reconstructions'),
            ('adc_limit_pads','Pads with an endpoint ADC code'),('events_with_invalid_pe','Events with invalid PE buckets'),
            ('invalid_pe_buckets','Invalid rebinned PE buckets'),('clipped_pe_buckets','Rebinned buckets containing ASIC clipping')])
        timings=''.join(f'<tr><td>{label}</td><td>{summary[key+"_mean"][()]:.4f}</td><td>{summary[key+"_median"][()]:.4f}</td></tr>' for key,label in [
            ('ion_preparation_s','Ion reconstruction preparation'),('ion_reconstruction_s','Ion reconstruction'),
            ('electron_reconstruction_s','PE correction and rebinning'),('detector_s','Detector transport'),('electronics_s','Electronics'),('controls_s','ZS controls')])
        options=''.join(f'<option value="{i}">{i}</option>' for i in ids)
        first=ids[0]
        validation=f['electron_inverse_validation']
        criterion=validation['criterion'][()].decode()
        validrows=[]
        for pattern in validation:
            g=validation[pattern]
            if not isinstance(g,h5py.Group):continue
            for j,expected in enumerate(g['expected_pe'][:]):
                validrows.append(f'<tr><td>{pattern}</td><td>{j+1}</td><td>{expected:.3g}</td><td>{g["bias_pe"][j]:.3g}</td><td>{g["standard_error_pe"][j]:.3g}</td><td>{g["tolerance_pe"][j]:.3g}</td></tr>')
        method=f["ion_reconstruction_method"][()].decode() if "ion_reconstruction_method" in f else "Linear inverse and iterative nonlinear fit (reference)"
        page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ITACA · Reconstruction and ZS controls</title>
<style>body{{font:17px/1.6 system-ui;background:#f3f6f8;color:#203343;margin:0}}main{{max-width:1250px;margin:auto;padding:28px;background:white}}h1,h2{{color:#185b83}}a{{color:#185b83}}img{{width:100%;height:auto}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{padding:8px;text-align:left;border-bottom:1px solid #d8e3e8}}select,button{{font:inherit;padding:6px;margin:6px}}.controls{{position:sticky;top:0;background:white;padding:10px;z-index:2}}summary{{cursor:pointer;font-weight:600}}</style></head><body><main>
<a href="../../../../notes/html/index.html">ITACA index</a> · <a href="../../../../notes/html/event_analysis/reco_signal.html">Reconstruction interfaces</a>
<h1>bb0nu · Reconstruction and ZS controls</h1>
<ul><li>TopMetal: 3 mm pitch, reconstructed ions per 20 ms bucket.</li><li>SiPMs: 4 mm pitch, corrected PE summed into 5 μs buckets. Both correspond to 3 mm in Z.</li><li>All channels are reconstructed. Signed values are retained; no ZS cut is applied.</li></ul>
<p>Ion reconstruction: {html.escape(method)}.</p><h2>Summary</h2><table><tr><th>Quantity</th><th>Count</th></tr>{rows}</table>
<img src="plots/summary.png" alt="Event charge residuals and zero suppression threshold scans">
<ul><li>Ion reference: collected ions in valid reconstructed cells.</li><li>PE reference: PDE × photons in valid reconstructed cells. Attachment remains included in both.</li><li>Threshold scans retain valid cells with reconstructed signal greater than the threshold.</li><li>Input charge below cut is the fraction of total input charge in valid cells at or below the threshold. Invalid cells are reported separately.</li></ul>
<h2>Events</h2><ul><li>Input signal &gt; 0: input ions or photons in that cell.</li><li>Zero input: no input in that time cell; the channel can have signal at other times.</li><li>Channels with no input in event: noise-only channels throughout the acquired frame.</li></ul>
<div class="controls"><button id="previous">Previous</button><label>Event number <select id="event">{options}</select></label><button id="next">Next</button><a id="png" href="plots/{event_labels[str(first)]}" download>Download PNG</a></div>
<img id="event-image" src="plots/{event_labels[str(first)]}" alt="Event {first} reconstruction distributions">
<h2>Noise check</h2><p>Variance of summed reconstructed ions in empty pads / predicted variance: <strong>{summary['quiet_pad_variance_ratio'][()]:.4f}</strong>. Mean empty-pad sum: {summary['quiet_pad_mean_total'][()]:.3g} ions. The prediction includes baseline subtraction, voltage correlations and ADC quantisation.</p>
<h2>Timing and memory</h2><p>One-time electronics preparation: {f['setup_s'][()]:.3f} s. Event loop, including the first event: {f['total_s'][()]:.2f} s. Means and medians below exclude the first event.</p>
<table><tr><th>Calculation</th><th>Mean (s/event)</th><th>Median (s/event)</th></tr>{timings}</table>
<p>Mean reconstruction allocations: {summary['reconstruction_allocated_bytes_mean'][()]/2**20:.2f} MiB/event. Maximum prepared reconstruction arrays: {summary['maximum_prepared_MiB'][()]:.2f} MiB.</p>
<h2>PE inverse check</h2><p>{validation['trials'][()]} repetitions per illumination pattern. Criterion: {html.escape(criterion)}. Tested up to 100,000 photons per original ASIC bucket; pilot maximum: {summary['maximum_photons_per_bucket'][()]:.0f} photons.</p>
<details><summary>Bias by pattern and bucket</summary><table><tr><th>Pattern</th><th>Bucket</th><th>Expected PE</th><th>Bias (PE)</th><th>Standard error (PE)</th><th>Tolerance (PE)</th></tr>{''.join(validrows)}</table></details>
<details><summary>Example pad vectors and voltage fits</summary>{''.join(examples)}</details>
<p>Full event vectors remain in memory. This report uses saved numerical controls and selected examples.</p>
</main><script>const images={json.dumps(event_labels)};const select=document.getElementById('event');function show(){{const value=select.value;document.getElementById('event-image').src='plots/'+images[value];document.getElementById('event-image').alt='Event '+value+' reconstruction distributions';document.getElementById('png').href='plots/'+images[value];}}select.addEventListener('change',show);document.getElementById('previous').onclick=()=>{{select.selectedIndex=Math.max(0,select.selectedIndex-1);show();}};document.getElementById('next').onclick=()=>{{select.selectedIndex=Math.min(select.options.length-1,select.selectedIndex+1);show();}};</script></body></html>'''
        (output/'index.html').write_text(page)
    print(output/'index.html')
if __name__=='__main__':render(Path(sys.argv[1]) if len(sys.argv)>1 else DEFAULT)
