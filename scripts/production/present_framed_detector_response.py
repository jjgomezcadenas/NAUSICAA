"""Make detector control plots from Julia HDF5 summaries; no transport calculation."""
from pathlib import Path
import sys, html, tomllib
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'results/readout/bb0nu_full_ion3mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42/detector_response/controls.h5'

def main(path):
    out=path.parent;plots=out/'plots';plots.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'figure.dpi':130})
    with h5py.File(path) as f:
        card=tomllib.loads(f['card'].asstr()[()])
        sample_label='Single electrons' if card.get('sample')=='single_electrons' else 'bb0nu'
        sections=[]
        def save(fig,name):
            fig.tight_layout();fig.savefig(plots/name);plt.close(fig)
            sections.append(f'<img src="plots/{name}" loading="lazy" alt="{html.escape(name)}">')
        def hist(ax,group,title,xlabel,source=None):
            h=(f if source is None else source)[group];ax.stairs(h['counts'][:],h['edges'][:],fill=True,alpha=.65)
            ax.set_title(title);ax.set_xlabel(xlabel);ax.set_ylabel('Events')
            ax.text(.98,.98,f"N = {h['entries'][()]:,}\nMean = {h['mean'][()]:.4g}\nRMS = {h['rms'][()]:.4g}",
                    ha='right',va='top',transform=ax.transAxes,fontsize=9,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
            ax.grid(alpha=.15)
        n=int(f['selection/accepted_events'][()])
        table='<table><tr><th>Selection</th><th>Events</th></tr>'
        inp=int(f['selection/input_events'][()]);outside=int(f['selection/outside_detector_events'][()]);v=int(f['selection/outside_vertex_events'][()])
        for title,count in [('Input',inp),('Inside detector',inp-outside),('Vertex within 50–1450 mm',inp-outside-v),('After gamma veto',n)]:
            table+=f'<tr><td>{title}</td><td>{count:,}</td></tr>'
        table+='</table><p>Efficiencies use the detector-contained sample. All following distributions use events passing the gamma veto.</p>'
        sections.append('<h2>Selection</h2>'+table)
        energy_path=out/'energy_controls.h5'
        if energy_path.exists():
            with h5py.File(energy_path) as ef:
                sections.append('<h2>Signal inside, outside and total</h2><p>Total = inside + outside the acquisition. Each entry is one accepted event.</p>')
                fig,axs=plt.subplots(2,3,figsize=(16,8))
                for row,(arm,label) in enumerate([('ions','Ions'),('photons','SiPM photons, attachment corrected')]):
                    for col,key in enumerate(['inside','outside','total']):
                        hist(axs[row,col],f'{arm}/corrected_counts/{key}',f'{label}: {key}',label,source=ef)
                save(fig,'signal_accounting.png')
                sections.append(f'<h2>Energy inside, outside and total</h2><p>Ions: N × W. Photons: N × W / (EL yield × optical acceptance). Photon energies below are corrected for attachment using exp(drift time / lifetime), with the simulated event drift time and {card["electron_lifetime_ms"]:g} ms lifetime.</p>')
                fig,axs=plt.subplots(2,3,figsize=(16,8))
                for row,(arm,label) in enumerate([('ions','Ions'),('photons','Photons, attachment corrected')]):
                    for col,key in enumerate(['inside','outside','total']):
                        ax=axs[row,col]
                        ax.set_title(f'{label}: {key}');ax.set_xlabel('Energy (keV)');ax.set_ylabel('Events')
                        if key=='outside':
                            h=ef[f'{arm}/outside_display']
                            ax.stairs(h['counts'][:],h['edges'][:],fill=True,alpha=.6);ax.set_xlim(0,10);ax.set_yscale('log')
                            text=f"N = {h['entries'][()]:,}\nMean = {h['mean'][()]:.3f} keV\nAbove 10 keV: {h['above_range'][()]}"
                        else:
                            h=ef[f'{arm}/peak_fits/{key}_keV']
                            ax.stairs(h['plot_counts'][:],h['plot_edges'][:],fill=True,alpha=.5)
                            ax.plot(h['curve_x'][:],h['curve_y'][:],color='darkred',lw=2)
                            ax.axvline(h['lower'][()],ls=':',color='gray');ax.axvline(h['upper'][()],ls=':',color='gray')
                            text=f"σ/μ = {h['relative_sigma_percent'][()]:.3f}%\nμ = {h['mu'][()]:.2f} keV\nσ = {h['sigma'][()]:.2f} ± {h['sigma_error'][()]:.2f} keV\nFit: {h['lower'][()]:.1f}–{h['upper'][()]:.1f} keV\nN fit = {h['entries'][()]:,}"
                        ax.text(.97,.97,text,ha='right',va='top',transform=ax.transAxes,fontsize=9,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
                        ax.grid(alpha=.15)
                save(fig,'energy_accounting.png')
                sections.append('<p>Gaussian peak fits use a fixed interval: median ± 2.5 × robust width (1.4826 × MAD). The likelihood includes the finite fit interval. Outside-energy plots show 0–10 keV; the mean uses all events.</p>')
                table='<table><tr><th>Readout</th><th>Energy</th><th>σ/μ, ±2 widths (%)</th><th>σ/μ, ±2.5 widths (%)</th><th>σ/μ, ±3 widths (%)</th></tr>'
                for arm,label in [('ions','Ions'),('photons','Photons, attachment corrected')]:
                    for key in ['inside','total']:
                        h=ef[f'{arm}/peak_fits/{key}_keV']
                        table+=f'<tr><td>{label}</td><td>{key}</td><td>{h["narrow/relative_sigma_percent"][()]:.3f}</td><td>{h["relative_sigma_percent"][()]:.3f}</td><td>{h["wide/relative_sigma_percent"][()]:.3f}</td></tr>'
                sections.append('<h3>Fit-interval check</h3>'+table+'</table>')
                sections.append('<h2>Photon total and drift distance</h2><p>The total includes photons outside the acquisition. Left: before attachment correction. Right: after correction.</p>')
                fig,axs=plt.subplots(1,2,figsize=(12,4))
                for ax,kind,title in zip(axs,['raw','corrected'],['Before attachment correction','After attachment correction']):
                    ax.scatter(ef['electron_drift_cm'][:],ef[f'photons/events_{kind}/total_keV'][:],s=3,alpha=.3)
                    ax.set_title(title);ax.set_xlabel('Electron drift distance (cm)');ax.set_ylabel('Total photon energy estimate (keV)')
                    h=ef[f'photons/{kind}_energy/total_keV']
                    ax.text(.98,.98,f"Mean = {h['mean'][()]:.2f} keV\nSD = {h['std'][()]:.2f} keV",transform=ax.transAxes,ha='right',va='top',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
                save(fig,'attachment_control.png')
        sections.append('<h2>Frame losses</h2><p>Whole channels are acquired. XY losses reach channels outside the acquired set. Time losses fall outside the acquisition interval. The total loss counts their overlap once.</p>')
        table='<table><tr><th>Readout</th><th>Loss</th><th>Affected events</th><th>Events (%)</th><th>Signal lost (%)</th></tr>'
        for arm,label in [('ions','Ions'),('photons','SiPM photons')]:
            for key,name in [('xy','Outside channels'),('time','Outside time interval'),('any','Either')]:
                affected=int(f[f'{arm}/summary/affected_{key}'][()]);frac=f[f'{arm}/summary/fraction_{key}'][()]
                table+=f'<tr><td>{label}</td><td>{name}</td><td>{affected:,}</td><td>{100*affected/n:.3f}</td><td>{100*frac:.6f}</td></tr>'
        sections.append(table+'</table>')
        sections.append('<h2>Collected signal</h2><p>Each event contributes one mean and one maximum over positive pad/bucket or SiPM/bucket cells. Photon counts are corrected for attachment and are before PDE.</p>')
        fig,axs=plt.subplots(2,2,figsize=(12,8))
        with h5py.File(energy_path) as ef:
            for row,(arm,label) in enumerate([('ions','Ions'),('photons','Photons, attachment corrected')]):
                for col,(key,name) in enumerate([('mean','Event mean per channel/bucket'),('maximum','Event maximum per channel/bucket')]):
                    hist(axs[row,col],f'{arm}/corrected_signal/{key}',f'{label}: {name.lower()}',label,source=ef)
        save(fig,'signal.png')
        sections.append('<h2>Acquisition size</h2>')
        fig,axs=plt.subplots(2,4,figsize=(18,8))
        for row,(arm,label,unit) in enumerate([('ions','Pads','ms'),('photons','SiPMs','μs')]):
            for col,(key,title,xlabel) in enumerate([('acquired_channels','Acquired channels','Channels'),('active_fraction','Channels receiving signal','Fraction of acquired channels'),('duration','Acquisition duration',f'Duration ({unit})'),('buckets','Buckets per channel','Buckets')]):
                hist(axs[row,col],f'{arm}/histograms/{key}',f'{label}: {title.lower()}',xlabel)
        save(fig,'acquisition.png')
        for group,heading in [('examples','First five accepted events'),('largest_losses','Largest fractional losses'),('frame_checks','Frame check: event 604')]:
            if group not in f:
                continue
            sections.append(f'<h2>{heading}</h2>')
            sections.append(f'<p>Photons are attachment corrected. Electron XZ displays sum {card.get("electron_display_rebin",1)} time buckets; original response vectors retain {card["electron_bucket_us"]:g} μs sampling. Dashed black: event frame. Blue: acquired XY channel boundaries. Red crosses: excluded signal. Positions are relative to the MC vertex.</p>')
            for id in f[group]:
                g=f[f'{group}/{id}'];lo=g['frame_low'][:];hi=g['frame_high'][:]
                fig,axs=plt.subplots(2,2,figsize=(12,9))
                for row,(arm,label) in enumerate([('ions','Ions'),('photons','Photons')]):
                    data=g[arm]
                    for col,proj in enumerate(['xy','xz']):
                        ax=axs[row,col];x=data[f'{proj}_x'][:];y=data[f'{proj}_y'][:];c=data[f'{proj}_count'][:]
                        if arm=='photons':
                            with h5py.File(energy_path) as ef:
                                c=ef[f'photons/example_counts/{group}/{id}/{proj}_count'][:]
                        if len(c):
                            sc=ax.scatter(x,y,c=c,s=14,marker='s',norm=LogNorm(vmin=max(1,int(c.min())),vmax=max(2,int(c.max()))),cmap='viridis')
                            fig.colorbar(sc,ax=ax,label=label)
                        ox=data[f'outside_{proj}_x'][:];oy=data[f'outside_{proj}_y'][:]
                        ax.scatter(ox,oy,c='red',marker='x',s=20,linewidths=.8)
                        d=1 if proj=='xy' else 2
                        ax.add_patch(Rectangle((lo[0],lo[d]),hi[0]-lo[0],hi[d]-lo[d],fill=False,ls='--',ec='black'))
                        if proj=='xy':
                            b=data['acquired_xy_bounds'][:];ax.add_patch(Rectangle((b[0],b[2]),b[1]-b[0],b[3]-b[2],fill=False,ec='royalblue'))
                        ax.set_xlabel('X (mm)');ax.set_ylabel(('Y' if d==1 else 'Z')+' (mm)');ax.set_aspect('equal',adjustable='datalim')
                        width=float(data['display_z_bin_mm'][()]) if 'display_z_bin_mm' in data else None
                        suffix=f' · ΔZ = {width:g} mm' if proj=='xz' and width is not None else ''
                        ax.set_title(f'Event {id} · {label} · {proj.upper()}{suffix}')
                save(fig,f'{group}_{id}.png')
        timing=f['timing']
        sections.append('<h2>Execution</h2><table><tr><th>Quantity</th><th>Time</th></tr>'+''.join(f'<tr><td>{label}</td><td>{timing[key][()]:.4f} s</td></tr>' for key,label in [('response_s','Detector response, all accepted events'),('response_mean_excluding_first_s','Mean per accepted event, excluding first compilation'),('response_median_s','Median per accepted event'),('selection_s','Reading and selection'),('total_s','Complete Julia run, before final file publication')])+'</table>')
        config=f'''<ul><li>{sample_label} · Xe–15% He · {card['pressure_bar']:g} bar · {card['drift_field_V_cm']:g} V/cm.</li>
<li>TopMetal pitch {card['ion_pitch_mm']:g} mm; time-buckets {card['ion_bucket_ms']:g} ms; drift velocity {card['ion_drift_velocity_mm_ms']:g} mm/ms.</li>
<li>SiPMs pitch {card['EL_pitch_mm']:g} mm; time-buckets {card['electron_bucket_us']:g} μs; electron drift velocity {card['electron_drift_velocity_mm_us']:g} mm/μs.</li>
<li>Electron lifetime {card['electron_lifetime_ms']:g} ms; E/P = {card['EL_reduced_field_kV_cm_bar']:g} kV/cm/bar.</li>
<li>Frame margin {card['frame_margin_mm']:g} mm on each side; gamma veto above {card['gamma_energy_cut_keV']:g} keV.</li></ul>'''
        page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>bb0nu detector response · ITACA</title><style>body{font:16px/1.6 system-ui;color:#203343;max-width:1500px;margin:30px auto;padding:20px}h1,h2{color:#185b83}img{width:100%;border:1px solid #ddd;margin:12px 0}table{border-collapse:collapse}td,th{padding:8px 18px;border:1px solid #ddd}a{color:#185b83}</style><a href="../../../../notes/html/index.html">ITACA index</a><h1>bb0nu: detector response</h1>'''+config+''.join(sections)+'</html>'
        page=page.replace('bb0nu detector response',sample_label+' detector response').replace('bb0nu: detector response',sample_label+': detector response')
        (out/'index.html').write_text(page)
        (out.parent/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>ITACA run</title><h1>'+html.escape(card['run_label'])+'</h1><p><a href="detector_response/index.html">Detector-response controls</a></p><p><a href="../../../notes/html/index.html">ITACA index</a></p></html>')
    print(out/'index.html')
if __name__=='__main__':
    main(Path(sys.argv[1]) if len(sys.argv)>1 else DEFAULT)
