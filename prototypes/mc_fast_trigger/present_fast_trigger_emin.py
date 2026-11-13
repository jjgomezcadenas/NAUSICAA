"""Plot saved energy-dependent S/N; numerical calculations remain in Julia."""
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path(__file__).resolve().parent;O=B/'results';P=O/'fast_trigger_emin';P.mkdir(exist_ok=True)
with h5py.File(O/'fast_trigger_emin.h5') as f:
    E=f['energy_keV'][:];L=f['drift_mm'][:];K=f['summed_cells'][:]
    fig,axes=plt.subplots(1,3,figsize=(15,4.5))
    for ax,drift in zip(axes,np.unique(L)):
        for k in [1,8,32,128]:
            m=(L==drift)&(K==k);y=f['sn_noise'][:][m];x=E[m]
            ax.plot(x,y,'o-',label=f'{k} SiPM–buckets')
            for xx,yy in zip(x,y):ax.annotate(f'{yy:.1f}',(xx,yy),xytext=(0,5),textcoords='offset points',ha='center',fontsize=8)
        ax.set_title(f'Electron drift {drift:g} mm');ax.set_xlabel('Deposited energy (keV)');ax.set_ylabel('S/N');ax.set_xticks([1,3,5,10]);ax.set_yscale('log');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(P/'signal_noise.png',dpi=140);plt.close(fig)
    worst=(L==1500)&(K==8)
    values=', '.join(f'{e:g} keV → {s:.1f}' for e,s in zip(E[worst],f['sn_noise'][:][worst]))
    page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Fast trigger: minimum energy</title><style>body{{font:16px system-ui;max-width:1500px;margin:35px auto;padding:0 25px;color:#203047}}p{{max-width:1100px;line-height:1.6}}img{{width:100%}}</style><a href="index.html">Fast-trigger studies</a><h1>Fast trigger: S/N at 1, 3, 5 and 10 keV</h1><p>Current geometry: 4 mm pitch, 2.5 mm holes, 3 × 3 mm² SiPMs; Xe–15% He, E/P = 2, PDE 0.32. ASIC noise is {f['asic_noise_rms_pe'][()]:.4g} PE RMS per SiPM–250 ns bucket. Dark noise adds 0.175 PE² per bucket.</p><p>S/N = mean detected PE / sqrt(K × (ASIC noise² + 0.175)), where K is the number of acquired SiPM–buckets summed for the deposit. The mean dark signal is subtracted. K is a sensitivity parameter, not a measured occupancy. Signal fluctuations are stored separately in the HDF5.</p><img src="fast_trigger_emin/signal_noise.png"><p>At 1500 mm drift, summing eight SiPM–buckets: <strong>{values}</strong>.</p><p>The calculation uses the production EL yield, optical acceptance, 1 mm/µs electron drift and 15 ms lifetime. It assumes linear recovered PE and independent noise. No clipping is expected for these mean signals; concentrating the full 10 keV signal instantaneously on one sensor gives an estimated 4.8% occupancy loss, reduced by recovery and signal sharing. Noise from baseline estimation and searching many possible clusters is not included.</p><p>Cards, mean photon/PE counts and all Fano, attachment, EL, collection, PDE, dark and ASIC variance components are saved in fast_trigger_emin.h5. Calculation: {f['elapsed_seconds'][()]:.2f} s.</p></html>'''
    with h5py.File(O/'fast_trigger_ion_emin.h5') as ions:
        ei=ions['energy_keV'][:];ki=ions['summed_cells'][:];arr=ions['arrangement'].asstr()[:]
        fig,axes=plt.subplots(1,2,figsize=(12,4.5))
        for ax,mode,title in zip(axes,['independent_pads','consecutive_buckets'],['One bucket per independent pad','Consecutive buckets on one pad']):
            for k in [1,8,32,128]:
                mask=(arr==mode)&(ki==k);x=ei[mask];y=ions['sn_noise'][:][mask]
                ax.plot(x,y,'o-',label=f'{k} pad–buckets')
                for xx,yy in zip(x,y):ax.annotate(f'{yy:.2f}',(xx,yy),xytext=(0,5),textcoords='offset points',ha='center',fontsize=8)
            ax.set_title(title);ax.set_xlabel('Deposited energy (keV)');ax.set_ylabel('S/N');ax.set_xticks([1,3,5,10]);ax.set_yscale('log');ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.tight_layout();fig.savefig(P/'ion_signal_noise.png',dpi=140);plt.close(fig)
        iontext=f'<h2>Ions: linear CSA, 4 mm pitch / 40 ms buckets</h2><p>Mean ions = E / 21.9 eV; no attachment. Fano variance = 0.15 × mean ions. The calibrated single-bucket noise is {ions["single_bucket_enc_ions"][()]:.2f} ions RMS. S/N below uses additive CSA noise; the HDF5 also saves signal/total RMS including Fano fluctuations.</p><img src="fast_trigger_emin/ion_signal_noise.png"><p>Independent pads: summed variance = K × single-bucket variance. Consecutive buckets on one pad: summed variance = sum of every entry of the K × K temporal covariance matrix. Correlations allow noise to cancel in a consecutive-bucket sum. These arrangements are sensitivity cases; they are not predictions of the deposit geometry.</p><p>The exact saved linear reconstructed-noise covariance is taken from the matched 4 mm production, with 5 ms voltage averaging and the three-frame inverse. Central windows of a 160-bucket covariance matrix avoid acquisition boundaries. Assume full charge collection and containment in the summed cells. No SPICE or nonlinear response is run; no ADC quantization is added. Calculation: {ions["elapsed_seconds"][()]:.2f} s.</p>'
        page=page.replace('</html>',iontext+'</html>')
    (O/'fast_trigger_emin.html').write_text(page)
print(O/'fast_trigger_emin.html')
