#!/usr/bin/env python3
"""Topmetal-style short-injection shaping and averaged 10 ms ion-frame study.

No circuit optimization or pole-zero correction. The published 180 us is an
older-chip shaping time; G=T is our explicit main flat-top assumption, not a
reported 2024 setting. Full waveform and noise use the identical linear filter.
"""
from pathlib import Path
import argparse,json,hashlib,os
import numpy as np
from build_ltspice_study import run
from ion_signal.physics import E_CHARGE
from ion_signal.shaping import (box_average,trapezoid,trapezoid_transfer,
    averaged_frame_transfer,filtered_rms)
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--simulator',type=Path,required=True)
    ap.add_argument('--include-ion-frames',action='store_true',help='Optional later stage; injection alone is the default')
    ap.add_argument('--rise-us',type=float,nargs='+',default=[180,200,250,300,350])
    args=ap.parse_args();sim=args.simulator.resolve()
    if any(t<=0 or not np.isfinite(t) for t in args.rise_us):ap.error('Shaping durations must be finite and positive')
    out=ROOT/'results/ltspice/shaping';decks=out/'decks';decks.mkdir(parents=True,exist_ok=True)
    source=ROOT/'results/ltspice/decks/uniform_N1_r3.cir'
    text=source.read_text();text='\n'.join(l for l in text.splitlines() if not l.lower().startswith(('.tran','.save','.end','is1 ')))+'\n'
    model=ROOT/'LTSpice/180nm_professional.inc'
    (decks/model.name).write_text('\n'.join(l.rstrip() for l in model.read_text().splitlines()).strip()+'\n')
    frames=json.loads((ROOT/'results/ltspice/study.json').read_text())
    report={'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (source,model,Path(__file__).resolve(),ROOT/'src/ion_signal/shaping.py',ROOT/'results/ltspice/noise_400.csv',ROOT/'results/ltspice/noise_200.csv',ROOT/'results/ltspice/study.json',ROOT/'results/ltspice/decks/uniform_r2.pwl',ROOT/'results/ltspice/uniform_N1_r1.csv.gz',ROOT/'results/ltspice/uniform_N1_r2.csv.gz')},
      'simulator':str(sim),'simulator_version':'LTspice 17.2.4 native Mac',
      'conditions':{'external_pad_F':5e-12,'feedback_F':3e-15,'temperature_K':293.15,
       'noise_band_Hz':[.1,1e6],'frame_s':.01,'estimator_nominal_tau_s':.3,
       'pulse_start_s':.005,'pulse_plateau_s':1e-6,'pulse_edge_s':1e-8,
       'main_flat_top':'G=T; scan 0,T,2T','pole_zero_correction':False,
       'averaging':'continuous analog waveform; ADC sampling/quantization excluded'},
      'limitations':['Original transistor dimensions/biases unchanged, including M5.',
        'Supplied model has zero flicker noise; stationary bias-point PSD.',
        '180 us belongs to an older uncalibrated chip result, not the corrected 2024 measurement.',
        'Fixed noiseless-template amplitude time; no maximization of noisy traces.',
        'Centered frame-boundary averages require latency T/2 and mix neighbouring frames.']}
    waves={};checks=[]
    for refine,step in ((1,1e-6),(2,.5e-6)):
        for n in (0,1,1000,16000):
            # PULSE area: amplitude * (plateau + (rise+fall)/2).
            amplitude=n*E_CHARGE/(1e-6+1e-8)
            deck=decks/f'injection_N{n}_r{refine}.cir'
            deck.write_text(text+f'Is1 0 in PULSE(0 {amplitude:.17e} .005 10n 10n 1u 1)\n.tran 0 .04 0 {step:.17e}\n.save V(in) V(out) I(Is1)\n.end\n')
            d,h=run(sim,deck);t=d['time'];v=d['V(out)'];i=d['I(Is1)']
            if n==0:base=(t,v);delta=v-v
            else:delta=v-np.interp(t,*base)
            q=np.trapezoid(i,t)
            if n and abs(q/(n*E_CHARGE)-1)>2e-5:raise ArithmeticError('Pulse charge error')
            waves[(n,refine)]=(t,delta)
            np.savetxt(out/f'injection_N{n}_r{refine}.csv.gz',np.c_[t,i,v,delta],delimiter=',',header='time_s,input_A,output_V,baseline_subtracted_output_V',comments='')
            checks.append({'ions':n,'refinement':refine,'charge_C':float(q),'peak_delta_V':float(-delta.min())})
    report['injection_transients']=checks
    spectra={}
    for density in (200,400):
        a=np.loadtxt(ROOT/f'results/ltspice/noise_{density}.csv',delimiter=',',skiprows=1)
        names=frames[f'noise_{density}_columns'];spectra[density]=(a[:,0],a[:,names.index('V(onoise)')])
    rises=[t*1e-6 for t in args.rise_us]
    report['conditions']['rise_scan_s']=rises
    short=[];ion=[];alpha=np.exp(-.01/.3)
    report['ion_frames_enabled']=args.include_ion_frames
    uniform={r:np.loadtxt(ROOT/f'results/ltspice/uniform_N1_r{r}.csv.gz',delimiter=',',skiprows=1) for r in (1,2)}
    finite_ion_waves={}
    raw_pwl=np.loadtxt(ROOT/'results/ltspice/decks/uniform_r2.pwl')
    def finite_ion_wave(n,refine):
        key=(n,refine)
        if key in finite_ion_waves:return finite_ion_waves[key]
        if n:base=finite_ion_wave(0,refine)
        step=5e-6/refine
        pwl=decks/f'ions_N{n}_r{refine}.pwl'
        np.savetxt(pwl,np.c_[raw_pwl[:,0],n*raw_pwl[:,1]],fmt='%.17e')
        deck=decks/f'ions_N{n}_r{refine}.cir'
        deck.write_text(text+f'Is1 0 in PWL FILE="{pwl.name}"\n.tran 0 .032 0 {step:.17e}\n.save V(in) V(out) I(Is1)\n.end\n')
        d,h=run(sim,deck);t=d['time']-.002;v=d['V(out)']
        delta=v-np.interp(t,*base) if n else v-v
        # Cache the raw baseline for subtraction, signal departures otherwise.
        finite_ion_waves[key]=(t,v if n==0 else delta)
        np.savetxt(out/f'ions_N{n}_r{refine}.csv.gz',np.c_[t,d['I(Is1)'],v,delta],delimiter=',',header='physical_time_s,current_A,output_V,baseline_subtracted_output_V',comments='')
        return finite_ion_waves[key]
    for T in rises:
        for flat_factor in (1,0,2):
            G=flat_factor*T
            H=lambda f:trapezoid_transfer(f,T,G)
            rms=filtered_rms(*spectra[400],H,2*T+G)
            rms_check=filtered_rms(*spectra[400],H,2*T+G,order=16)
            rms_coarse=filtered_rms(*spectra[200],H,2*T+G)
            # All filter support lies within the simulation; locate a deterministic
            # template peak on a grid much finer than the shaping time.
            times=np.arange(.005,.012,.25e-6)
            gains=[]
            for refine in (1,2):
                t,v=waves[(1,refine)];y=trapezoid(t,v,times,T,G)
                j=np.argmin(y);gains.append(float(-y[j]))
            tpeak=float(times[j]);g=gains[-1]
            measured={}
            for n in (1000,16000):
                t,v=waves[(n,2)]
                measured[n]=float(-trapezoid(t,v,tpeak,T,G)/n/g)
            row={'rise_s':T,'flat_s':G,'output_noise_rms_V':rms,'signal_V_per_e':g,
              'enc_e':rms/g,'amplitude_time_s':tpeak,'gain_ratio_at_injected_counts':measured,
              'gain_refinement_relative':abs(gains[0]/g-1),
              'noise_grid_relative':abs(rms_coarse/rms-1),'noise_quadrature_relative':abs(rms_check/rms-1)}
            short.append(row)
            if flat_factor==1:print('SHORT',json.dumps(row),flush=True)
        if not args.include_ion_frames: continue
        # A short trapezoid applied directly to a 10 ms ramp measures local current,
        # not total frame charge. Instead average the two frame-boundary voltages.
        Hframe=lambda f:averaged_frame_transfer(f,T,.01,alpha)
        noise=filtered_rms(*spectra[400],Hframe,.01+T)
        noise_check=filtered_rms(*spectra[200],Hframe,.01+T)
        noise_order_check=filtered_rms(*spectra[400],Hframe,.01+T,order=16)
        gain=[]
        for refine in (1,2):
            a=uniform[refine];avg=box_average(a[:,0],a[:,-1],np.array([0,.01])+T/2,T)
            gain.append(float(-(avg[1]-alpha*avg[0])))
        g=gain[-1];threshold=5*noise/g
        a=uniform[2]
        point=-(np.interp(.01,a[:,0],a[:,-1])-alpha*np.interp(0.,a[:,0],a[:,-1]))
        actual=[]
        candidate=int(np.floor(threshold))+1
        for n in range(max(1,candidate-1),candidate+2):
            responses=[]
            for refine in (1,2):
                t,v=finite_ion_wave(n,refine)
                av=box_average(t,v,np.array([0.,.01])+T/2,T)
                responses.append(float(-(av[1]-alpha*av[0])))
            actual.append({'ions':n,'snr':responses[-1]/noise,
                'calibrated_recovered_ions':responses[-1]/g,
                'gain_relative_to_single_ion':responses[-1]/(n*g),
                'timestep_refinement_relative':abs(responses[0]/responses[-1]-1)})

        # Shift the same packet by +- one frame; include baseline before the
        # modeled induction start and measured residual through +30 ms.
        a=uniform[2];t=a[:,0];v=a[:,-1]
        t=np.r_[-.03,t];v=np.r_[0.,v]
        cross={}
        for shift in (-.01,0,.01):
            av=box_average(t+shift,v,np.array([0,.01])+T/2,T)
            cross[str(shift)]=float(-(av[1]-alpha*av[0])/g)
        row={'averaging_s':T,'output_noise_rms_V':noise,'signal_V_per_ion':g,'enc_equivalent_ions':noise/g,
            'N_for_snr5_linear':threshold,'minimum_integer_N_linear':int(np.floor(threshold))+1,
            'gain_refinement_relative':abs(gain[0]/g-1),'noise_grid_relative':abs(noise_check/noise-1),
            'neighbour_packet_response_relative':cross,'signal_retained_vs_point_samples':g/point,
            'noise_quadrature_relative':abs(noise_order_check/noise-1),'finite_count_checks':actual}
        ion.append(row);print('FRAME',json.dumps(row),flush=True)
    report['short_injection']=short;report['ion_frames']=ion
    for row in ion:
        passing=sorted(c['ions'] for c in row['finite_count_checks'] if c['snr']>5)
        row['verified_minimum_integer_N']=passing[0] if passing and passing[0]-1 in [c['ions'] for c in row['finite_count_checks']] else None
    if ion and any(r['verified_minimum_integer_N'] is None for r in ion):
        raise ArithmeticError('Expand finite-count bracket to verify the threshold')
    if ion and max(c['timestep_refinement_relative'] for r in ion for c in r['finite_count_checks'])>1e-3:
        raise ArithmeticError('Finite-count transistor response not converged')
    if max(r['gain_refinement_relative'] for r in short+ion)>1e-3:raise ArithmeticError('Signal not converged')
    if max(r['noise_grid_relative'] for r in short+ion)>1e-3:raise ArithmeticError('Noise grid not converged')
    if max(r['noise_quadrature_relative'] for r in short+ion)>1e-5:raise ArithmeticError('Quadrature not converged')
    # Independent time-domain covariance and transfer checks are in test_shaping.py.
    (out/'study.json').write_text(json.dumps(report,indent=2)+'\n')
    with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
        fig,axs=plt.subplots(1,3 if ion else 2,figsize=(11 if ion else 7.5,3.3),layout='constrained')
        for factor in (0,1,2):
            rows=sorted([r for r in short if np.isclose(r['flat_s'],factor*r['rise_s'],rtol=1e-8,atol=1e-12)],key=lambda r:r['rise_s'])
            axs[0].plot([r['rise_s']*1e6 for r in rows],[r['enc_e'] for r in rows],'o-',label=f'G={factor}T')
        axs[0].axhspan(120,140,color='0.7',alpha=.3);axs[0].axhline(200,c='0.5',ls='--')
        axs[0].set(xlabel='Shaping rise T (us)',ylabel='Short-injection ENC (electrons RMS)');axs[0].legend(frameon=False)
        t,v=waves[(1000,2)];times=np.arange(.004,.007,.5e-6)
        for T in (180e-6,200e-6):axs[1].plot((times-.005)*1e3,-trapezoid(t,v,times,T,T)*1e3,label=f'{T*1e6:g} us')
        axs[1].set(xlabel='Time since injection (ms)',ylabel='Shaped signal (mV)',title='1,000 injected electrons');axs[1].legend(frameon=False)
        if ion:
            rows=sorted(ion,key=lambda r:r['averaging_s'])
            axs[2].plot([r['averaging_s']*1e6 for r in rows],[r['enc_equivalent_ions'] for r in rows],'o-')
            axs[2].axhline(315,color='0.5',ls='--',label='Point samples');axs[2].legend(frameon=False)
            axs[2].set(xlabel='Boundary average duration (us)',ylabel='10 ms frame ENC (ions RMS)')
        fig.savefig(out/'shaping_study.pdf',metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(ROOT/'notes/physics/work/shaping/study.png',dpi=160)
    if ion:
        table=[]
        for r in ion:
            shot=next(s for s in short if np.isclose(s['rise_s'],r['averaging_s']) and np.isclose(s['flat_s'],s['rise_s']))
            cross=r['neighbour_packet_response_relative']
            table.append([r['averaging_s']*1e6,shot['enc_e'],r['enc_equivalent_ions'],r['verified_minimum_integer_N'],
                100*r['signal_retained_vs_point_samples'],100*cross['-0.01'],100*cross['0.01']])
        a=np.array(table)
        np.savetxt(out/'shaping_scan.csv',a,delimiter=',',header='averaging_us,injection_ENC_e_RMS,ion_frame_ENC_ions_RMS,minimum_N_snr5,mean_signal_retained_percent,previous_frame_response_percent,next_frame_response_percent',comments='')
        with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
            fig,axs=plt.subplots(1,2,figsize=(8,3.4),layout='constrained')
            axs[0].plot(a[:,0],a[:,1],'o-',label='Short injection, G=T')
            axs[0].plot(a[:,0],a[:,2],'s-',label='10 ms ion frame')
            axs[0].set(xlabel='Averaging duration T (us)',ylabel='Input-referred ENC\n(elementary charges RMS)')
            axs[1].plot(a[:,0],a[:,5],'o-',label='Previous frame')
            axs[1].plot(a[:,0],a[:,6],'s-',label='Next frame')
            axs[1].set(xlabel='Averaging duration T (us)',ylabel='Equal-charge neighbour response (%)')
            for ax in axs:ax.legend(frameon=False);ax.grid(alpha=.2)
            fig.savefig(out/'shaping_scan.pdf',metadata={'CreationDate':None,'ModDate':None})
            fig.savefig(ROOT/'notes/physics/work/shaping/scan.png',dpi=150)
    rows=[r for r in short if np.isclose(r['flat_s'],r['rise_s'])]
    summary='''# First Topmetal-style shaping study

Generated by scripts/build_shaping_study.py. Circuit dimensions/biases unchanged;
5 pF external pad plus intrinsic capacitances; 20 C; stationary output noise
integrated over 0.1 Hz–1 MHz. Supplied models include no flicker noise. These
results are not an exact replication of the 2024 Topmetal measurement.

The trapezoid is a difference of two moving averages of duration T, separated
by T+G. Main choice G=T is a study assumption. No pole-zero correction is used.
Pulse injection is 1 us wide with 10 ns rise/fall; its integral is checked.
ENC uses the actual filtered signal gain, not 1/Cf. The amplitude is read at a
fixed noiseless-template peak time, not maximized over a noisy waveform.

| T (us) | G (us) | Short-injection ENC (electrons RMS) | 10 ms frame ENC (ions RMS) | Verified minimum N for S/N > 5 |
|---:|---:|---:|---:|---:|
'''
    for r in rows:
        f=next((x for x in ion if x['averaging_s']==r['rise_s']),None)
        if f is None:
            summary+=f"| {r['rise_s']*1e6:g} | {r['flat_s']*1e6:g} | {r['enc_e']:.3f} | Deferred | Deferred |\n"
            continue
        summary+=f"| {r['rise_s']*1e6:g} | {r['flat_s']*1e6:g} | {r['enc_e']:.3f} | {f['enc_equivalent_ions']:.3f} | {f['verified_minimum_integer_N']} |\n"
    summary+=('\nThe ion-frame stage was run with `--include-ion-frames`. Finite-count SPICE checks bracket S/N=5 at every scanned duration; see study.json.\n' if ion else '\nIon frames are deferred in this injection-only run.\n')
    summary+='''

For 10 ms arrivals, the frame estimator uses centered moving averages at the two
frame boundaries, followed by the same nominal 0.3 s discharge correction.
Its normalization is calibrated with the uniform finite-gap ion signal. A short
trapezoid pulse height must not be interpreted as total charge of a 10 ms ramp.
The continuous-average result assumes sufficiently frequent voltage access;
one point read every 10 ms cannot implement these averages digitally. A circuit
aperture, repeated reads near the boundary, or faster continuous sampling is needed.

The extra latency is T/2. Neighbour-frame response is in study.json; it includes
induction before collection and physical feedback history, not only smoothing.
Pulse-amplitude linearity is checked at 1,000 and 16,000 electrons against the
one-electron gain. Noise is evaluated about the zero-signal operating point.
Finite sampling, ADC noise/quantisation and calibrated transistor flicker noise
remain to be included. Manuscript text has not been changed.
'''
    if ion:
        summary+='\n## Direct finite-count verification\n\n'
        for row in ion:
            if row['verified_minimum_integer_N'] is None:continue
            n=row['verified_minimum_integer_N'];checks={r['ions']:r for r in row['finite_count_checks']}
            summary+=f"At {row['averaging_s']*1e6:g} us, N={n-1} gives S/N={checks[n-1]['snr']:.6f}; N={n} gives S/N={checks[n]['snr']:.6f}. Averaging retains {100*row['signal_retained_vs_point_samples']:.4f}% of the point-sample mean signal.\n\n"
            cross=row['neighbour_packet_response_relative']
            summary+=f"Equal-size packets shifted by one frame contribute {100*cross['-0.01']:.3f}% (preceding packet) and {100*cross['0.01']:.3f}% (following packet) of the target-packet response. These are separate small-signal response tests, not a fully simulated stochastic track.\n\n"
    (out/'README.md').write_text(summary.rstrip()+'\n')

if __name__=='__main__':main()
