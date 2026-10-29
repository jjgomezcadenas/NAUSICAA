#!/usr/bin/env python3
"""Reproduce explicit-arrival finite-gain CSA, noise and SPICE benchmarks."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
from dataclasses import asdict
import hashlib
import json
import tomllib
import platform
import scipy
import numpy as np
from ion_signal.physics import Detector,E_CHARGE,square_phi
from ion_signal.readout import csa_response
from ion_signal.realistic_csa import (CSA,Noise,response,ions_waveform,
    sampled_voltage,recovered_charge,frame_charge,noise_metrics,pwl_current)

ROOT=Path(__file__).resolve().parents[2]


def spice_comparator(path,pwl_name,model,stop_s,maxstep_s):
    """LTspice-style comparator; generated but not executed/validated in SPICE.

    The amplifier pole uses a noiseless VCCS and capacitor, avoiding an
    unintended resistor thermal-noise contribution in the amplifier model.
    The real Rf contributes its passive thermal noise in a noise analysis.
    """
    text=f'''* Illustrative linear CSA comparator; NOT a transistor design or PDK model.
* LTspice-style PWL file syntax. Verify in the chosen simulator before use.
* All node voltages are small-signal deviations; input current flows 0 -> in.
* Do not add Cin again if a transistor model already includes its components.
.param A0={model.open_loop_gain:.15g} wp={2*np.pi*model.amplifier_pole_Hz:.15g}
.param Cin={model.input_capacitance_F:.15g} Cf={model.feedback_capacitance_F:.15g} Rf={model.feedback_resistance_ohm:.15g}
Iion 0 in PWL FILE="{pwl_name}"
Cinput in 0 {{Cin}}
Cfeedback in out {{Cf}}
Rfeedback in out {{Rf}}
* vo=A(s)*(en-vi), en=0 for signal benchmark.
Vseries en 0 0
Eamp drive 0 en in {{A0}}
Gpole pole 0 pole drive 1
Cpole pole 0 {{1/wp}}
Ebuffer out 0 pole 0 1
.options reltol=1e-6 abstol=1e-22 vntol=1e-10 numdgt=15
.tran 0 {stop_s:.15g} 0 {maxstep_s:.15g}
.save V(in) V(out) I(Iion)
* Halve maxstep until sampled voltages and recovered charges converge.
* The PWL input is interpolated: refine its generation step independently.
* No invented phase margin or stochastic waveform is supplied.
.end
'''
    path.write_text(text)


def configured_arrivals(case):
    kind=case['arrival_kind']
    if kind=='simultaneous':
        return np.full(case['ions'],case['collection_time_s'])
    if kind=='evenly_spaced':
        return case['arrival_start_s']+(np.arange(case['ions'])+.5)*case['arrival_duration_s']/case['ions']
    if kind=='groups':
        if len(case['ions_per_group'])!=len(case['collection_times_s']):
            raise ValueError('Group counts and collection times differ')
        return np.concatenate([np.full(n,t) for n,t in zip(case['ions_per_group'],case['collection_times_s'])])
    raise ValueError('Unknown configured arrival kind')


def main():
    path=ROOT/'config/realistic_csa.toml'
    config=tomllib.loads(path.read_text())
    model=CSA(**config['electronics']); noise=Noise(**config['noise'])
    d=Detector.from_toml(ROOT/'config/detector.toml')
    readout=config['readout']; num=config['numerics']
    samples=np.asarray(readout['sample_times_s'])
    interval=readout['sample_interval_s']; aperture=readout['sampling_aperture_s']
    if not np.allclose(np.diff(samples),interval,rtol=1e-12,atol=0):
        raise ValueError('Configured samples must be equally spaced')
    metrics=noise_metrics(model,noise,interval,aperture,num['noise_relative_tolerance'])
    tighter=noise_metrics(model,noise,interval,aperture,num['noise_relative_tolerance']/10)
    convergence=abs(tighter['recovered_charge_rms_C']/metrics['recovered_charge_rms_C']-1)
    if convergence>1e-5:
        raise RuntimeError('Noise quadrature has not converged')
    out=ROOT/'results/realistic_csa'; spice=out/'spice'; spice.mkdir(parents=True,exist_ok=True)
    gen=ROOT/'notes/physics/generated'; gen.mkdir(parents=True,exist_ok=True)
    cases=config['cases']
    qg=E_CHARGE*(1-square_phi(d.gap_m,d.pad_width_m/2))
    target=int(np.argmin(abs(samples[1:]-readout['target_interval_end_s'])))
    if not np.isclose(samples[target+1],readout['target_interval_end_s']):
        raise ValueError('Target interval end must be a sample time')
    records=[]; rows=[]
    for case in cases:
        name=case['id']; label=case['label']; arrivals=configured_arrivals(case)
        wave=ions_waveform(d,arrivals,num['start_s'],num['end_s'],num['max_step_s'])
        state=response(wave,model)
        ideal=csa_response(wave,model.feedback_capacitance_F,model.nominal_discharge_time_s)
        voltage=sampled_voltage(wave.edges_s,state[:,1],samples,aperture)
        ideal_sample=sampled_voltage(wave.edges_s,ideal,samples,aperture)
        qhat=recovered_charge(voltage,model,interval)
        qi=frame_charge(wave,samples)
        # This calibration scales the ENTIRE fixed arrival template, including
        # previous-frame arrivals. It is not a count estimate for arbitrary tracks.
        unit_template=qhat[target]/len(arrivals)
        snr=abs(qhat[target])/metrics['recovered_charge_rms_C']
        threshold=readout['snr_threshold']*metrics['recovered_charge_rms_C']/abs(unit_template)
        prior_arrivals=arrivals[arrivals<=samples[target]]
        target_arrivals=arrivals[(arrivals>samples[target])&(arrivals<=samples[target+1])]
        if len(target_arrivals)==0 or np.any(arrivals>samples[target+1]):
            raise ValueError('This benchmark requires target arrivals and no future arrivals')
        history_leakage=0.; history_true=0.
        if len(prior_arrivals):
            history=ions_waveform(d,prior_arrivals,num['start_s'],num['end_s'],num['max_step_s'])
            history_v=response(history,model)[:,1]
            history_q=recovered_charge(sampled_voltage(history.edges_s,history_v,samples,aperture),model,interval)
            history_leakage=float(history_q[target])
            history_true=float(frame_charge(history,samples)[target])
        corrected=float(qhat[target]-history_leakage)
        current_frame_coefficient=corrected/len(target_arrivals)
        current_frame_threshold=readout['snr_threshold']*metrics['recovered_charge_rms_C']/abs(current_frame_coefficient)
        np.savetxt(out/(name+'_waveform.csv'),np.c_[wave.edges_s,
                   np.r_[wave.current_A,0.],state,ideal],delimiter=',',comments='',
                   header='time_s,current_A_left_bin_average,input_node_V,output_V,ideal_output_V')
        np.savetxt(out/(name+'_samples.csv'),np.c_[samples[1:],voltage[1:],ideal_sample[1:],qi,qhat],
                   delimiter=',',comments='',header='sample_time_s,output_V,ideal_output_V,true_interval_charge_C,recovered_charge_C')
        arrival_file=out/(name+'_arrivals.csv')
        np.savetxt(arrival_file,arrivals,delimiter=',',comments='',header='collection_time_s')
        t,current=pwl_current(wave)
        pwl=spice/(name+'.pwl'); np.savetxt(pwl,np.c_[t,current],fmt='%.15e')
        exported=np.loadtxt(pwl); exported_charge=float(np.trapezoid(exported[:,1],exported[:,0]))
        error=abs(exported_charge/wave.total_charge_C-1)
        if error>1e-10:
            raise RuntimeError('PWL serialization failed charge check')
        dt=float(wave.edges_s[1]-wave.edges_s[0])
        pole_time=1/max(abs(model.poles_per_s()))
        suggested_step=min(dt/2,pole_time/10)
        spice_comparator(spice/(name+'.cir'),pwl.name,model,t[-1],suggested_step)
        # Independent physical-bin refinement at identical sample/aperture times.
        fine=ions_waveform(d,arrivals,num['start_s'],num['end_s'],num['max_step_s']/2)
        fine_state=response(fine,model)
        fine_q=recovered_charge(sampled_voltage(fine.edges_s,fine_state[:,1],samples,aperture),model,interval)
        sample_error=float(np.max(abs(fine_q-qhat))/(len(arrivals)*qg))
        if sample_error>2e-4:
            raise RuntimeError('Sampled-charge convergence failed')
        record=dict(id=name,label=label,ions=len(arrivals),arrival_configuration=case,
            arrivals_file=str(arrival_file.relative_to(ROOT)),
            arrivals_sha256=hashlib.sha256(arrival_file.read_bytes()).hexdigest(),
            input_charge_C=wave.total_charge_C,
            expected_finite_gap_charge_C=len(arrivals)*qg,peak_output_magnitude_V=float(-np.min(state[:,1])),
            peak_input_node_magnitude_V=float(np.max(abs(state[:,0]))),
            target_interval_s=[float(samples[target]),float(samples[target+1])],
            true_target_interval_charge_C=float(qi[target]),recovered_target_charge_C=float(qhat[target]),
            target_bias_C=float(qhat[target]-qi[target]),target_snr=float(snr),
            prior_arrival_ions=len(prior_arrivals),target_frame_arrival_ions=len(target_arrivals),
            history_target_recovered_charge_C=history_leakage,
            history_true_target_interval_charge_C=history_true,
            history_subtracted_target_charge_C=corrected,
            history_subtracted_target_bias_C=corrected-float(qi[target]-history_true),
            target_frame_template_response_C_per_ion=float(current_frame_coefficient),
            target_frame_ion_count_at_snr_threshold=float(current_frame_threshold),
            history_subtracted_target_snr=abs(corrected)/metrics['recovered_charge_rms_C'],
            history_subtraction='Subtract expected response of known prior-arrival template; uncertainty in estimating that history is not included',
            unit_arrival_template_response_C=float(unit_template),
            template_ion_count_at_snr_threshold=float(threshold),
            calibration='Scale the entire fixed arrival template, including any previous-frame history; not an arbitrary-track count estimator',
            sampled_charge_step_refinement_error_fraction_total_input=sample_error,
            physical_time_offset_s=float(wave.edges_s[0]),numerical_step_s=dt,
            spice_maxstep_starting_s=float(suggested_step),pwl_charge_relative_error=error,
            pwl_sha256=hashlib.sha256(pwl.read_bytes()).hexdigest(),
            spice_status='Generated behavioral comparator, not executed in SPICE; timestep and PWL refinement remain required')
        records.append(record)
        rows.append(f"{name.replace('_', ' ')} & {len(arrivals)} & {record['peak_output_magnitude_V']*1e3:.4g} & {qhat[target]/E_CHARGE:.3g} & {snr:.3g}"+r' \\')
    f=np.geomspace(noise.low_frequency_Hz,noise.high_frequency_Hz,10001)
    hi,he=model.transfer(f)
    np.savetxt(out/'noise_spectrum.csv',np.c_[f,noise.output_psd(model,f),noise.output_psd(model,f,aperture),abs(hi),abs(he)],
               delimiter=',',comments='',header='frequency_Hz,output_PSD_V2_per_Hz,aperture_filtered_PSD_V2_per_Hz,current_transfer_ohm,voltage_noise_gain')
    payload=dict(model='Finite-gain one-pole CSA; identical finite-gap pad trajectories',
       config_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
       detector_config_sha256=hashlib.sha256((ROOT/'config/detector.toml').read_bytes()).hexdigest(),
       detector=asdict(d),numerics=num,
       environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
       electronics=asdict(model),noise=asdict(noise),readout=readout,provenance=config['provenance'],
       aperture_convention='Trailing rectangular average ending at each stated sample time',
       poles_per_s=[dict(real=float(x.real),imag=float(x.imag)) for x in model.poles_per_s()],
       charge_gain_reduction_fast_amplifier=model.charge_gain_reduction(),
       noise_metrics=metrics,noise_quadrature_relative_change=convergence,cases=records,
       omitted='Transistor nonlinearities, saturation, extra amplifier poles, source correlations, event tracks and foundry-calibrated parameters')
    (out/'benchmarks.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    values={'RealisticGain':f'{model.open_loop_gain:g}',
            'RealisticInputCapacitance':f'{model.input_capacitance_F*1e15:g}',
            'RealisticPoleFrequency':f'{model.amplifier_pole_Hz:g}',
            'RealisticDischargeTime':f'{model.nominal_discharge_time_s:g}',
            'RealisticGainReduction':f'{model.charge_gain_reduction():.4f}',
            'RealisticSampleNoise':f"{metrics['sample_voltage_rms_V']*1e6:.3g}",
            'RealisticChargeNoise':f"{metrics['induced_charge_ENC_e']:.3g}"}
    (gen/'realistic_csa_values.tex').write_text('% Generated illustrative assumptions, NOT validated circuit performance.\n'+''.join(fr'\newcommand{{\{k}}}{{{v}}}'+'\n' for k,v in values.items()))
    (gen/'realistic_csa_table.tex').write_text('\n'.join([r'\begin{tabular}{lrrrr}',r'\toprule',
        r'Arrival case & Ions & Peak $|V|$ (mV) & $\widehat Q/e$ & S/N \\',r'\midrule',*rows,r'\bottomrule',r'\end{tabular}','']))
    print('Noise:',json.dumps(metrics))
    for r in records:
        print(r['id'],'peak_V=',r['peak_output_magnitude_V'],'target_SNR=',r['target_snr'],
              'refinement_error=',r['sampled_charge_step_refinement_error_fraction_total_input'])
    print('Created',out/'benchmarks.json')

if __name__=='__main__':
    main()
