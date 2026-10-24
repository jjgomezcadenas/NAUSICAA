#!/usr/bin/env python3
"""Uniform-arrival mean response per ion; simultaneous mid-frame cross-check.

N is left free. The finite-gap current is convolved with a uniform probability
of collection over a frame. Multiplication by N gives the mean signal for N
ions; individual stochastic arrival fluctuations are not electronic noise.
"""
from pathlib import Path
from dataclasses import asdict
import hashlib
import json
import os
import tomllib
import numpy as np
from ion_signal.physics import Detector, E_CHARGE, square_phi
from ion_signal.readout import csa_response
from ion_signal.realistic_csa import (CSA, uniform_ion_waveform, ions_waveform,
    response, sampled_voltage, recovered_charge, frame_charge)

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'work/matplotlib-cache'))
os.environ.setdefault('XDG_CACHE_HOME', str(ROOT/'work/cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    config_path=ROOT/'config/realistic_csa.toml'
    detector_path=ROOT/'config/detector.toml'
    cfg=tomllib.loads(config_path.read_text())
    model=CSA(**cfg['electronics']); detector=Detector.from_toml(detector_path)
    dt=cfg['readout']['sample_interval_s']; samples=np.array([0.,dt])
    num=cfg['numerics']
    out=ROOT/'results/realistic_csa'; out.mkdir(parents=True,exist_ok=True)
    gen=ROOT/'notes/physics/generated'; gen.mkdir(parents=True,exist_ok=True)
    records=[]; uniform=None
    qg=E_CHARGE*(1-square_phi(detector.gap_m,detector.pad_width_m/2))
    for name in ('uniform_per_ion','simultaneous_per_ion'):
        estimates=[]
        for divisor in (1,2):
            step=num['max_step_s']/divisor
            if name=='uniform_per_ion':
                wave=uniform_ion_waveform(detector,0.,dt,ions=1.,start_s=num['start_s'],end_s=num['end_s'],max_step_s=step)
            else:
                wave=ions_waveform(detector,[dt/2],num['start_s'],num['end_s'],step)
            finite=response(wave,model)[:,1]
            ideal=csa_response(wave,model.feedback_capacitance_F,model.nominal_discharge_time_s)
            fv=sampled_voltage(wave.edges_s,finite,samples,0.)
            iv=sampled_voltage(wave.edges_s,ideal,samples,0.)
            qhat=float(recovered_charge(fv,model,dt)[0])
            estimates.append(qhat)
        error=abs(estimates[0]/estimates[1]-1)
        if error>1e-5 or abs(wave.total_charge_C/qg-1)>1e-9:
            raise ArithmeticError(f'{name}: refinement or charge-conservation failure')
        truth=float(frame_charge(wave,samples)[0])
        ideal_qhat=float(recovered_charge(iv,model,dt)[0])
        path=out/(name+'_waveform.csv')
        np.savetxt(path,np.column_stack([wave.edges_s,np.r_[wave.current_A,0.],finite,ideal]),
                   delimiter=',',header='time_s,current_A_per_ion_left_bin_average,output_V_per_ion,ideal_output_V_per_ion',comments='')
        records.append(dict(id=name,normalization_ions=1.,sample_times_s=samples.tolist(),
                            input_charge_per_ion_C=truth,finite_gain_estimated_charge_per_ion_C=qhat,
                            ideal_estimated_charge_per_ion_C=ideal_qhat,
                            finite_gain_recovered_fraction=qhat/truth,ideal_recovered_fraction=ideal_qhat/truth,
                            finite_gain_sample_voltage_per_ion_V=fv.tolist(),ideal_sample_voltage_per_ion_V=iv.tolist(),
                            timestep_refinement_relative_difference=error,
                            finite_path_charge_relative_error=abs(wave.total_charge_C/qg-1),
                            waveform_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        if name=='uniform_per_ion':
            uniform=(wave,finite,ideal,fv,iv)
    (out/'arrival_cases_point_samples.json').write_text(json.dumps(dict(
        electronics=asdict(model),sample_interval_s=dt,sampling_aperture_s=0.,noise_included=False,
        capacitance_reference=cfg['capacitance_reference'],
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (config_path,detector_path)},
        assumptions=['N is free; waveforms are normalized to mean response per ion.',
                     'Uniform collection probability on [0, Delta t]; simultaneous mid-frame case is a cross-check.',
                     'Central square pad in grounded plane; constant gap speed; finite-path charge not renormalized.',
                     'Point voltage samples at frame boundaries; zero initial signal; no random-arrival variance.'],
        cases=records),indent=2)+'\n')
    values={
        'ArrivalFrameMs':f'{dt*1e3:g}',
        'ArrivalCfFf':f'{model.feedback_capacitance_F*1e15:g}',
        'ArrivalCinPf':f'{model.input_capacitance_F*1e12:g}',
        'ArrivalGain':f'{model.open_loop_gain:g}',
        'ArrivalPoleKhz':f'{model.amplifier_pole_Hz/1e3:g}',
        'ArrivalRfOhm':f'{model.feedback_resistance_ohm:.0e}',
        'ArrivalTauSeconds':f'{model.nominal_discharge_time_s:g}',
        'TopmetalElectrodeDiameterMm':f"{cfg['capacitance_reference']['electrode_diameter_m']*1e3:g}",
        'ArrivalUniformIdealPercent':f"{100*records[0]['ideal_recovered_fraction']:.2f}",
        'ArrivalUniformRecoveredPercent':f"{100*records[0]['finite_gain_recovered_fraction']:.0f}",
    }
    (gen/'arrival_cases_values.tex').write_text('% Generated by scripts/build_arrival_cases.py; mean per-ion response.\n'+
        ''.join(fr'\newcommand{{\{k}}}{{{v}}}'+'\n' for k,v in values.items()))
    r=records[0]
    rows=[f"Uniform & {r['input_charge_per_ion_C']/E_CHARGE:.4f} & {r['ideal_estimated_charge_per_ion_C']/E_CHARGE:.4f} & {r['finite_gain_estimated_charge_per_ion_C']/E_CHARGE:.4f}"+r' \\']
    (gen/'arrival_cases_table.tex').write_text('\n'.join([r'\begin{tabular}{lrrr}',r'\toprule',
        r'Arrivals & $Q_n/(Ne)$ & Ideal $\widehat Q_n/(Ne)$ & Finite gain $\widehat Q_n/(Ne)$ \\',
        r'\midrule',*rows,r'\bottomrule',r'\end{tabular}','']))
    wave,finite,ideal,fv,iv=uniform
    with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42}):
        fig,axes=plt.subplots(1,2,figsize=(6.6,2.9),constrained_layout=True)
        axes[0].stairs(wave.current_A*1e18,wave.edges_s*1e3,color='#17877E',linewidth=1.4)
        axes[0].set(ylabel=r'$i_N/N$ (aA)',title='Uniform arrivals')
        axes[1].plot(wave.edges_s*1e3,-ideal*1e6,'--',color='#355F85',label='Ideal CSA',linewidth=1.4)
        axes[1].plot(wave.edges_s*1e3,-finite*1e6,color='#B67B24',label='Finite-gain CSA',linewidth=1.5)
        axes[1].plot(samples*1e3,-iv*1e6,'s',color='#355F85',markersize=4)
        axes[1].plot(samples*1e3,-fv*1e6,'o',color='#B67B24',markersize=4)
        axes[1].set(ylabel=r'$-V_N/N$ ($\mu$V)',title='Mean voltage per ion')
        for ax in axes:
            ax.set(xlabel='Time (ms)',xlim=(-.2*dt*1e3,1.5*dt*1e3))
            ax.axvspan(0,dt*1e3,color='0.5',alpha=.055)
            for t in samples:
                ax.axvline(t*1e3,color='0.45',linestyle=':',linewidth=.8)
            ax.grid(axis='y',alpha=.15)
        handles,labels=axes[1].get_legend_handles_labels()
        fig.legend(handles,labels,loc='outside lower center',ncols=2,frameon=False)
        fig.savefig(gen/'arrival_cases.pdf',metadata={'CreationDate':None,'ModDate':None})
        work=ROOT/'notes/physics/work';work.mkdir(parents=True,exist_ok=True)
        fig.savefig(work/'uniform_arrivals.png',dpi=180)
        plt.close(fig)
    print('Generated uniform mean per-ion signals and simultaneous cross-check')
    for r in records:
        print(r['id'], 'g/e=',r['finite_gain_estimated_charge_per_ion_C']/E_CHARGE,
              'recovery=',r['finite_gain_recovered_fraction'])


if __name__=='__main__':
    main()
