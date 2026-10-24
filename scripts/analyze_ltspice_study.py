#!/usr/bin/env python3
"""Integrate circuit spectra with frame estimator; summarize convergence/results."""
from pathlib import Path
import json, os, tomllib
import numpy as np
from ion_signal.physics import E_CHARGE
from ion_signal.spice import charge_noise
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    out=ROOT/'results/ltspice'; report=json.loads((out/'study.json').read_text())
    cfg=tomllib.loads((ROOT/'config/realistic_csa.toml').read_text())
    dt=cfg['readout']['sample_interval_s'];cf=cfg['electronics']['feedback_capacitance_F'];tau=cf*cfg['electronics']['feedback_resistance_ohm']
    report['measurement']={'frame_interval_s':dt,'sampling_aperture_s':0.,'feedback_capacitance_F':cf,'nominal_estimator_tau_s':tau,'external_pad_capacitance_F':5e-12,'temperature_K':293.15,'noise_band_Hz':[.1,1e6]}
    noise=[]
    for density in (200,400):
        a=np.loadtxt(out/f'noise_{density}.csv',delimiter=',',skiprows=1)
        names=report[f'noise_{density}_columns']; f=a[:,0]
        values={name:a[:,j] for j,name in enumerate(names)}
        sigma=charge_noise(f,values['V(onoise)'],dt,cf,tau)
        check=charge_noise(f,values['V(onoise)'],dt,cf,tau,order=16)
        components={name:charge_noise(f,asd,dt,cf,tau)/E_CHARGE for name,asd in values.items()
                    if name.startswith('V(') and name not in ('V(onoise)',) and not any(s in name for s in ('.',))}
        flicker={name:float(np.max(asd)) for name,asd in values.items() if '1overf' in name}
        noise.append(dict(points_per_decade=density,charge_rms_e=sigma/E_CHARGE,
            quadrature_relative_difference=abs(check/sigma-1),components_rms_e=components,
            flicker_max_output_asd_V_sqrtHz=flicker))
    report['noise']=noise
    for row in noise:
        if abs(sum(v*v for v in row['components_rms_e'].values())/row['charge_rms_e']**2-1)>1e-6:
            raise ArithmeticError('Device noise variances do not sum to total')
    for row in noise:
        row['component_variance_percent']={name:100*(v/row['charge_rms_e'])**2 for name,v in row['components_rms_e'].items()}
    r={x['ions']:x for x in report['transients'] if x['refinement']==2}
    g=r[1]['gain_e_per_ion']; sigma=noise[-1]['charge_rms_e'];threshold=5*sigma/g
    report['threshold']={'linearized_ions_for_snr5':threshold,'linearized_minimum_integer_ions':int(np.floor(threshold))+1,
        'charge_noise_rms_e_in_estimator_units':sigma,'equivalent_input_ions_rms':sigma/g,'mean_estimator_e_per_ion':g,
        'noise_grid_relative_difference':abs(noise[0]['charge_rms_e']/sigma-1),
        'maximum_timestep_relative_difference':max(abs(x['estimated_charge_C']/r[x['ions']]['estimated_charge_C']-1) for x in report['transients'] if x['refinement']==1 and x['ions']),
        'snr_at_simulated_counts':{n:x['estimated_charge_C']/E_CHARGE/sigma for n,x in r.items() if n}}
    tight={x['ions']:x for x in report['transients'] if x['refinement']==3}
    report['threshold']['tolerance_relative_difference']=max(abs(x['estimated_charge_C']/r[n]['estimated_charge_C']-1) for n,x in tight.items() if n)
    verified=sorted(n for n,x in tight.items() if x['estimated_charge_C']/E_CHARGE/sigma>5)
    report['threshold']['verified_minimum_ions']=verified[0] if verified and verified[0]-1 in tight else None
    report['threshold']['tight_snr']={n:x['estimated_charge_C']/E_CHARGE/sigma for n,x in tight.items() if n}
    if report['threshold']['tolerance_relative_difference']>1e-3: raise ArithmeticError('Tolerance not converged')
    if report['threshold']['noise_grid_relative_difference']>1e-3: raise ArithmeticError('Noise grid not converged')
    if max(x['quadrature_relative_difference'] for x in noise)>1e-5: raise ArithmeticError('Noise quadrature not converged')
    if report['threshold']['maximum_timestep_relative_difference']>1e-3: raise ArithmeticError('Transient not converged')
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
        fig,axs=plt.subplots(1,3,figsize=(11,3.3),layout='constrained')
        for n in (68,300,1000):
            a=np.loadtxt(out/f'uniform_N{n}_r2.csv.gz',delimiter=',',skiprows=1)
            axs[0].plot(a[:,0]*1e3,-a[:,-1]*1e3,label=f'{n} ions')
        axs[0].set(xlabel='Time (ms)',ylabel='Output change magnitude (mV)',title='5 pF external pad')
        axs[0].legend(frameon=False)
        a=np.loadtxt(out/'noise_400.csv',delimiter=',',skiprows=1);names=report['noise_400_columns']
        for name in ('V(onoise)','V(m1)','V(mf)'):
            if name in names: axs[1].loglog(a[:,0],a[:,names.index(name)],label=name)
        axs[1].set(xlabel='Frequency (Hz)',ylabel=r'Output noise (V/$\sqrt{\mathrm{Hz}}$)',title='Model-default noise')
        axs[1].legend(frameon=False)
        ns=[n for n in r if n];snrs=[report['threshold']['snr_at_simulated_counts'][n] for n in ns]
        axs[2].loglog(ns,snrs,'o-');axs[2].axhline(5,color='0.5',ls='--')
        axs[2].set(xlabel='Ions per frame',ylabel='Signal / noise',title='10 ms frame')
        fig.savefig(out/'circuit_study.pdf',metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(ROOT/'notes/physics/work/ltspice-review/circuit_study.png',dpi=150)
    text=f'''# LTspice CSA study

Generated by `scripts/build_ltspice_study.py` and `scripts/analyze_ltspice_study.py`.
Native LTspice 17.2.4; original source files preserved. See study.json for hashes,
assumptions and all simulated counts. External pad capacitance is 5 pF in addition
to transistor capacitances. Temperature is 20 C; original reproduction uses 27 C.
Noise is integrated over 0.1 Hz to 1 MHz with point samples at the frame boundaries.

- Uniform mean arrivals during 10 ms; central grounded-plane pad; finite-gap charge.
- One-ion mean recovered charge: {g:.6g} e in the paper's uncalibrated estimator.
- Noise: {sigma:.6g} e RMS in those same estimator units.
- Direct transistor simulation, with tighter tolerances, brackets S/N=5 at {report['threshold']['verified_minimum_ions']} ions (see adjacent counts in analysis.json).
- Linearized S/N > 5 threshold: N > {threshold:.6g}; minimum integer {int(np.floor(threshold))+1}.
- Maximum frame-charge change on halving PWL spacing and timestep: {100*report['threshold']['maximum_timestep_relative_difference']:.5g}%.
- Noise change on doubling spectral density: {100*report['threshold']['noise_grid_relative_difference']:.5g}%.

The threshold uses stationary small-signal device noise about the DC bias point.
The `1overf` columns report the noise actually included by the supplied models;
zero values must not be taken as evidence that a fabricated circuit has no flicker noise.
The generated `circuit_study.pdf` is a diagnostic, not a manuscript figure.

Cf=3 fF and nominal tau=0.3 s remain the paper estimator parameters, not fitted
circuit parameters. The recovered fraction is not charge-collection efficiency.
This is a model calculation, not foundry or experimental validation.
'''
    (out/'README.md').write_text(text)
    print(json.dumps(report['threshold'],indent=2));print(json.dumps(noise[-1],indent=2))
if __name__=='__main__':main()
