#!/usr/bin/env python3
"""Compare 5 pF benchmark with an area-scaled 1.25 pF external load.

Only Cd1 changes. The finite-gap input already describes a 0.5 mm square pad.
1.25 pF is a design estimate, not a geometrical capacitance calculation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import numpy as np
from build_ltspice_study import run
from ion_signal.physics import E_CHARGE
from ion_signal.shaping import box_average, trapezoid, trapezoid_transfer, averaged_frame_transfer, filtered_rms

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--simulator', required=True, type=Path)
    ap.add_argument('--pad-pf', type=float, default=1.25)
    ap.add_argument('--averaging-us', type=float, default=300.)
    args = ap.parse_args()
    if not (np.isfinite(args.pad_pf) and args.pad_pf > 0 and np.isfinite(args.averaging_us) and args.averaging_us > 0):
        ap.error('Capacitance and averaging duration must be finite and positive')
    T = args.averaging_us*1e-6
    alpha = np.exp(-.01/.3)
    source = ROOT/'results/ltspice/decks/uniform_N1_r3.cir'
    model = ROOT/'LTSpice/180nm_professional.inc'
    pwl_source = ROOT/'results/ltspice/decks/uniform_r2.pwl'
    baseline_source = ROOT/'results/ltspice/shaping/study.json'
    out = ROOT/f'results/ltspice/capacitance/pad_{args.pad_pf:g}pF'
    decks = out/'decks'
    decks.mkdir(parents=True, exist_ok=True)
    text = source.read_text()
    assert text.count('Cd1 in 0 5p') == 1
    text = text.replace('Cd1 in 0 5p', f'Cd1 in 0 {args.pad_pf*1e-12:.17e}')
    body = '\n'.join(l for l in text.splitlines() if not l.lower().startswith(('.tran', '.save', '.end', 'is1 ')))+'\n'
    (decks/model.name).write_text('\n'.join(l.rstrip() for l in model.read_text().splitlines()).strip()+'\n')
    raw_pwl = np.loadtxt(pwl_source)
    report = {'external_pad_F': args.pad_pf*1e-12, 'averaging_s': T,
        'conditions': {'temperature_K':293.15, 'noise_band_Hz':[.1,1e6], 'feedback_F':3e-15,
            'frame_s':.01, 'estimator_nominal_tau_s':.3, 'pad_width_m':.0005,
            'weighting_geometry':'square pad in grounded plane; finite 5 mm path, no charge renormalization',
            'capacitance_estimate':'5 pF times (0.5/1)^2; area-scaling hypothesis only',
            'intrinsic_transistor_capacitances':'unchanged; additional to external load',
            'noise':'stationary supplied transistor models, zero flicker; ADC noise excluded'},
        'simulator':str(args.simulator.resolve()), 'simulator_version':'LTspice 17.2.4',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (source,model,pwl_source,baseline_source,Path(__file__).resolve(),ROOT/'src/ion_signal/shaping.py')},
        'transient_charge_checks':[]}
    spectra = {}
    for density in (200,400):
        deck = decks/f'noise_{density}.cir'
        deck.write_text(body+f'Is1 0 in 0 AC 1\n.noise V(out) Is1 dec {density} .1 1Meg\n.end\n')
        d,h = run(args.simulator.resolve(),deck)
        names = list(d)
        a = np.column_stack([d[x].real for x in names])
        np.savetxt(out/f'noise_{density}.csv',a,delimiter=',',header=','.join(names),comments='')
        spectra[density] = (a[:,0],d['V(onoise)'].real)
        report[f'noise_{density}_columns'] = names
        report['flicker_max_asd_V_sqrtHz'] = {x:float(np.max(d[x].real)) for x in names if '1overf' in x}
    waves = {}
    def wave(kind,n,refine):
        key = kind,n,refine
        if key in waves:return waves[key]
        if n:base = wave(kind,0,refine)
        if kind == 'ions':
            pwl = decks/f'ions_N{n}_r{refine}.pwl'
            np.savetxt(pwl,np.c_[raw_pwl[:,0],n*raw_pwl[:,1]],fmt='%.17e')
            source_line = f'Is1 0 in PWL FILE="{pwl.name}"'
            stop,step,offset = .032,5e-6/refine,.002
        else:
            source_line = f'Is1 0 in PULSE(0 {n*E_CHARGE/1.01e-6:.17e} .005 10n 10n 1u 1)'
            stop,step,offset = .04,1e-6/refine,0.
        deck = decks/f'{kind}_N{n}_r{refine}.cir'
        deck.write_text(body+f'{source_line}\n.tran 0 {stop} 0 {step:.17e}\n.save V(in) V(out) I(Is1)\n.end\n')
        d,h = run(args.simulator.resolve(),deck)
        t,v,i = d['time']-offset,d['V(out)'],d['I(Is1)']
        delta = v-np.interp(t,*base) if n else v-v
        waves[key] = (t,v if n==0 else delta)
        np.savetxt(out/f'{kind}_N{n}_r{refine}.csv.gz',np.c_[t,i,v,delta],delimiter=',',
            header='physical_time_s,current_A,output_V,baseline_subtracted_output_V',comments='')
        if n:
            expected = n*(np.trapezoid(raw_pwl[:,1],raw_pwl[:,0]) if kind=='ions' else E_CHARGE)
            error = abs(np.trapezoid(i,t)/expected-1)
            report['transient_charge_checks'].append(dict(kind=kind,ions=n,refinement=refine,relative_error=float(error)))
            assert error < 2e-5, ('charge',kind,n,error)
        return waves[key]
    def frame(t,v):
        av = box_average(t,v,np.array([0.,.01])+T/2,T)
        return float(-(av[1]-alpha*av[0]))
    for kind in ('injection','ions'):
        H = (lambda f:trapezoid_transfer(f,T,T)) if kind=='injection' else (lambda f:averaged_frame_transfer(f,T,.01,alpha))
        span = 3*T if kind=='injection' else .01+T
        noise = filtered_rms(*spectra[400],H,span)
        gains = []
        for refine in (1,2):
            t,v = wave(kind,1,refine)
            gains.append(float(-np.min(trapezoid(t,v,np.arange(.005,.012,.25e-6),T,T))) if kind=='injection' else frame(t,v))
        g = gains[-1]
        row = dict(output_noise_rms_V=noise,signal_V_per_ion=g,enc_e=noise/g,
            gain_refinement_relative=abs(gains[0]/g-1),
            noise_grid_relative=abs(filtered_rms(*spectra[200],H,span)/noise-1),
            noise_quadrature_relative=abs(filtered_rms(*spectra[400],H,span,order=16)/noise-1))
        assert row['gain_refinement_relative'] < 1e-3
        assert row['noise_grid_relative'] < 1e-3
        assert row['noise_quadrature_relative'] < 1e-5
        if kind=='ions':
            candidate = int(np.floor(5*noise/g))+1
            checks = []
            for n in range(max(1,candidate-2),candidate+3):
                responses = [frame(*wave(kind,n,r)) for r in (1,2)]
                check = dict(ions=n,snr=responses[-1]/noise,gain_relative_to_single_ion=responses[-1]/(n*g),
                    timestep_refinement_relative=abs(responses[0]/responses[-1]-1))
                assert check['timestep_refinement_relative'] < 1e-3
                checks.append(check)
            passing = [c['ions'] for c in checks if c['snr']>5]
            minimum = min(passing)
            assert any(c['ions']==minimum-1 and c['snr']<=5 for c in checks)
            row.update(verified_minimum_integer_N=minimum,finite_count_checks=checks)
            t,v = wave(kind,1,2)
            row['signal_retained_vs_point_samples'] = g/(-(np.interp(.01,t,v)-alpha*np.interp(0.,t,v)))
            t,v = np.r_[-.03,t],np.r_[0.,v]
            row['neighbour_packet_response_relative'] = {str(s):frame(t+s,v)/g for s in (-.01,.01)}
        report[kind] = row
        print(kind,json.dumps(row),flush=True)
    baseline = json.loads(baseline_source.read_text())
    bshot = next(r for r in baseline['short_injection'] if np.isclose(r['rise_s'],T) and np.isclose(r['flat_s'],T))
    bion = next(r for r in baseline['ion_frames'] if np.isclose(r['averaging_s'],T))
    rows = [[5.,bshot['enc_e'],bion['enc_equivalent_ions'],bion['verified_minimum_integer_N'],bion['signal_V_per_ion'],bion['output_noise_rms_V']],
        [args.pad_pf,report['injection']['enc_e'],report['ions']['enc_e'],report['ions']['verified_minimum_integer_N'],report['ions']['signal_V_per_ion'],report['ions']['output_noise_rms_V']]]
    np.savetxt(out/'comparison.csv',rows,delimiter=',',header='external_pad_pF,injection_ENC_e_RMS,ion_frame_ENC_e_RMS,minimum_N_snr5,ion_signal_V_per_ion,ion_output_noise_rms_V',comments='')
    (out/'study.json').write_text(json.dumps(report,indent=2)+'\n')
    fig,ax = plt.subplots(figsize=(5.6,3.6),layout='constrained')
    a = np.array(rows)
    for j,label in ((1,'Short injection, trapezoid G=T'),(2,'Uniform 10 ms ion frame')):
        ax.plot(a[:,0],a[:,j],'o-',label=label)
    ax.set(xlabel='External pad capacitance (pF)',ylabel='Input-referred ENC (electrons RMS)',title=f'{args.averaging_us:g} us averaging')
    ax.legend(frameon=False);ax.grid(alpha=.2)
    fig.savefig(out/'comparison.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(ROOT/'notes/physics/work/shaping/capacitance.png',dpi=150)
    summary = f'''# External pad capacitance comparison

Generated by `scripts/build_capacitance_study.py --simulator PATH --pad-pf {args.pad_pf:g} --averaging-us {args.averaging_us:g}`.

The {args.pad_pf:g} pF load is a design estimate from area scaling, 5 pF times
(0.5/1)^2 for the default run. The 5 pF benchmark is not an independently
measured 1 mm square pad capacitor. Fringing, substrate stack, guard geometry
and fixed wiring parasitics require an electrostatic/layout estimate. The input
current already uses a 0.5 mm square pad in a grounded plane, a finite 5 mm gap
and uniform mean arrivals over 10 ms. It is unchanged, with no renormalization.

Only Cd1 changes; original transistor dimensions, biases, Cf=3 fF and intrinsic
capacitances remain. Temperature 20 C; noise band 0.1 Hz–1 MHz; supplied model
has zero flicker noise. Short injection uses a trapezoid with rise and flat top
{args.averaging_us:g} us. Ion frames use centered boundary averages of that duration,
nominal tau=0.3 s discharge correction and actual circuit signal calibration.
ADC noise and sampling are not included. The 5 pF results are read from the
committed shaping study; both transient response and noise PSD are freshly
simulated for the new capacitance. Source hashes are in study.json.

| External pad (pF) | Injection ENC (e RMS) | Ion-frame ENC (e RMS) | Minimum N for S/N > 5 |
|---:|---:|---:|---:|
'''
    for r in rows:summary += f'| {r[0]:g} | {r[1]:.3f} | {r[2]:.3f} | {r[3]} |\n'
    summary += '\nAdjacent integer counts directly bracket S/N=5 in LTspice. Halved timestep, doubled spectral density, noise quadrature and injected charge checks pass. See study.json for individual checks and neighbouring-frame response.\n'
    (out/'README.md').write_text(summary)

if __name__=='__main__':main()
