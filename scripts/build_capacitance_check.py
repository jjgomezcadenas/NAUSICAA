#!/usr/bin/env python3
"""Historical diagnostic: compare the former 250 fF example with the published Topmetal-S 5 pF value.

Run `make capacitance-check`. Only Cin changes. Noise spectra, feedback and
amplifier parameters remain the illustrative assumptions in realistic_csa.toml.
The literature reports electrode capacitance to ground; this check substitutes
5 pF for total Cin (excluding Cf), without adding an unknown transistor/layout
capacitance. It is a benchmark, not an extracted NAUSICA capacitance.
"""
from pathlib import Path
from dataclasses import asdict, replace
import hashlib
import json
import platform
import tomllib
import numpy as np
import scipy
from ion_signal.physics import Detector, E_CHARGE
from ion_signal.realistic_csa import (CSA, Noise, noise_metrics, ions_waveform,
                                    response, sampled_voltage, recovered_charge, frame_charge)
from ion_signal.noise_study import direct_charge_variance
from build_realistic_csa import configured_arrivals

ROOT = Path(__file__).resolve().parents[1]


def main():
    circuit_path = ROOT/'config/realistic_csa.toml'
    detector_path = ROOT/'config/detector.toml'
    cfg = tomllib.loads(circuit_path.read_text())
    detector = Detector.from_toml(detector_path)
    baseline = CSA(**cfg['electronics'])
    noise = Noise(**cfg['noise'])
    dt = cfg['readout']['sample_interval_s']
    stop = cfg['readout']['target_interval_end_s']
    samples = np.array([stop-dt, stop])
    rho = cfg['readout']['snr_threshold']
    num = cfg['numerics']
    records = []
    cases = [c for c in cfg['cases'] if c['id'] in ('simultaneous_1000', 'spread_1000')]
    if len(cases) != 2:
        raise ValueError('Both paper arrival cases required')
    for label, capacitance in [('historical_250fF_example', 250e-15), ('topmetal_s', 5e-12)]:
        model = replace(baseline, input_capacitance_F=capacitance)
        total = noise_metrics(model, noise, dt, aperture_s=0., relative_tolerance=1e-9)
        direct = direct_charge_variance(model, noise, dt, order=16)
        direct_error = abs(direct/total['recovered_charge_rms_C']**2-1)
        if direct_error > 1e-6:
            raise ArithmeticError('Independent noise integrals disagree')
        blank = replace(noise, voltage_white_V_per_sqrtHz=0., voltage_flicker_V2=0.,
                        current_white_A_per_sqrtHz=0., passive_feedback_thermal=False)
        sources = dict(white_voltage=replace(blank, voltage_white_V_per_sqrtHz=noise.voltage_white_V_per_sqrtHz),
                       flicker_voltage=replace(blank, voltage_flicker_V2=noise.voltage_flicker_V2),
                       feedback=replace(blank, passive_feedback_thermal=noise.passive_feedback_thermal),
                       additional_current=replace(blank, current_white_A_per_sqrtHz=noise.current_white_A_per_sqrtHz))
        components = {name: noise_metrics(model, source, dt, 0., 1e-9) for name, source in sources.items()}
        sum_variance = sum(c['recovered_charge_rms_C']**2 for c in components.values())
        addition_error = abs(sum_variance/total['recovered_charge_rms_C']**2-1)
        if addition_error > 1e-6:
            raise ArithmeticError('Noise variances do not add')
        signals = []
        for case in cases:
            arrivals = configured_arrivals(case)
            estimates = []
            for divisor in (1, 2):
                wave = ions_waveform(detector, arrivals, num['start_s'], num['end_s'], num['max_step_s']/divisor)
                volts = response(wave, model)[:, 1]
                readings = sampled_voltage(wave.edges_s, volts, samples, aperture_s=0.)
                qhat = float(recovered_charge(readings, model, dt)[0])
                estimates.append(qhat)
            refinement = abs(estimates[0]/estimates[1]-1)
            if refinement > 1e-4:
                raise ArithmeticError('Sampled signal has not converged')
            truth = float(frame_charge(wave, samples)[0])
            snr = abs(qhat)/total['recovered_charge_rms_C']
            threshold = rho*len(arrivals)/snr
            signals.append(dict(id=case['id'], arrival_configuration=case, ions=len(arrivals),
                                input_charge_e=truth/E_CHARGE, recovered_charge_e=qhat/E_CHARGE,
                                recovered_fraction=qhat/truth, sample_voltage_V=readings.tolist(),
                                signal_to_noise=snr, target_signal_to_noise=rho,
                                template_ion_count_at_target=threshold,
                                minimum_integer_template_ion_count=int(np.ceil(threshold)),
                                timestep_refinement_relative_difference=refinement))
        records.append(dict(id=label, electronics=asdict(model), noise=total, components=components,
                            charge_gain_reduction=model.charge_gain_reduction(), cases=signals,
                            checks=dict(direct_noise_relative_difference=direct_error,
                                        variance_addition_relative_difference=addition_error)))
    payload = dict(literature_source=dict(url='https://arxiv.org/pdf/2010.09226', section='2.3',
                                         electrode_diameter_m=.001, electrode_capacitance_F=5e-12),
                   assumptions=[__doc__, 'Central square pad in a grounded plane; constant mobility; finite-path charge not renormalized.',
                                'Point voltage samples; stationary independent noise over configured band.',
                                'Thresholds scale each full arrival template; no random sparse-arrival simulation.'],
                   noise_inputs=asdict(noise), sample_interval_s=dt, sampling_aperture_s=0.,
                   sample_times_s=samples.tolist(), numerics=num,
                   environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                   source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in (circuit_path, detector_path)}, results=records)
    out = ROOT/'results/realistic_csa'
    out.mkdir(parents=True, exist_ok=True)
    (out/'capacitance_check.json').write_text(json.dumps(payload, indent=2)+'\n')
    selected = records[1]
    components = selected['components']
    snr_rounded = {f"{5*round(c['signal_to_noise']/5):.0f}" for c in selected['cases']}
    recovery_rounded = {f"{100*c['recovered_fraction']:.0f}" for c in selected['cases']}
    if len(snr_rounded) != 1 or len(recovery_rounded) != 1:
        raise ValueError('Two cases no longer share the rounded values quoted in the paper')
    voltage_variance = sum(components[name]['recovered_charge_rms_C']**2
                           for name in ('white_voltage', 'flicker_voltage'))
    values = {
        'TopmetalCapacitancePf': f"{selected['electronics']['input_capacitance_F']*1e12:g}",
        'TopmetalElectrodeDiameterMm': f"{payload['literature_source']['electrode_diameter_m']*1e3:g}",
        'TopmetalNoiseTotalElectrons': f"{selected['noise']['induced_charge_ENC_e']:.2f}",
        'TopmetalNoiseWhiteElectrons': f"{components['white_voltage']['induced_charge_ENC_e']:.2f}",
        'TopmetalNoiseFlickerElectrons': f"{components['flicker_voltage']['induced_charge_ENC_e']:.2f}",
        'TopmetalNoiseFeedbackElectrons': f"{components['feedback']['induced_charge_ENC_e']:.2f}",
        'TopmetalVoltageNoisePercent': f"{100*voltage_variance/selected['noise']['recovered_charge_rms_C']**2:.0f}",
        'TopmetalRecoveredPercent': recovery_rounded.pop(),
        'TopmetalApproxSnr': snr_rounded.pop(),
        'TopmetalExampleIons': str(selected['cases'][0]['ions']),
    }
    generated = ROOT/'notes/physics/generated'
    generated.mkdir(parents=True, exist_ok=True)
    (generated/'capacitance_check_values.tex').write_text(
        '% Generated by scripts/build_capacitance_check.py; 5 pF benchmark.\n'+
        ''.join(fr'\newcommand{{\{key}}}{{{value}}}'+'\n' for key, value in values.items()))
    lines = ['# Topmetal-S capacitance check', '', 'Reproduce: `make capacitance-check`.', '',
             'Use the published 5 pF electrode capacitance as Cin, excluding Cf, with no additional transistor capacitance. All other circuit and noise inputs remain fixed.',
             'Source: https://arxiv.org/pdf/2010.09226, Section 2.3 (1 mm electrode). This is a benchmark substitution, not a geometrical extraction for NAUSICA.', '',
             f'Point readings separated by {dt:g} s; assumed noise band {noise.low_frequency_Hz:g}–{noise.high_frequency_Hz:g} Hz; passive-resistor feedback noise.', '',
             '| Cin (pF) | White voltage noise (e RMS) | Flicker (e RMS) | Feedback (e RMS) | Total (e RMS) |',
             '| ---: | ---: | ---: | ---: | ---: |']
    for r in records:
        c = r['components']
        lines.append(f"| {r['electronics']['input_capacitance_F']*1e12:g} | {c['white_voltage']['induced_charge_ENC_e']:.3f} | {c['flicker_voltage']['induced_charge_ENC_e']:.3f} | {c['feedback']['induced_charge_ENC_e']:.3f} | {r['noise']['induced_charge_ENC_e']:.3f} |")
    lines += ['', '| Cin (pF) | Arrivals | Recovered charge / input (%) | S/N for 1000 ions | Ions for S/N = 5 (rounded up) |',
              '| ---: | --- | ---: | ---: | ---: |']
    for r in records:
        for c in r['cases']:
            lines.append(f"| {r['electronics']['input_capacitance_F']*1e12:g} | {c['id']} | {100*c['recovered_fraction']:.2f} | {c['signal_to_noise']:.2f} | {c['minimum_integer_template_ion_count']} |")
    lines += ['', 'Checks: halve the current time step; compare covariance and independent direct noise integration; verify addition of independent source variances. Detailed values and input hashes are in capacitance_check.json.', '']
    (out/'capacitance_check.md').write_text('\n'.join(lines))
    print('\n'.join(lines))
    print('Maximum signal refinement:', max(c['timestep_refinement_relative_difference'] for r in records for c in r['cases']))
    print('Maximum noise cross-check:', max(r['checks']['direct_noise_relative_difference'] for r in records))


if __name__ == '__main__':
    main()
