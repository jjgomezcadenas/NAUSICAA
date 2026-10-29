#!/usr/bin/env python3
"""Reproduce paper noise estimates with point samples at frame boundaries.

Run `make noise-study`. Configuration and previously generated signal records
are hashed in the JSON output. No manuscript edits or LaTeX compilation.
Independent sources add in variance. Ion thresholds use the mean signal per ion for each arrival
pattern; they are S/N criteria, not detection probabilities.
"""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
from dataclasses import asdict, replace
import csv
import hashlib
import json
import platform
import tomllib
import numpy as np
import scipy
from ion_signal.physics import E_CHARGE
from ion_signal.realistic_csa import CSA, Noise, noise_metrics
from ion_signal.noise_study import direct_charge_variance

ROOT = Path(__file__).resolve().parents[2]


def main():
    config_path = ROOT/'config/realistic_csa.toml'
    signal_path = ROOT/'results/realistic_csa/arrival_cases_point_samples.json'
    cfg = tomllib.loads(config_path.read_text())
    signal = json.loads(signal_path.read_text())
    model, noise = CSA(**cfg['electronics']), Noise(**cfg['noise'])
    dt = cfg['readout']['sample_interval_s']
    if (signal['electronics'] != asdict(model) or
            signal['sample_interval_s'] != dt or signal['sampling_aperture_s'] != 0.):
        raise ValueError('Regenerate point-sample signals with the current circuit')
    blank = replace(noise, voltage_white_V_per_sqrtHz=0., voltage_flicker_V2=0.,
                    current_white_A_per_sqrtHz=0., passive_feedback_thermal=False)
    sources = {
        'white_voltage': replace(blank, voltage_white_V_per_sqrtHz=noise.voltage_white_V_per_sqrtHz),
        'flicker_voltage': replace(blank, voltage_flicker_V2=noise.voltage_flicker_V2),
        'feedback': replace(blank, passive_feedback_thermal=noise.passive_feedback_thermal),
        'additional_current': replace(blank, current_white_A_per_sqrtHz=noise.current_white_A_per_sqrtHz),
        'total': noise,
    }
    results = {}
    for name, source in sources.items():
        coarse = noise_metrics(model, source, dt, aperture_s=0., relative_tolerance=1e-7)
        fine = noise_metrics(model, source, dt, aperture_s=0., relative_tolerance=1e-9)
        direct8 = direct_charge_variance(model, source, dt, order=8)
        direct16 = direct_charge_variance(model, source, dt, order=16)
        variance = fine['recovered_charge_rms_C']**2
        scale = variance if variance else 1.
        checks = dict(covariance_tolerance_relative_change=abs(coarse['recovered_charge_rms_C']**2-variance)/scale,
                      direct_quadrature_relative_change=abs(direct16-direct8)/scale,
                      direct_vs_covariance_relative_difference=abs(direct16-variance)/scale)
        if max(checks.values()) > 1e-6:
            raise ArithmeticError(f'{name}: unconverged integral: {checks}')
        results[name] = dict(**fine, variance_C2=variance,
                             direct_variance_C2=direct16, checks=checks)
    total = results['total']
    sum_variance = sum(row['variance_C2'] for name, row in results.items() if name != 'total')
    addition_error = abs(sum_variance/total['variance_C2']-1)
    if addition_error > 1e-6:
        raise ArithmeticError('Independent source variances do not add')
    cases = []
    rho = cfg['readout']['snr_threshold']
    for case in signal['cases']:
        g = case['finite_gain_estimated_charge_per_ion_C']
        boundary = rho*total['recovered_charge_rms_C']/abs(g)
        # Strict S/N > rho: floor + 1 also handles an exactly integral boundary.
        minimum = int(np.floor(boundary))+1
        cases.append(dict(id=case['id'], mean_estimated_charge_per_ion_C=g,
                          mean_estimated_charge_per_ion_e=g/E_CHARGE,
                          signal_to_noise_per_ion=abs(g)/total['recovered_charge_rms_C'],
                          target_signal_to_noise=rho, ion_count_boundary=boundary,
                          criterion='S/N > target', minimum_integer_ion_count=minimum,
                          snr_at_minimum=minimum*abs(g)/total['recovered_charge_rms_C'],
                          snr_at_one_fewer=(minimum-1)*abs(g)/total['recovered_charge_rms_C']))
        if not cases[-1]['snr_at_one_fewer'] <= rho < cases[-1]['snr_at_minimum']:
            raise ArithmeticError('Integer threshold does not bracket requested S/N')
    out = ROOT/'results/realistic_csa'
    payload = dict(electronics=asdict(model), noise=asdict(noise), sample_interval_s=dt,
                   sampling_aperture_s=0., source_configuration_sampling_aperture_s=cfg['readout']['sampling_aperture_s'],
                   provenance=cfg['provenance'], signal_assumptions=signal['assumptions'],
                   assumptions=['Stationary independent equivalent input noise sources; one-sided PSDs per Hz.',
                                'Point voltage readings at frame boundaries; no acquisition averaging.',
                                'Feedback noise is the passive-resistor benchmark, not a MOS noise prediction.',
                                'S/N threshold scales the specified arrival template, including its timing.'],
                   source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in (config_path, signal_path)},
                   environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                   components=results, variance_addition_relative_error=addition_error, cases=cases)
    (out/'noise_point_samples.json').write_text(json.dumps(payload, indent=2)+'\n')
    with (out/'noise_point_samples.csv').open('w', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(['source', 'sample_voltage_rms_V', 'sample_correlation', 'charge_rms_C', 'charge_rms_e', 'variance_fraction'])
        for name, row in results.items():
            writer.writerow([name, row['sample_voltage_rms_V'], row['adjacent_sample_correlation'],
                             row['recovered_charge_rms_C'], row['induced_charge_ENC_e'], row['variance_C2']/total['variance_C2']])
    # Values used by the approved numerical example; preserve one source of truth.
    if noise.current_white_A_per_sqrtHz != 0 or not noise.passive_feedback_thermal:
        raise ValueError('Paper example assumes no additional current noise and passive feedback noise')
    generated = ROOT/'notes/physics/generated'
    generated.mkdir(parents=True, exist_ok=True)
    uniform = next(c for c in cases if c['id']=='uniform_per_ion')
    simultaneous = next(c for c in cases if c['id']=='simultaneous_per_ion')
    values = {
        'NoiseWhiteNv': f'{noise.voltage_white_V_per_sqrtHz*1e9:g}',
        'NoiseFlickerVTwo': f'{noise.voltage_flicker_V2:g}',
        'NoiseTemperatureK': f'{noise.temperature_K:g}',
        'NoiseLowHz': f'{noise.low_frequency_Hz:g}',
        'NoiseHighMHz': f'{noise.high_frequency_Hz/1e6:g}',
        'NoiseTotalElectrons': f"{total['induced_charge_ENC_e']:.2f}",
        'NoiseWhiteElectrons': f"{results['white_voltage']['induced_charge_ENC_e']:.2f}",
        'NoiseFlickerElectrons': f"{results['flicker_voltage']['induced_charge_ENC_e']:.2f}",
        'NoiseFeedbackElectrons': f"{results['feedback']['induced_charge_ENC_e']:.2f}",
        'NoiseFeedbackPercent': f"{100*results['feedback']['variance_C2']/total['variance_C2']:.0f}",
        'NoiseVoltagePercent': f"{100*(results['white_voltage']['variance_C2']+results['flicker_voltage']['variance_C2'])/total['variance_C2']:.0f}",
        'NoiseMeanChargePerIon': f"{uniform['mean_estimated_charge_per_ion_e']:.4f}",
        'NoiseSnrPerIon': f"{uniform['signal_to_noise_per_ion']:.4f}",
        'NoiseUniformBoundary': f"{uniform['ion_count_boundary']:.2f}",
        'NoiseUniformMinimumIons': str(uniform['minimum_integer_ion_count']),
        'NoiseSimultaneousMinimumIons': str(simultaneous['minimum_integer_ion_count']),
        'NoiseTargetSnr': f'{rho:g}',
    }
    (generated/'noise_values.tex').write_text(
        '% Generated by scripts/benchmarks/build_noise_study.py; point samples.\n'+
        ''.join(fr'\newcommand{{\{key}}}{{{value}}}'+'\n' for key, value in values.items()))
    (generated/'noise_inputs_table.tex').write_text('\n'.join([
        '% Generated by scripts/benchmarks/build_noise_study.py from config/realistic_csa.toml.',
        r'\begin{tabular}{lll}', r'\toprule',
        r'Parameter & Symbol & Value \\', r'\midrule',
        r'White input voltage noise & $e_w$ & $\NoiseWhiteNv\,\mathrm{nV}/\sqrt{\mathrm{Hz}}$ \\',
        r'Flicker coefficient & $a_f$ & $\num{\NoiseFlickerVTwo}\,\mathrm{V^2}$ \\',
        r'Additional white input current noise & $i_w$ & $'+f'{noise.current_white_A_per_sqrtHz:g}'+r'\,\mathrm{A}/\sqrt{\mathrm{Hz}}$ \\',
        r'Temperature & $T$ & \SI{\NoiseTemperatureK}{K} \\',
        r'Lower integration limit & $f_{\min}$ & \SI{\NoiseLowHz}{Hz} \\',
        r'Upper integration limit & $f_{\max}$ & \SI{\NoiseHighMHz}{MHz} \\',
        r'\bottomrule', r'\end{tabular}', '']))
    (generated/'noise_sensitivity_table.tex').write_text('\n'.join([
        '% Generated by scripts/benchmarks/build_noise_study.py; uniform mean response.',
        r'\begin{tabular}{lrr}', r'\toprule',
        r'Arrivals & S/N per ion & Minimum $N$ for S/N $>\NoiseTargetSnr$ \\',
        r'\midrule',
        r'Uniform & $\NoiseSnrPerIon$ & $\NoiseUniformMinimumIons$ \\',
        r'\bottomrule', r'\end{tabular}', '']))
    lines = ['# Point-sample noise calculation', '',
             'Reproduce from the repository root: `make noise-study`.', '',
             'Inputs: `config/realistic_csa.toml`; signal: `arrival_cases_point_samples.json`.',
             'The JSON records all input values, assumptions, hashes and numerical checks.', '',
             f'Frame interval: {dt:g} s. Sampling aperture: 0 s. Noise band: {noise.low_frequency_Hz:g}–{noise.high_frequency_Hz:g} Hz.',
             'The sources are study assumptions; feedback uses the passive-resistor thermal spectrum.', '',
             '| Source | Charge noise (e RMS) | Fraction of variance |', '| --- | ---: | ---: |']
    for name, row in results.items():
        lines.append(f"| {name} | {row['induced_charge_ENC_e']:.6f} | {row['variance_C2']/total['variance_C2']:.6f} |")
    lines += ['', '| Arrivals | Mean estimated charge per ion (e) | S/N per ion | N boundary for S/N = 5 | Minimum integer N for S/N > 5 |',
              '| --- | ---: | ---: | ---: | ---: |']
    for case in cases:
        lines.append(f"| {case['id']} | {case['mean_estimated_charge_per_ion_e']:.6f} | {case['signal_to_noise_per_ion']:.6f} | {case['ion_count_boundary']:.3f} | {case['minimum_integer_ion_count']} |")
    lines += ['', 'Checks: covariance quadrature tolerance 1e-7 versus 1e-9; independent positive spectral integral with Gauss–Legendre orders 8 and 16; variance addition. All relative discrepancies must be below 1e-6.',
              'The ion thresholds scale each fixed arrival pattern. A different arrival history requires its own signal calculation.', '']
    (out/'noise_point_samples.md').write_text('\n'.join(lines))
    print('\n'.join(lines))
    print('Maximum relative numerical discrepancy:', max(max(r['checks'].values()) for r in results.values()))


if __name__ == '__main__':
    main()
