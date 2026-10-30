#!/usr/bin/env python3
"""Freeze the nominal 20 s CSA, its runnable LTspice inputs and saved results."""
from pathlib import Path
import csv
import hashlib
import json
import re
import subprocess
import sys
import tomllib

import h5py

ROOT = Path(__file__).resolve().parents[3]
CARD = ROOT / 'config/csa_design_v3.toml'
DEST = ROOT / 'LTSpice/v3'
sys.path.insert(0, str(ROOT / 'src'))
from ion_signal.spice import schematic_netlist


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plain(value):
    if isinstance(value, bytes):
        return value.decode()
    return value.item() if hasattr(value, 'item') else value


def number(value, parameters):
    if value.startswith('{'):
        return parameters[value[1:-1].lower()]
    factors = {'u': 1e-6, 'p': 1e-12, 'f': 1e-15}
    return float(value[:-1]) * factors[value[-1].lower()] if value[-1].lower() in factors else float(value)


def main():
    c = tomllib.loads(CARD.read_text())
    source = ROOT / c['study_directory']
    selected = c['selected_variant']
    decks = source / 'decks' / selected
    DEST.mkdir(parents=True, exist_ok=True)
    records = []
    copied_inputs = set()

    def save(relative, contents, origin=None):
        target = DEST / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = contents.encode() if isinstance(contents, str) else contents
        target.write_bytes(payload)
        records.append({'file': str(relative), 'sha256': sha(target), 'bytes': len(payload),
                        **({'source': str(origin.relative_to(ROOT)), 'source_sha256': sha(origin)} if origin else {})})

    def copy(origin, relative):
        save(relative, origin.read_bytes(), origin)

    # Localize the archived decks' input paths; preserve every electrical value.
    for deck in sorted(decks.glob('*.cir')):
        body = deck.read_text().replace('../../inputs/', '../inputs/')
        save(Path('netlists') / deck.name, body, deck)
        for dep in re.findall(r'PWL FILE="([^"]+)"', deck.read_text()):
            target = Path('inputs') / Path(dep).name
            if target not in copied_inputs:
                copy((deck.parent / dep).resolve(), target)
                copied_inputs.add(target)
        log = deck.with_suffix('.run.txt')
        if log.exists():
            copy(log, Path('netlists') / log.name)
    copy(decks / 'corrected.inc', Path('netlists/corrected.inc'))

    # Reuse Rogelio's editable drawing, with electrical equivalence checked.
    drawing = ROOT / 'LTSpice/Spice2/v3/CSAfolded180nm-2_Astronio.asc'
    template = '\n'.join(line for line in drawing.read_text(encoding='latin1').splitlines()
                         if not line.startswith(('TEXT ', 'DATAFLAG '))) + '\n'
    cases = {'NAUSICA_20ms_120ions': '20ms_single_N120_r2',
             'NAUSICA_40ms_400ions': '40ms_single_N400_r2',
             'NAUSICA_20ms_calibration': '20ms_single_N100_r2',
             'NAUSICA_40ms_calibration': '40ms_single_N100_r2',
             'NAUSICA_20ms_baseline': '20ms_single_N0_r2',
             'NAUSICA_40ms_baseline': '40ms_single_N0_r2',
             'NAUSICA_20ms_discharge': '20ms_tail_single_N100_r2',
             'NAUSICA_noise': 'noise_400'}
    copy(decks / 'corrected.inc', Path('schematic/180nm_professional.inc'))
    for name, reference in cases.items():
        deck = DEST / 'netlists' / (reference + '.cir')
        lines = deck.read_text().splitlines()
        elements = {x.split()[0]: x for x in lines if x and not x.startswith(('*', '.'))}

        def replace_symbol(match):
            block = match.group()
            instance = re.search(r'^SYMATTR InstName (\S+)', block, re.M)[1]
            if instance in ('Vion1', 'B1'):
                return ''
            if instance in ('Is1', 'Vnload', 'Vp1', 'Vf1', 'Cd1'):
                block = re.sub(r'^SYMATTR Value .*$', 'SYMATTR Value ' + elements[instance].split(maxsplit=3)[3], block, flags=re.M)
                block = re.sub(r'^SYMATTR Value2 .*\n?', '', block, flags=re.M)
            return block

        text = re.sub(r'^SYMBOL .*?(?=^SYMBOL |\Z)', replace_symbol, template, flags=re.M | re.S)
        text = text.replace('SHEET 1 1220 900\n', 'SHEET 1 1220 900\nWIRE -688 48 -688 -32\n')
        directives = [x.replace('.lib corrected.inc', '.lib 180nm_professional.inc')
                      for x in lines if x.startswith('.') and x != '.end']
        text += ''.join(f'TEXT -888 {-456 + 24*i} Left 2 !{line}\n' for i, line in enumerate(directives))
        relative = Path('schematic') / (name + '.asc')
        save(relative, text, drawing)
        actual, commands = schematic_netlist(DEST / relative)
        normalized = lambda xs: {x.split()[0]: ' '.join(x.lower().split()) for x in xs}
        assert normalized(actual) == normalized(elements.values()), name
        assert commands == directives, name

    text = (decks / '20ms_single_N100_r2.cir').read_text()
    parameters = {key.lower(): number(val, {}) for key, val in
                  re.findall(r'^\.(?:param|parameter)\s+(\w+)=(\S+)', text, re.M | re.I)}
    devices = {}
    biases = {}
    capacitors = {}
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0].startswith('M'):
            attrs = dict(x.lower().split('=', 1) for x in parts[6:])
            devices[parts[0]] = {'model': parts[5], 'terminals_D_G_S_B': parts[1:5],
                                'length_m': number(attrs['l'], parameters),
                                'width_m': number(attrs['w'], parameters),
                                'multiplicity': number(attrs['m'], parameters)}
        elif parts[0].startswith('V'):
            biases[parts[0]] = {'positive_node': parts[1], 'negative_node': parts[2], 'voltage_V': float(parts[3])}
        elif parts[0].startswith('C'):
            capacitors[parts[0]] = number(parts[3], parameters)

    result = {'version': c['version'], 'nominal_discharge_time_s': c['nominal_discharge_time_s'],
              'selected_variant': selected, 'averaging_s': c['averaging_s'],
              'external_input_capacitance_F': capacitors['Cd1'],
              'feedback_capacitance_F': capacitors['Cf'],
              'temperature_K': float(re.search(r'^\.temp (\S+)', text, re.M)[1]) + 273.15,
              'transistors': devices, 'voltage_sources': biases, 'frames': {}}
    assert devices['Mf']['length_m'] == c['feedback']['length_m']
    assert devices['Mf']['width_m'] == c['feedback']['width_m']
    assert biases['Vf1']['voltage_V'] == c['feedback']['gate_V']
    assert biases['Vp1']['voltage_V'] == c['bias']['m4_gate_V']
    with h5py.File(source / 'analysis.h5') as h:
        g = h[selected]
        result['circuit'] = {k: plain(v) for k, v in g.attrs.items()}
        for key in ('20ms', '40ms'):
            f = g[key]
            result['frames'][key] = {k: plain(v) for k, v in f.attrs.items()}
            result['frames'][key].update({k: float(f[k][0]) for k in ('gain_V_per_ion', 'noise_rms_V', 'ENC_e_rms')})
        assert result['circuit']['model_sha256'] == sha(decks / 'corrected.inc')
    with h5py.File(source / 'variants' / (selected + '.h5')) as h:
        result['operating_point'] = {key: float(h['op'][key][0]) for key in h['op']}
        # Every included netlist matches the saved acquisition after path localization.
        for deck in decks.glob('*.cir'):
            key = deck.stem
            if key.startswith(('20ms_', '40ms_')):
                frame, rest = key.split('_', 1)
                key = frame + '/' + rest.replace('tail_', 'tail/', 1)
            assert h[key].attrs['deck_sha256'] == sha(deck), deck
    with h5py.File(source / 'inputs.h5') as h:
        result['physical_input'] = {k: plain(v) for k, v in h.attrs.items()}
    save(Path('parameters.json'), json.dumps(result, indent=2) + '\n')
    save(Path('operating_point.csv'), 'quantity,value\n' + ''.join(f'{k},{v:.17g}\n' for k, v in result['operating_point'].items()))
    for name in ('calibration.csv', 'noise_budget.csv', 'compression.csv'):
        origin = ROOT / 'results/ltspice/retention_comparison' / name
        with origin.open() as stream:
            rows = list(csv.reader(stream))
        rows = [rows[0]] + [r for r in rows[1:] if float(r[0]) == 1e6 * devices['Mf']['length_m']]
        target = DEST / name
        with target.open('w', newline='') as stream:
            csv.writer(stream, lineterminator='\n').writerows(rows)
        records.append({'file': name, 'sha256': sha(target), 'source': str(origin.relative_to(ROOT)), 'source_sha256': sha(origin), 'bytes': target.stat().st_size})
    for path in (CARD, ROOT / 'config/retention_scan.toml', ROOT / 'config/detector.toml',
                 source / 'cards/length_1200um_readout.toml', source / 'cards/length_1200um_viewer.toml'):
        copy(path, Path('reference') / path.name)
    for name in ('report.tex', 'retention_comparison.pdf', 'retention_comparison.png', 'sources.txt'):
        copy(ROOT / 'latex_generated/retention_comparison' / name, Path('report') / name)
    copy(ROOT / 'notes/physics/figures/Nausica_Image.jpg', Path('Nausica_Image.jpg'))
    for deck in [*DEST.rglob('*.cir'), *DEST.rglob('*.asc')]:
        for dep in re.findall(r'(?:\.lib\s+|PWL FILE=)"?([^\s"\n]+)', deck.read_text(), re.I):
            assert (deck.parent / dep).is_file(), (deck, dep)
    f20, f40 = result['frames']['20ms'], result['frames']['40ms']
    readme = f'''# NAUSICA CSA v3 — nominal 20 s discharge

Frozen provisional design for Vicente Herrero and Rogelio Palomo.
The selected circuit is `{selected}`. Its fitted small-signal discharge time
is {f20['tau_s']:.2f} s; the design label is 20 s. The saved alpha factors use
the fitted value. The discharge comes from Mf in the transistor circuit.

## Open and run

Open an `.asc` in `schematic/` in LTspice and press Run. Keep the complete
v3 folder together: currents are in `inputs/`, with time in seconds and current
in amperes. Plot V(out), V(in) and I(Is1).

- `NAUSICA_20ms_120ions.asc`: illustrative 120-ion packet over 20 ms.
- `NAUSICA_40ms_400ions.asc`: illustrative 400-ion packet over 40 ms.
- `NAUSICA_20ms_calibration.asc` and `NAUSICA_40ms_calibration.asc`: 100 ions.
- `NAUSICA_20ms_baseline.asc` and `NAUSICA_40ms_baseline.asc`: zero current.
- `NAUSICA_20ms_discharge.asc`: 100 ions and a six-second discharge tail.
- `NAUSICA_noise.asc`: output noise from 0.0001 Hz to 1 MHz.

The exact model is `schematic/180nm_professional.inc`, also provided as
`netlists/corrected.inc`. All eight schematics are electrically checked against
the archived netlists. `netlists/` also includes both time resolutions, the
10/100-ion tail checks, all tested packet sizes and readable simulator logs.
The circuit runs used native macOS LTspice 17.2.4.

## Circuit and readout

Supply 1.8 V, temperature 20 C, external input capacitance 1.25 pF plus
intrinsic model capacitances, feedback capacitance 3 fF. Mf: W=0.36 um,
L=1200 um, gate=0 V. M4 gate=0.6621 V. `parameters.json` lists the actual
W, L and multiplicity of every transistor, every bias and the operating point.
The schematic and netlists preserve the exact values used in the saved runs.

| Quantity | 20 ms frame | 40 ms frame |
|---|---:|---:|
| Fitted discharge time (s) | {f20['tau_s']:.4f} | {f40['tau_s']:.4f} |
| alpha | {f20['alpha']:.9f} | {f40['alpha']:.9f} |
| Gain (uV / collected ion) | {1e6*f20['gain_V_per_ion']:.3f} | {1e6*f40['gain_V_per_ion']:.3f} |
| Filtered output noise (mV RMS) | {1e3*f20['noise_rms_V']:.3f} | {1e3*f40['noise_rms_V']:.3f} |
| ENC (electrons RMS) | {f20['ENC_e_rms']:.2f} | {f40['ENC_e_rms']:.2f} |

Subtract the zero-input baseline, average the output in centred 5 ms windows
at successive frame boundaries, and calculate
`D_n = -(Vbar_n - alpha*Vbar_(n-1))`, with `alpha=exp(-frame/tau)`.
The recovered ion count is `D_n/g`, where g is the calibrated gain above.
The same averaging and frame-difference filter is applied to the one-sided
SPICE output-noise spectrum. These operations are performed offline in Julia.
The external capacitance and averaging width are common to both frame choices.

The physical input is uniform mean collection over a frame, with an on-axis
constant-speed 3 m/s trajectory through the final 5 mm above a 0.5 mm square
pad. The pad weighting potential is one; the surrounding infinite grounded
plane and infinity are zero, with homogeneous dielectric. The integrated
current per ion is e[1-phi_w(g)]. Collection begins at simulation time 10 ms;
induction begins before collection. The waveform keeps the finite-path charge.

## Saved results and checks

`parameters.json`: complete circuit, readout, physical input, calibration and
convergence results, with source hashes. `operating_point.csv`: DC voltages and
currents. `calibration.csv`, `noise_budget.csv` and `compression.csv`: selected
results from the saved Julia analysis. The latter uses 8000 double-beta and
8000 electron events for each geometry, including empty buckets. Its event
fractions count at least one locally compressed bucket.

`report/report.tex` and its figure describe the complete retention scan and
validation. Numerical checks cover charge conservation, transient resolution,
tail fits, noise integration, gain linearity, and nonlinear-history refinement.
`reference/` saves the design and study cards. `manifest.json` records file
checksums and their repository sources. The HDF5 numerical data remain in
the repository's local results directories.

## Reproduction in the repository

    python scripts/electronics/feedback/run_retention_comparison.py --simulator /path/to/LTspice
    python notes/physics/tools/package_csa_v3.py
'''
    save(Path('README.md'), readme)
    save(Path('.gitignore'), '/inputs/\n*.raw\n*.log\n*.op.raw\n*.h5\n')
    manifest = {'version': c['version'], 'nominal_discharge_time_s': c['nominal_discharge_time_s'],
                'selected_variant': selected, 'generator': str(Path(__file__).relative_to(ROOT)),
                'generator_sha256': sha(Path(__file__)),
                'repository_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'analysis_sha256': sha(source / 'analysis.h5'), 'files': records,
                'validation': {'schematics_electrically_matched': len(cases), 'all_dependencies_local': True}}
    (DEST / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'{DEST}: {len(records)} files; {len(cases)} schematics checked; all dependencies local.')


if __name__ == '__main__':
    main()
