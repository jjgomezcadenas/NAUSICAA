"""VH1: the longer-load circuit (M3 at 8 um, present cascode bias) with the
production feedback capacitors.

vh1_cascode_bias.py measured the L3_8u_casc0.66 variant at 3 fF. This script
takes its tuned biases (the operating point does not depend on Cf) and, at
20 fF and 30 fF, measures the open-loop gain, the charge fraction, the gain
per ion and the ENC for 20 ms frames, saves the noise spectra for the VH1
cards, and scans the output window (gain and fraction against the output
level) as in the cascode study.

Decks go to LTSpice/designs/VH1/cascode_bias/, results to
results/ltspice/design/VH1/cascode_bias/ (cf_scan_L3_8u_casc0.66.json,
noise_L3_8u_casc0.66_cf20f.csv, noise_L3_8u_casc0.66_cf30f.csv).

    PYTHONPATH=src work/event-viewer-venv/bin/python scripts/electronics/design/vh1_cascode_cf.py
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ion_signal.spice_runtime import resolve_simulator  # noqa: E402
import vh1_reset_csa as base  # noqa: E402
import vh1_improvements as imp  # noqa: E402
import vh1_cascode_bias as cb  # noqa: E402  (patches imp.core with the cascode gate and sets the folders)

VARIANT = 'L3_8u_casc0.66'
CAPS = {'20f': 20e-15, '30f': 30e-15}
LEVELS = [0.886, 0.80, 0.70, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30]
E = base.E


def measure_cf(r, p, ns, label, cf_label, cf):
    d = r(f'ac_{label}_cf{cf_label}', imp.ac_deck(p, ns, cf=cf_label))
    f = d['frequency'].real; A = np.abs(-d['V(out)'] / d['V(in)']); zin = np.abs(d['V(x)']); tf = np.abs(d['V(out)'])
    i = int(np.argmin(np.abs(f - 1.0)))
    A0 = float(A[0]); fraction = float(tf[i] * 2 * np.pi * f[i] * cf)
    c_total = 1 / (2 * np.pi * f[i] * zin[i]) - (1 + A0) * cf
    dn = r(f'noise_{label}_cf{cf_label}', imp.noise_deck(p, ns, cf=cf_label))
    fn = dn['frequency'].real; asd = dn['V(onoise)'].real
    gain = fraction * E / cf
    enc = {k: base.enc_weighted(fn, asd, w, gain) for k, w in base.WEIGHTS.items()}
    csv = cb.OUT / f'noise_{label}_cf{cf_label}.csv'
    np.savetxt(csv, np.column_stack([fn, asd]), delimiter=',', header='frequency_Hz,output_V_per_sqrtHz', comments='')
    return dict(Cf_F=cf, A0=A0, charge_fraction=fraction, gain_uV_per_ion=gain * 1e6, C_total_pF=c_total * 1e12,
                enc=enc, noise_csv=str(csv.relative_to(ROOT)))


def window(r, p, ns, label, cf_label, cf):
    rows = []
    gate = cb.p_gate(ns)
    for v in LEVELS:
        shift = v - gate
        d = r(f'win_{label}_cf{cf_label}_{v:.2f}'.replace('.', 'p'),
              imp.ac_deck(p, ns, cf=cf_label).replace(f'Vshift out r {imp.VSHIFT:.6f}', f'Vshift out r {shift:.6f}'), keep=False)
        f = d['frequency'].real; A = np.abs(-d['V(out)'] / d['V(in)']); tf = np.abs(d['V(out)'])
        i = int(np.argmin(np.abs(f - 1.0)))
        rows.append(dict(out_V=v, A0=float(A[0]), charge_fraction=float(tf[i] * 2 * np.pi * f[i] * cf)))
    return rows


def main():
    sim = resolve_simulator()
    r = imp.Runner(sim)
    prev = json.loads((cb.OUT / 'summary.json').read_text())
    p = dict(prev[VARIANT]['params'])
    op = imp.operating_point(r, p, f'op_{VARIANT}')
    ns = imp.nodeset(op)
    results = dict(variant=VARIANT, params=p, operating_point=op, at_3f=prev[VARIANT])
    for cf_label, cf in CAPS.items():
        results[cf_label] = measure_cf(r, p, ns, VARIANT, cf_label, cf)
        results[cf_label]['window'] = window(r, p, ns, VARIANT, cf_label, cf)
        v = results[cf_label]
        print(cf_label, {k: round(v[k], 4) for k in ('A0', 'charge_fraction', 'gain_uV_per_ion', 'C_total_pF')},
              {k: round(x, 1) for k, x in v['enc'].items()}, flush=True)
        print('  window', [(w['out_V'], round(w['A0']), round(w['charge_fraction'], 3)) for w in v['window']], flush=True)
    results['runs'] = r.n
    (cb.OUT / f'cf_scan_{VARIANT}.json').write_text(json.dumps(results, indent=1))
    report(results)


def report(res):
    key5 = 'frame difference of 5 ms averages'; key20 = 'frame difference of 20 ms averages'
    L = [f'# VH1 {VARIANT}: production feedback capacitors', '',
         'The longer-load circuit (M1 110 um, M2 2 um, M3 200/8 um, M4 23 um, cascode gate 0.66 V) with the biases tuned at 3 fF; '
         'open-loop gain, charge fraction and ENC for 20 ms frames at each feedback capacitor.', '',
         '| Cf | A0 | charge fraction | uV/ion | ENC 5 ms | ENC 20 ms |', '|---|---|---|---|---|---|']
    a = res['at_3f']
    L.append(f"| 3 fF | {a['A0']:.0f} | {a['charge_fraction']:.3f} | {a['gain_uV_per_ion']:.1f} | {a['enc'][key5]:.1f} | {a['enc'][key20]:.1f} |")
    for cf_label in CAPS:
        v = res[cf_label]
        L.append(f"| {cf_label[:-1]} fF | {v['A0']:.0f} | {v['charge_fraction']:.4f} | {v['gain_uV_per_ion']:.3f} | {v['enc'][key5]:.1f} | {v['enc'][key20]:.1f} |")
    L += ['', '## Output window: open-loop gain and charge fraction against the output level', '']
    for cf_label in CAPS:
        L += [f'### {cf_label[:-1]} fF', '', '| out (V) | A0 | charge fraction |', '|---|---|---|']
        L += [f"| {w['out_V']:.3f} | {w['A0']:.0f} | {w['charge_fraction']:.4f} |" for w in res[cf_label]['window']]
        L.append('')
    L += [f"Simulator runs: {res['runs']}.", '']
    (cb.OUT / f'cf_scan_{VARIANT}.md').write_text('\n'.join(L)); print('\n'.join(L))


if __name__ == '__main__':
    main()
