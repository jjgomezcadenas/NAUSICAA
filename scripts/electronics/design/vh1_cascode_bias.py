"""VH1: cascode bias and load length on the improved circuit.

Section 13.5 of notes/physics/nausicaa_VH.tex argues that the cascode gate
bias limits the length of the load M3. This script closes the point on the
combined circuit with the 23 um cascode (combined_W4_23u of
vh1_improvements.py):

1. M3 at 4, 8 and 16 um with the cascode gate at 0.66 V (present) and
   raised by 0.1 and 0.2 V; every variant re-tuned with the loop closed as in
   vh1_improvements.py; open-loop gain, node med, input capacitance, charge
   fraction and ENC for 20 ms frames.
2. The lower edge of the output window at each cascode bias: the level shift
   of the closed reset loop is lowered step by step, which moves the output
   down at a fixed gate voltage as accumulated charge does, and the small-
   signal gain and charge fraction are read at each level.

Decks go to LTSpice/designs/VH1/cascode_bias/, results to
results/ltspice/design/VH1/cascode_bias/ (summary.json, summary.md, noise
spectra as CSV).

    PYTHONPATH=src work/event-viewer-venv/bin/python scripts/electronics/design/vh1_cascode_bias.py
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ion_signal.spice_runtime import resolve_simulator  # noqa: E402
import vh1_reset_csa as base  # noqa: E402
import vh1_improvements as imp  # noqa: E402

DESIGN = ROOT / 'LTSpice/designs/VH1/cascode_bias'
OUT = ROOT / 'results/ltspice/design/VH1/cascode_bias'
imp.DESIGN = DESIGN; imp.OUT = OUT; imp.RUN = OUT / 'run'
VNCASC0 = 0.6621


def core(p, voff=0.0):
    """The VH1 core of vh1_improvements.py with the cascode gate as a parameter."""
    text = imp_core_original(p, voff)
    return text.replace(f'Vp1 Vncasc 0 {VNCASC0}', f"Vp1 Vncasc 0 {p.get('Vncasc', VNCASC0):.4f}")


imp_core_original = imp.core
imp.core = core
E = base.E; CF = base.CF


def window_scan(r, p, ns, label, levels):
    """Small-signal gain and charge fraction against the output level, gate fixed."""
    rows = []
    for v in levels:
        shift = v - p_gate(ns)
        d = r(f'win_{label}_{v:.2f}'.replace('.', 'p'), imp.ac_deck(p, ns).replace(
            f'Vshift out r {imp.VSHIFT:.6f}', f'Vshift out r {shift:.6f}'), keep=False)
        f = d['frequency'].real; A = np.abs(-d['V(out)'] / d['V(in)']); tf = np.abs(d['V(out)'])
        i = int(np.argmin(np.abs(f - 1.0)))
        rows.append(dict(out_V=v, A0=float(A[0]), charge_fraction=float(tf[i] * 2 * np.pi * f[i] * CF)))
    return rows


def p_gate(ns):
    for tok in ns.split():
        if tok.startswith('V(in)='):
            return float(tok[6:])
    raise ValueError(ns)


def main():
    sim = resolve_simulator()
    for d in (DESIGN, imp.RUN):
        d.mkdir(parents=True, exist_ok=True)
    shutil.copy(base.MODEL, imp.RUN / 'corrected.inc'); shutil.copy(base.MODEL, DESIGN / 'corrected.inc')
    r = imp.Runner(sim)
    prev = json.loads((ROOT / 'results/ltspice/design/VH1/improvements/summary.json').read_text())
    combo = dict(prev['combined_W4_23u']['params'])
    variants = [('L3_4u_casc0.66', dict(L3=4e-6, Vncasc=0.6621)),
                ('L3_8u_casc0.66', dict(L3=8e-6, Vncasc=0.6621)),
                ('L3_4u_casc0.76', dict(L3=4e-6, Vncasc=0.7621)),
                ('L3_8u_casc0.76', dict(L3=8e-6, Vncasc=0.7621)),
                ('L3_16u_casc0.76', dict(L3=16e-6, Vncasc=0.7621)),
                ('L3_8u_casc0.86', dict(L3=8e-6, Vncasc=0.8621)),
                ('L3_16u_casc0.86', dict(L3=16e-6, Vncasc=0.8621))]
    results = {}
    for label, change in variants:
        p = dict(combo, **change)
        try:
            res = imp.measure(r, p, label)
            op = imp.operating_point(r, dict(res['params']), f'op_{label}')
            res['med_V'] = op['med']; res['s_V'] = op['s']
            results[label] = res
        except Exception as exc:
            results[label] = dict(error=str(exc))
        v = results[label]
        print(label, {k: (round(v[k], 4) if isinstance(v.get(k), float) else v.get(k)) for k in ('Vnmos', 'Vnload', 'A0', 'med_V', 'C_amp_pF', 'charge_fraction')},
              {k: round(x, 1) for k, x in v.get('enc', {}).items()}, flush=True)

    # output window: three representative circuits
    levels = [0.886, 0.80, 0.70, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30]
    windows = {}
    for label in ('L3_4u_casc0.66', 'L3_8u_casc0.76', 'L3_8u_casc0.86'):
        if 'enc' not in results[label]:
            continue
        p = dict(results[label]['params'])
        op = imp.operating_point(r, p, f'op_{label}')
        windows[label] = window_scan(r, p, imp.nodeset(op), label, levels)
        print('window', label, [(w['out_V'], round(w['A0']), round(w['charge_fraction'], 3)) for w in windows[label]], flush=True)
    results['windows'] = windows
    results['runs'] = r.n
    (OUT / 'summary.json').write_text(json.dumps(results, indent=1))
    report(results, [l for l, _ in variants])


def report(res, order):
    key5 = 'frame difference of 5 ms averages'; key20 = 'frame difference of 20 ms averages'
    L = ['# VH1 cascode bias and load length', '',
         'Combined circuit with the 23 um cascode (M1 110 um, M2 2 um, M4 23 um); M3 length and cascode gate varied.',
         'Every variant re-tuned with the reset loop closed: Vnload (M3) for I(M1) = 210.3 uA, Vnmos (M2) for V(out) = 0.886 V.',
         'ENC in ions for 20 ms frames, frame difference of 5 ms averages and of 20 ms averages.', '',
         '| variant | Vncasc (V) | Vnmos (V) | Vnload (V) | V(med) (V) | A0 | C_amp (pF) | charge fraction | uV/ion | ENC 5 ms | ENC 20 ms |',
         '|---|---|---|---|---|---|---|---|---|---|---|']
    for l in order:
        v = res.get(l, {})
        if 'enc' not in v:
            L.append(f"| {l} | failed: {v.get('error','')[:60]} | | | | | | | | | |"); continue
        L.append(f"| {l} | {v['params']['Vncasc']:.2f} | {v['Vnmos']:.3f} | {v['Vnload']:.3f} | {v['med_V']:.3f} | {v['A0']:.0f} | {v['C_amp_pF']:.2f} | "
                 f"{v['charge_fraction']:.3f} | {v['gain_uV_per_ion']:.1f} | {v['enc'][key5]:.1f} | {v['enc'][key20]:.1f} |")
    L += ['', '## Output window: open-loop gain and charge fraction against the output level (3 fF, gate at its bias)', '']
    for l, rows in res.get('windows', {}).items():
        L += [f'### {l}', '', '| out (V) | A0 | charge fraction |', '|---|---|---|']
        L += [f"| {w['out_V']:.3f} | {w['A0']:.0f} | {w['charge_fraction']:.3f} |" for w in rows]
        L.append('')
    L += [f"Simulator runs: {res['runs']}.", '']
    (OUT / 'summary.md').write_text('\n'.join(L)); print('\n'.join(L))


if __name__ == '__main__':
    main()
