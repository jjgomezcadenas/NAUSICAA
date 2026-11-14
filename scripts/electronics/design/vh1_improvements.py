"""VH1 improvements: device sweeps on the reset CSA prototype.

Three studies on the VH1 core (see vh1_reset_csa.py):
1. Noise of the current source M2 and of the load M3: length sweeps with the
   gate biases re-tuned to keep the M1 current (210 uA) and the output level
   (0.886 V).
2. Input capacitance: width sweeps of the cascode M4 and of the input M1.
3. A buildable reset level shift: resistor and current source instead of the
   ideal voltage source, tested in transient.

For every variant the script tunes Vnmos (M1 current) and Vnload (output DC)
with closed-loop operating points, then runs the AC deck (open-loop gain,
input capacitance, charge fraction) and the noise deck (ENC for 20 ms frames).

    PYTHONPATH=src work/event-viewer-venv/bin/python scripts/electronics/design/vh1_improvements.py
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
from ion_signal.spice_acquisition import run  # noqa: E402
import vh1_reset_csa as base  # noqa: E402

DESIGN = ROOT / 'LTSpice/designs/VH1/improvements'
OUT = ROOT / 'results/ltspice/design/VH1/improvements'
RUN = OUT / 'run'
E = base.E; CF = base.CF; T = base.T_FRAME
TARGET_I = 2.1031e-4
TARGET_OUT = base.BASELINE_V
VSHIFT = base.VSHIFT

BASE = dict(W1=220e-6, L1=0.5e-6, W2=27e-6, L2=0.5e-6, M2=3, W3=200e-6, L3=2e-6,
            W4=2.88e-6, L4=1.44e-6, Vnmos=1.6, Vnload=0.538477478921)


def core(p, voff=0.0):
    return f"""* VH1 variant
M1 s in med Vdda 180nmP l={p['L1']:g} w={p['W1']:g} m=1
Vdd Vdda 0 1.8
M2 Vdda Vnmos s 0 180nmN l={p['L2']:g} w={p['W2']:g} m={p['M2']}
Vnmos Vnmos 0 {p['Vnmos']:.6f}
M3 med Vnload 0 0 180nmN l={p['L3']:g} w={p['W3']:g} m=1
Vnload Vnload 0 {p['Vnload']:.6f}
M4 med Vncasc out 0 180nmN l={p['L4']:g} w={p['W4']:g} m=1
M5 vd6 Vp5 out Vdda 180nmP l=0.5u w=220u m=1
M6 Vdda Vp6 vd6 Vdda 180nmP l=6u w=1u m=1
Vp6 Vp6 0 1.1
Vp5 Vp5 0 0.8
Vp1 Vncasc 0 0.6621
Cd1 x 0 {p.get('Cd', '1.25p')}
Voff in x {voff}
.lib corrected.inc
.temp 20
.options reltol=1e-6 abstol=1e-19 vntol=1e-10 numdgt=15 plotwinsize=0
.options gmin=1e-18
"""


def nodeset(op):
    keys = ['in', 'x', 'out', 'med', 's', 'vd6']
    return '.nodeset ' + ' '.join(f'V({k})={op[k]:.6f}' for k in keys if k in op) + '\n'


def op_deck(p, ns=''):
    return core(p) + ns + f"""Cf out x 3f
Vshift out r {VSHIFT:.6f}
Rrs r in 1k
.op
.end
"""


def ac_deck(p, ns, cd=None):
    q = dict(p)
    if cd is not None:
        q['Cd'] = cd
    return core(q) + ns + f"""Cf out x 3f
Vshift out r {VSHIFT:.6f}
Rdc r x 1e15 noiseless
Is1 0 x 0 AC 1
.ac dec 40 1e-3 1e9
.save V(x) V(in) V(out)
.end
"""


def noise_deck(p, ns):
    return core(p) + ns + f"""Cf out x 3f
Vshift out r {VSHIFT:.6f}
Rdc r x 1e15 noiseless
Is1 0 x 0 AC 1
.noise V(out) Is1 dec 100 1e-4 1e6
.end
"""


class Runner:
    def __init__(self, sim):
        self.sim = sim; self.n = 0

    def __call__(self, name, text, keep=True):
        self.n += 1
        deck = RUN / f'{name}.cir'; deck.write_text(text)
        if keep:
            (DESIGN / f'{name}.cir').write_text(text)
        data, _ = run(self.sim, deck, timeout_s=900)
        return data


def operating_point(r, p, name, ns=''):
    d = r(name, op_deck(p, ns), keep=False)
    v = {k[2:-1]: float(np.real(d[k][0])) for k in d if k.startswith('V(')}
    return dict(I1=abs(float(np.real(d['Id(M1)'][0]))), **v)


def secant(f, x0, x1, target, tol, maxit=8):
    y0 = f(x0) - target; y1 = f(x1) - target
    for _ in range(maxit):
        if abs(y1) < tol:
            return x1
        x2 = x1 - y1 * (x1 - x0) / (y1 - y0)
        x0, y0, x1 = x1, y1, x2
        y1 = f(x1) - target
    return x1


def tune(r, p, label):
    """Re-tune Vnmos for the M1 current and Vnload for the output level."""
    p = dict(p); log = []
    for round_ in range(3):
        # With the reset loop closed, Vnload (M3) sets the M1 current and Vnmos (M2) the output level.
        def current(v):
            q = dict(p, Vnload=v); return operating_point(r, q, f'op_{label}')['I1']
        p['Vnload'] = secant(current, p['Vnload'], p['Vnload'] + 0.005, TARGET_I, 1e-7)
        def out(v):
            q = dict(p, Vnmos=v); return operating_point(r, q, f'op_{label}')['out']
        p['Vnmos'] = secant(out, p['Vnmos'], p['Vnmos'] + 0.02, TARGET_OUT, 1e-3)
        op = operating_point(r, p, f'op_{label}')
        log.append((p['Vnmos'], p['Vnload'], op['I1'], op['out']))
        if abs(op['I1'] - TARGET_I) < 1e-7 and abs(op['out'] - TARGET_OUT) < 1e-3:
            break
    return p, op, log


def measure(r, p, label):
    p, op, log = tune(r, p, label)
    ns = nodeset(op)
    d = r(f'ac_{label}', ac_deck(p, ns))
    f = d['frequency'].real; A = np.abs(-d['V(out)'] / d['V(in)']); zin = np.abs(d['V(x)']); tf = np.abs(d['V(out)'])
    i = int(np.argmin(np.abs(f - 1.0)))
    A0 = float(A[0]); c_total = 1 / (2 * np.pi * f[i] * zin[i]) - (1 + A0) * CF
    fraction = float(tf[i] * 2 * np.pi * f[i] * CF)
    d2 = r(f'ac_{label}_nopad', ac_deck(p, ns, cd='1e-18'), keep=False)
    zin2 = np.abs(d2['V(x)']); A02 = float(np.abs(-d2['V(out)'] / d2['V(in)'])[0])
    c_amp = 1 / (2 * np.pi * f[i] * zin2[i]) - (1 + A02) * CF
    f0 = float(np.interp(A0 / np.sqrt(2), A[::-1], f[::-1]))
    unity = float(np.interp(1.0, A[::-1], f[::-1]))
    dn = r(f'noise_{label}', noise_deck(p, ns))
    fn = dn['frequency'].real; asd = dn['V(onoise)'].real
    gain = fraction * E / CF
    enc = {k: base.enc_weighted(fn, asd, w, gain) for k, w in base.WEIGHTS.items()}
    np.savetxt(OUT / f'noise_{label}.csv', np.column_stack([fn, asd]), delimiter=',',
               header='frequency_Hz,output_V_per_sqrtHz', comments='')
    return dict(params={k: (float(v) if isinstance(v, (int, float)) else v) for k, v in p.items()},
                Vnmos=p['Vnmos'], Vnload=p['Vnload'], I1_uA=op['I1'] * 1e6, out_V=op['out'], in_V=op['in'],
                A0=A0, f_pole_Hz=f0, f_unity_Hz=unity, C_total_pF=c_total * 1e12, C_amp_pF=c_amp * 1e12,
                charge_fraction=fraction, gain_uV_per_ion=gain * 1e6, enc=enc, tuning=log)


def practical_reset_deck(p, ns, rshift, ishift):
    """Reset through a resistor and a current sink instead of the ideal source."""
    reset_end = 1e-3; ends = [reset_end + k * T for k in (1, 2, 3)]
    src = base.rect_source([(reset_end + 0.005, reset_end + 0.015, 100), (reset_end + 0.025, reset_end + 0.035, 200)])
    ph = base.phases(reset_end, ends)
    return core(p) + ns + f"""Cf1 x a1 3f
Cf2 x a2 3f
S1r a1 vref p1r 0 SWM
S1o a1 out p1o 0 SWM
S2r a2 vref p2r 0 SWM
S2o a2 out p2o 0 SWM
Rshift out r {rshift:g}
Ishift r 0 {ishift:g}
Srs r in prs 0 SWM
Vref vref 0 {TARGET_OUT:.6f}
.model SWM SW(Ron=1k Roff=1e18 Vt=0.5 Vh=0.05)
Vprs prs 0 {base.pwl01(ph['prs'])}
Vp1o p1o 0 {base.pwl01(ph['p1o'])}
Vp1r p1r 0 {base.pwl01(ph['p1r'])}
Vp2o p2o 0 {base.pwl01(ph['p2o'])}
Vp2r p2r 0 {base.pwl01(ph['p2r'])}
{src}
.tran 0 {ends[-1]:.9g} 0 1e-5
.save V(x) V(in) V(out) V(r)
.end
""", ends


def main():
    sim = resolve_simulator()
    for d in (DESIGN, RUN):
        d.mkdir(parents=True, exist_ok=True)
    shutil.copy(base.MODEL, RUN / 'corrected.inc'); shutil.copy(base.MODEL, DESIGN / 'corrected.inc')
    r = Runner(sim)
    results = {}
    results['base'] = measure(r, BASE, 'base')
    sweeps = {
        'M2 length': [('L2_1u', dict(L2=1e-6)), ('L2_2u', dict(L2=2e-6)), ('L2_4u', dict(L2=4e-6))],
        'M3 length': [('L3_4u', dict(L3=4e-6)), ('L3_8u', dict(L3=8e-6)), ('L3_16u', dict(L3=16e-6))],
        'M4 width': [('W4_5.76u', dict(W4=5.76e-6)), ('W4_11.5u', dict(W4=11.52e-6)), ('W4_23u', dict(W4=23.04e-6)), ('W4_46u', dict(W4=46.08e-6))],
        'M1 width': [('W1_110u', dict(W1=110e-6)), ('W1_55u', dict(W1=55e-6))],
    }
    for group, items in sweeps.items():
        for label, change in items:
            try:
                results[label] = measure(r, dict(BASE, **change), label)
            except Exception as exc:  # keep going; report the failure
                results[label] = dict(error=str(exc))
            print(group, label, {k: results[label].get(k) for k in ('Vnmos', 'Vnload', 'charge_fraction', 'C_amp_pF')},
                  results[label].get('enc'), flush=True)

    # combination: best feasible of each sweep by ENC (5 ms averages), gate biases within supply
    key = 'frame difference of 5 ms averages'
    def best(labels):
        """Lowest ENC among feasible variants; a change must beat the base by 2 %."""
        ok = [l for l in labels if 'enc' in results[l] and results[l]['Vnmos'] < 1.78]
        b = min(ok, key=lambda l: results[l]['enc'][key]) if ok else None
        if b is not None and results[b]['enc'][key] > 0.98 * results['base']['enc'][key]:
            b = 'base'
        return b
    choice = {}
    for group, items in sweeps.items():
        b = best(['base'] + [l for l, _ in items])
        choice[group] = b
    combo = dict(BASE)
    for group, items in sweeps.items():
        for l, change in items:
            if l == choice[group]:
                combo.update(change)
    results['combined'] = measure(r, combo, 'combined')
    results['combined']['choice'] = choice
    # the wider cascode restores the open-loop gain that the longer M3 costs
    results['combined_W4_23u'] = measure(r, dict(combo, W4=23.04e-6), 'combined_W4_23u')

    # practical reset on the combined circuit: 8.85 Mohm and 0.1 uA (the output stage carries 0.38 uA)
    p, op, _ = tune(r, combo, 'combined_reset'); ns = nodeset(op)
    deck, ends = practical_reset_deck(p, ns, 8.85e6, 0.1e-6)
    d = r('swap_practical_reset', deck)
    t = d['time']; vo = d['V(out)']
    vend = base.frame_values(t, vo, ends); drop = TARGET_OUT - vend
    results['practical_reset'] = dict(frame_end_V=vend.tolist(), drop_mV=(drop * 1e3).tolist(),
                                      uV_per_ion=(drop / np.array([100, 300, 300]) * 1e6).tolist(),
                                      frame2_uV_per_ion=float((vend[0] - vend[1]) / 200 * 1e6),
                                      out_during_reset_V=float(np.interp(0.9e-3, t, vo)), in_during_reset_V=float(np.interp(0.9e-3, t, d['V(in)'])))
    np.savetxt(OUT / 'swap_practical_reset.csv', np.column_stack([t, vo, d['V(x)'], d['V(r)']]), delimiter=',',
               header='time_s,V(out),V(x),V(r)', comments='')
    results['runs'] = r.n
    (OUT / 'summary.json').write_text(json.dumps(results, indent=1))
    report(results, sweeps, key)


def report(res, sweeps, key):
    L = ['# VH1 improvements', '',
         'Every variant re-tuned with the reset loop closed: Vnload (M3) for I(M1) = 210.3 uA, Vnmos (M2) for V(out) = 0.886 V.',
         'ENC for 20 ms frames, frame difference of 5 ms averages (as the present calibration) and of 20 ms averages.', '',
         '| variant | Vnmos (V) | Vnload (V) | A0 | C_amp (pF) | charge fraction | uV/ion | ENC 5 ms avg | ENC 20 ms avg |',
         '|---|---|---|---|---|---|---|---|---|']
    order = ['base'] + [l for items in sweeps.values() for l, _ in items] + ['combined', 'combined_W4_23u']
    for l in order:
        v = res.get(l, {})
        if 'enc' not in v:
            L.append(f"| {l} | failed: {v.get('error','')[:60]} | | | | | | | |"); continue
        L.append(f"| {l} | {v['Vnmos']:.3f} | {v['Vnload']:.3f} | {v['A0']:.0f} | {v['C_amp_pF']:.2f} | {v['charge_fraction']:.3f} | {v['gain_uV_per_ion']:.1f} | {v['enc'][key]:.1f} | {v['enc']['frame difference of 20 ms averages']:.1f} |")
    L += ['', f"Combined variant choices: {res['combined'].get('choice')}", '',
          '## Practical reset (8.85 Mohm and 0.1 uA level shift) on the combined circuit', '',
          f"Output during reset {res['practical_reset']['out_during_reset_V']:.4f} V, gate {res['practical_reset']['in_during_reset_V']:.4f} V.",
          f"Frame-end drops: {', '.join(f'{x:.3f}' for x in res['practical_reset']['drop_mV'])} mV for 100, 300, 300 accumulated ions: "
          f"{', '.join(f'{x:.2f}' for x in res['practical_reset']['uV_per_ion'])} uV per ion; "
          f"frame 2 alone {res['practical_reset']['frame2_uV_per_ion']:.2f} uV per ion.", '',
          f"Simulator runs: {res['runs']}.", '']
    (OUT / 'summary.md').write_text('\n'.join(L)); print('\n'.join(L))


if __name__ == '__main__':
    main()
