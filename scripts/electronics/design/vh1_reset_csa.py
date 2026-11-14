"""VH1: reset and two-capacitor CSA prototype built on the csa_v3 amplifier core.

The amplifier M1-M6 with its biases, the 1.25 pF pad capacitance and the
3 fF feedback value are those of LTSpice/releases/csa_v3. The discharge
transistor Mf is removed. Two feedback capacitors alternate between the
feedback position and a parking position at V_REF through ideal switches,
following notes/physics/nausicaa_VH.tex, Section 7.

Decks are written to LTSpice/designs/VH1 and run in
results/ltspice/design/VH1/run. Results (CSV, JSON, PNG, summary.md) go to
results/ltspice/design/VH1.

    PYTHONPATH=src work/event-viewer-venv/bin/python scripts/electronics/design/vh1_reset_csa.py
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
from ion_signal.spice_runtime import resolve_simulator  # noqa: E402
from ion_signal.spice_acquisition import run  # noqa: E402
from ion_signal.spice import charge_noise  # noqa: E402

DESIGN = ROOT / 'LTSpice/designs/VH1'
OUT = ROOT / 'results/ltspice/design/VH1'
RUN = OUT / 'run'
RELEASE = ROOT / 'LTSpice/releases/csa_v3'
MODEL = RELEASE / 'netlists/corrected.inc'
E = 1.602176634e-19
CF = 3e-15
T_FRAME = 0.020
BASELINE_V = 0.8863611214728842   # csa_v3 operating point, V(out)
INPUT_BIAS_V = 0.0013494938197126996  # csa_v3 operating point, V(in): the gate of M1 sits near 0 V
VSHIFT = BASELINE_V - INPUT_BIAS_V    # level shift in the reset path, out = in + VSHIFT
NODESET = '.nodeset V(in)=0.00135 V(x)=0.00135 V(out)=0.886 V(med)=0.184 V(s)=0.836 V(vd6)=1.25\n'
PRESENT_GAIN_V_PER_ION = 40.157e-6  # csa_v3 20 ms frame-difference calibration
PRESENT_TAU_S = 19.5441

CORE = """* VH1 reset CSA prototype: csa_v3 amplifier core (M1-M6 and biases) without Mf
M1 s in med Vdda 180nmP l=0.5u w=220u m=1
Vdd Vdda 0 1.8
M2 Vdda Vnmos s 0 180nmN l=0.5u w=27u m=3
Vnmos Vnmos 0 1.6
M3 med Vnload 0 0 180nmN l=2u w=200u m=1
Vnload Vnload 0 0.538477478921
M4 med Vncasc out 0 180nmN l=1.44u w=2.88u m=1
M5 vd6 Vp5 out Vdda 180nmP l=0.5u w=220u m=1
M6 Vdda Vp6 vd6 Vdda 180nmP l=6u w=1u m=1
Vp6 Vp6 0 1.1
Vp5 Vp5 0 0.8
Vp1 Vncasc 0 0.6621
Cd1 x 0 1.25p
Voff in x {voff}
.lib corrected.inc
.temp 20
.options reltol=1e-6 abstol=1e-19 vntol=1e-10 numdgt=15 plotwinsize=0
* gmin=1e-12 would put 1 Tohm across every off switch: 3 ms leak on 3 fF
.options gmin=1e-18
"""


def ac_deck(cf):
    return CORE.format(voff=0) + NODESET + f"""Cf out x {cf}
Vshift out r {VSHIFT:.6f}
Rdc r x 1e15 noiseless
Is1 0 x 0 AC 1
.ac dec 40 1e-3 1e9
.save V(x) V(in) V(out)
.end
"""


def noise_deck(cf):
    return CORE.format(voff=0) + NODESET + f"""Cf out x {cf}
Vshift out r {VSHIFT:.6f}
Rdc r x 1e15 noiseless
Is1 0 x 0 AC 1
.noise V(out) Is1 dec 100 1e-4 1e6
.end
"""


def pwl01(on_intervals, ramp=1e-7):
    """0/1 control waveform, 1 inside the intervals, with short ramps."""
    pts = [(0.0, 0.0)]
    for a, b in on_intervals:
        if a <= 0:
            pts = [(0.0, 1.0)]
        else:
            pts += [(a, 0.0), (a + ramp, 1.0)]
        pts += [(b, 1.0), (b + ramp, 0.0)]
    return 'PWL(' + ' '.join(f'{t:.9g} {v:g}' for t, v in pts) + ')'


def phases(reset_end, frame_ends, gap=1e-6):
    """Switch schedules: reset until reset_end, then frames ending at frame_ends.
    Odd frames integrate on Cf1, even frames on Cf2. Break before make."""
    starts = [reset_end] + list(frame_ends[:-1])
    fb1 = [(s + gap, e) for k, (s, e) in enumerate(zip(starts, frame_ends)) if k % 2 == 0]
    fb2 = [(s + gap, e) for k, (s, e) in enumerate(zip(starts, frame_ends)) if k % 2 == 1]
    tstop = frame_ends[-1]

    def parked(fb):
        on, t = [], 0.0
        for a, b in fb:
            on.append((t, a - gap))
            t = b + gap
        if t < tstop:
            on.append((t, tstop))
        return on
    return dict(prs=[(0.0, reset_end - gap)], p1o=fb1, p1r=parked(fb1), p2o=fb2, p2r=parked(fb2))


def tran_deck(vref, voff, source_line, reset_end, frame_ends, step=1e-5, cf=CF):
    p = phases(reset_end, frame_ends)
    tstop = frame_ends[-1]
    return CORE.format(voff=voff) + NODESET + f"""Cf1 x a1 {cf}
Cf2 x a2 {cf}
S1r a1 vref p1r 0 SWM
S1o a1 out p1o 0 SWM
S2r a2 vref p2r 0 SWM
S2o a2 out p2o 0 SWM
Vshift out r {VSHIFT:.6f}
Srs r in prs 0 SWM
Vref vref 0 {vref}
.model SWM SW(Ron=1k Roff=1e18 Vt=0.5 Vh=0.05)
Vprs prs 0 {pwl01(p['prs'])}
Vp1o p1o 0 {pwl01(p['p1o'])}
Vp1r p1r 0 {pwl01(p['p1r'])}
Vp2o p2o 0 {pwl01(p['p2o'])}
Vp2r p2r 0 {pwl01(p['p2r'])}
{source_line}
.tran 0 {tstop:.9g} 0 {step:g}
.save V(x) V(in) V(out) V(a1) V(a2) I(Is1)
.end
"""


def rect_source(frames):
    """Rectangular current pulses: list of (t_start, t_end, ions)."""
    pts = [(0.0, 0.0)]
    for a, b, n in frames:
        i = n * E / (b - a)
        pts += [(a, 0.0), (a + 1e-6, i), (b - 1e-6, i), (b, 0.0)]
    return 'Is1 0 x PWL(' + ' '.join(f'{t:.9g} {v:.6g}' for t, v in pts) + ')'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_and_run(sim, name, text, timeout_s=900):
    (DESIGN / f'{name}.cir').write_text(text)
    deck = RUN / f'{name}.cir'
    deck.write_text(text)
    data, _ = run(sim, deck, timeout_s=timeout_s)
    return data


def frame_values(time, v, frame_ends, before=5e-5):
    return np.array([np.interp(t - before, time, v) for t in frame_ends])


def enc_weighted(f, asd, weight, gain_v_per_ion, order=8, half_period=T_FRAME):
    """ENC in ions from an output voltage spectral density and a weighting |W|^2."""
    f = np.asarray(f); asd = np.asarray(asd)
    boundaries = np.unique(np.r_[f, np.arange(f[0], f[-1], 0.5 / half_period)])
    x, w = np.polynomial.legendre.leggauss(order)
    half = np.diff(boundaries) / 2; mid = (boundaries[:-1] + boundaries[1:]) / 2
    nodes = mid[:, None] + half[:, None] * x
    spectrum = np.exp(np.interp(np.log(nodes), np.log(f), 2 * np.log(np.maximum(asd, 1e-150))))
    variance = np.sum(half[:, None] * w * spectrum * weight(nodes))
    return float(np.sqrt(variance) / gain_v_per_ion)


def sinc2(x):
    return np.sinc(x / np.pi) ** 2


WEIGHTS = {
    'frame difference of point samples': lambda f: 4 * np.sin(np.pi * f * T_FRAME) ** 2,
    'frame difference of 5 ms averages': lambda f: 4 * np.sin(np.pi * f * T_FRAME) ** 2 * sinc2(np.pi * f * 0.005),
    'frame difference of 20 ms averages': lambda f: 4 * np.sin(np.pi * f * T_FRAME) ** 2 * sinc2(np.pi * f * T_FRAME),
}


def main():
    sim = resolve_simulator()
    for d in (DESIGN, RUN):
        d.mkdir(parents=True, exist_ok=True)
    shutil.copy(MODEL, DESIGN / 'corrected.inc')
    shutil.copy(MODEL, RUN / 'corrected.inc')
    (OUT / '.gitignore').write_text('run/\n')
    summary = {'model_sha256': sha(MODEL), 'simulator': str(sim)}

    # ---- AC: open-loop gain from the closed loop, A = -V(out)/V(in) ----
    ac = {}
    for cf, label in ((3e-15, '3f'), (30e-15, '30f')):
        d = write_and_run(sim, f'ac_cf{label}', ac_deck(f'{cf:g}'))
        f = d['frequency'].real
        A = -d['V(out)'] / d['V(in)']
        zin = d['V(x)']           # current source is 1 A: V(x) is Z_in
        tf = d['V(out)']          # V/A: transimpedance
        np.savetxt(OUT / f'ac_cf{label}.csv',
                   np.column_stack([f, np.abs(A), np.degrees(np.angle(A)), np.abs(zin), np.abs(tf)]),
                   delimiter=',', header='frequency_Hz,|A|,phase_A_deg,|Zin|_ohm,|Vout/Iin|_ohm', comments='')
        A0 = float(np.abs(A[0]))
        f0 = float(np.interp(A0 / np.sqrt(2), np.abs(A)[::-1], f[::-1]))
        unity = float(np.interp(1.0, np.abs(A)[::-1], f[::-1])) if np.abs(A).min() < 1 else None
        # closed-loop bandwidth: where |Vout/Iin| falls to 1/sqrt2 of 1/(w Cf)
        ideal = 1 / (2 * np.pi * f * cf)
        ratio = np.abs(tf) / ideal
        hi = f > 1e3
        fcl = float(np.interp(1 / np.sqrt(2), ratio[hi][::-1], f[hi][::-1])) if ratio[hi].min() < 0.7 else None
        plateau = float(np.min(np.abs(zin)[(f > 1e3) & (f < 1e8)])) if label == '3f' else None
        ac[label] = dict(A0=A0, f_pole_Hz=f0, f_unity_Hz=unity, GBP_Hz=A0 * f0,
                         closed_loop_bandwidth_Hz=fcl, Zin_minimum_1kHz_100MHz_ohm=plateau)
    summary['ac'] = ac
    # amplifier input capacitance from the input impedance without the pad
    d = write_and_run(sim, 'ac_cf3f_nopad', ac_deck('3f').replace('Cd1 x 0 1.25p', 'Cd1 x 0 1e-18'))
    f = d['frequency'].real; zin = np.abs(d['V(x)']); A0 = float(np.abs(-d['V(out)'] / d['V(in)'])[0])
    i = int(np.argmin(np.abs(f - 1.0)))
    c_total = 1 / (2 * np.pi * f[i] * zin[i])
    summary['amplifier_input_capacitance_F'] = c_total - (1 + A0) * 3e-15
    d = write_and_run(sim, 'ac_cf3f_pad', ac_deck('3f'))
    zin = np.abs(d['V(x)']); A0 = float(np.abs(-d['V(out)'] / d['V(in)'])[0])
    c_total = 1 / (2 * np.pi * f[i] * zin[i])
    summary['total_input_capacitance_with_pad_F'] = c_total - (1 + A0) * 3e-15
    summary['charge_fraction_predicted'] = (1 + A0) * 3e-15 / c_total
    tf = np.abs(d['V(out)']); summary['charge_fraction_ac_1Hz'] = float(tf[i] * 2 * np.pi * f[i] * 3e-15)

    # ---- Noise: VH1 core in feedback, and the present circuit with Mf ----
    d = write_and_run(sim, 'noise_cf3f', noise_deck('3f'))
    f_vh1 = d['frequency'].real; asd_vh1 = d['V(onoise)'].real
    np.savetxt(OUT / 'noise_vh1.csv', np.column_stack([f_vh1, asd_vh1]), delimiter=',',
               header='frequency_Hz,output_V_per_sqrtHz', comments='')
    present = RUN / 'noise_present.cir'
    present.write_text((RELEASE / 'netlists/noise_200.cir').read_text())
    dp, _ = run(sim, present, timeout_s=900)
    f_p = dp['frequency'].real; asd_p = dp['V(onoise)'].real
    np.savetxt(OUT / 'noise_present.csv', np.column_stack([f_p, asd_p]), delimiter=',',
               header='frequency_Hz,output_V_per_sqrtHz', comments='')

    # ---- Transients ----
    reset_end = 1e-3
    ends3 = [reset_end + k * T_FRAME for k in (1, 2, 3)]
    rect3 = rect_source([(reset_end + 0.005, reset_end + 0.015, 100),
                         (reset_end + 0.025, reset_end + 0.035, 200)])
    runs = {
        'swap_vref0886_voff0': tran_deck(BASELINE_V, 0.0, rect3, reset_end, ends3),
        'swap_vref0886_voff50mV': tran_deck(BASELINE_V, 0.05, rect3, reset_end, ends3),
        'swap_vref1000_voff0': tran_deck(1.0, 0.0, rect3, reset_end, ends3),
    }
    expected = {k: [100, 300, 300] for k in runs}
    # realistic 100-ion packet from the release: collection starts at 10 ms
    pwl = RELEASE / 'inputs/20ms_N100_r1.pwl'
    shutil.copy(pwl, RUN / pwl.name)
    ends_real = [0.009 + k * T_FRAME for k in (1, 2, 3)]
    runs['swap_packet100'] = tran_deck(BASELINE_V, 0.0, f'Is1 0 x PWL FILE="{pwl.name}"', 0.009, ends_real)
    expected['swap_packet100'] = [100, 100, 100]
    # range: 2000 ions per frame, six frames
    ends6 = [reset_end + k * T_FRAME for k in range(1, 7)]
    runs['swap_range2000'] = tran_deck(BASELINE_V, 0.0, rect_source(
        [(reset_end + k * T_FRAME + 0.005, reset_end + k * T_FRAME + 0.015, 2000) for k in range(6)]),
        reset_end, ends6)
    expected['swap_range2000'] = [2000 * (k + 1) for k in range(6)]

    tran = {}
    waves = {}
    for name, text in runs.items():
        d = write_and_run(sim, name, text)
        t = d['time']; vo = d['V(out)']; vx = d['V(x)']
        fe = ends6 if name == 'swap_range2000' else (ends_real if name == 'swap_packet100' else ends3)
        vref = 1.0 if 'vref1000' in name else BASELINE_V
        vend = frame_values(t, vo, fe)
        drop = vref - vend
        exp_ions = np.array(expected[name], float)
        tran[name] = dict(frame_end_V=vend.tolist(), drop_V=drop.tolist(),
                          accumulated_ions_expected=exp_ions.tolist(),
                          uV_per_ion=(drop / exp_ions * 1e6).tolist(),
                          vref_V=vref, x_during_reset_V=float(np.interp(reset_end * 0.9 if name != 'swap_packet100' else 0.008, t, vx)),
                          x_end_V=float(vx[-1]))
        waves[name] = (t, vo, vx)
        np.savetxt(OUT / f'{name}.csv', np.column_stack([t, vo, vx, d['V(a1)'], d['V(a2)']]), delimiter=',',
                   header='time_s,V(out),V(x),V(a1),V(a2)', comments='')
    summary['transient'] = tran

    # gain per ion from the first frame of the rectangular test, and ENC
    g_vh1 = tran['swap_vref0886_voff0']['uV_per_ion'][0] * 1e-6
    summary['vh1_gain_V_per_ion'] = g_vh1
    summary['ideal_gain_V_per_ion'] = E / CF
    enc = {}
    for label, w in WEIGHTS.items():
        enc[label] = dict(vh1=enc_weighted(f_vh1, asd_vh1, w, g_vh1),
                          present=enc_weighted(f_p, asd_p, w, PRESENT_GAIN_V_PER_ION))
    enc['present calibrated estimator (point samples, alpha)'] = dict(
        present=charge_noise(f_p, asd_p, T_FRAME, CF, PRESENT_TAU_S) / E)
    summary['enc_ions'] = enc

    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1))
    report(summary, f_vh1, asd_vh1, f_p, asd_p, waves, ends3, ends6)


def report(s, f1, a1, fp, ap, waves, ends3, ends6):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ac = np.loadtxt(OUT / 'ac_cf3f.csv', delimiter=',', skiprows=1)
    ax[0].loglog(ac[:, 0], ac[:, 1], label='|A|, Cf=3 fF')
    ac30 = np.loadtxt(OUT / 'ac_cf30f.csv', delimiter=',', skiprows=1)
    ax[0].loglog(ac30[:, 0], ac30[:, 1], '--', label='|A|, Cf=30 fF')
    ax[0].set_xlabel('frequency (Hz)'); ax[0].set_ylabel('open-loop gain'); ax[0].legend(); ax[0].set_title('VH1 core: open-loop gain')
    ax[1].loglog(fp, ap, label='present circuit (with Mf)')
    ax[1].loglog(f1, a1, label='VH1 core (Mf removed)')
    ax[1].set_xlabel('frequency (Hz)'); ax[1].set_ylabel('output noise (V/$\\sqrt{\\mathrm{Hz}}$)'); ax[1].legend()
    ax[1].set_title('Output noise spectral density')
    fig.tight_layout(); fig.savefig(OUT / 'ac_and_noise.png', dpi=130); plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    t, vo, vx = waves['swap_vref0886_voff0']
    ax[0].plot(t * 1e3, vo * 1e3, label='V(out)')
    for e in ends3:
        ax[0].axvline(e * 1e3, color='gray', ls=':')
    ax[0].set_xlabel('time (ms)'); ax[0].set_ylabel('mV'); ax[0].set_title('Reset, 100 ions in frame 1, 200 in frame 2, 0 in frame 3')
    ax[0].legend()
    t, vo, vx = waves['swap_range2000']
    ax[1].plot(t * 1e3, vo, label='V(out)')
    ax[1].plot(t * 1e3, vx, label='V(x) input node')
    for e in ends6:
        ax[1].axvline(e * 1e3, color='gray', ls=':')
    ax[1].set_xlabel('time (ms)'); ax[1].set_ylabel('V'); ax[1].set_title('2000 ions per frame, six frames'); ax[1].legend()
    fig.tight_layout(); fig.savefig(OUT / 'transients.png', dpi=130); plt.close(fig)

    L = ['# VH1 reset CSA prototype: results', '',
         f"Core: csa_v3 amplifier (M1-M6, biases, Cd = 1.25 pF), Mf removed; model `{s['model_sha256'][:12]}`.", '',
         '## Open-loop response (A = -V(out)/V(in), closed loop with Cf)', '',
         '| Cf | A0 | pole (Hz) | unity-gain (Hz) | GBP (Hz) | closed-loop bandwidth (Hz) |', '|---|---|---|---|---|---|']
    for k, v in s['ac'].items():
        L.append(f"| {k} | {v['A0']:.0f} | {v['f_pole_Hz']:.3g} | {v['f_unity_Hz'] if v['f_unity_Hz'] is None else format(v['f_unity_Hz'],'.3g')} | {v['GBP_Hz']:.3g} | {v['closed_loop_bandwidth_Hz'] if v['closed_loop_bandwidth_Hz'] is None else format(v['closed_loop_bandwidth_Hz'],'.3g')} |")
    L += ['', f"Input impedance minimum between 1 kHz and 100 MHz (Cf = 3 fF): {s['ac']['3f']['Zin_minimum_1kHz_100MHz_ohm']:.3g} ohm.", '',
          '## Transients (output read 50 us before each frame end)', '',
          '| run | Vref (V) | frame-end V(out) | drop from Vref (mV) | expected accumulated ions | uV per ion |', '|---|---|---|---|---|---|']
    for k, v in s['transient'].items():
        L.append(f"| {k} | {v['vref_V']:.3f} | {', '.join(f'{x:.4f}' for x in v['frame_end_V'])} | {', '.join(f'{x*1e3:.3f}' for x in v['drop_V'])} | {', '.join(str(int(x)) for x in v['accumulated_ions_expected'])} | {', '.join(f'{x:.2f}' for x in v['uV_per_ion'])} |")
    L += ['', f"Ideal gain e/Cf = {s['ideal_gain_V_per_ion']*1e6:.2f} uV per ion; measured {s['vh1_gain_V_per_ion']*1e6:.2f} uV per ion "
          f"(fraction {s['vh1_gain_V_per_ion']/s['ideal_gain_V_per_ion']:.3f}; AC transimpedance at 1 Hz gives {s['charge_fraction_ac_1Hz']:.3f}).",
          f"Input capacitance seen at the gate: amplifier alone {s['amplifier_input_capacitance_F']*1e12:.2f} pF, with the 1.25 pF pad {s['total_input_capacitance_with_pad_F']*1e12:.2f} pF; "
          f"predicted charge fraction (1+A)Cf/[C_in+(1+A)Cf] = {s['charge_fraction_predicted']:.3f}.", '',
          '## ENC (ions RMS, 20 ms frames)', '', '| weighting | VH1 | present circuit |', '|---|---|---|']
    for k, v in s['enc_ions'].items():
        L.append(f"| {k} | {v.get('vh1', float('nan')):.1f} | {v['present']:.1f} |")
    L += ['', 'VH1 uses the measured gain per ion of the first test frame; the present circuit uses its calibrated 40.16 uV per ion.',
          '', '![AC and noise](ac_and_noise.png)', '', '![Transients](transients.png)', '']
    (OUT / 'summary.md').write_text('\n'.join(L))
    print('\n'.join(L))


if __name__ == '__main__':
    main()
