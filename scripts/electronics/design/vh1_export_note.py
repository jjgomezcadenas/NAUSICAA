"""Export VH1 results to the physics note: LaTeX macros, a table body and PDF figures.

Reads results/ltspice/design/VH1/summary.json, improvements/summary.json and
the saved spectra and waveforms; writes
notes/physics/generated/vh1_numbers.tex, vh1_sweeps.tex, vh1_buckets.tex and
notes/physics/figures/vh1_ac_noise.pdf, vh1_transients.pdf. No simulation.

    PYTHONPATH=src work/event-viewer-venv/bin/python scripts/electronics/design/vh1_export_note.py
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import vh1_reset_csa as base  # noqa: E402

RES = ROOT / 'results/ltspice/design/VH1'
GEN = ROOT / 'notes/physics/generated'
FIG = ROOT / 'notes/physics/figures'
FIVE = 'frame difference of 5 ms averages'
FULL = 'frame difference of 20 ms averages'
POINT = 'frame difference of point samples'


def macro(name, value):
    return f'\\newcommand{{\\{name}}}{{{value}}}\n'


def main():
    GEN.mkdir(exist_ok=True)
    s = json.loads((RES / 'summary.json').read_text())
    imp = json.loads((RES / 'improvements/summary.json').read_text())
    ac3 = s['ac']['3f']; ac30 = s['ac']['30f']
    tr = s['transient']
    m = ''
    m += macro('vhAzero', f"{ac3['A0']:.0f}")
    m += macro('vhPolekHz', f"{ac3['f_pole_Hz']/1e3:.0f}")
    m += macro('vhUnityMHz', f"{ac3['f_unity_Hz']/1e6:.0f}")
    m += macro('vhClBwkHz', f"{ac3['closed_loop_bandwidth_Hz']/1e3:.0f}")
    m += macro('vhClBwThirtyMHz', f"{ac30['closed_loop_bandwidth_Hz']/1e6:.1f}")
    m += macro('vhZinOhm', f"{ac3['Zin_minimum_1kHz_100MHz_ohm']:.0f}")
    m += macro('vhGmmS', f"{1e3/ac3['Zin_minimum_1kHz_100MHz_ohm']:.2f}")
    m += macro('vhCamp', f"{s['amplifier_input_capacitance_F']*1e12:.2f}")
    m += macro('vhCtot', f"{s['total_input_capacitance_with_pad_F']*1e12:.2f}")
    m += macro('vhFraction', f"{s['charge_fraction_ac_1Hz']:.3f}")
    m += macro('vhGain', f"{s['vh1_gain_V_per_ion']*1e6:.1f}")
    m += macro('vhIdealGain', f"{s['ideal_gain_V_per_ion']*1e6:.1f}")
    t = tr['swap_vref0886_voff0']
    m += macro('vhDropOne', f"{t['drop_V'][0]*1e3:.2f}")
    m += macro('vhDropTwo', f"{t['drop_V'][1]*1e3:.2f}")
    m += macro('vhDropThree', f"{t['drop_V'][2]*1e3:.2f}")
    rg = tr['swap_range2000']
    m += macro('vhRangeDropV', f"{rg['drop_V'][-1]:.2f}")
    m += macro('vhRangeOutV', f"{rg['frame_end_V'][-1]:.2f}")
    m += macro('vhRangeDev', f"{100*(1-rg['uV_per_ion'][-1]/rg['uV_per_ion'][0]):.0f}")
    e = s['enc_ions']
    m += macro('vhEncPresFive', f"{e[FIVE]['present']:.1f}")
    m += macro('vhEncPresFull', f"{e[FULL]['present']:.1f}")
    m += macro('vhEncPresPoint', f"{e[POINT]['present']:.1f}")
    m += macro('vhEncBaseFive', f"{e[FIVE]['vh1']:.1f}")
    m += macro('vhEncBaseFull', f"{e[FULL]['vh1']:.1f}")
    m += macro('vhEncBasePoint', f"{e[POINT]['vh1']:.1f}")
    c = imp['combined_W4_23u']
    m += macro('vhEncImpFive', f"{c['enc'][FIVE]:.1f}")
    m += macro('vhEncImpFull', f"{c['enc'][FULL]:.1f}")
    m += macro('vhImpFraction', f"{c['charge_fraction']:.3f}")
    m += macro('vhImpCamp', f"{c['C_amp_pF']:.2f}")
    m += macro('vhImpAzero', f"{c['A0']:.0f}")
    pr = imp['practical_reset']
    m += macro('vhResetOutV', f"{pr['out_during_reset_V']:.3f}")
    m += macro('vhResetGatemV', f"{pr['in_during_reset_V']*1e3:.1f}")
    m += macro('vhResetGain', f"{pr['frame2_uV_per_ion']:.1f}")
    (GEN / 'vh1_numbers.tex').write_text(m)

    # sweep table rows
    names = {'base': 'VH1 as built', 'L2_1u': 'M2 length 1 \\si{\\um}', 'L2_2u': 'M2 length 2 \\si{\\um}',
             'L2_4u': 'M2 length 4 \\si{\\um}', 'L3_4u': 'M3 length 4 \\si{\\um}', 'L3_8u': 'M3 length 8 \\si{\\um}',
             'L3_16u': 'M3 length 16 \\si{\\um}', 'W4_5.76u': 'M4 width 5.8 \\si{\\um}', 'W4_11.5u': 'M4 width 11.5 \\si{\\um}',
             'W4_23u': 'M4 width 23 \\si{\\um}', 'W4_46u': 'M4 width 46 \\si{\\um}', 'W1_110u': 'M1 width 110 \\si{\\um}',
             'W1_55u': 'M1 width 55 \\si{\\um}', 'combined': 'M2 2, M3 4, M1 110 \\si{\\um}',
             'combined_W4_23u': 'M2 2, M3 4, M1 110, M4 23 \\si{\\um}'}
    rows = []
    for k, label in names.items():
        v = imp[k]
        rows.append(f"{label} & {v['Vnmos']:.2f} & {v['Vnload']:.3f} & {v['A0']:.0f} & {v['C_amp_pF']:.1f} & "
                    f"{v['charge_fraction']:.2f} & {v['enc'][FIVE]:.1f} & {v['enc'][FULL]:.1f} \\\\")
    (GEN / 'vh1_sweeps.tex').write_text('\n'.join(rows) + '\n\\bottomrule%\n')  # the rule travels with the rows: \\input cannot precede \\bottomrule

    # bucket-length table from the spectra
    def sinc2(x):
        return np.sinc(x / np.pi) ** 2

    def enc(f, a, T, avg, gain):
        return base.enc_weighted(f, a, lambda x: 4 * np.sin(np.pi * x * T) ** 2 * sinc2(np.pi * x * avg), gain, half_period=T)
    spectra = [('Present circuit', RES / 'noise_present.csv', base.PRESENT_GAIN_V_PER_ION),
               ('VH1 as built', RES / 'noise_vh1.csv', s['vh1_gain_V_per_ion']),
               ('VH1 with device changes', RES / 'improvements/noise_combined_W4_23u.csv', c['gain_uV_per_ion'] * 1e-6)]
    rows = []
    for label, path, g in spectra:
        d = np.loadtxt(path, delimiter=',', skiprows=1)
        full = [enc(d[:, 0], d[:, 1], T, T, g) for T in (0.01, 0.02, 0.04, 0.08)]
        brief = [enc(d[:, 0], d[:, 1], T, 0.005, g) for T in (0.01, 0.02, 0.04, 0.08)]
        rows.append(label + ' & ' + ' & '.join(f'{x:.1f}' for x in full) + ' & ' + ' & '.join(f'{x:.1f}' for x in brief) + ' \\\\')
    (GEN / 'vh1_buckets.tex').write_text('\n'.join(rows) + '\n\\bottomrule%\n')

    # figures
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif'})
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7))
    a3 = np.loadtxt(RES / 'ac_cf3f.csv', delimiter=',', skiprows=1)
    a30 = np.loadtxt(RES / 'ac_cf30f.csv', delimiter=',', skiprows=1)
    ax[0].loglog(a3[:, 0], a3[:, 1], 'k', label='$C_f=3$ fF')
    ax[0].loglog(a30[:, 0], a30[:, 1], '--', color='tab:orange', label='$C_f=30$ fF')
    ax[0].set_xlabel('frequency (Hz)'); ax[0].set_ylabel('$|A|$'); ax[0].legend(frameon=False); ax[0].set_title('Open-loop gain')
    p = np.loadtxt(RES / 'noise_present.csv', delimiter=',', skiprows=1)
    v = np.loadtxt(RES / 'noise_vh1.csv', delimiter=',', skiprows=1)
    w = np.loadtxt(RES / 'improvements/noise_combined_W4_23u.csv', delimiter=',', skiprows=1)
    ax[1].loglog(p[:, 0], p[:, 1], 'k', label='present, with Mf')
    ax[1].loglog(v[:, 0], v[:, 1], color='tab:orange', label='VH1 as built')
    ax[1].loglog(w[:, 0], w[:, 1], color='tab:purple', label='VH1, device changes')
    ax[1].set_xlabel('frequency (Hz)'); ax[1].set_ylabel('output noise (V/$\\sqrt{\\mathrm{Hz}}$)')
    ax[1].legend(frameon=False); ax[1].set_title('Output noise')
    fig.tight_layout(); fig.savefig(FIG / 'vh1_ac_noise.pdf'); plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7))
    d = np.loadtxt(RES / 'swap_vref0886_voff0.csv', delimiter=',', skiprows=1)
    ax[0].plot(d[:, 0] * 1e3, d[:, 1] * 1e3, 'k')
    for k in (1, 2, 3):
        ax[0].axvline(1 + 20 * k, color='gray', ls=':', lw=0.8)
    ax[0].set_xlabel('time (ms)'); ax[0].set_ylabel('$v_o$ (mV)'); ax[0].set_title('100, 200 and 0 ions in three frames')
    d = np.loadtxt(RES / 'swap_range2000.csv', delimiter=',', skiprows=1)
    ax[1].plot(d[:, 0] * 1e3, d[:, 1], 'k', label='$v_o$')
    ax[1].plot(d[:, 0] * 1e3, d[:, 2], color='tab:orange', label='gate node')
    for k in range(1, 7):
        ax[1].axvline(1 + 20 * k, color='gray', ls=':', lw=0.8)
    ax[1].set_xlabel('time (ms)'); ax[1].set_ylabel('V'); ax[1].set_title('2000 ions per frame'); ax[1].legend(frameon=False)
    fig.tight_layout(); fig.savefig(FIG / 'vh1_transients.pdf'); plt.close(fig)
    print('exported', GEN, FIG)


if __name__ == '__main__':
    main()
