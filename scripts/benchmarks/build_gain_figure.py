#!/usr/bin/env python3
"""Normalized open-loop amplifier gain for the paper; no device prediction."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import os
import tomllib
import numpy as np
from ion_signal.realistic_csa import CSA

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'work/matplotlib-cache'))
os.environ.setdefault('XDG_CACHE_HOME', str(ROOT/'work/cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    config = tomllib.loads((ROOT/'config/realistic_csa.toml').read_text())
    model = CSA(**config['electronics'])
    ratio = np.geomspace(.01, 100., 1001)
    magnitude = abs(model.amplifier_gain(ratio*model.amplifier_pole_Hz))/model.open_loop_gain
    at_pole = abs(model.amplifier_gain(model.amplifier_pole_Hz))/model.open_loop_gain
    gen = ROOT/'notes/physics/generated'
    preview = ROOT/'notes/physics/work'
    gen.mkdir(parents=True, exist_ok=True)
    preview.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42}):
        fig, ax = plt.subplots(figsize=(5.7, 3.0), constrained_layout=True)
        ax.semilogx(ratio, magnitude, color='#17877E', linewidth=2)
        ax.axhline(1, color='0.6', linewidth=.8, linestyle=':')
        ax.vlines(1, 0, at_pole, color='0.45', linewidth=.9, linestyle='--')
        ax.hlines(at_pole, ratio[0], 1, color='0.45', linewidth=.9, linestyle='--')
        ax.plot(1, at_pole, 'o', color='#17877E', markersize=5)
        ax.annotate(r'$|A|/A_0=1/\sqrt{2}$', xy=(1, at_pole), xytext=(2.3, .84),
                    arrowprops={'arrowstyle': '-', 'color': '0.4'}, fontsize=10)
        ax.set(xlim=(ratio[0], ratio[-1]), ylim=(0, 1.08),
               xlabel=r'Frequency $f/f_a$',
               ylabel=r'Normalized gain $|A(f)|/A_0$')
        ax.set_xticks([.01, .1, 1, 10, 100], labels=['0.01', '0.1', '1', '10', '100'])
        ax.grid(axis='y', alpha=.15)
        fig.savefig(gen/'amplifier_gain.pdf', metadata={'CreationDate': None, 'ModDate': None})
        fig.savefig(preview/'amplifier_gain.png', dpi=180)
        plt.close(fig)
    print('Created', gen/'amplifier_gain.pdf')


if __name__ == '__main__':
    main()
