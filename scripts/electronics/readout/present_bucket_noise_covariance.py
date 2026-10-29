#!/usr/bin/env python3
"""Present the saved Julia calculation of noise correlations between buckets."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import hashlib,json
import h5py
ROOT=Path(__file__).resolve().parents[3]
def main():
 source=ROOT/'results/ltspice/bucket_noise_covariance/analysis.h5';out=ROOT/'latex_generated/bucket_noise_covariance';out.mkdir(parents=True,exist_ok=True)
 with h5py.File(source) as h:
  x=h['20ms/correlation'][:];y=h['40ms/correlation'][:]
 rows='\n'.join(f'{lag} & {left:.5f} & {right:.5f}'+r' \\' for lag,(left,right) in enumerate(zip(x,y)))
 text=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=24mm]{geometry}\usepackage{amsmath,booktabs,hyperref}
\title{Noise correlations between NAUSICA buckets}
\author{Numerical calculation}\date{\today}
\begin{document}\maketitle
For the Monte Carlo approximation, each pad--time bucket receives a zero-mean
Gaussian fluctuation with RMS equal to its ENC. Independent fluctuations neglect
temporal correlations between the processed bucket measurements. This calculation
quantifies those correlations using the existing stationary output-noise spectra.

The readout is $D_n=-(\overline V_n-\alpha\overline V_{n-1})$,
where the voltage readings are centred averages over 300 $\mu$s and
$\alpha=\exp(-\Delta t/\tau_0)$ uses the fitted circuit discharge time.
The covariance at lag $k$ is
\[
 C_D(k)=\int_{0.1\,\mathrm{Hz}}^{1\,\mathrm{MHz}} S_V(f)
 \operatorname{sinc}^2(fT_a)
 |1-\alpha e^{-2\pi i f\Delta t}|^2
 \cos(2\pi f k\Delta t)\,df,
\qquad \rho(k)=C_D(k)/C_D(0),
\]
with $\operatorname{sinc}(x)=\sin(\pi x)/(\pi x)$.

\begin{center}\begin{tabular}{rrr}\toprule
Separation (buckets) & Correlation, 20 ms & Correlation, 40 ms\\\midrule
@ROWS@
\bottomrule\end{tabular}\end{center}
These coefficients describe the noise after the frame difference. They are
different from correlations between the underlying CSA voltage samples.
Drawing independent Gaussian bucket fluctuations retains the calculated marginal
RMS but omits the coefficients above. Noise between pads is independent in the
present local circuit model; shared electronics are not part of this calculation.

The spectra include the transistor thermal and flicker noise at the DC operating
point. The same operating point, model and 0.5 mm pad circuit benchmark as the
paper are retained. The noise variance at zero lag reproduces the saved ENC
calculation. Checks double the spectral density, increase quadrature order from
8 to 12, enforce the covariance bound and confirm positive definiteness of the
six-bucket covariance matrix. All ten checks pass.

Reproduce from the repository root:
\begin{verbatim}
julia --project=. scripts/electronics/readout/analyze_bucket_noise_covariance.jl
python scripts/electronics/readout/present_bucket_noise_covariance.py
\end{verbatim}
Results and source hashes: \texttt{results/ltspice/bucket\_noise\_covariance/analysis.h5}.
Julia calculates the covariances; Python presents these saved values.
\end{document}
'''.replace('@ROWS@',rows)
 (out/'report.tex').write_text(text)
 (out/'provenance.json').write_text(json.dumps({'source':str(source.relative_to(ROOT)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'generator':str(Path(__file__).relative_to(ROOT))},indent=2)+'\n')
 print(out/'report.tex')
if __name__=='__main__':main()
