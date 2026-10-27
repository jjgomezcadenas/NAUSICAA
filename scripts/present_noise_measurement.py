"""Generate a standalone report from saved noise integrals."""
from pathlib import Path
import hashlib,json
import h5py

ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'results/ltspice/noise_measurement/analysis.h5'
out=ROOT/'latex_generated/noise_measurement';out.mkdir(parents=True,exist_ok=True)
with h5py.File(source) as h:
    a=dict(h.attrs);p=dict(h['paper_rounded/400'].attrs);v=dict(h['v3/400'].attrs)
    rows=[]
    for label,key in [('Unaveraged output','raw_rms_V'),('One averaged reading','average_rms_V'),('Corrected frame difference','difference_rms_V')]:
        rows.append(f'{label} & {p[key]*1e3:.6f} & {v[key]*1e3:.6f} '+r'\\')
    density=max(abs(h[name+'/200'].attrs[key]/h[name+'/400'].attrs[key]-1) for name in h for key in ['raw_rms_V','average_rms_V','difference_rms_V'])
    quadrature=max(h[name+'/400'].attrs[key] for name in h for key in ['raw_quadrature_relative_error','average_quadrature_relative_error','difference_quadrature_relative_error'])
    text=r'''\documentclass[11pt]{article}
\usepackage[a4paper,margin=25mm]{geometry}
\usepackage{amsmath,booktabs,siunitx}
\title{Output noise and the noise of the frame measurement}
\author{NAUSICAA calculation note}
\date{\today}
\begin{document}\maketitle
Rogelio's integrated output noise is reproduced with the supplied v3 circuit.
LTspice reports $\sigma_V=@@LT@@$ mV over 0.1 Hz--1 MHz.
Independent integration of the saved spectrum gives $@@RAW@@$ mV.
The LTspice noise measurement integrates the squared spectral amplitude and
returns its square root. The reported noise is therefore consistent with the
circuit used in the paper.

The paper estimates charge from the change between two averaged voltage readings,
\[
 D_n=-(\overline V_n-\alpha\overline V_{n-1}),\qquad
 \alpha=e^{-\Delta t/\tau},
\]
with $\Delta t=@@FRAME@@$ ms, an averaging window $T_a=@@WIDTH@@$ $\mu$s,
and $\tau=@@TAU@@$ s. The same coefficient is applied to both spectra below.
For the one-sided output power spectral density $S_V(f)$,
\begin{align*}
 \sigma_V^2 &=\int S_V(f)\,df,\\
 \sigma_{\overline V}^2 &=\int S_V(f)\operatorname{sinc}^2(fT_a)\,df,\\
 \sigma_D^2 &=\int S_V(f)\operatorname{sinc}^2(fT_a)
 [1+\alpha^2-2\alpha\cos(2\pi f\Delta t)]\,df,
\end{align*}
where $\operatorname{sinc}(x)=\sin(\pi x)/(\pi x)$ and all integrals run
from 0.1 Hz to 1 MHz. Each averaged reading uses a rectangular window.

\begin{center}\begin{tabular}{lrr}\toprule
Measurement & Rounded paper circuit (mV RMS) & v3 (mV RMS)\\\midrule
@@ROWS@@
\bottomrule\end{tabular}\end{center}

The averaging alone changes the noise only slightly. Most of the reduction
comes from subtracting the slowly varying noise shared by the two readings.
Their correlation coefficient is $@@CORR@@$ for v3. Equivalently,
\[
 \sigma_D^2=(1+\alpha^2)\sigma_{\overline V}^2
 -2\alpha\operatorname{Cov}(\overline V_n,\overline V_{n-1}).
\]
Treating the two samples as independent would omit this covariance.

For the rounded paper circuit the calibrated frame response is
$g=@@GAIN@@$ $\mu$V per ion, giving $\mathrm{ENC}=\sigma_D/g=@@ENC@@$
elementary charges RMS. This gain is the response of the same frame estimator
to ions collected uniformly during the frame. It includes the 300 $\mu$s windows.
Rogelio's ratio of the unaveraged output noise to the reported single-ion
peak, $@@RAW@@\,\mathrm{mV}/27.8\,\mu\mathrm V=@@PEAKENC@@$, describes
a different measurement. Both signal and noise must be evaluated with the
chosen measurement procedure. The short comparison applies the reference
filter to v3; it does not recalibrate v3's frame gain.

The stated exponential source of 0.096 fA and 1.7 ms integrates to
$@@CHARGE@@\,e$ with negligible rise time. This small charge-normalization
difference does not explain the difference between the two noise estimates.

The two spectra use $C_d=1.25$ pF, $C_f=3$ fF, 20 $^\circ$C and the supplied
180 nm transistor models, including their thermal and flicker noise.
The rounded paper biases are $V_{\rm ncasc}=0.5568$ V and
$V_{\rm nload}=0.5385$ V; v3 uses 0.556803 V and 0.53847 V.
The reference gain uses a central square pad of side 0.5 mm in an infinite
grounded plane, a 5 mm approach at 3 m/s and the finite-path induced charge.
Uniform ion collection times span 20 ms.

The frequency sampling was doubled from 200 to 400 points per decade;
the largest change in the integrated RMS values was $@@DENSITY@@$ relative.
Increasing Gaussian quadrature order from eight to twelve changed the
values by at most $@@QUAD@@$ relative. The supplied v3 calculation was run
with LTspice 17.2.4 for macOS. Numerical results and provenance are saved in
\texttt{results/ltspice/noise\_measurement/analysis.h5}; spectra and simulator
logs are in the adjacent \texttt{raw.h5}. The pipeline and parameters are
documented in its README.
\end{document}
'''
    values={'LT':f"{v['ltspice_measured_rms_V']*1e3:.4f}",'RAW':f"{v['raw_rms_V']*1e3:.6f}",'FRAME':f"{a['frame_s']*1e3:g}",'WIDTH':f"{a['average_s']*1e6:g}",'TAU':f"{a['tau_s']:.6f}",'ROWS':'\n'.join(rows),'CORR':f"{v['averaged_noise_correlation']:.5f}",'GAIN':f"{a['gain_V_per_ion']*1e6:.5f}",'ENC':f"{p['noise_over_reference_gain_e']:.3f}",'PEAKENC':f"{v['raw_over_reported_peak_e']:.1f}",'CHARGE':f"{a['reported_injected_charge_e']:.6f}",'DENSITY':r'\num{'+f'{density:.3e}'+'}','QUAD':r'\num{'+f'{quadrature:.3e}'+'}'}
    for k,val in values.items():text=text.replace('@@'+k+'@@',val)
    assert '@@' not in text
    (out/'report.tex').write_text(text)
    (out/'provenance.json').write_text(json.dumps({'analysis_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'presentation_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
