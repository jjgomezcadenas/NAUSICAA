#!/usr/bin/env python3
"""Generate the ADC range/resolution estimate from saved Julia results."""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import hashlib,json
import h5py
ROOT=Path(__file__).resolve().parents[3]
def main():
 source=ROOT/'results/ltspice/adc_estimate/analysis.h5';out=ROOT/'latex_generated/adc_estimate';out.mkdir(parents=True,exist_ok=True)
 with h5py.File(source) as h:
  a=dict(h.attrs);p=h['populations'];rows=[];adcrows=[]
  for label,n,v,snr in zip(p['labels'].asstr()[:],p['ions_per_frame'][:],p['estimated_frame_signal_V'][:],p['analog_snr'][:]):
   rows.append(f'{label.capitalize()} & {n:g} & {v*1e3:.2f} & {snr:.2f}'+r' \\')
  for key in ['span_1V','span_2V']:
   g=h[key]
   for j,bits in enumerate(g['bits'][:]):
    adcrows.append(f"{g.attrs['span_V']:g} & {bits} & {g['LSB_V'][j]*1e3:.4f} & {g['ADC_ENC_e_rms'][j]:.2f} & {g['total_ENC_e_rms'][j]:.2f} & {g['total_ENC_relative_increase'][j]*100:.3f}"+r' \\')
  text=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}
\usepackage[margin=22mm]{geometry}\usepackage{amsmath,booktabs,hyperref}
\title{Preliminary ADC range and resolution for NAUSICA}
\author{NAUSICA calculation record}\date{\today}
\begin{document}\maketitle
A nominal 12-bit ADC over a 1--2 V input span gives a small quantization
contribution relative to the present CSA noise. At 2 V, the estimated ADC
contribution is @ADCENC@ electrons RMS and increases the total ENC by
@INCREASE@\%. The voltage span must accommodate the CSA baseline and
charge retained from previous frames.

\section*{Inputs and frame signals}
The provisional populations supplied for this estimate are 50, 300 and
2000 ions per pad in a 10 ms frame. Constant-rate scaling to 20 ms gives
the populations below. They represent illustrative minimum, average and
maximum cases; they are inputs to this estimate rather than measured
distribution bounds. Track calculations will establish their distribution
at the selected frame duration.

The rounded paper circuit uses M4 and M3 gate voltages of 0.5568 and
0.5385 V, external input capacitance 1.25 pF, feedback capacitance 3 fF,
and a temperature of $20\,^{\circ}$C. The saved calculation gives
$g=@GAIN@\,\mu$V/ion and $\sigma_{D,\rm CSA}=@NOISE@$ mV RMS for
20 ms frames and 300 $\mu$s averaging. Its ENC is @ENC@ electrons RMS.
The linear estimate of the corrected frame signal is $D=Ng$.
The gain source uses a central square pad of 0.5 mm side in an infinite
grounded plane, homogeneous gas and weighting potential zero at infinity.
Ions approach through the final 5 mm at 3 m/s, with finite-path charge retained.

\begin{center}\begin{tabular}{lrrr}\toprule
Case & Ions per 20 ms & $D$ (mV) & Analogue $S/N$\\\midrule
@ROWS@
\bottomrule\end{tabular}\end{center}
The 4000-ion value is a linear extrapolation of the saved response;
linearity at this population has not been simulated in this estimate.
The low-population example remains below unit analogue $S/N$.
ADC resolution preserves the available analogue information; additional
voltage gain does not increase the CSA's signal-to-noise ratio.

\section*{Quantization budget}
For full-scale span $V_{\rm FS}$ and resolution $b$, the quantization step is
$q=V_{\rm FS}/2^b$. Independent uniform quantization errors on the two
readings give
\begin{equation}
 \sigma_{D,\rm ADC}=\frac{q}{\sqrt{12}}\sqrt{1+\alpha^2},
 \qquad \alpha=\exp(-\Delta t/\tau),\qquad \tau=@TAU@\,\mathrm{s}.
\end{equation}
This estimate assigns no additional noise reduction to averaging the ADC
samples. Uniform independent quantization is a statistical model; correlated
errors and inadequate dither can change the result. The combined noise is
\begin{equation}
 \mathrm{ENC}_{\rm total}=
 \frac{\sqrt{\sigma_{D,\rm CSA}^2+\sigma_{D,\rm ADC}^2}}{g}.
\end{equation}

\begin{center}\begin{tabular}{rrrrrr}\toprule
Span (V) & Bits & LSB (mV) & ADC ENC ($e$) & Total ENC ($e$) & Increase (\%)\\\midrule
@ADCROWS@
\bottomrule\end{tabular}\end{center}

Requiring the ADC noise RMS to be below 10\% of the CSA noise RMS
limits its variance to 1\% of the CSA variance and increases total ENC
by less than approximately 0.5\%. This budget requires
@BITS1@ effective quantization bits at 1 V or @BITS2@ at 2 V, corresponding
to @INT1@ and @INT2@ integer bits. A nominal 12-bit choice therefore leaves
room for a 2 V span. Actual ADC effective resolution and driver noise
should be specified at the intended sampling rate and included in the
same noise budget.

\section*{Voltage range}
The ADC samples the absolute CSA output, whose baseline is @BASE@ V.
That voltage includes charge retained from earlier frames. The @MAXSIGNAL@ mV
maximum illustrative frame signal does not determine the required full
input span. Range selection needs the total charge accumulated on each pad
over the event, the discharge history, and the usable linear output swing.
A 1 V input span can cover the present nominal baseline; a 2 V span can
cover the 1.8 V supply range. Buffer common-mode limits and the ADC input
range must be matched. The full dynamic-range requirement remains tied to
the cumulative track signal.

\section*{Reproduction}
The user-supplied populations and ADC choices are recorded in
\texttt{config/adc\_estimate.toml}. Julia reads the saved rounded-circuit
result and writes \texttt{results/ltspice/adc\_estimate/analysis.h5}.
Python generates this report from those values. Commands are in the
result directory's README. No new SPICE transient or manuscript edit is
part of this estimate.
\end{document}
'''
  g=h['span_2V'];j=list(g['bits'][:]).index(12)
  values={'ADCENC':f"{g['ADC_ENC_e_rms'][j]:.2f}",'INCREASE':f"{g['total_ENC_relative_increase'][j]*100:.3f}",
          'GAIN':f"{a['gain_V_per_ion']*1e6:.3f}",'NOISE':f"{a['CSA_noise_rms_V']*1e3:.3f}",'ENC':f"{a['CSA_ENC_e_rms']:.2f}",
          'TAU':f"{a['tau_s']:.6f}",'BASE':f"{a['baseline_V']:.4f}",'MAXSIGNAL':f"{p['estimated_frame_signal_V'][-1]*1e3:.1f}",'ROWS':'\n'.join(rows),'ADCROWS':'\n'.join(adcrows),
          'BITS1':f"{h['span_1V'].attrs['effective_bits_for_ADC_noise_budget']:.2f}",'BITS2':f"{h['span_2V'].attrs['effective_bits_for_ADC_noise_budget']:.2f}",
          'INT1':str(h['span_1V'].attrs['integer_bits_for_ADC_noise_budget']),'INT2':str(h['span_2V'].attrs['integer_bits_for_ADC_noise_budget'])}
  for k,v in values.items():text=text.replace('@'+k+'@',v)
  (out/'report.tex').write_text(text)
 (out/'provenance.json').write_text(json.dumps({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,Path(__file__))},indent=2)+'\n')
if __name__=='__main__':main()
