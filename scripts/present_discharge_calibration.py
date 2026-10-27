#!/usr/bin/env python3
"""Generate the discharge-calibration report and figure from saved Julia results."""
from pathlib import Path
import json,tomllib
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def main():
    cfg=tomllib.loads((ROOT/'config/discharge_calibration.toml').read_text())
    result=ROOT/cfg['output'];out=ROOT/'latex_generated/discharge_calibration';out.mkdir(parents=True,exist_ok=True)
    with h5py.File(result/'analysis.h5') as h:
        a=dict(h.attrs);old=dict(h['old'].attrs);new=dict(h['calibrated'].attrs)
        ts=h['tail_time_s'][:];y=h['tail_magnitude_V'][:];fit=h['tail_fit_V'][:]
        truth=h['calibrated/train_true_counts'][:];qo=h['old/train_estimated_counts'][:];qn=h['calibrated/train_estimated_counts'][:]
        counts=h['calibrated/counts'][:];snrs=h['calibrated/snr'][:]
        local=h['local_tau_s'][:];localenc=h['calibrated/local_fit_ENC_e'][:]
        localstart=h['fit_window_start_s'][:];localstop=h['fit_window_stop_s'][:]
        counttaus=h['calibrated/tail_tau_by_count_s'][:]
        transienterrors=h['calibrated/refinement_error'][:]
    with h5py.File(ROOT/cfg['historical']/'analysis.h5') as h:oldn=int(h.attrs['verified_N_snr5'])
    threshold=int(new['verified_N_snr5'])
    fig,axes=plt.subplots(1,2,figsize=(9,3.6),layout='constrained')
    ax=axes[0];ax.plot(ts[::20]*1000,y[::20]*1000,'o',ms=3,color='#285f8d',label='LTspice')
    ax.plot(ts*1000,fit*1000,color='#bb5a28',lw=1.5,label='Exponential fit')
    ax.set(xlabel='Time since first collection (ms)',ylabel='Output signal magnitude (mV)',title=f"100 ions; fitted decay {a['tau_s']:.3f} s")
    ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False)
    ax=axes[1];x=np.arange(1,len(truth)+1)
    ax.bar(x-.18,qo-truth,width=.36,color='#285f8d',label='0.3 s correction')
    ax.bar(x+.18,qn-truth,width=.36,color='#bb5a28',label='Calibrated correction')
    ax.axhline(0,color='.4',lw=.8);ax.set(xlabel='Bucket',ylabel='Recovered minus input ions',xticks=x)
    ax.legend(frameon=False,fontsize=9);ax.spines[['top','right']].set_visible(False)
    fig.savefig(out/'discharge_calibration.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(out/'discharge_calibration.png',dpi=160);plt.close(fig)
    rows=[]
    for label,key,factor,fmt in [('Correction time (s)','tau_s',1,'.4f'),(r'$\alpha$','alpha',1,'.6f'),(r'Gain ($\mu$V/ion)','gain_V_per_ion',1e6,'.3f'),('Noise in voltage difference (mV)','noise_rms_V',1e3,'.3f'),('ENC (elementary charges RMS)','ENC_e_rms',1,'.2f'),(r'Previous-packet residual (\%)','previous_frame_residual_fraction',100,'.3f'),(r'Following-packet contribution (\%)','following_frame_fraction',100,'.3f')]:
        rows.append(label+' & '+format(old[key]*factor,fmt)+' & '+format(new[key]*factor,fmt)+r' \\')
    rows.append(f'Minimum ions for $S/N>5$ & {oldn} & {threshold}'+r' \\')
    sequence='\n'.join(f'{k} & {n} & {o:.2f} & {v:.2f}'+r' \\' for k,(n,o,v) in enumerate(zip(truth,qo,qn),1))
    windows='\n'.join(f'{start*1e3:g}--{stop*1e3:g} & {tau:.6f} & {enc:.5f}'+r' \\' for start,stop,tau,enc in zip(localstart,localstop,local,localenc))
    brackets='\n'.join(f'{n} & {snrs[list(counts).index(n)]:.6f}'+r' \\' for n in (threshold-1,threshold))
    text=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=22mm]{geometry}
\usepackage{amsmath,booktabs,graphicx,siunitx,hyperref}
\title{NAUSICAA: discharge correction calibrated from the CSA}
\author{Numerical study}\date{\today}
\begin{document}\maketitle
\section{Calibration}
The final 180 nm CSA uses M3 with $W=\SI{200}{\micro\metre}$,
$L=\SI{2}{\micro\metre}$ and gate bias \SI{0.538477478921}{V};
M4 has $W=\SI{2.88}{\micro\metre}$, $L=\SI{1.44}{\micro\metre}$
and gate bias \SI{0.556803417206}{V}. The external input capacitance is
\SI{1.25}{pF}, the feedback capacitor is \SI{3}{fF}, the supply is
\SI{1.8}{V}, and the temperature is \SI{20}{\celsius}. The supplied
NOIMOD=2 model includes intrinsic MOS capacitances and transistor noise.

The input is the finite-gap current for uniform ion collection over
\SI{20}{ms}. The benchmark is a central trajectory at \SI{3}{m/s} through
a \SI{5}{mm} gap to a square pad of side \SI{0.5}{mm} in an infinite grounded
plane. The gas is a homogeneous half-space. The weighting potential is one
on the pad, zero on the surrounding plane and tends to zero at infinity.
The integral per ion remains $e[1-\phi_w(g)]$.

After collection ends, the input current is zero and the CSA signal decays.
A matched zero-input run supplies the baseline. We fit the output change
from a 100-ion packet between 21 and 120 ms with
\[
 V(t)=B\exp[-(t-t_*)/\tau_0].
\]
The fit is linear in $\log[-V(t)]$, evaluated on a uniform \SI{200}{\micro\second}
grid, with free amplitude $B$. The fit starts \SI{1}{ms} after the last ion
is collected. It gives
\[
 \tau_0=\SI{@TAU@}{s},\qquad
 \alpha=\exp[-\SI{20}{ms}/\tau_0]=@ALPHA@.
\]
The fitted time describes the decay of the complete CSA response. It includes
the transistor feedback path and the finite amplifier response.

The maximum relative fit residual is \num{@FITERR@}. A 10-ion calibration
changes the fitted time by a relative \num{@AMPERR@}; halving the numerical
step changes it by \num{@TIMEERR@}. Fits to separate parts of the tail give:
\begin{center}\begin{tabular}{rrr}\toprule
Fit interval (ms) & $\tau_0$ (s) & Recalculated ENC\\\midrule
@WINDOWS@
\bottomrule\end{tabular}\end{center}
The fitted times over the tested 10--2000 ion range lie between
\SI{@TAUMIN@}{s} and \SI{@TAUMAX@}{s}.

\section{Readout and noise}
The readings are centred averages over \SI{300}{\micro\second}. The
measurement and calibrated ion estimate are
\[
 D_n=-[\overline V_n-\alpha\overline V_{n-1}],\qquad
 \widehat N_n=D_n/g.
\]
Both the gain $g$ and the noise are recalculated with the new $\alpha$.
The gain is obtained from the uniform mean response per ion. The same saved
stationary noise spectrum of the final circuit is integrated from
\SI{0.1}{Hz} to \SI{1}{MHz} with the averaging and difference filter:
\[
 \sigma_D^2=\int S_V(f)\operatorname{sinc}^2(fT_a)
 \left[(1-\alpha)^2+4\alpha\sin^2(\pi f\Delta t)\right]\,\mathrm{d}f,
 \qquad \mathrm{ENC}=\sigma_D/g.
\]
Here $\operatorname{sinc}(x)=\sin(\pi x)/(\pi x)$.

\begin{center}\begin{tabular}{lrr}\toprule
Quantity & Previous correction & Calibrated correction\\\midrule
@ROWS@
\bottomrule\end{tabular}\end{center}
The previous-packet entry is the signal in the first empty bucket following
an isolated packet, divided by its signal in its own bucket. The following-packet
entry arises from induction before collection and averaging across the boundary.
The correction removes most of the previous-packet contribution. The remaining
boundary response is included in the multi-bucket test below.

Direct finite-charge LTspice transients bracket the new threshold:
\begin{center}\begin{tabular}{rr}\toprule
Ions & $S/N$\\\midrule
@BRACKETS@
\bottomrule\end{tabular}\end{center}
The noise is stationary about the quiescent operating point. Signal amplitudes
come from the nonlinear transistor transients.

\section{Consecutive and empty buckets}
An independent transient contains the following sequence of uniform packets.
Each estimate uses the same gain and correction calibrated above.
\begin{center}\begin{tabular}{rrrr}\toprule
Bucket & Input ions & Previous estimate & Calibrated estimate\\\midrule
@SEQUENCE@
\bottomrule\end{tabular}\end{center}
The largest absolute error changes from @OLDMAX@ to @NEWMAX@ ions over this
sequence, compared with an ENC of @ENC@ elementary charges RMS. These
estimates use the deterministic waveforms; random noise has not been added
to the displayed sequence. The noise level is calculated separately from
the output spectrum. The test establishes the accuracy over this sequence
of charge levels and arrival times.

\begin{center}\includegraphics[width=\linewidth]{discharge_calibration.pdf}\end{center}

\section{Validation and reproduction}
The input-bin width and maximum transient timestep are reduced together from
5 to \SI{2.5}{\micro\second}. Integrated input charge is checked against the
finite-path value. Noise calculations compare 200 and 400 spectral points per
decade and quadrature orders 8 and 16. The largest relative change in tested
finite-charge signals is \num{@TRANSERR@}; the noise-density check changes
RMS noise by \num{@SPECERR@}. Both threshold brackets have numerical errors
less than one tenth of their distance from $S/N=5$.

The original 10 and 20 ms studies remain unchanged. Configuration:
\texttt{config/discharge\_calibration.toml}. Inputs, raw and reduced data,
model copies, PWL sources and LTspice logs are stored under
\texttt{results/ltspice/discharge\_calibration/20ms/}.
The scripts record SHA-256 hashes for the configuration, currents, circuit,
model, noise source, simulator and numerical code. Julia performs the physical
calculation; Python acquires LTspice outputs and generates this report from HDF5.
The README gives reproduction commands.
\end{document}
'''
    subst={'TAU':f"{a['tau_s']:.6f}",'ALPHA':f"{a['alpha']:.8f}",'FITERR':f"{a['tail_fit_max_relative_error']:.3g}",'AMPERR':f"{abs(a['tau_10_vs_100_relative_change']):.3g}",'TIMEERR':f"{a['tau_timestep_relative_error']:.3g}",'WINDOWS':windows,'TAUMIN':f'{min(counttaus):.6f}','TAUMAX':f'{max(counttaus):.6f}','ROWS':'\n'.join(rows),'BRACKETS':brackets,'SEQUENCE':sequence,'OLDMAX':f'{max(abs(qo-truth)):.2f}','NEWMAX':f'{max(abs(qn-truth)):.2f}','ENC':f"{new['ENC_e_rms']:.2f}",'TRANSERR':f'{max(transienterrors):.3g}','SPECERR':f"{new['spectral_error']:.3g}"}
    for key,value in subst.items():text=text.replace('@'+key+'@',value)
    (out/'report.tex').write_text(text)
    values={'SpiceCorrectionSeconds':f"{a['tau_s']:.3f}",'SpiceGainUv':f"{new['gain_V_per_ion']*1e6:.3f}",'SpiceNoiseMv':f"{new['noise_rms_V']*1e3:.3f}",'SpiceEncIons':f"{new['ENC_e_rms']:.2f}",'SpiceMinimumIons':str(threshold),'SpiceBelowIons':str(threshold-1),'SpiceBelowSnr':f'{snrs[list(counts).index(threshold-1)]:.4f}','SpiceAboveSnr':f'{snrs[list(counts).index(threshold)]:.4f}','SpicePreviousPercent':f"{new['previous_frame_residual_fraction']*100:.3f}",'SpiceFollowingPercent':f"{new['following_frame_fraction']*100:.3f}"}
    values={k:(v.rstrip('0').rstrip('.') if '.' in v else v) for k,v in values.items()}
    (out/'values.tex').write_text('% Calibrated 20 ms result; proposed manuscript values, not yet loaded.\n'+''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k,v in values.items()))
    (out/'values.json').write_text(json.dumps(values,indent=2)+'\n')
    print(out/'report.tex')
if __name__=='__main__':main()
