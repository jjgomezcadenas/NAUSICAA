# Finite-gain CSA: shared paper and Python specification

Agreed scope, 24 October 2026. This is a linear small-signal model for one pad,
not a validated X-FAB transistor model. Rogelio's supplied HTML remains unchanged
as provenance. `specs_to_spice.md` records the separate circuit handoff task.

## Circuit and conventions

- Positive current enters the input node; signal output is negative relative to
  the DC baseline. Ion trajectories and finite-gap weighting-current model are
  unchanged. The injected charge per ion is Qg, not silently rescaled to e.
- Cin is the capacitance from the input node to ground, excluding Cf. It includes
  the detector and an effective amplifier/stray contribution. A real feedback
  capacitance must not also be included in Cin. Output loading is represented
  through the effective amplifier pole in this first model.
- A(s) = A0 / (1 + s/omega_a), vo = -A(s) vi, omega_a = 2 pi f_a.
  A0 is dimensionless, f_a is the open-loop pole in Hz, not the closed-loop
  bandwidth. The approximate gain-bandwidth product is A0 f_a when A0 >> 1.
- Yf(s) = Gf + s Cf, Gf = 1/Rf. Rf is an effective linear discharge resistance;
  there is no artificial upper cap. Parameter values are labelled illustrative
  until supported by Rogelio's device model or measurements.
- H_I(s) = -A(s) / [s Cin + (1+A(s)) Yf(s)] (V/A).
- For equivalent input voltage noise defined by vo=A(s)(e_n-vi),
  H_e(s) = A(s)[s Cin+Yf(s)] / [s Cin+(1+A(s))Yf(s)] (V/V).

The independent time-domain form is

    dot(vo) = -omega_a (vo + A0 vi)
    (Cin+Cf) dot(vi) = i + Cf dot(vo) - Gf (vi-vo).

For static A0 and negligible discharge, the settled charge gain is

    vo/Q = -1 / [Cf + (Cin+Cf)/A0].

For finite Rf the DC transimpedance is -A0 Rf/(1+A0). The ideal CSA limit
requires large gain over the signal frequency range, not just large DC gain.
All source spectra, frequency limits and sampling conventions use SI units.

## Signals and frame estimator

Compute a one-ion waveform and explicit simultaneous/distributed arrival
sequences, including previous-frame arrivals. A uniform deterministic schedule
or its expected density is a controlled input, not a simulated track. Preserve
and report the difference between k collected ions and the induced input charge
integrated within a given frame: induction starts before collection.

Retain circuit state across frames. Use exact zero-order-hold state propagation
for numerical-bin average currents; refine the physical-current bins to check
convergence. Time bins and readout frames are independent choices.

Use the nominal two-sample estimator for comparison with the ideal model:

    tau = Rf Cf; alpha = exp(-Delta/tau)
    K = Cf (Delta/tau)/(1-alpha)
    Qhat_n = -K [Vbar_n - alpha Vbar_(n-1)].

Vbar is an instantaneous sample or a specified aperture average; use the same
aperture for signal and noise. The estimator is exact only for the ideal CSA,
constant current within the frame and instantaneous endpoints. Otherwise report
bias and history dependence, not automatic recovery. Calibrate template response
explicitly if converting this estimator into an ion-count sensitivity.

## Noise

Use one-sided PSDs in Hz, independent stationary sources:

    Se(f) = ew^2 + af/f      [V^2/Hz], af in V^2
    Si(f) = iw^2 + Sfb(f)   [A^2/Hz]
    Sout(f) = |H_I|^2 Si + |H_e|^2 Se.

Feedback noise may be supplied independently. The optional Sfb=4 kB T/Rf is
only a passive-resistor benchmark; an effective MOS resistance does not by itself
establish this noise level. Declare low/high frequency cutoffs and test numerical
convergence within the declared band. Low-frequency cutoff is a modelling
assumption, not a replacement for a baseline-removal transfer function.

For a rectangular aperture Ta, |P(f)|^2=sinc(pi f Ta)^2. Define

    C(lag) = integral Sout(f) |P(f)|^2 cos(2 pi f lag) df
    var(Qhat) = K^2 [(1+alpha^2) C(0) - 2 alpha C(Delta)].

Use a nonnegative spectral expression as an independent check of the covariance
subtraction. Sigma_Q/e is the estimator noise expressed in elementary-charge
units; it is not an ion detection threshold unless the template response is
calibrated. For per-ion template response g in coulombs per ion,
SNR(k)=k |g|/sigma_Q and k_min=target_SNR sigma_Q/|g|. State the template,
phase and prior history. Distinguish mean ion-count sensitivity from fluctuations
of an actual random arrival process. No transistor noise prediction is implied.

## Validation and outputs

Recover ideal and finite-gain analytic limits, charge conservation and linear
scaling. Compare independent ODE and transfer-function calculations, refine time
and frequency integration, test nonzero history and noise covariance/scaling.
Generate parameter/provenance JSON, waveform and sampled-estimator CSVs, and
charge-checked PWL files. A behavioral SPICE comparator tests this model; it does
not establish that a transistor circuit meets the required response.

## Rewrite audit and validation status

| Original HTML feature | Python replacement |
| --- | --- |
| Artificial maximum feedback resistance | Explicit positive Rf with no upper cap |
| Prescribed difference-of-exponentials pulse | Kirchhoff state equations driven by the physical current |
| Transistor controls not consistently connected to results | Explicit effective circuit parameters; no unsupported device sizing claim |
| Fixed phase-margin display | Model poles and complex response; no transistor stability claim |
| Inconsistent noise formulas and unrelated shaping interval | One-sided equivalent-input spectra propagated through the circuit and actual sampling estimator |
| Instantaneous-charge example alone | Single-ion, simultaneous, distributed, late and previous-frame cases |

`make test physics realistic-readout` passes 37 tests (15 new). They include
independent ODE/frequency response, high-gain and charge-sharing limits,
underdamped poles, retained state, physical-current refinement, white-noise
Lyapunov comparison and direct spectral-estimator integration. The parent also
checked the configured full-band noise result with an independently refined
positive spectral integral. These checks validate the stated model and code,
not an X-FAB device implementation. The supplied HTML remains unchanged.

The generated history benchmark reports both scaling of the complete two-frame
template and sensitivity to target-frame ions after subtracting the known prior
signal. The latter subtraction assumes a deterministic known history; uncertainty
from an additional noisy measurement is not included in its quoted noise.
