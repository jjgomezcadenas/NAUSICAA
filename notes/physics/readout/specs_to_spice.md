# Specs to SPICE: one-pad CSA comparison

Task agreed with J. J. Gómez-Cadenas on 24 October 2026.

**Purpose:** give Rogelio reproducible physical input currents and a comparison
procedure, then obtain circuit-supported response and noise parameters for the
paper. A behavioral reproduction of the analytical model is a consistency
check; the transistor circuit is a separate validation step.

**Status (LTSpice branch):** Rogelio's transistor schematic, BSIM3 models and
reports are saved unchanged under `../../../LTSpice/`.
Native LTspice 17.2.4 has now run the converted circuit, the Python behavioral
comparator, uniform ion currents with 5 pF external pad capacitance, AC and noise.
The [study results](../../../results/ltspice/README.md) and
[analysis](../../../results/ltspice/analysis.json) document the calculations.
All 46 tests pass. This supersedes the earlier "not run/not supplied" status.

The model gives S/N > 5 at approximately 1,600 ions per 10 ms frame (1578 in the
adjacent-count numerical check), with the paper's current estimator. The supplied
model produces zero flicker-noise contribution. M4 and M6 dominate the included
noise. These results have not been introduced into the manuscript. Foundry/model
provenance and calibrated flicker parameters still need confirmation from Rogelio;
nothing has been sent to him.

## 1. Package and fixed inputs

The model is defined in [realistic_csa_spec.md](realistic_csa_spec.md), with
reproducible inputs in `../../config/realistic_csa.toml`. Use generated metadata
as the numerical authority; do not copy numbers from the figure labels.

The active paper now fixes Cin = 5 pF from the Topmetal-S electrode benchmark
(Mei et al., arXiv:2010.09226, Section 2.3). Its main signal is the uniform
mean current per ion, scaled by a free N, with point samples at frame boundaries.
`make noise-study` generates the paper signals and threshold (68 ions for
S/N > 5 under the illustrative spectra). `make realistic-readout` retains the
explicit diagnostic pulses and its configured 10 microsecond aperture. Apply
the paper's zero aperture when comparing its threshold; do not mix these two
sampling conventions. Use `uniform_per_ion_waveform.csv` for the normalized
paper waveform and scale its current by the selected N for a transient check.

For each benchmark, supply:

- The PWL file: time in seconds, current in amperes, positive into the CSA input.
- Physical time offset: SPICE starts at zero, but frame/arrival times refer to
  the original physical clock. Apply the recorded offset to every sample time.
- Collection times and ion count, detector configuration, and total induced
  input charge. Every complete finite-gap trajectory contributes Qg, not exactly
  e. Check the PWL integral without renormalizing it.
- The ideal and finite-gain voltage predictions and the actual induced charge
  per frame. State that the two-sample recovered charge can be biased.
- Cf, input-to-ground Cin excluding Cf, Rf, A0 and amplifier pole. These are
  effective small-signal comparator parameters, not prescribed transistor sizes.
- Frame interval, phase, acquisition-aperture convention and duration, initial
  state, noise frequency band, noise coefficients and desired S/N criterion.

The minimum input set is one ion; k simultaneous arrivals; k arrivals distributed
over a frame; and a packet with preceding-frame history. Use multiple k values
to test scaling. Distributed arrivals are diagnostics, not a Geant4 event.

The generated package is under `../../results/realistic_csa/`:

| File | Contents |
| --- | --- |
| `benchmarks.json` | Circuit, noise and sampling assumptions; charge checks and comparison metrics |
| `<case>_arrivals.csv` | Physical ion collection times in seconds |
| `<case>_waveform.csv` | Physical time, current, input-node voltage, output voltage and ideal output |
| `<case>_samples.csv` | Sampled voltage, true frame charge and estimated charge |
| `noise_spectrum.csv` | One-sided output spectra before/after aperture filtering and response magnitudes |
| `spice/<case>.pwl` | SPICE time and current; physical-time offset is in the JSON |
| `spice/<case>.cir` | Behavioral comparator using the configured gain, pole and feedback |

Run `make PYTHON=/Users/jjgomezcadenas/miniconda3/bin/python3 realistic-readout`
from the repository root to reproduce the package. The generated comparator
uses LTspice-style `PWL FILE` syntax; this syntax must be checked in the selected
simulator. It has not been run here.

## 2. Circuit connection and operating point

Connect the external current source from the reference node into the detector
input. The simplified comparator is linear about zero signal baseline. The
transistor circuit must instead start at its physical bias point; record that
baseline and subtract it from the output before comparison.

For the transistor netlist, add only the detector/external capacitance that is
not already represented by transistor and wiring capacitances. The analytical
Cin includes those internal contributions: adding its entire value to the
transistor input would double-count them. Determine the effective input
capacitance and feedback capacitance at the chosen operating point, and report
how they were extracted. Avoid double-counting output loading when matching the
amplifier pole.

Extract the actual DC gain and frequency response, and the local discharge
behaviour over the signal range. Map these to A0, the effective pole and Rf for
a first comparison. If a single pole or linear discharge law does not describe
the circuit, return the measured response rather than forcing a fit.

## 3. Transient comparison

1. Run the same PWL current through the behavioral comparator and transistor
   circuit. Retain the state between frames. Include the full induction interval
   before the earliest collection time and the specified post-arrival window.
2. Record input-node voltage and output voltage relative to baseline. Use enough
   waveform precision to resolve the single-ion signal.
3. Halve the transient maximum timestep and tighten solver tolerances until the
   peak voltage and recovered frame charge are stable to the stated comparison
   tolerance. Refinement of the PWL current grid is a separate check. The numerical
   integration step is not the readout frame interval. Resolve the circuit's fast
   response if it affects the quantities being compared.
4. Apply the same sampling aperture and two-sample estimator offline to both
   waveforms. Compare frame charge, fractional bias, peak voltage, peak time and
   residual voltage after the packet. For history cases, distinguish raw recovered
   charge from subtraction of a known previous-frame template; a noisy measured
   history subtraction would require additional covariance terms. Report absolute error for empty frames,
   where fractional charge error is undefined.
5. Repeat with increasing k and with a previous-frame packet. Report the onset
   of amplitude-dependent gain, discharge or clipping. A small single-ion result
   alone does not establish the usable dynamic range.

## 4. Noise comparison

Obtain the output noise PSD about the same DC operating point, in V²/Hz with a
one-sided Hz convention. If SPICE supplies V/sqrt(Hz), square it before
integration. Use the same external capacitance and feedback network as for the
transient calculation. Request separate source contributions where supported.
Do not add a passive-resistor Johnson source on top of transistor feedback noise
already represented by the device model.

Integrate the output spectrum using the specified acquisition aperture and
frame-difference filter. Compare sample variance, adjacent-sample covariance,
and recovered-charge RMS, not just a broadband output RMS or an unrelated ENC.
Report the frequency bounds and frequency-grid convergence. Low-frequency
cutoffs must match; a proposed baseline-removal method needs its own filter.

The stationary operating-point calculation is the initial comparison. If the
circuit is periodically reset/switched, or its noise changes substantially with
signal level, a stationary small-signal noise analysis will not describe all of
its behaviour. Record that difference for a subsequent circuit calculation.

## 5. Requested return from Rogelio

- Simulator/version, transistor netlist, device-model identity/corner,
  temperature, biases, supply and DC operating point.
- Detector/load capacitances and extracted effective gain, response and discharge.
- Per-case transient CSVs with units and timing, plus timestep/tolerance checks.
- Output noise spectra and integration bounds; distinguish transistor feedback
  noise from the optional passive-resistor analytical benchmark.
- Output range and amplitude-dependent behaviour, bias current/power, and any
  additional poles relevant to stability. These quantities must be calculated,
  not copied from the original HTML's fixed labels.

## Completion criteria

- [x] Numerical handoff generated and its current integrals verified.
- [x] Behavioral comparator checked against the Python waveform (state clearly
      whether a SPICE executable was actually run).
- [x] Supplied transistor circuit and operating conditions identified; foundry provenance remains unconfirmed.
- [x] Matched uniform-arrival transient and sampled-charge comparison completed.
- [x] Matched noise/covariance comparison completed for the supplied model defaults (zero flicker contribution).
- [ ] Circuit-supported parameters transferred to the paper and track model.

Sending this package is a separate action and requires the user's instruction.
