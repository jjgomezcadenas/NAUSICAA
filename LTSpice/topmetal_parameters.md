# Topmetal ENC benchmark and NAUSICA shaping study

Agreed scope, 24 October 2026. Branch: `LTSpice`.
Step 1 is complete. A first shaping study has been authorized, starting at 180 microseconds. This note does not change the
reviewed manuscript or its illustrative numerical results.

## Aim

Determine whether Rogelio's circuit can reach the measured Topmetal-S ENC,
then determine the noise and charge recovery for ions arriving over 10 ms.
Use 120–140 electrons RMS as the published benchmark, and approximately
200 electrons RMS as the engineering target discussed with the author.
Do not impose either number on the simulation.

The existing result is approximately 315 electrons RMS referred to the input
for our uniform-arrival, two-point frame estimator, with 5 pF external pad
capacitance. It is specific to that signal and processing method. See
[the existing LTspice study](../results/ltspice/README.md).

## Work plan

1. **Establish the measurement benchmark — complete.** Read the corrected
   Topmetal-S measurement and supporting primary sources. Record injection,
   calibration, capacitances, feedback, sampling and filtering. Distinguish
   reported settings from missing information and earlier chip generations.
2. **Repeat a comparable charge-injection measurement — pending.** Keep
   Rogelio's circuit unchanged initially. Inject calibrated short pulses, retain
   the 5 pF external pad benchmark and apply the same defined trapezoidal filter
   to the signal and output-noise spectrum. Scan shaping duration and flat-top
   duration; check pulse-amplitude linearity and the finite CSA rise time.
   Repeat around the published 200 ms feedback decay if the circuit can be set
   there, extracting its actual decay rather than imposing an ideal RfCf.
   Report both the as-supplied bias point and any deliberate feedback adjustment.
   Unreported Topmetal settings are scan parameters, not claimed replications.
3. **Identify the limiting circuit elements — pending.** If the shaped ENC
   remains above the benchmark, examine transistor noise contributions, gain,
   capacitances and bias. Confirm the intended M5 dimensions, foundry/model
   provenance and flicker-noise parameters. Change one circuit parameter at a
   time, retaining the original comparator and accounting for power/bandwidth.
4. **Return to the ion signal — pending.** Apply suitable processing to uniform
   mean arrivals over 10 ms and to neighbouring frames. Recalibrate the filtered
   signal; calculate ENC, charge recovery and the minimum N for S/N > 5. Measure
   charge mixing between frames and the effect on longitudinal resolution.

All calculations belong in tracked Python scripts and reproducible SPICE decks.
Check calibration, charge conservation, filter normalization, transient and
spectral convergence, and agreement between time- and frequency-domain noise
calculations. Report ENC versus shaping time and compare short charge injection,
uniform ion arrivals and the existing two-sample estimator. Do not transfer
results into the paper before review. Nothing is to be sent to Rogelio without
an instruction to send it.

## Step 1: corrected benchmark

Primary source: Tianyu Liang et al., *Performance of a novel charge sensor on
the ion detection for the development of a high-pressure avalancheless ion TPC*,
JINST 19 (2024) C04004, Sections 2.1–2.2, Figure 2 and Section 4.
[DOI](https://doi.org/10.1088/1748-0221/19/04/C04004).
Bibliography key: `Liang2024TopmetalS`.
Local source: `/Users/jjgomezcadenas/Dropbox-Personal/2026/ITACA/Liang_2024_J._Inst._19_C04004.pdf`.

| Quantity | Published information | Consequence for our comparison |
| --- | --- | --- |
| Chip | Second-generation Topmetal-S, 130 nm CMOS; one relatively low-noise variant of four | Distinct from Rogelio's 180 nm circuit and the older 350 nm Topmetal-S |
| Electrode | Hexagonal, 1 mm diameter, surrounded by a guard ring | Relevant electrode scale; exact input capacitance of this measured variant is not tabulated |
| Charge injection | Pulse on guard ring, coupled capacitively to collection electrode | Calibrate injected charge independently of CSA voltage gain |
| Injection capacitance | 255 ± 10 fF, calibrated against external 200, 300 and 500 fF capacitors | This is the guard-ring injection capacitance, NOT the total input/pad capacitance |
| Main noise measurement | Approximately 16,000 electrons injected, repeated 100 times | Figure 2(c) is a distribution of filtered pulse amplitudes |
| Processing | Trapezoidal filter | Filtering is explicit, but numerical filter parameters are not supplied |
| Main feedback decay | Approximately 200 ms | Associated with the Figure 2(c) noise result |
| Main ENC | 121 ± 9 electrons RMS | Primary calibrated benchmark |
| Feedback dependence | Approximately 140 electrons at 10 ms decay, versus 121 at 200 ms | Both settings are useful comparisons once actual circuit decay is extracted |
| Injection response | Approximately 3 ms CSA waveform rise time | This is a signal rise time, NOT a quoted trapezoidal shaping time |
| Charge dependence | ENC reported independent of injected charge up to approximately 40,000 electrons | Check our signal linearity; do not infer a pulse-size-independent ENC without checking it |
| Adjustable decay | Approximately 3.5 ms to 1.5 s | Broader characterization range, not all settings used for the 121-electron result |

The paper's 320 mV/fC gain, 50-electron ENC, 49.5 mW power and 35,000-electron
linear range are explicitly **design specifications**. They should not replace
measured values or be treated as the measured operating point. The exact feedback
capacitor and six bias voltages of the selected variant are not listed.

For an injected charge Qinj and filtered amplitude A, the comparison uses

    ENC [electrons RMS] = (Qinj / e) × sigma(A) / mean(A).

The signal and noise must undergo the same filter, with a defined amplitude
measurement. For small-signal simulation, the equivalent form is filtered RMS
noise divided by the filtered response per injected elementary charge.

## Missing settings for an exact reproduction

The 2024 paper does not specify:

- trapezoid rise/averaging duration and flat-top duration;
- numerical filter implementation, pole-zero correction and baseline treatment;
- amplitude extraction rule (fixed sample, plateau average or peak search);
- actual sampling frequency, record length and analog bandwidth in the ENC run;
- injection pulse edge/width and repetition period;
- exact input and feedback capacitances and bias voltages of the tested variant;
- noise-run temperature and whether external calibration capacitors remained
  connected during that run.

The R&S RTM3004 and approximately 297 K are stated for the **ion-detection
setup**, not explicitly as settings for the Figure 2(c) ENC run. The 2022
prototype poster describes a board capable of up to 125 MHz sampling; this
is a capability, not the actual rate used for the corrected ENC measurement.

These details would permit a direct replication. Until available, step 2 can
scan explicit filter settings and compare the achieved ENC with the published
range. It must be called a comparable benchmark, not an exact reproduction.

## Earlier Topmetal-S result: keep generations separate

Chaosong Gao et al., *A Low-Noise Charge-Sensitive Amplifier for Gainless Charge
Readout in High-Pressure Gas TPC*, PoS(TWEPP2018)083 (2019), Sections 3–5.
[Primary proceedings PDF](https://pos.sissa.it/343/083/pdf).
Bibliography key: `Gao2019TopmetalS`.

This older 350 nm implementation explicitly uses approximately 5 pF detector
capacitance and approximately 3 fF feedback capacitance. Section 5 reports an
optimum shaping time of approximately **180 microseconds** for an off-chip
digital trapezoidal shaper. It also states that the injection and feedback
capacitances had not been calibrated, and that the intended open-loop gain was
not reached during measurement because of stability issues.

The 2024 paper explicitly identifies injection-capacitance calibration as the
reason earlier Topmetal-S ENC values were underestimated. Therefore the older
28.7-electron result is not our calibrated target. Nor may its 180-microsecond
setting be assigned to the 2024 121-electron measurement. It is useful as a
historical scan point only. The 5 pF load remains our explicit benchmark; the
2024 article does not independently establish its exact value for that variant.

## Supporting material checked

- [2022 prototype poster](https://indico.ihep.ac.cn/event/16065/contributions/43627/attachments/62123/71663/Poster.pdf):
  corrected guard-ring capacitance near 250 fF, approximately 120-electron ENC,
  trapezoidal processing and maximum board sampling capability. It specifies
  8 mV injection through approximately 250 fF for that poster measurement. It
  reports approximately 600-electron ENC for the first chip version and 120
  for the second; this supports the importance of circuit development rather
  than a universal ENC bound for any Topmetal-style circuit. No complete
  numerical filter prescription is given. Do not identify this poster run
  with the later 16,000-electron journal run.
- [2023 conference attachment](https://indico.cern.ch/event/1255624/contributions/5445276/attachments/2663135/4614176/Attach.pdf):
  pulse and alpha-track waveforms, with fast injection described at the ns scale;
  no trapezoid parameters. This does not establish the exact injection edge in
  the later journal noise measurement.
- L. Lang et al., *Design and Demonstration of Digital Readout Chain in NνDEx
  Experiment*, IEEE TNS 70 (2023), 2499: discusses 2 kHz sampling for the intended
  readout architecture, not the sampling of the corrected ENC measurement.
  The digital-chain demonstrator excludes the CSA and ADC. Local PDF is in the
  same ITACA Dropbox directory. Do not use 2 kHz as the measured ENC-run rate.

## Step 1 conclusion

The corrected 120–140-electron benchmark is supported by a calibrated experiment.
Its processing method is established, but its full numerical recipe is not
published in the material checked. Step 2 should use a calibrated short-pulse
response and a documented shaping scan, with feedback comparisons near 200 ms
and 10 ms where achievable. The resulting settings must then be tested with
our 10 ms ion signal and inter-frame charge history.


## First injection study

The author requested injection alone be understood before the 10 ms ion case.
`make shaping-study` now defaults to short injections only. The optional
`--include-ion-frames` is now enabled in the saved run following the author's next instruction to run the 10 ms ion signal.

See [injection results](../results/ltspice/shaping/README.md). First settings:
180 and 200 microseconds shaping rise, with flat top G=T as an explicit study
choice; flat-top alternatives G=0 and G=2T are also checked. Circuit dimensions,
biases and feedback are unchanged. The noise model still has zero flicker.
The 5 pF external pad remains in addition to intrinsic transistor capacitance.
These are model results, not a demonstrated improvement over measured Topmetal.


## First 10 ms ion study

The author has now authorized the ion case. `make ion-shaping-study` reproduces
the complete saved run. It uses short centered averages of the boundary voltages,
with the existing discharge correction and calibration from the uniform ion
signal. This is a frame-charge measurement; it does not interpret the peak of a
short-pulse trapezoid as the total charge of a long ion signal.

The generated report gives the ENC, direct neighbouring-count S/N checks at
180 and 200 microseconds, signal change relative to point samples, and the
response to packets shifted by one frame. Access to the output during the
averaging window is required; a single instantaneous voltage read per 10 ms
frame is insufficient. Flicker noise and finite ADC/readout effects remain
outside the supplied simulation. The manuscript is unchanged.


## Requested shaping scan

The current saved scan covers 180, 200, 250, 300 and 350 microseconds, with
5 pF external pad capacitance held fixed. Each ion threshold is bracketed by
direct transistor simulations at adjacent integer counts. See
[scan table](../results/ltspice/shaping/shaping_scan.csv) and
[scan plot](../results/ltspice/shaping/shaping_scan.pdf).

Next: vary the external pad capacitance while leaving transistor dimensions,
intrinsic capacitances and biases unchanged. Recompute both signal response and
noise spectrum at each capacitance; do not scale the existing ENC by capacitance.
Keep a common averaging duration for the first comparison so capacitance and
shaping changes can be distinguished. The smaller NAUSICA pad capacitance remains
a design estimate until geometry/parasitics are established.
