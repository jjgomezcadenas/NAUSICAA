# AXEL-inspired MPPC specifications and response-model scope

Date: 2026-10-30. Status: first model approved and four-event pilot implemented;
full-sample response awaits pilot review. MPPC is Hamamatsu's name for its
SiPM devices. This file retains the requested filename `mpcc_axel.md`.

## Adopted NAUSICAA inputs

These are user design inputs, not published AXEL measurements.

| Input | Adopted value or responsibility |
| --- | --- |
| Response input | Incident VUV photons per SiPM and time bucket, separately for pure Xe and Xe–15% He |
| New photon and electronics bucket | **0.25 µs = 250 ns**; electronics cannot resolve a shorter output interval |
| Within-bucket photon timing | Uniform incident-photon arrival times within each 250 ns bucket; user-adopted first approximation, 2026-10-30 |
| Electronics | NAUSICAA ASIC; adopted bucket-level response below |
| Electronics noise | Zero-mean Gaussian, **0.5 PE RMS per 250 ns bucket**, adopted as the provisional ASIC noise level |
| Electronics dynamic range | Upper saturation at **1 mA bucket-average current**; supersedes the earlier approximately 300-PE limit |
| New photon production | Four 250 ns photon histories are available; paths and provenance are in the handoff |
| S/N assessment | Calculations with the other agent will assess the provisional noise model |
| Existing PDE baseline | 0.30, adopted for the first model |
| Recovery parameter | **73.4 ns**, adopted by the user from the AXEL reference |
| Dark-count area scaling | Rate proportional to active area at fixed temperature, overvoltage and sensor technology; adopted by the user |
| Existing optical geometry | A 5 mm diameter exit-covering sensor in the photon-production model |

The user clarified that the range refers to the SiPM's response, expressed in
photoelectrons, rather than incident photons. We will express charge in
photoelectron equivalents using the calibrated charge of one fully recovered
avalanche. Gain variation, partial recovery and correlated avalanches can make
that quantity differ from the integer number of primary photodetections.
The later ASIC specification replaces the approximately 300-PE limit with
1 mA, interpreted in the established bucket-average-current convention.
These are electronics parameters for our ASIC; AXEL supplies sensor-model
lineage and historical readout comparisons.

The four current incident-photon HDF5 files contain 0.25 µs buckets and their
production cards and provenance. Their paths are recorded in
[`sipm_electronics_handoff.md`](sipm_electronics_handoff.md). Older 2 µs files
remain historical references. The parent-electron-arrival timing convention
remains an explicit upstream approximation.

## Approved first implementation

The complete model and source discussion are in
[`sipm_electronics_model.tex`](sipm_electronics_model.tex), with the generated
[model and pilot PDF](../latex_generated/sipm_response/sipm_electronics_model.pdf).
The configurable defaults are in
[`sipm_v1.toml`](../config/specs/response/sipm/sipm_v1.toml) and
[`readout_v1.toml`](../config/specs/response/sipm_electronics/readout_v1.toml).
The effective 5 mm sensor has 7854 microcells, gain 2.55 million and a
73.4 ns exponential recovery time. PDE stays constant during recovery.
Primary dark pulses follow an area-scaled Poisson process, using 1 MHz per
9 mm² at the reference conditions. Crosstalk, afterpulsing and gain spread
are deferred by agreement.

The pilot acquires all channels illuminated somewhere in each event, across
the full photon-bucket span plus four quiet buckets at each end. Microcell
initial ages follow stationary dark noise. Additional dark-only channels and
the empty-event acquisition interval are TOML-configurable. Avalanche charge
is assigned to its trigger bucket; the ASIC integrates that charge, adds
independent 0.5 PE RMS Gaussian noise, then clips average current above 1 mA.

The saved SiPM stage contains detected PE, dark pulses and optical/dark/total
charge. The independently saved ASIC stage adds noise, unclipped/clipped
current and saturation flags. See [`sipm_response.md`](sipm_response.md) for
schemas, reproducibility commands and the pilot boundary. The research scope
below records the broader effects considered; the approved first model above
sets the current implementation.

## Published sensor information

Primary source: Ardavan Ghassemi, *Latest MPPC development from Hamamatsu*,
Hamamatsu Photonics USA, 19 November 2024,
[local presentation](../papers/ELCC/axel-mppc.pdf), slides 14–18.
Bibliography key: `Ghassemi2024MPPC` in
[`references.bib`](../notes/physics/references.bib).

| Property | Source statement and scope |
| --- | --- |
| AXEL custom sensor | VUV-sensitive chip based on the S13370 series; slide 18 |
| Microcell pitch | 50 µm within each sensor; slide 18 |
| Interconnection | Through-silicon vias (TSV); slide 18 |
| Array | 64 readout channels; slide 18 |
| Packaging | Flexible polyimide substrate and PEEK stiffener; slide 18 |
| PDE curve conditions | Overvoltage 4 V, temperature −25 °C; wavelengths below 230 nm measured in vacuum; slide 18 |
| PDE near 175 nm | Approximately 24%, a graphical reading of the curve displayed on slide 18; not a tabulated calibration or a measurement of our device |
| Active area and microcell count | Not specified for the custom device in this presentation |
| Recovery time and avalanche gain | Not specified for the custom device in this presentation |

Slide 16 separately gives a VUV4 reference PDE of 30% at 175 nm and a VUV5
target of 35%, explicitly in a comparison for 75 µm microcells. These entries
do not calibrate the 50 µm AXEL device. The presentation's detailed
visible-light-device gain, dark-count and crosstalk table describes a different
device family. Its VUV dark-count development results concern cryogenic
temperatures and require temperature-specific interpretation.

The AXEL analysis paper, Akiyama et al., PTEP 2025 053H03, names the
S13370-3050CN for the measured prototype. Its specifications must be
distinguished from the custom 64-channel assembly above.
See [local paper](../papers/ELCC/elcc_analysis.pdf), sections 2.1 and 5;
bibliography key `Akiyama:2025dia`.

The assumed 5 mm optical diameter in our photon files is a detector design
input. A physical sensor choice that changes collection area or alignment
needs coordination with photon production.

## Web specification search: standard S13370-3050CN

Search date: 2026-10-30. The search index exposes the manufacturer's
2017-03-31 document `KSX-I50133-E`, *VUV-MPPC (VUV4)*, including its
specification table. Its [distributor-hosted PDF URL](https://hamamatsu.su/files/uploads/pdf/3_mppc/s13370_vuv4-mppc_b_%281%29.pdf)
currently redirects to a product page and that page's download link returns
404. The following values are transcribed from the indexed manufacturer
document, pending inspection of a recoverable original PDF. Bibliography key:
`Hamamatsu2017VUV4Indexed`. They are reference specifications for the standard
device, not newly adopted values for the custom AXEL sensor.

The table states **25 °C and 4.0 V overvoltage**, unless otherwise noted.

| Standard S13370-3050CN parameter | Indexed manufacturer value |
| --- | --- |
| Active area | 3 × 3 mm² |
| Microcell pitch / count | 50 µm / 3600 |
| Geometrical fill factor | 60% |
| Dark-count rate, typical | 1.0 Mcps per device |
| Dark-count rate, maximum | 3.0 Mcps per device |
| Crosstalk probability | 3% |
| Gain | 2.55 × 10⁶ |
| Terminal capacitance | 320 pF |
| Breakdown voltage | 53 ± 5 V |
| Recommended overvoltage | 4 V |
| Operating-voltage temperature coefficient near room temperature | 54 mV/°C |
| Peak PDE | 35% at the 500 nm peak; distinct from VUV PDE |

The table states that PDE excludes crosstalk and afterpulsing. No usable
room-temperature afterpulse probability or delay distribution for this exact
device was found in this search. The recovered table text also does not define
enough of the dark-count measurement procedure to identify it unambiguously
with the independent primary-dark-avalanche rate; that interpretation must be
explicit before adding separately sampled afterpulses.

Primary VUV4 characterization papers located include
[Arneodo et al.](https://arxiv.org/abs/1707.08004), on S13370-3050CN at about
175 K, and [Gallina et al.](https://arxiv.org/abs/1903.03663), on S13370-6152
at 163–233 K. These are useful characterization references, but their
cryogenic correlated-noise measurements do not fix our room-temperature
parameters. The operating temperature and overvoltage remain to be selected.
The user adopts area-proportional dark-count scaling for the custom sensor at
the same temperature and overvoltage. For the assumed 5 mm diameter active
disk, use `R_dark = R_reference * (pi * (2.5 mm)^2) / (9 mm^2)`, with the
standard 3 × 3 mm² device as reference. This uses geometrical active area;
the microcell fill factor is already embodied in the reference rate and must
not be applied a second time. The reference rate and operating point remain
separate from this adopted scaling law.

## Adopted provisional ASIC response

For bucket duration `dt = 250e-9 s`, let `Q_sensor_C` include all sensor
avalanches, including dark and correlated avalanches when enabled. Set
`I_sensor_A = Q_sensor_C/dt`, with `q0_C = e_C * gain` defining one fully
recovered PE equivalent. The provisional electronics model is

```text
epsilon_PE ~ Normal(0, 0.5^2)
I_unclipped_A = I_sensor_A + epsilon_PE * q0_C/dt
I_output_A = min(I_unclipped_A, 1e-3)
saturated = I_unclipped_A > 1e-3
```

The noise is added before the upper range limit. Preserve sensor charge,
sampled noise, unclipped current and the saturation flag separately.
This prescription specifies the per-bucket Gaussian distribution; temporal
and interchannel correlations remain an explicit implementation-review choice.
The current minimal model has a known zero baseline and upper clipping only;
an ASIC shaping function, baseline drift and ADC quantization can be added
when their specifications are supplied. The 0.5 PE term is electronics noise
in addition to explicitly simulated sensor dark pulses. It is a user-adopted
design approximation, not a decomposition of AXEL's measured baseline.

The upper current limit remains 1 mA if the sensor gain changes. Its PE
equivalent is `1e-3 * dt/q0_C`; likewise the noise in A follows the chosen
gain. The gain of 2.55 × 10⁶ used in the discussion remains a provisional
standard-device reference, separately from the adopted ASIC specifications.

## AXEL electronics-noise reference

Akiyama et al. (2025), section 4 and figures 11–12, report a mean baseline
standard deviation of 0.46 ADC counts during operation and 0.45 without high
voltage. Their sampling interval is 200 ns and one ADC count is 0.488 mV.
These are measured baseline statistics, not a specification of independent
Gaussian noise or of its decomposition into sensor, analog and ADC terms.

Yoshida et al. (2024), section 4.1, use a hit threshold of 3.5 ADC counts,
equivalent to 3.8 detected-photon units, and calibrate gains with dark pulses.
They integrate from five samples before to five samples after the threshold
crossings. Section 5.2.3 treats baseline-offset uncertainty within one ADC
count separately, with an estimated contribution to energy resolution of at
most 0.09% FWHM at 1836 keV. Such an offset persists over the integrated
waveform and differs from independently fluctuating sample noise.

The relevant AxFEB electronics reference is Nakamura et al.,
[arXiv:2001.02104](https://arxiv.org/abs/2001.02104), bibliography key
`Nakamura2020AxFEB`. Its section 2.3 describes a 17 ns first-order stage and
an approximately 400 ns second-order Sallen–Key shaping stage for the 5 MS/s
physics readout. The 40 MS/s auxiliary readout measures single-PE dark pulses
for gain calibration. The local `papers/axel_electronics.pdf` is an earlier
2017 design and has different sampling and filtering parameters.

The 2024 ADC-to-PE conversion gives an order-of-magnitude reference when
combined with the 2025 baseline RMS, but it is not a same-run calibration.
These measurements provide context for the adopted ASIC approximation above.

## Broader characterization scope

The forward model will map incident light to sensor charge, then expose that
charge to the separately configured electronics. It will quantify mean
response, fluctuations, saturation, recovery and the accuracy of a calibrated
linearity correction. Physical calculations and response types belong in
`src/ion_signal`; Julia writes numerical results and provenance to HDF5.

| Component | Required behaviour | Parameters or choices still needed |
| --- | --- | --- |
| Primary detection | Conditional binomial photon detection in the fully recovered, low-occupancy limit; preserve integer primary-detection ledger | PDE versus wavelength, overvoltage and temperature; selected operating point |
| Avalanche response | Convert fired microcells into charge and a time response; preserve the integrated single-avalanche charge | Mean gain, gain dispersion, pulse shape and parameter provenance |
| Saturation and recovery | Account for finite microcells and their changing availability and charge after firing; carry state across buckets | Active microcell count, recovery law/time, illumination distribution, recovery dependence of detection probability |
| Dark counts and dark current | Generate spontaneous avalanches through the same recovery and charge model; account separately for any leakage-current contribution | Dark-count rate at operating temperature and bias; rate normalization; leakage current and baseline treatment |
| Optical crosstalk | Describe additional correlated avalanches and their contribution to mean charge and variance | Probability, multiplicity law, spatial/timing approximation |
| Afterpulsing | Describe delayed correlated avalanches, including partial recovery | Probability, delay distribution and recovery dependence |
| Linearity and calibration | Compare incident light, primary detections and output charge; test correction bias and residual variance | Calibration prescription and accepted operating range |
| Electronics boundary | Convert bucket charge to average current, add 0.5 PE RMS Gaussian noise and clip above 1 mA; record clipping | Gain calibration and noise correlations; later ASIC transfer function if needed |

Each sensor effect needs a measured, published or explicitly adopted parameter
set before a physical production run. Unknown effects remain marked pending;
an ideal reference may explicitly disable them for comparisons after review.

### Recovery-time evidence

Yoshida et al., PTEP 2024 013H01, section 4.1.1, report a mean measured
recovery time of **73.4 ns** for their 3600-microcell MPPCs. They measured
high-intensity LED response using a reference MPPC behind a 5% neutral-density
filter. Section 4.3.5 then calibrates small run-dependent shifts of the effective
recovery time. Source: <https://doi.org/10.1093/ptep/ptad146>, bibliography key
`Yoshida2024AXELResponse`.

The earlier Ban et al., PTEP 2020 033H01, section 4.1, used approximately
**120 ns**, from their linearity measurement and adopted commonly across
channels. Source: <https://doi.org/10.1093/ptep/ptaa030>, bibliography key
`Ban2020AXELPrototype`.

These establish a recovery scale of order 0.1 µs for the studied AXEL devices.
They are parameters in a calibrated nonlinearity model, not yet a measurement
of the custom sensor proposed here. The user adopts 73.4 ns. The 3600-cell
count applies to the measured prototype device.

Appendix A of Ban et al. explicitly derives the correction from exponential
gain recovery: `g(t) = g0 * (1 - exp(-t/tau))`. Thus tau is the time to recover
to `1 - 1/e` of the full gain, rather than a hard dead time. For stationary
Poisson arrivals at rate lambda per microcell, averaging this gain over the
waiting-time density `lambda * exp(-lambda*t)` gives
`mean(g)/g0 = 1/(1 + lambda*tau)`. The publisher's HTML shows a positive
exponent in the waiting-time density; normalization and the stated integral
result require the negative sign used here.

Writing `n` for the ideal detected count, `m` for observed charge in fully
recovered PE equivalents, `M` for microcell count and `dt` for sample duration,
the resulting mean response is `m = n/(1 + n*tau/(M*dt))`. AXEL applies its
inverse, `n_corrected = m/(1 - m*tau/(M*dt))`, to measured waveform samples.
The 2024 analysis uses `dt = 200 ns`, `M = 3600` and individually measured
recovery times with mean 73.4 ns; it subsequently refines their calibration
using energy peaks. This count is after photon detection, not the incident
VUV count before PDE.

The derivation assumes equal illumination of microcells and a stationary
arrival process. It models the mean loss of avalanche gain. It does not by
itself specify recovery-dependent PDE, correlated noise, fluctuations of the
response or transients between changing illumination levels. These are the
remaining choices for a stochastic forward model with uniform arrival times
conditional on each bucket's photon count.

### Timing and computational representation

The 250 ns output interval fixes what electronics records. A recovery model
may need finer internal time resolution or an analytic treatment within that
interval. Neither method supplies missing physical arrival-time information:
the within-bucket illumination law must be stated and its influence tested.
The user adopts independent uniform photon times conditional on the count in
each bucket. The physical motivation is approximately uniform EL emission
along the electron path in the 5 mm diameter channel at E/P = 3 kV/cm/bar
and P = 15 bar. This is an adopted approximation for incident light; a more
detailed calculation would also weight emission by position-dependent optical
acceptance. Recovery state continues across bucket boundaries. This timing
choice preserves each input bucket's photon count and its upstream timing
convention.

Start by comparing a calibrated effective recovery model with a finite-cell
reference on small synthetic illumination cases. Select an efficient
count-level or microcell-state representation based on its accuracy in our
intensity and timing regime. Avoid individual-photon loops for production.
Treat finite-cell saturation and the electronics' 1 mA current limit as
separate effects; preserve the sensor response before electronics clipping.

AXEL's analysis corrects recovery-driven nonlinearity and refines its effective
recovery calibration using the concentration of light across waveform samples
and channels (Akiyama et al., section 5.3.3). Our forward response and any
inverse correction will be validated separately: a correction can remove a
mean bias while leaving additional fluctuations or lost information.

### Dark response and acquisition window

Dark avalanches contribute both average current and fluctuations. A quoted
DC dark current alone does not determine the avalanche rate or noise
correlations. Establish whether a source's current includes correlated
avalanches and leakage before converting it into model parameters.

Define the acquisition window, channel set, prehistory or initial recovery
state, and post-signal tail explicitly. Include zero-photon buckets and the
chosen channels with no incident light so dark counts and electronics noise
have a defined exposure. Preserve signed event-relative bucket indices.

## Interfaces, records and validation scope

The existing `EventAnalysis.read_photon_history(group)` supplies owned
`PhotonSiPMHistory` data with signed Int64 cell indices and incident counts.
The response input boundary must also read `pair_checksum` and file provenance,
which are stored outside that returned history object. Preserve sample, event,
gas, placement, geometry, drift length and shared-pair identity.

Propose separate records for primary detections, sensor response and electronics
response. Save the count ledgers, sensor charge in C, time in s, calibrated PE
equivalents, noise contributions, clipping indicators, exact cards, upstream
content identities, seeds and code versions. Cache identities must permit
electronics changes to reuse the same sensor realization. Exact types, TOML
contracts and versioned HDF5 schemas will be presented in the HTML review guide
before implementing the interfaces.

Validation will cover:

- Binomial limits, means and variances at low occupancy; primary detections
  bounded by incident photons before correlated sensor effects.
- Single-avalanche integrated charge, gain dispersion and temporal response.
- Finite-cell occupancy under short illumination, recovery under separated
  pulses, steady illumination and convergence with internal time resolution.
- Dark-only response versus exposure and temperature/bias assumptions;
  correlated-noise means, variances and delays.
- Mean response and residual nonlinearity versus intensity and pulse shape;
  calibration bias and variance after correction.
- Electronics clipping near 1 mA, independently of sensor saturation;
  Gaussian noise moments at 0.5 PE RMS, gain-to-current conversion and S/N.
- Signed bucket alignment, zero-photon gaps, empty events, acquisition tails,
  deterministic replay across event order/workers, gas identity and HDF5 round trips.

Use controlled synthetic cases first, then a small pilot on the new 0.25 µs
photon histories. Save Julia results in `results/`; Python will generate figures,
tables and a self-contained LaTeX/PDF report under `latex_generated/` from those
results. Full-sample response awaits agreement on sensor and noise assumptions.

## Next review

Review the saved pilot's sensor losses and substantial ASIC clipping before
full-sample production. The effective sensor geometry, reference gain and
dark rate remain assumptions to replace with device characterization. Later
calibration must address the difference between stationary AXEL correction
and transient illumination, and the information lost in clipped samples.
