# SiPM and ASIC response

The user has authorized co-worker audit and full-sample production. See
[`sipm_production_handoff.md`](sipm_production_handoff.md) for the handoff and
production-runner adaptations; the commands below still reproduce the pilot.

The approved first model is implemented in Julia. The physical note is
[`sipm_electronics_model.tex`](sipm_electronics_model.tex); the compiled
[model and pilot report](../latex_generated/electron_readout/sipm_response/sipm_electronics_model.pdf)
includes parameter provenance, assumptions, signal recovery and validation.
The [short pilot report](../latex_generated/electron_readout/sipm_response/report.pdf) is generated
entirely from saved numerical results.

## Configuration and commands

- Sensor: `config/specs/response/sipm/sipm_v1.toml`.
- ASIC: `config/specs/response/sipm_electronics/readout_v1.toml`.
- Pilot: `config/specs/runs/sipm_response_pilot.toml`, event 0 in each of four
  existing 250 ns incident-photon files. The pilot runner permits at most
  16 explicit event IDs. Full-sample production is authorized after audit and
  requires the runner adaptations described in the handoff.

From the repository root:

```sh
make sipm-pilot PYTHON=python3.13
```

This runs the Julia tests and pilot, generates tables from HDF5, and compiles
both reports. `make sipm-report PYTHON=python3.13` regenerates the presentation
from saved results. `make sipm-test` runs the focused physical/interface tests.
The standard repository checks are `make test physics PYTHON=python3.13` and
`make note PYTHON=python3.13`.

## Interfaces

`SiPMResponse.Sensor` owns scalar physical parameters. `simulate_channel`
consumes successive integer photon counts, bucket duration and a caller-owned
RNG, returning owned detected-PE, dark-pulse and charge vectors. It preserves
the input, samples binomial PDE once per bucket and loops over detected
avalanches. Their uniform times and uniform microcell indices are random;
the recovered charge conditional on each interval is deterministic. Cell
state carries across buckets. Local Float64 times and charge sums retain
precision independently of the source's signed bucket origin.

`SiPMElectronics.Electronics` holds duration, Gaussian noise RMS in PE and
current limit in A. `apply_electronics` consumes saved sensor charge and the
charge calibration of a recovered avalanche. It returns independent Gaussian
noise, current before/after clipping and flags, preserving negative noise
estimates and the input charge. The limit applies to average current over
the bucket. The first model assigns all avalanche charge to its trigger
bucket; a future peak-current or shaped-pulse model needs an ASIC pulse shape.

`scripts/production/sipm_pipeline.jl` resolves cards, uses the established
`EventAnalysis.read_photon_history` interface, adds zero-photon buckets and
persists each stage. The acquisition includes channels illuminated anywhere
in the event, across its global photon span with four quiet buckets at each
end. Additional channel indices and an empty-event interval are configurable.
There is no charge threshold: an isolated distant photon adds its receiving
SiPM, and intervening spatial channels are not filled. The user has explicitly
confirmed this selection policy.
This channel selection defines the dark/noise exposure; it is not a
detector-wide trigger/noise simulation. Each microcell starts with an
independent stationary dark-process age; at zero dark rate it starts recovered.

## HDF5 contract and identity

Each sample/gas has two separate files under its existing realisation tree:

```text
response/electrons/<photon_output_id>/<sensor_id>/sipm_response.h5
response/electrons/<photon_output_id>/<electronics_id>/electron_response.h5
```

The schemas are `nausicaa.sipm_file.v1` and
`nausicaa.electron_response_file.v1`. Each contains `metadata` and `events`.
Metadata records the exact card, canonical card hash, code hash and source
snapshot, Julia version, upstream file SHA-256, event selection, gas/sample,
bucket duration and charge per PE. Sensor metadata also embeds the upstream
photon metadata. Files are written through temporary paths; completion is
recorded after event validation, followed by read-back checks and atomic rename.

Event schemas are `nausicaa.sipm_history.v1` and
`nausicaa.electron_response_history.v1`. Rows are sorted by SiPM then bucket.
All selected channels have every acquisition bucket, including zero-photon
rows. Array columns have equal lengths:

| Dataset | Meaning/unit |
| --- | --- |
| `pixel_ix`, `pixel_iy`, `bucket` | Signed Int64 physical SiPM and bucket indices |
| `incident_photons` | Int64 photons from the detector stage |
| `detected_pe` | Int64 primary optical detections before charge loss |
| `dark_pulses` | Int64 primary dark pulses |
| `optical_charge_C`, `dark_charge_C`, `charge_C` | Float64 optical, dark and total sensor charge |
| `noise_pe` | ASIC Gaussian noise in recovered avalanche equivalents |
| `current_unclipped_A`, `current_A` | ASIC average currents before/after upper clipping |
| `saturated` | Boolean flag, unclipped current strictly above the limit |

The last three rows belong to the ASIC file. Each event also stores its ID,
gas, pair checksum, acquisition limits, channel indices and stage-specific
UInt64 channel seeds. `photon_input` retains the versioned upstream history,
including placement, drift distances, grid, primary pairs and optical counts.
Sensor charge divided by metadata `charge_per_pe_C` yields calibrated PE;
charge divided by `bucket_width_s` yields noiseless average current.

The sensor model identity uses only its card, source code and Julia environment.
Its stream identity adds photon-file contents. Per-event/channel seeds derive
from that identity and the configured seed, independently of event ordering
or pilot selection. The cache ID additionally records the selected event IDs.
The electronics identity adds its own card and code. Changing the ASIC card
reuses exactly the saved sensor file; cached electronics verify its SHA-256.
Individual avalanche times/microcell indices are reproducible from seeds and
the preserved environment rather than stored as large event tables.

`results/sipm_response_pilot/latest.toml` locates the pilot manifest and numerical
summary. `results/sipm_response_validation/validation.h5` contains independent
numerical checks. Python only reads these saved values for tables/LaTeX.

## Recovery and pilot interpretation

`axel_mean` and `axel_inverse` are stationary uniform-illumination diagnostics.
The inverse returns a validity flag and rejects clipping or a nonpositive
denominator. It retains signed finite noise estimates. The pilot saves the
forward response; event-energy reconstruction and calibration remain later
analysis. Mean dark subtraction must account for the nonlinear response and
illumination-dependent recovery. Clipped buckets require a separate treatment
and cannot be restored by the scalar AXEL inverse.

The first four events show substantial clipping at the adopted 1 mA limit.
Their saved results demonstrate the need for a range review; they do not
establish an ensemble saturation fraction. The report compares the stationary
mean with independent Poisson illumination and a transient first bucket,
exposing finite-bin inverse bias as well as variance.
