# VH1 readout: handoff for production

Date: 2026-11-14. Branch: `VH1` (everyone works here). Author of the module: the
auditing agent; physics in `notes/physics/nausicaa_VH.tex`, Sections 7, 13
and 14, and `notes/html/event_analysis/sim_vh1_readout.html`.

## What exists

| Item | Path |
|---|---|
| Module (electronics and reconstruction) | `src/vh1/VH1Readout.jl` |
| Cards | `config/specs/response/ion_csa/vh1_device_changes.toml` (default), `vh1_base.toml` |
| Tests | `tests/vh1/test_vh1.jl` (run as the other Julia tests, with the `work/julia-depot` depot) |
| SPICE inputs the cards point to | `results/ltspice/design/VH1/` (committed summaries and spectra) |
| Design page | `notes/html/event_analysis/sim_vh1_readout.html` |

The module is an independent realisation beside `ElectronicsResponse` and
`SignalReconstruction`. It imports their current builder, noise sampler, ADC
quantiser, averaged covariance and output types, and changes none of them.
Do not edit `src/electronics`, `src/reconstruction` or `src/vh1`; ask for
changes to `src/vh1` and keep the cards' hashes as they are unless the SPICE
study is rerun (see "Cards").

## Physics in one paragraph

The amplifier has no discharge element. Two 3 fF capacitors alternate between
feedback and a parking node every 20 ms frame; the output is V_REF minus the
gain per ion (46.9 uV for the default circuit) times the ions accumulated
since the reset, linear until the output leaves the 0.4 to 0.95 V window. A
sigma-delta converter integrates the output over two contiguous 10 ms windows
per frame and delivers one 16-bit value per window. Noise is the SPICE output
spectrum of the VH1 circuit, Gaussian and stationary, independent between
pads. The reconstruction solves readings = response matrix x ion counts by
covariance-weighted least squares; the matrix holds the spread of the 1.67 ms
ion pulse and of the window integration over neighbouring frames. Nothing is
nonlinear; a pad whose noiseless output leaves the window is invalid.

Noise after reconstruction, 20 ms buckets, default circuit: 31 ions per
bucket (present chain: 50). With `conversions_per_frame = 1` it is 43 at a
5 ms window and 59 at a full-frame window; the second reading per frame is
what makes the response matrix well conditioned.

## Interfaces

Load after `DetectorResponse`, `ElectronicsResponse` and `SignalReconstruction`:

```julia
include("src/vh1/VH1Readout.jl"); using .VH1Readout
s=load_vh1_settings(ROOT,joinpath(ROOT,card["vh1_card"]);minimum_samples=32768)
```

`minimum_samples` sizes the noise covariance period: it must be at least the
internal sample count of the longest event, 200 per frame (100 us sampling
over 20 ms). 32768 covers 163 frames and costs a few MB. Load once per run.

Per event, with `ions::IonDetectorResponse` from `sim_ion_detector_response`:

```julia
r=sim_vh1_electronics(erng("vh1_noise"),erng("vh1_adc"),ions,s)      # VH1FrameResponse
plan=prepare_vh1_reconstruction(r.time,ions.ions.time,s)             # VH1ReconstructionPlan
ri=reco_vh1_ions(r,plan)                                             # ReconstructedIons
```

- `VH1FrameResponse`: `codes::Matrix{UInt16}` with one row per reading
  (`conversions_per_frame` per frame, time order) and one column per pad;
  `time::VH1FrameTime` (start_ms from event t0, frame_ms, frame_count,
  reference_frames); `saturated::BitVector` per pad; ADC out-of-range counts.
  The first frame is a reference frame with no induced current.
- `ReconstructedIons` is the existing type: rows are collection buckets,
  columns pads, signed ions, `valid`, and `checks` (baseline is the reference
  reading, `nonlinear` always false, `adc_at_limit`, `converged` false for
  saturated or clipped pads, RMS residual of the readings). Hits, energies
  and tracks downstream need no change.
- `frame_voltages(r,s)` converts codes to volts (interval centres) for plots.
- `plan.ion_covariance` is the predicted noise covariance of the bucket
  counts; `sqrt.(diag(...))` is the per-bucket noise to compare with quiet pads.
  The plan depends on the event only through its bucket count and can be
  cached by `bucket_count`.

Seeds: derive the two RNGs from the electronics seed with new tags
(`"vh1_noise"`, `"vh1_adc"`) through the existing `seedfor`, so the present
chain's streams stay unchanged when both run on the same events.

## Production loop

In `run_framed_reconstruction.jl` the three calls
`electronics_settings` / `sim_ion_electronics_prepared` /
`prepare_ion_reconstruction` + `reconstruct_ions` are replaced by the three
above. The electron arm (SiPM, ASIC, PE reconstruction) is untouched. For the
pilot, run both chains on the same detector response in the same loop so the
comparison is event by event.

Record per event: ions input and reconstructed (valid cells), saturated pads,
pads at ADC limit, quiet-pad statistics (sum over buckets of reconstructed
ions in empty pads; its variance against `sum(plan.ion_covariance)`), the
per-bucket predicted noise, timing of the three calls, and the same
quantities from the present chain. Keep a few example pads: input counts,
noiseless readings, codes, reconstructed counts.

Results go under
`results/readout/<run_label>/vh1_pilot_100/` with `controls.h5`, a Python
presentation (`index.html`, plots) in the style of
`present_framed_electronics.py`, and the results column of the VH1 row in
`notes/html/index.html`.

## Policies

- Saturated pads are invalid (NaN), as nonlinear failures are today; report
  the count and the maximum noiseless drop. The 3 fF card saturates 1.5 % of
  the pads at 4 mm pitch (12.7 % of the ions); production uses the larger
  capacitors below, which left no saturated pad in the 100-event pilots.
- Reference frame: the module places it so that it ends where the earliest
  induced current of the first bucket begins. No baseline estimate is needed.
- Out-of-range converter codes make the pad invalid (same rule as today).
- No SPICE is ever launched; the cards carry the spectrum and the gain.

## Cards

Production cards (2026-11-14 pilots, zero saturated pads, about 33 ions per
bucket): `vh1_device_changes_cf30f.toml` for 4 mm pitch (30 fF, 5.3 uV per
ion, range 92 300 ions) and `vh1_device_changes_cf20f.toml` for 3 mm pitch
(20 fF, 7.9 uV per ion, range 61 900 ions); both span the converter over
0.3 to 1.0 V. `vh1_device_changes.toml` (3 fF) is the SPICE-study baseline
and saturates busy pads; `vh1_device_changes_cf10f.toml` is intermediate.

Keys, with the values of `vh1_device_changes.toml`: `gain_V_per_ion`
(4.686e-5), `charge_fraction` (0.8775), `feedback_capacitance_F` (3e-15),
`v_ref_V` (0.886), `output_window_V` ([0.4, 0.95]), `frame_ms` (20),
`conversions_per_frame` (2), `conversion_window_ms` (10), `reference_frames`
(1), `adc_bits` (16), `adc_minimum_V`/`adc_maximum_V` (0/1.8),
`adc_input_noise_rms_V` (0), `noise_sample_us` (100), the single-ion current
prescription (same as the present calibration), and the SPICE sources with
their SHA-256. `load_vh1_settings` refuses a card whose recorded hashes do not
match the files. If the SPICE study is rerun, regenerate the hashes:

```sh
shasum -a 256 results/ltspice/design/VH1/improvements/noise_combined_W4_23u.csv results/ltspice/design/VH1/improvements/summary.json
```

`vh1_base.toml` is the same architecture on the unchanged amplifier (42.9 uV
per ion, about 45 ions per bucket) for comparison runs.

## Validation before production

1. `tests/vh1/test_vh1.jl` passes (exact recovery, spill, window flag,
   noise covariance within 8 %).
2. On the 100-event pilot: quiet-pad variance agrees with the predicted
   covariance; reconstructed totals equal input totals within the predicted
   noise; no pad invalid for reasons other than saturation; timing per
   event reported.
3. The comparison with the present chain on the same events: per-bucket
   residual RMS, totals, invalid pads.

## Open items (not blocking)

- Switch charge injection, reset noise and clock feedthrough are not modelled.
- Whether the signal-to-noise per pad-bucket at 3 mm pitch (6.8 mean) is
  sufficient for tracking is to be settled with the hit-cut studies.
