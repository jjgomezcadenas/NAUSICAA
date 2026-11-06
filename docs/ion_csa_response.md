# Ion CSA response and reconstructed matrices

> The pilot paths and completion statements below describe the original CSA
> workflow. The current completed 8,000-event matched 4 mm production is
> documented in [elm_matched_production.md](elm_matched_production.md).
> New production uses the validated fast nonlinear CSA by default, with
> `mode = "production"`. Set `mode = "validation"` to repeat SPICE controls. See the
> [model and interfaces](ion_csa_nonlinear_model.html) and
> [validation report](../latex_generated/ion_readout/ion_csa_nonlinear/index.html).
> The SPICE/linear workflow below remains selectable as `hybrid_spice`.

The ion response uses the adopted v3 circuit: nominal 20 s retention, fitted
19.544 s discharge, 5 ms voltage averages and the three-frame charge
reconstruction. The matched pilots contain 100 evenly spaced event IDs from
each of the double-beta and single-electron samples at 2, 3 and 4 mm pitch,
with 20, 30 and 40 ms buckets respectively. S/N measures signal quality; the
response applies no S/N cut or clipping of negative estimates.

## Production entry point

From the repository root, using the Julia scientific depot:

```sh
JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia" \
  julia --startup-file=no --project=. scripts/production/run_ion_csa.jl \
  config/specs/runs/ion_csa_v3_pilot_all.toml
```

This card selects the six completed corrected pilot files. The Julia stage
replays their source ion histories, checks the circuit/reconstruction cards
and matched event identities, and publishes versioned production HDF5 nodes.
The numerical matrices are copied exactly. Each node adds explicit channel
selection and parent-response metadata. A run manifest references the six
outputs, source files, selected event IDs and hashes. An ion-stage manifest
also appears under each source realisation's `runs/<graph_id>/` directory for
later assembly with the SiPM response branch.

`config/specs/runs/ion_csa_v3_production.toml` selects all 8,000 existing events
per sample at 3 and 4 mm. Running this card starts a new full simulation.
It has been prepared for an explicit production request; the current work
ran the 100-event pilots and a separate one-event integration test.

Each `[[branches]]` entry selects a pitch. A `corrected_study` path imports a
saved corrected study. With that key omitted, the stage runs the calibrated
nonlinear response, domain fallback, noise and three-frame
reconstruction. Set `mode = "validation"` and `response_model = "hybrid_spice"` at the run-card top level
to use the SPICE reference below the guard. The nonlinear calibration HDF5 is
resolved and hash-pinned in the execution's `resolved_response.toml`.
Python launches LTspice; physical response, numerical reduction and HDF5
writing use Julia. The selected circuit calibration must already exist; generate it with
`config/csa_v3_response_calibration.toml` before production. Its content hash
enters the execution identity and must match imported corrected studies.
The API accepts `python`, `simulator` and `workers` keyword
arguments for local tool paths and acquisition concurrency. The numerical
backend restores the configured Julia depot in its child processes.

## Channel and bucket selection

The saved rows are a compact list of **pixels with input charge**. Each pixel
appears once, with its signed `(pixel_ix, pixel_iy)` key. A distant isolated
deposit adds its occupied pixels without filling the intervening spatial
rectangle. Every saved pixel has the complete event bucket interval, including
zero-input buckets, two leading quiet frames and the configured recovery tail.

Absent pixels received no input charge; their electronic noise was not
simulated. This input-based channel selection is recorded as
`pixels_with_positive_input_charge`. Measured zero suppression can be studied
as a separate analysis. The primary response retains weak deposits and signed
noise fluctuations.

## Julia analysis interface

`src/ion_signal/ion_csa_pipeline.jl` defines `IonCSAProduction`:

- `settings(root, card)` validates the production card and resolves paths and
  configuration identities.
- `run(root, card)` executes or imports each response branch, validates it and
  returns the complete ion-stage manifest path.
- `read_event(path, event_id)` returns an owned `ReconstructedIonEvent`.
- `read_covariance(path, nframes)` returns the exact reconstructed temporal
  noise covariance in elementary-charge squared.
- `validate_response(path)` checks all selected events against the source ion
  file and returns its count, identity and signed-noise ledger. An optional
  `source_path` supports relocation with unchanged source content.

Example downstream analysis:

```julia
include("src/ion_signal/ion_csa_pipeline.jl")
const CSA = IonCSAProduction
trace = CSA.read_event(response_path, 3790)
counts = trace.reconstructed_ions  # Julia axes: frame, pad
input = trace.input_ions
buckets = trace.first_bucket .+ (0:size(counts, 1)-1)
pixel_keys = collect(zip(trace.pixel_ix, trace.pixel_iy))
covariance = CSA.read_covariance(response_path, size(counts, 1))
```

`ReconstructedIonEvent` owns Int64 input counts, Float32 noiseless/noise/final
matrices, Float64 time edges in seconds and per-frame noise sigma in ion
charges. It preserves the Float32 reference/phase coordinates in mm, bucket
width in mm and complementary drift lengths in cm. Pad and bucket indices are
signed. Mutating returned arrays leaves the HDF5 file unchanged. The reader
rejects incompatible units, unsupported reconstruction, duplicate pixel keys,
invalid shapes, missing event IDs, nonfinite samples, inconsistent time axes,
zero-input-only pixels and broken charge/additivity ledgers.

The HDF5/Python matrix axes are **pad, frame**. Julia reads these HDF5 arrays
as **frame, pad**. Metadata records both conventions. The voltage arrays use
frame boundaries and contain one extra time sample.

## HDF5 outputs

The production schema is `nausicaa.ion_csa_response.v1`. Its numerical payload
is the validated `nausicaa.ion_csa_response.pilot.v2` corrected response.
Versioned files live under:

```text
results/simulations/production_pilot/<sample>/<realisation_id>/
  response/ions/<ion_output_id>/<response_id>/ion_response.h5
```

The run manifest gives exact paths and SHA-256 hashes. Large HDF5, raw SPICE
and PWL files stay outside Git; source, cards, small summaries, source snapshots
and report artifacts are tracked.

| Event dataset | Meaning |
|---|---|
| `pixel_ix`, `pixel_iy` | Compact signed Int64 pixel keys, in matrix order. |
| `first_bucket`, `edge_time_s` | First signed bucket and complete boundary times in seconds. |
| `counts` | Exact Int64 input ions, including restored zero buckets. |
| `linear_estimated_ions` | Reconstructed calibrated linear reference. |
| `noiseless_estimated_ions` | Reconstructed mean response, including transistor compression. |
| `noise_ions` | Sampled temporally correlated reconstructed noise. |
| `recovered_ions` | Final signed reconstructed ions: mean plus noise. |
| `noise_sigma_ions` | Per-frame standard deviation after reconstruction. |
| `output_at_edges_V`, `average_at_edges_V` | Mean output voltage and 5 ms average, in V. |
| `model_used` | Per pixel: 0 calibrated linear propagation, 1 direct SPICE. |
| `noise_seed` | Stable per-pixel UInt64 seed for replay. |
| `compression_level` | Local response-loss diagnostic on positive input cells. |

The event retains ID, pair checksum, geometry, placement and complementary
drift lengths. `covariance_by_frames/<N>/covariance_e2` gives the exact
frame covariance for a history of length N. Metadata snapshots the CSA and
reconstruction cards, physical/code identities, source path/hash, acquisition
policy and parent corrected response. Publication is atomic; completeness is
validated before the final file is renamed. Existing outputs are verified
against their saved node checksums before reuse.

## Numerical response

The existing backend is `src/ion_signal/ion_csa_response.jl`; its stage entry
point is `scripts/production/ion_csa_pilot.jl`. It consumes
`EventAnalysis.read_ion_history`, preserving the original Int64 sparse keys and
exact sampled ion counts. Counts are converted explicitly to SI charge for
the finite-gap current. Uniform mean collection times within each bucket and
the central square-pad weighting benchmark give an integral of
`sum(N) e [1 - phi_w(g)]`.

The calibrated SPICE kernel propagates every complete pad history. Direct
transistor transients replace the linear prediction for candidates and
controls across accumulated-voltage strata. The direct range expands until
remaining control errors meet the configured ENC bound. Charge quadrature
failures trigger finer transient acquisition; the original charge-integral
bound is retained. Independent timestep refinement checks direct histories.

`src/ion_signal/ion_csa_correction.jl` solves the three-frame response system
for mean and noise. The same linear transform propagates the temporal noise
covariance. Quiet/recovery buckets provide zero-collected-charge boundary
conditions. The original scalar-response files remain available as parent
artifacts.

The one-sided quiescent SPICE spectrum supplies full-band frame variance and
temporal covariance after 5 ms averaging and discharge-corrected differencing.
Pad channels are independent. Noise retains the quiescent spectrum for all
mean histories, including compressed signals. Cholesky sampling uses stable
source/model/event/pixel seeds. The three-frame inverse transforms the sampled
noise without drawing a new realisation.

Noise mean, second moment and adjacent-frame products are checked against
expectations. Their validation limits use the larger of the configured
sampling tolerance and six standard errors derived from the actual Gaussian
temporal covariance and number of independent pads. This accounts for small
runs while retaining the existing numerical integration and inverse bounds.
Float64 physical calculations precede validated Float32 storage.

## Diagnostics and presentation

`src/ion_signal/ion_csa_summary.jl` saves numerical histograms, per-event S/N,
compression flags, charge bias and representative histories. S/N is noiseless
reconstructed signal divided by the frame ENC, over every positive input
pixel-frame cell. Cell-weighted and event-weighted averages are distinct.
Compression flags compare mean response with the linear reference; they do
not discard events. Signed total charge includes quiet and recovery buckets.

Python presentation reads saved summaries:

```sh
work/event-viewer-venv/bin/python scripts/production/present_ion_csa_pilot.py \
  --results CORRECTED_STUDY_DIRECTORY
```

The three-pilot design note is `notes/physics/csa_design.tex`. Its numerical
reduction and presentation commands are:

```sh
julia --project=. scripts/production/reduce_csa_design.jl
work/event-viewer-venv/bin/python scripts/production/present_csa_design.py
```

Generated values, tables, figures and standalone source are in
`latex_generated/ion_readout/csa_design/`. Compilation remains a separate explicit step.

## Validation commands

```sh
julia --project=. tests/test_ion_csa_pipeline.jl
julia --project=. tests/test_ion_csa_response.jl
julia --project=. tests/test_ion_csa_correction.jl
python tests/production/check_ion_csa_pipeline.py RUN_MANIFEST
python tests/production/check_ion_csa_response.py FILTERED_STUDY_DIRECTORY
python tests/production/check_ion_csa_correction.py CORRECTED_STUDY_DIRECTORY
python tests/production/check_ion_csa_html.py HTML_FILE
make test physics
```

The production validator independently compares every published event,
calibration and covariance dataset with its corrected parent, replays every
upstream sparse cell and checks shared event/pair identities across pitches.
The Julia tests include distant pixels, signed/empty buckets, ownership,
changed-source rejection, publication/cache checks and analytic correlated
Gaussian moment identities. The integration smoke card exercises acquisition,
noise, correction and publication for one event; it is separate from the
matched 100-event comparison.

## Acquired-resource identity (2026-11-06)

`src/ion_signal/csa_resources.jl` verifies the actual production deck, transistor
model and linear waveform against their saved acquisition hashes before reuse.
New response identities and metadata include these digests. The analysis pins
the waveform, which pins the acquired deck and model. The detailed audit
follow-up is in [provenance fixes](reviews/provenance_fixes_2026-11-06.md).
