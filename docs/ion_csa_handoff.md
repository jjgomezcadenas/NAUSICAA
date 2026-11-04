# Ion CSA agent handoff

## Audit handoff — completed CSA work, 2026-10-30

This section records the completed CSA work for the agent assembling the final
full-statistics runs. The input-file contract and earlier implementation plan
remain below. Repository root: `/Users/jjgomezcadenas/Projects/NAUSICAA`.
The branch observed when this handoff was written is `prod`. Other agents have
concurrent changes in the working tree; preserve their work and Git staging.

### Decisions to preserve

- Adopt circuit v3: nominal 20 s retention, fitted approximately 19.544 s,
  5 ms centred voltage averages, external input capacitance 1.25 pF,
  feedback capacitance 3 fF and 0.5 mm square pads. The feedback transistor
  has W = 0.36 micrometres, L = 1200 micrometres and gate at 0 V. Exact circuit
  parameters are in the cards, calibration and netlists.
- Compare 3 mm / 30 ms and 4 mm / 40 ms for the design decision. The user has
  set aside 2 mm / 20 ms because of its lower S/N; its pilot remains the
  comparison baseline. S/N is a quality metric. Apply no S/N acceptance cut.
- Retain the three-frame inverse. A scalar frame calibration leaves positive
  total-charge bias because adjacent frames receive part of the response.
  Transform both the signal and its temporal noise covariance.
- Save a compact list of pixels with positive input charge. Every listed
  pixel retains the entire event bucket interval, including zero-input,
  leading quiet and recovery buckets. A remote deposit adds its own pixels
  without filling the intervening spatial rectangle. Absent pixels received
  no input charge and their electronic noise was not simulated.
- Keep signed reconstructed counts, including negative noise fluctuations.
  Measured zero suppression is a later analysis. No ADC quantisation enters
  these response files.
- Keep the finite-path current integral, uniform mean arrival times within
  each bucket, and the same square-pad weighting benchmark across geometries.
  One Xe_2^+ carries +e. No extra collection impulse or charge renormalisation.
- Noise is stationary Gaussian noise from the quiescent one-sided SPICE PSD,
  with full temporal covariance and independent pad channels. Each event
  starts at the quiescent circuit state. These assumptions are recorded in
  the output; pad histories include accumulation within each event.

### Completed outputs and reference results

All three pilots use the same 100 deterministic, evenly spaced event IDs per
sample: 200 events per geometry, 600 response events over six files.

| Pitch / frame | ENC, elementary charges RMS | Cell mean S/N, bb0nu | Cell mean S/N, single electron | Mean noiseless charge bias, bb0nu / electron |
|---|---:|---:|---:|---:|
| 2 mm / 20 ms | 51.263 | 1.265 | 1.077 | -0.0893% / -0.0701% |
| 3 mm / 30 ms | 55.322 | 2.959 | 2.501 | -0.2897% / -0.2366% |
| 4 mm / 40 ms | 58.952 | 5.207 | 4.363 | -0.6766% / -0.5807% |

This rounded table is an audit reference. Generated note tables obtain their
values from the saved HDF5 summaries. At the >=50% local-loss diagnostic,
affected-event fractions are 0% for both samples at 2 and 3 mm, and 6% / 7%
at 4 mm. An affected event contains a locally compressed pixel-frame; the
response retains that event. Pitch and frame width change together.

Corrected studies, each containing `summary.h5`, CSV summaries, `outputs.json`
and a standalone `ion_csa.html`:

```text
results/ion_csa_v3_corrected_2mm/621ee0d8943bf79c
results/ion_csa_v3_corrected/bfe6dd3720d0e9b9
results/ion_csa_v3_corrected_4mm/1a907c5dff36993d
```

The final published pilot manifest is:

```text
results/ion_csa_v3_production_pilot/d3ec668fa3aa637f/manifest.toml
results/ion_csa_v3_production_pilot/latest.toml
```

Use its six `responses[].path` and `sha256` entries to locate the final
`nausicaa.ion_csa_response.v1` files. These contain the corrected numerical
payload unchanged, plus explicit channel-selection and parent-response
metadata. The exact input ion counts, noiseless counts, sampled noise, signed
reconstructed counts, voltage averages, calibration and covariance are saved.
HDF5/Python matrix axes are `(pad, frame)`; Julia axes are `(frame, pad)`.
The source event and pair identities match across pitches within each sample.
The two samples have separate primary realisations.

### Code and documentation map

| File | Purpose |
|---|---|
| `src/ion_signal/ion_csa_pipeline.jl` | Production card resolution, source replay, publication, run manifests and owned event/covariance readers. |
| `scripts/production/run_ion_csa.jl` | Julia production entry point. |
| `src/ion_signal/ion_csa_response.jl` | Input adaptation, finite-gap current, covariance, linear propagation, SPICE reduction and filtered response writer. |
| `src/ion_signal/ion_csa_correction.jl` | Three-frame inverse, transformed covariance, corrected writer and noise-statistic checks. |
| `src/ion_signal/ion_csa_summary.jl` | Numerical histograms, compression, charge and quality reductions. |
| `scripts/production/run_ion_csa_pilot.py` | Parallel LTspice acquisition and Julia stage orchestration, including automatic quadrature retries. |
| `scripts/production/run_ion_csa_correction.py` | Corrected reconstruction and presentation orchestration. |
| `scripts/production/present_ion_csa_pilot.py` | HTML, scientific figures and report source from saved summaries. |
| `config/specs/response/ion_csa/csa_v1.toml` | Circuit, input, noise, averaging, seed and numerical controls. |
| `config/specs/response/ion_csa/reconstruction_v1.toml` | Three-frame inverse, covariance and signed-estimate policy. |
| `config/specs/runs/ion_csa_v3_pilot_all.toml` | Registers the six completed corrected pilot files. |
| `config/specs/runs/ion_csa_v3_production.toml` | New full-statistics calculation: 8,000 events per sample, pitches 3 and 4 mm. |
| `docs/ion_csa_response.md` | Current API, output schema, assumptions and commands. |
| `docs/ion_csa_pilot_review.html` | Reviewed input types and integration contracts. |

`IonCSAProduction.read_event(path, event_id)` returns an owned
`ReconstructedIonEvent`. It exposes input, noiseless, noise and reconstructed
matrices, signed pixel/bucket keys, time edges in seconds, reference/phase in
mm, drift lengths, frame sigma, model flags and seeds. `read_covariance` gives
the exact reconstructed covariance for the saved frame count.

The technical note is `notes/physics/csa_design.tex`. Its standalone source,
local bibliography, circuit image, generated values/tables and four scientific
figure pairs are in `latex_generated/csa_design/`. Numerical reduction is
`scripts/production/reduce_csa_design.jl`; presentation is
`scripts/production/present_csa_design.py`; the card is
`config/specs/runs/csa_design_note.toml`; saved numerical results and hashes are
in `results/csa_design/`. Julia reduction passed 38 matched-input/summary
checks. The LaTeX was checked statically, without compilation. The manuscript
`notes/physics/nausica.tex` was unchanged. Compile only when the user requests
it; visually inspect a PDF only on an explicit request.

### Validation already performed

- 2 mm pilot: 1,609 independent filtered-response checks, 2,807 independent
  corrected-response checks, 13 HTML-control checks and 410 cross-pitch/hash
  checks. All 22 frozen source hashes remained unchanged during that run.
- Production publication: `independent_validation.json` beside the final
  manifest records 626,305 checks across six files / 600 events. Every event,
  calibration and covariance dataset exactly matches its corrected parent.
  Independent source replay checks every sparse ion cell and paired event.
- Julia tests: response 18; correction 25 plus 8 finite-sample Gaussian-moment
  checks; production reader/publication 30. `make test physics` passed the
  57 Python tests and physical benchmark.
- A separate one-event end-to-end execution exercised new acquisition,
  noise, inverse reconstruction and publication. Its saved manifest is
  `results/ion_csa_v3_integration_smoke/103b01baa236c4e1/manifest.toml`;
  independent verification passed 434 checks. This predates the final
  provenance-hardening edits; the final production tests separately cover
  saved-summary checksums and changed-calibration rejection.

Relevant small fixes made while integrating:

1. Event selection handles a one-event request by choosing the first sorted
   source ID. The previous `range(...; length=1)` failed for a larger source.
2. The noise ensemble checks account for finite sample size. Mean, second
   moment and adjacent-product limits are the larger of the configured
   tolerance and six covariance-derived standard errors. Gaussian moment
   identities independently test those variances. This accommodates small
   event subsets; the physical integration and inverse bounds are unchanged.
3. Saved-study import verifies the response hash recorded in its summary.
   Pipeline identities include calibration contents and numerical backend
   sources; imported physical model IDs must match the selected calibration.

### Audit focus before full statistics

1. Check the production interfaces, selected calibration and exact numerical
   equivalence to the frozen pilot outputs. Audit the complete circuit/cache
   dependency set, including transient-calibration HDF5, effective netlist
   biases and transistor include files. The existing acquisition reuse checks
   PWL-current and transistor-model hashes; verify that any proposed change to
   circuit resources invalidates all affected caches.
2. Confirm the independent SiPM branch interfaces with its owner. The ion
   stage emits reusable nodes and per-realisation `ion_manifest.toml`; final
   detector-graph assembly with SiPM outputs remains the integrating agent's
   task. Join by original sample/event ID, shared realisation and pair checksum.
   Preserve the different ion/photon bucket units and event-relative clocks.
3. Run and time a modest subset through the final combined driver before
   full production. Runtime depends on pad histories, their lengths and the
   number requiring direct SPICE. The 2 mm 200-event pilot took about 191 s
   after preparation, including 1,533 direct histories, solver controls and
   three automatic quadrature retries. Treat this as one measured case.
4. Execute the full 3 and 4 mm ion studies when the final run is launched;
   no 8,000-event CSA response run has been performed. Verify all selected
   IDs, source/pair ledgers, covariance and retained signed/zero buckets.
5. Generate full-statistics S/N distributions, event means, compression
   severity/affected-event rates and signed charge-bias distributions. Retain
   the pilot summaries and regenerate new reports from the full HDF5 outputs.
   The 3-versus-4 mm decision requires track continuity and blob observables.
6. Keep channel selection explicit. Noise-only spatial occupancy requires a
   defined acquisition region; the present files simulate noise on input-
   charged pixels and preserve all time buckets on those channels.

### Commands for the next agent

Run from the repository root. The known local environments are:

```sh
export JULIA_DEPOT_PATH="$PWD/work/julia-depot:$PWD/notes/physics/work/julia-depot:$HOME/.julia"
TASK_JULIA=/Users/jjgomezcadenas/.juliaup/bin/julia
TASK_PYTHON=/Users/jjgomezcadenas/miniconda3/bin/python3

"$TASK_JULIA" --startup-file=no --project=. tests/test_ion_csa_response.jl
"$TASK_JULIA" --startup-file=no --project=. tests/test_ion_csa_correction.jl
"$TASK_JULIA" --startup-file=no --project=. tests/test_ion_csa_pipeline.jl
"$TASK_PYTHON" tests/production/check_ion_csa_pipeline.py \
  results/ion_csa_v3_production_pilot/d3ec668fa3aa637f/manifest.toml
make test physics PYTHON="$TASK_PYTHON"
```

To validate/register the existing pilots through the final interface:

```sh
"$TASK_JULIA" --startup-file=no --project=. scripts/production/run_ion_csa.jl \
  config/specs/runs/ion_csa_v3_pilot_all.toml
```

To start the full ion run once the audit/combined-run preparation is complete:

```sh
"$TASK_JULIA" --startup-file=no --project=. scripts/production/run_ion_csa.jl \
  config/specs/runs/ion_csa_v3_production.toml
```

The selected circuit calibration must exist under
`results/ltspice/csa_v3_response_calibration/`; its content hash is checked.
The acquisition executable defaults to the local LTspice installation used by
the pilots. The Julia API accepts explicit Python/simulator paths and worker
count if needed. Presentation uses `work/event-viewer-venv/bin/python`.
The full-run manifest path is printed and saved to the result root's
`latest.toml`; use it rather than assuming the pilot graph ID survives edits.

Relevant source/cards/reports are staged, with no commit or push performed by
this CSA task. Preserve the concurrent agent's unstaged work, including its
SiPM Makefile block. HDF5, raw SPICE, PWL and large embedded-Plotly HTML files
remain outside Git. Never configure Git LFS without explicit user approval.

---


The CSA agent owns the Julia response implementation, its tests, the CSA
cards, the pilot runner and the numerical HDF5 summaries. Physical response
code lives in `src/ion_signal/ion_csa_response.jl` and
`src/ion_signal/ion_csa_summary.jl`; the pilot entry point is
`scripts/production/ion_csa_pilot.jl`. Python presents the saved results.
Read `docs/ion_csa_response.md` and `docs/ion_csa_pilot_review.html` for the
current production interfaces. The corrected 100-event-per-sample pilots now
cover 2, 3 and 4 mm pitch (20, 30 and 40 ms buckets).

## Production entry point and saved reconstruction

Run `julia --project=. scripts/production/run_ion_csa.jl
config/specs/runs/ion_csa_v3_pilot_all.toml` with the scientific Julia depot.
This registers the six corrected pilot outputs as versioned production response
nodes and writes a complete ion-stage manifest. Each source realisation also
receives `runs/<graph_id>/ion_manifest.toml` for final detector-graph assembly.
The manifest supplies output/source paths, content hashes, selected event IDs,
frame widths, geometry and shared primary-pair identities. A branch without
`corrected_study` executes acquisition, noise and reconstruction. The full-sample
card is `config/specs/runs/ion_csa_v3_production.toml`; execute it when a full
production run is requested.

`src/ion_signal/ion_csa_pipeline.jl` exports the owned
`ReconstructedIonEvent` reader and exact temporal covariance. The production
schema is `nausicaa.ion_csa_response.v1`. Its matrices copy the corrected v2
payload exactly, including signed `recovered_ions`, input counts, mean, noise,
calibration and covariance. Matrix axes are `(pad, frame)` in HDF5/Python and
`(frame, pad)` in Julia.

The agreed acquisition policy saves a compact list of **pixels with input
charge**, with the complete bucket interval for each listed pixel, including
zeros and quiet/recovery frames. Far deposits add their own pixels without
filling a spatial rectangle. Absent pixels received no input charge and their
noise was not simulated. Record `channel_selection =
"pixels_with_positive_input_charge"`, `absent_channel_policy =
"no_input_charge_noise_not_simulated"`, measured zero suppression false and
negative-estimate clipping false. Measured zero suppression belongs to later
analysis. See the response documentation for datasets, API and controls.

## Files to consume

The local production run contains 8,000 events from each 1 mm Geant4 hit
sample. Its ion outputs are under
`results/simulations/production_pilot/`, outside Git because the six HDF5
files total about 777 MB. Paths in the table are relative to that directory.

| Sample | Pad pitch | Frame | Ion file |
|---|---:|---:|---|
| `bb0nu` | 2 mm | 20 ms | `bb0nu/c25084b9e721e407c0688bfdba837f8f7fe5e453d8dad73034681d2d994ad95c/detector/ions/17bb8250a71cf1c20a3e6efe80fc0726464623a10f087946322f29ad971743fe/ions.h5` |
| `bb0nu` | 3 mm | 30 ms | `bb0nu/c25084b9e721e407c0688bfdba837f8f7fe5e453d8dad73034681d2d994ad95c/detector/ions/d978512060a41657396c2741c45a5cc4925c5f22daab718af0831cfac1d153c0/ions.h5` |
| `bb0nu` | 4 mm | 40 ms | `bb0nu/c25084b9e721e407c0688bfdba837f8f7fe5e453d8dad73034681d2d994ad95c/detector/ions/14d59c9436d9cdfa37465d0e788feefa4b84f153c118c1bc1ed5c937e3cd06ac/ions.h5` |
| `single_electrons` | 2 mm | 20 ms | `single_electrons/a82036783cb95a9651bf8ca7bc87056fa74f54cd8869d9a236614d41b6f590c7/detector/ions/7906df4ce60797d995dba62e6aeaa8c6a6652c7ef4e3311774226942e08be85e/ions.h5` |
| `single_electrons` | 3 mm | 30 ms | `single_electrons/a82036783cb95a9651bf8ca7bc87056fa74f54cd8869d9a236614d41b6f590c7/detector/ions/0a1582d36d85eb59dafe112515fe9632579faa63703124e5deb8caf263612e8d/ions.h5` |
| `single_electrons` | 4 mm | 40 ms | `single_electrons/a82036783cb95a9651bf8ca7bc87056fa74f54cd8869d9a236614d41b6f590c7/detector/ions/efb4c32699fa7b814ce163d7a8f01ad5b927f93100c59438013ed5918de9d956/ions.h5` |

The machine-readable index
`results/simulations/production_pilot/grid_scan_summary_238ddc5a75e0fd9b.h5`
also gives each absolute path in `/<sample>/ions/<pitch>/path`, the frame
width, count totals, and the resolved ion cards. The three pitches for a
sample share the same event placement and primary-pair realisation. The two
samples are separate realisations. Start with the 3 mm/30 ms files if one
configuration is needed for the first CSA integration.

## HDF5 contract

The file root has `schema = "nausicaa.ion_pad_file.pilot.v1"`. Check
`metadata/complete`, `metadata/event_count`, `metadata/sample`,
`metadata/signal_unit = "ions"`, `metadata/time_unit = "ms"`, and the embedded
`metadata/branch_card_toml`. The card gives `pad_pitch_mm`, `bucket_width_ms`,
`drift_speed_m_s`, temperature, and drift field. The file also records the
source SHA-256, seed, card hash, code hash, and realisation card.

Each `/events/<event_id>` group has schema
`nausicaa.ion_pad_history.v1` and contains:

| Dataset | Meaning |
|---|---|
| `pixel_ix`, `pixel_iy`, `bucket` | Signed `Int64` pad and frame indices, one entry per occupied pad–frame cell. The four sparse columns are sorted by `(pixel_ix, pixel_iy, bucket)`. |
| `count` | Positive `Int64` ion count in that cell. One ion carries `+e`; convert explicitly to coulombs for the CSA. |
| `initial_pairs` | Event's sampled primary-pair count. `sum(count) == initial_pairs`. |
| `pitch_mm`, `bucket_width_mm` | Pad pitch and longitudinal width; here `bucket_width_mm` equals 2, 3 or 4 mm. |
| `reference_mm`, `phase_mm` | Source-coordinate grid edge reference and transverse phase, in mm. The x/y edge for index zero is `reference_mm[0:2] + phase_mm`; the z edge is `reference_mm[2]`. |
| `drift_direction` | `+1` for ions. |
| `ion_drift_cm`, `electron_drift_cm` | Complementary event drift lengths; they sum to 150 cm. |
| `event_id`, `pair_checksum` | Original event ID within the sample and checksum of the shared primary-pair draw. |

Use the HDF5 event ID, rather than an assumed numeric position, to identify
an event. The `bucket` index is an event-relative longitudinal bin: for the
ion branch it follows the diffused source `z` from `reference_mm[2]`, with
width `bucket_width_mm`. With drift speed 0.1 m/s, one bin is the card's
20, 30 or 40 ms frame. Indices can be negative when diffusion crosses the
reference edge. Preserve gaps in the sparse bucket sequence as zero-ion
frames. The file does not encode an absolute acquisition clock or an
individual ion's time within its frame.

Julia readers can use `EventAnalysis.read_ion_history(group)` from
`src/event_analysis/production_readout.jl`. A direct Python/HDF5 reader can
discover the selected path and read an event as follows:

```python
from pathlib import Path
import h5py

root = Path("results/simulations/production_pilot")
with h5py.File(root / "grid_scan_summary_238ddc5a75e0fd9b.h5") as index:
    path = Path(index["bb0nu/ions/3/path"][()].decode())
with h5py.File(path) as source:
    assert source.attrs["schema"] == b"nausicaa.ion_pad_file.pilot.v1"
    assert bool(source["metadata/complete"][()])
    event = source["events/0"]
    ix, iy = event["pixel_ix"][:], event["pixel_iy"][:]
    frame, ions = event["bucket"][:], event["count"][:]
    assert int(ions.sum()) == int(event["initial_pairs"][()])
```

## CSA response task

### TOML graph and output location

The agreed stage-2 configuration layout is in
`docs/detector_response_pipeline_plan.md`:

```text
config/specs/response/ion_csa/csa_v1.toml
config/specs/response/sipm/sipm_v1.toml
config/specs/response/sipm_electronics/readout_v1.toml
config/specs/runs/<detector_and_response_name>.toml
```

The ion CSA agent owns `config/specs/response/ion_csa/csa_v1.toml`, its Julia
reader/runner, tests and response HDF5 files. The **ion CSA card owns both the
circuit response and its noise model**. Its v1 schema must define, with
explicit units and provenance:

| Section | Required content |
|---|---|
| Circuit | Selected design and variant; circuit/transfer-function source and content hash; charge and voltage calibration; discharge model. |
| Input and readout | Conversion from ion counts to the assumed current or per-ion response within a frame; frame width from the selected ion card; sampling aperture and frame estimator. |
| Noise | One-sided spectrum or covariance source and content hash; spectral bounds; channel independence/correlation assumption; frame covariance calculation and RNG seed/policy for a sampled realisation. |
| Provenance | Schema/model version, code identity and the physical assumptions used for timing and noise. |

The current card defines the concrete key names; its Julia parser validates
them before running. Distinct response identities are needed when circuit,
noise, readout or seed changes. The ion response identity also includes the selected
upstream `ions.h5` identity. SiPM and SiPM-electronics cards belong to the
other response branch. The final-run TOML will select one detector
configuration and all three response cards when those branches are ready;
the ion CSA agent should expose its card and output identity for that later
assembly. The SiPM agent creates the two SiPM response cards.

The implemented pilot selects variant `length_1200um`, nominal 20 s
discharge and 5 ms averaging through `config/csa_design_v3.toml` and
`config/specs/response/ion_csa/csa_v1.toml`. Its calibration card is
`config/csa_v3_response_calibration.toml`, with numerical output under
`results/ltspice/csa_v3_response_calibration/`. The pilot derives the 30 ms
calibration from the selected circuit and spectrum. The run card is
`config/specs/runs/ion_csa_v3_pilot.toml`; its response schema and validation
are documented in `docs/ion_csa_response.md`. The final detector-and-response
run card remains to be assembled with the SiPM branch.

Write one versioned ion response per sample, ion file and CSA/noise card under
the existing production root:

```text
results/simulations/production_pilot/<sample>/<realisation_id>/
  response/ions/<ion_output_id>/<csa_spec_id>/ion_response.h5
```

Here `<ion_output_id>` identifies the selected `ions.h5`; `<csa_spec_id>`
identifies the parsed CSA/noise card and implementation. The run manifest
under `runs/<final_graph_id>/manifest.toml` later references that response
file and the detector branch files. One CSA card applied to this scan yields
six ion-response files, one for each sample/pitch combination. Write each
response atomically, mark it complete after its ledgers pass, and store the
exact source path/hash, resolved cards, seed and code identity. Put numerical
analysis summaries in `results/`; generate any final plots and tables with
Python and a self-contained report in `latex_generated/`. The raw branch HDF5
files remain outside Git because of their size; scripts, cards, notes and
generated report artifacts are tracked.

### Response and validation

1. Treat each `(pixel_ix, pixel_iy)` as an independent pad channel and its
   ordered `bucket` counts as the ion input history. Keep all positive cells;
   the `n_min` setting used in some plots is a display selection, not detector
   charge loss. Include zero-ion frames and quiet channels where the noise or
   baseline model needs them.
2. Specify the assumed ion input-current or per-ion CSA response within a
   frame. The files give the integrated arriving-ion counts per frame, so a
   within-frame timing or current-shape model is needed for a time-domain CSA
   calculation. Earlier CSA benchmarks used uniform arrivals within a frame;
   assess whether that approximation serves this response and consult the
   user before adopting a different physical timing model. Use the existing
   `src/ion_signal` weighting-field/current and electronics work where
   applicable, and record the chosen timing model.
3. Add CSA noise, baseline memory and frame-to-frame noise covariance with the
   electronics model's own units and parameters. Preserve the full per-pad
   frame sequence through the response; integrated pad charge is a separate
   energy observable.
4. Save event and pad/frame keys, noiseless response, noise realisation or its
   covariance, and final sampled response in the versioned `ion_response.h5`.
5. Validate one-ion charge/current normalization, zero-input noise, frame
   alignment, carryover across adjacent frames, and conservation of the input
   ion ledger. Compare the 2, 3 and 4 mm configurations on their shared event
   IDs; each has 8,000 events. The photon/SiPM branch is separate and does not
   enter this CSA input.

The detector pipeline design and provenance rules are in
`docs/detector_response_pipeline_plan.md`. The numerical comparison is in
`results/readout_pulses_8000/readout_pulses.html` and its reduced
`readout_pulses.h5`; those are summaries, while the six `ions.h5` files are
the CSA inputs.

## Verified handoff checks

On 2026-10-30, all six local files had complete metadata and 8,000 event
groups. Every event passed `sum(count) == initial_pairs`, and every file's
summed count agreed with the grid-scan index. The total was 897,306,359 ions
for `bb0nu` and 897,244,198 for `single_electrons` at each pitch. Event-ID
sets matched across the three pitches within each sample.
