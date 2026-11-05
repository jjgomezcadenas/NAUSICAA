# SiPM response: audit and full-production handoff

## Authorization and task

On 2026-10-30, after reviewing the model, saved-channel policy and pilot
clipping fractions, the user authorized the co-worker to **audit the work
and run all events**. This supersedes the earlier instruction to stop after
the pilot. Use the approved sensor and ASIC cards, including the 1 mA limit;
the observed clipping is a result to quantify, not a reason to change the
range silently or request the same production authorization again.

The task is to audit the implementation and its provenance, make the
production orchestration ready, run the four complete photon histories,
and generate a self-contained report from saved numerical results. Full
production has not yet been run. The existing runner intentionally still
requires `pilot_only = true` and at most 16 explicit event IDs.

## Read first

1. `AGENTS.md` and `docs/sipm_electronics_handoff.md` for upstream input paths,
   units, photon-production assumptions and interfaces.
2. `docs/sipm_electronics_model.tex` and its compiled report
   `latex_generated/electron_readout/sipm_response/sipm_electronics_model.pdf` for the physics,
   sources, adopted approximations and later AXEL correction.
3. `docs/sipm_response.md` for output schemas and reproduction commands.
   `docs/mpcc_axel.md` records sensor specifications and the discussion.
4. The cards, implementation and tests listed below. The approved interface
   is also in the SiPM subsection of `docs/event_analysis_review/index.html`.

## Settled model and output policy

| Item | Adopted first model |
| --- | --- |
| Input | Incident photons, 250 ns buckets; optical acceptance and EL already applied |
| PDE | Binomial detection, probability 0.30; constant during recovery |
| Sensor | Effective 5 mm disk, 7854 microcells, 50 µm pitch assumption |
| Gain | 2.55 million electrons per fully recovered avalanche |
| Recovery | Charge fraction `1-exp(-elapsed/73.4 ns)`; each trigger resets its cell |
| Arrival/cell assignment | Uniform time within each bucket and uniform microcell index |
| Dark pulses | Independent primary Poisson process, 1 MHz per 9 mm² scaled by active area |
| Reference conditions | 25 °C and +4 V overvoltage; gain/dark values are reference assumptions |
| Initial state | Independent stationary dark-process age per microcell; recovered if dark rate is zero |
| Deferred effects | Crosstalk, afterpulsing, intrinsic gain spread, recovery-dependent PDE |
| ASIC | Bucket charge divided by duration; independent zero-mean 0.5 PE RMS Gaussian noise |
| Range | Upper clipping at 1 mA **bucket-average current**, after adding noise; retain negative estimates |
| Charge timing | All avalanche charge is assigned to the trigger bucket |

The user explicitly confirmed the following storage policy after the pilot:

- No minimum-charge or minimum-PE cut.
- Select every SiPM receiving any incident photon anywhere in the event,
  even if its sampled detected count is zero.
- Save all time buckets of those selected channels, including empty buckets,
  through the event's global photon span plus four buckets on each side.
- A photon far from the main track adds only its receiving SiPM. **Do not
  fill intervening SiPMs or a spatial bounding rectangle.**
- Dark pulses and ASIC noise populate the full time window of selected
  channels. Additional channels and the empty-event window are configurable.

The production readout has the same logical event–SiPM–bucket format as the
incident-photon input: signed `pixel_ix`, `pixel_iy` and `bucket` keys, with a
PE value for each saved cell. The sensor stage records the integer primary
`detected_pe` count. Recovery and dark pulses change the charge, so their
sensor output is also expressible as a real-valued PE equivalent
(`charge_C / charge_per_pe_C`). After ASIC noise and clipping, the final
readout is a signed, real-valued PE equivalent
(`current_A * bucket_width_s / charge_per_pe_C`), with a saturation flag.
Preserve the independent sensor HDF5 stage so an ASIC-card change can reuse
it. The selected channels form a sparse spatial set; the approved acquisition
window includes quiet buckets for dark pulses and noise.

## Files to audit

| Role | File |
| --- | --- |
| Sensor card | `config/specs/response/sipm/sipm_v1.toml` |
| ASIC card | `config/specs/response/sipm_electronics/readout_v1.toml` |
| Existing pilot card | `config/specs/runs/sipm_response_pilot.toml` |
| Sensor physics | `src/ion_signal/sipm_response.jl` |
| ASIC response | `src/ion_signal/sipm_electronics.jl` |
| Card resolution, acquisition, HDF5 and summaries | `scripts/production/sipm_pipeline.jl` |
| Pilot entry point | `scripts/production/run_sipm_pilot.jl` |
| Saved-result presentation | `scripts/production/present_sipm_pilot.py` |
| Focused tests | `tests/test_sipm_response.jl` |
| Upstream reader/types | `src/event_analysis/production_readout.jl`, `src/event_analysis/EventAnalysis.jl` |
| Canonical hashing/seeds | `scripts/production/specs.jl` |

The existing pilot v1 files save incident counts, primary detected PE, dark
pulses, optical/dark/total charge, channel seeds, a copy of the photon history,
and event identity. The ASIC pilot also saves noise, unclipped/clipped current
and saturation flags. These arrays support model audit; they do not define the
compact production data product. Keep the per-bucket PE outputs, required
flags and provenance; avoid copying the upstream photon history and all
intermediate arrays into both production files. See `docs/sipm_response.md`
for the current pilot schemas and paths.

## Inputs and current artifacts

Use this existing machine-readable index:

```text
results/simulations/production_pilot/electron_250ns_summary_464425eac32a2c34.h5
```

Its datasets `/<sample>/<gas>_path` resolve the four files. Samples are
`bb0nu` and `single_electrons`; gas alternatives are `pure_xe` and `xe_he15`.
Each file contains 8000 events: 32000 event/gas responses in total. The gas
alternatives share each sample's event/pair realization. Select actual sorted
event IDs from the source groups, and verify their count and gas-pair matches.
Do not regenerate upstream photons or use the older 2 µs histories.

Current validated pilot locator:

```text
results/sipm_response_pilot/latest.toml
```

It points to study `a237f33cb4c7a933`, whose `manifest.toml` identifies the
four independent sensor and four ASIC files. The pilot uses event 0 in each
input. Numerical results are `summary.h5`; focused statistical validation is
`results/sipm_response_validation/validation.h5`. The generated tables in
`latex_generated/electron_readout/sipm_response/pilot_results.tex` contain the reported
clipping counts and fractions. Their denominator is **all acquired
SiPM–bucket cells**, including quiet buckets; these are not charge-loss
fractions or ensemble estimates.

## Validation already completed

The final focused run passed 62 checks, covering binomial limits/moments,
independent two-trigger quadrature, stochastic recovery, stationary AXEL
mean, transient/inverse bias, area scaling and dark moments, ASIC noise and
clipping boundaries, signed gaps, empty events, HDF5 round trips, event-order
replay, independent sensor reuse and stale electronics regeneration.
Float32 charge-storage error is measured against Float64; current storage and
accumulators remain Float64. Existing production-readout tests passed 51
checks. `make test physics PYTHON=python3.13` and `make note PYTHON=python3.13`
passed. Both SiPM PDFs compiled with clean final logs, without visual PDF
inspection, following repository instructions.

Commands for the existing pilot and audit:

```sh
make sipm-test
make sipm-pilot PYTHON=python3.13
make test physics PYTHON=python3.13
make note PYTHON=python3.13
```

The Makefile uses `work/julia-depot:$HOME/.julia` and the project environment.
The complete pilot log is `work/sipm_complete_pipeline.log`; the existing
readout-test log is `work/sipm_existing_readout_tests.log`.

## Production work still required

1. **Audit first.** Check model equations, counted versus charge-equivalent
   PE, initial conditions, charge/current units, channel selection and source
   identity. Add an explicit distant-photon regression if useful: selected
   channels must equal the illuminated set, irrespective of their separation.
   Keep the adopted physics unless the audit finds a concrete defect.
2. **Add an explicit full-production mode/card and PE output schema.** `settings` currently checks
   `pilot_only` and limits the explicit ID list to 16. Keep the small pilot
   available and add a reproducible all-source-events selection. Write
   separate full-run results/report locations and publish the exact command.
   Save per-event SiPM–bucket PE values with the input's signed keys. Keep
   primary detected PE distinct from real-valued sensor and ASIC PE equivalents.
   Version the production HDF5 schema and preserve the current pilot for audit.
   `make sipm-pilot` currently remains a pilot command.
3. **Prepare scalable execution and recovery.** The code processes events
   individually but writes each complete stage to one temporary file before
   atomic rename; it currently has no per-event restart. Benchmark a modest
   batch for total runtime, peak memory and disk growth before scheduling the
   full run. Preserve final completion checks. If adding workers/shards,
   derive seeds from the same event/channel identities and avoid concurrent
   writers to a shared HDF5 file. Include provenance and completeness of any
   shards in the final manifest.
4. **Preserve the two stages and audit identities.** Changing only the ASIC
   card reuses the sensor stage. Per-event seeds do not depend on selection
   order, but `sipm_pipeline.jl` itself participates in the code hash: changing
   orchestration currently changes both model identities and sampled streams.
   Keep old pilot outputs for audit; use matching code/card/source identities
   when testing exact replay across pilot/full subsets. A recreated sensor
   file can differ in timing metadata and thus byte hash; stale electronics
   now rebuild automatically and bind to its current SHA-256.
5. **Extend summaries and presentation for the ensemble.** The existing
   presenter explicitly describes event 0 and the pilot; it is not a full-run
   analysis report. Calculate full-run statistics in Julia and save them to
   HDF5. Report clipped cell fraction with explicit denominator, event-wise
   distributions, fraction of events with clipping, sensor recovery losses,
   and an explicitly defined measure of integrated charge/current removed by
   electronics clipping. Distinguish all acquired buckets from photon-occupied
   buckets if reporting both. Do not treat clipped-cell fraction as charge
   loss. Keep sample and gas alternatives separate.
6. **Validate completion.** Check all 8000 IDs per file, source/pair/gas
   identities, photon preservation, detected PE bounds, charge bounds, signed
   alignment, zero-input behavior, clipping flags and stage provenance. Save
   numerical validation and a manifest, then generate the full TeX/PDF from
   saved results. The known substantial clipping should remain visible.

The AXEL inverse is currently a diagnostic helper and documented later
analysis. Full forward simulation does not require implementing energy
reconstruction, fitting recovery parameters or attempting to recover clipped
samples. Those remain separate analysis tasks.

## Workspace state

The SiPM work has not been committed. The shared checkout also contains
co-worker changes for ion CSA and photon-pulse analysis; leave those intact.
Large response HDF5 files are ignored by Git. There is no request to commit,
publish or push as part of this handoff.
