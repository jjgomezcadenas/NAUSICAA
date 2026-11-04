# Handoff: SiPM and electron-readout electronics

**Current authorization (2026-10-30):** the user has now authorized the
co-worker to audit this implementation and run all events. Start with
[`sipm_production_handoff.md`](sipm_production_handoff.md), which records the
approved output policy, current artifacts, audit checks and remaining
production-runner work. This supersedes the earlier pilot-only stopping
instruction below; the existing pilot command itself remains limited.

## Assignment and first review

The first SiPM/ASIC model has now been approved and implemented. The current
cards, independent saved SiPM stage, schemas and reproduction commands are in
[`sipm_response.md`](sipm_response.md); the physical model, sources, signal
recovery discussion and four-event pilot results are in
[`sipm_electronics_model.tex`](sipm_electronics_model.tex). The pilot uses event
0 from each input below. Full-sample response is now authorized after audit.
The scope and questions later in this handoff record the initial assignment;
the approved model note records the settled first-model choices.

Develop the response stage that consumes **incident VUV photons per SiPM and
0.25 µs (250 ns) bucket** and produces gas-specific SiPM/electronics histories. Work
directly with the user on the sensor and electronics model before a full
sample run. The first deliverable is a reviewable model: TOML cards, physical
assumptions, data types/interfaces, output schema, tests, and a pilot run.
The ion CSA branch has a separate handoff in `docs/ion_csa_handoff.md`.

Read, in order:

1. `AGENTS.md` for repository conventions and scientific boundaries.
2. This handoff, then `docs/detector_response_pipeline_plan.md` for the
   detector-to-response dependency graph.
3. `docs/production_formulas.tex` (or its compiled PDF) for the adopted
   electron diffusion, EL, optical collection and their provenance.
4. `config/specs/detector/electrons/elcc_250ns.toml` and the embedded cards
   in the photon HDF5 files for the actual production settings.
5. `src/event_analysis/production_readout.jl` for `PhotonSiPMHistory`,
   sparse-cell semantics and the HDF5 reader; `scripts/production/run_electron_250ns.jl`
   for the replay; `scripts/production/README.md` and
   `scripts/production/reduce_readout_pulses.jl` for the current 250 ns scan.
6. `docs/event_analysis_review/index.html`, Production detector model, for
   the reviewed detector-stage interfaces. The tracked numerical comparison
   is `results/readout_pulses_8000/readout_pulses.html` and its reduced HDF5.
7. `Project.toml` and `Makefile` for dependencies and test/report commands.
   `papers/ELCC/Axel-ELCC.pdf` gives ELCC/SiPM design context, and
   `papers/ELCC/NEXT-EL-Xe-He.pdf` supplies the EL-yield fit already used
   upstream. Cite primary sources through `notes/physics/references.bib`.
8. `src/event_analysis/EventAnalysis.jl` for the exported photon reader and
   `tests/event_analysis/test_production_readout.jl` for its HDF5 round trip,
   sparse-cell ordering and ledger checks. Start the response interface from
   `read_photon_history`, then preserve its event, gas and signed bucket keys.
9. `scripts/production/specs.jl` for canonical card hashing and event-keyed
   RNG seeds. `src/ion_signal/ion_csa_response.jl` and
   `scripts/production/run_ion_csa_pilot.py` show the current response-stage
   provenance and pilot orchestration. Their ion waveform physics belongs to
   the separate CSA branch.

`docs/decisions.md` records the development history. Its 2026-10-29
direct-PE proposal was superseded by the production detector/response split:
the current photon HDF5 files store incident photons; SiPM PDE and electronics
belong to this response stage. Use the current HDF5 schema and pipeline plan
when implementing the interface.

## Input files available locally

The four current 250 ns files are under
`results/simulations/production_pilot/` (large data, ignored by Git):

| Sample | Realisation and photon directory | Files |
|---|---|---|
| `bb0nu` | `bb0nu/c25084b9e721e407c0688bfdba837f8f7fe5e453d8dad73034681d2d994ad95c/detector/photons/56c68bc327d71dccc426312f4c157054d3b3bf368ad05be966ee93b1ee1777fe/` | `photons_pure_xe.h5`, `photons_xe_he15.h5` |
| `single_electrons` | `single_electrons/a82036783cb95a9651bf8ca7bc87056fa74f54cd8869d9a236614d41b6f590c7/detector/photons/a9db39ffb9369e26c005d1de51dbbdc76159b07860ca835700e6f7c87d318a41/` | `photons_pure_xe.h5`, `photons_xe_he15.h5` |

The local machine-readable index is
`results/simulations/production_pilot/electron_250ns_summary_464425eac32a2c34.h5`.
It resolves absolute paths at `/<sample>/pure_xe_path` and
`/<sample>/xe_he15_path`, and records
the event count and aggregate incident-photon total. The two gas files for a
sample use the same primary-pair realisation and event placement; they are
alternative gas models, not two simultaneous readout arms. The ion files
share each event's primary-pair draw, identified by `pair_checksum`.

| Sample | Gas | Events | Total incident photons | File size |
|---|---|---:|---:|---:|
| `bb0nu` | Pure Xe | 8,000 | 381,689,031,793 | 226,681,314 bytes |
| `bb0nu` | Xe–15% He | 8,000 | 375,318,926,831 | 166,781,897 bytes |
| `single_electrons` | Pure Xe | 8,000 | 381,740,368,430 | 246,802,508 bytes |
| `single_electrons` | Xe–15% He | 8,000 | 375,369,610,548 | 176,463,555 bytes |

## Current photon-pulse comparison

The updated [pulse comparison](../results/readout_pulses_8000/readout_pulses.html)
has a [numerical HDF5 source](../results/readout_pulses_8000/readout_pulses.h5)
with schema
`nausicaa.readout_pulses.v2`. Reproduce both, including the generated PDF,
with `python3.13 scripts/production/run_readout_pulses.py` from the repository
root (using the Julia depot specified in `scripts/production/README.md`). The
reducer combines the unchanged ion grid-scan index with the 250 ns electron
index above and records both index hashes.

Each positive sparse `count` is **one incident-photon count for one SiPM in one
250 ns bucket**. For each event, the plotted mean is the total count divided
by its number of occupied SiPM–bucket cells; the plotted maximum is its
largest single SiPM–bucket count. Each histogram contains one event value per
event, for all 8,000 events in that sample and gas model. There is no count
threshold and no integration over a SiPM's full event. The saved values are
incident photons before PDE, PE sampling and electronics.

Use the generated table and saved event arrays for occupancy tails and
electronics dynamic-range studies.

## HDF5 input contract

The file root has `schema = "nausicaa.photon_sipm_file.pilot.v1"`. Require
`metadata/complete = true`, `metadata/event_count = 8000`, the expected
`metadata/sample` and `metadata/branch` (`pure_xe` or `xe_he15`),
`metadata/signal_unit = "incident_photons"`, and
`metadata/time_unit = "us"`. The metadata also holds the exact realisation
and electron-detector TOML, source path/SHA-256, card and code hashes, seed and
stream identity. Identify events by their HDF5 group name and saved event ID.

Each `/events/<event_id>` group has schema
`nausicaa.photon_sipm_history.v1` and contains:

| Dataset | Meaning |
|---|---|
| `pixel_ix`, `pixel_iy`, `bucket` | Signed `Int64` SiPM grid and time-bucket indices. The four sparse vectors are sorted by the full `(pixel_ix, pixel_iy, bucket)` key; zero-count cells are omitted. |
| `count` | Positive `Int64` **VUV photons incident on that SiPM in that bucket**, before PDE and sensor/electronics response. |
| `gas_model` | `pure_xe` or `xe_he15`; keep this label in every response output. |
| `initial_pairs`, `surviving_electrons`, `emitted_photons` | Shared primary count, post-attachment electrons and pre-optical EL photons. Check `surviving_electrons ≤ initial_pairs` and `sum(count) ≤ emitted_photons`. |
| `event_id`, `pair_checksum` | Source event identity and checksum of its shared primary-pair draw. |
| `pitch_mm`, `bucket_width_mm`, `reference_mm`, `phase_mm`, `drift_direction` | Grid geometry and event phase in source-coordinate mm; current values are 10 mm pitch, 0.25 mm longitudinal width and direction `-1`. |
| `ion_drift_cm`, `electron_drift_cm` | Complementary lengths summing to 150 cm. Preserve the electron length for later attachment correction/energy analysis. |

Electron drift speed is 1 mm/µs, so each longitudinal bucket represents
0.25 µs. For electrons, `reference_mm[2]` is the maximum source-hit z; the
bucket index increases with `(reference_z - diffused_z)/0.25 mm`. Diffusion can
produce negative indices. The time axis is event-relative; the files do not
provide an absolute acquisition clock or photon emission times within the
5 mm EL hole. The current detector model assigns photons to the parent
electron's arrival bucket. Preserve signed indices and gaps as zero-photon
buckets when making a time sequence. `pixel_ix` and `pixel_iy` are SiPM
indices despite the generic sparse-cell names.

Use `EventAnalysis.read_photon_history(group)` in Julia or read the HDF5
datasets directly. A minimal Python check from the repository root is:

```python
from pathlib import Path
import h5py

index_path = Path("results/simulations/production_pilot/electron_250ns_summary_464425eac32a2c34.h5")
with h5py.File(index_path) as index:
    path = Path(index["bb0nu/pure_xe_path"][()].decode())
with h5py.File(path) as source:
    assert source.attrs["schema"] == b"nausicaa.photon_sipm_file.pilot.v1"
    assert bool(source["metadata/complete"][()])
    event = source["events/0"]
    photons = event["count"][:]
    assert int(photons.sum()) <= int(event["emitted_photons"][()])
```

## Settled detector inputs and physical boundary

The user-specified design is 15 bar, with pure Xe and Xe–15% He alternatives;
the electron drift is complementary to the ion drift in a 1.5 m chamber.
These files use the 1 mm Geant4 samples, 8,000 events each, common
`W = 21.9 eV`, `F = 0.15`, event seed 2026 and at least 10 cm of ion drift;
the resolved realisation card is embedded in each file. Electron attachment
uses a 15 ms lifetime in the detector stage. Keep the saved electron drift
length for a later energy estimator's attachment correction; the SiPM
response consumes the already attenuated photon history.
The EL grid has 5 mm diameter holes, 10 mm pitch and 5 mm depth, operated at
`E/P = 3 kV/cm/bar`. The source electron drift/attachment, gas-dependent EL
yield and `J = 0.01` at that field (a user-supplied approximate reading of
Oliveira et al., Figure 6), and analytic direct optical acceptance
have already been applied to the photon files. The optical model assumes
absorbing hole walls and a SiPM disk covering the exit. The SiPM response
starts with the saved incident-photon count; it must not resample EL or
geometric acceptance. The source for EL yield is
`papers/ELCC/NEXT-EL-Xe-He.pdf` (Table 1), with primary citations in
`notes/physics/references.bib`; the formulas and approximations are written
in `docs/production_formulas.tex`. The SiPM PDE baseline of **0.30** is a
user design input, recorded in `config/electron_transport.toml`; the AXEL
ELCC paper establishes a VUV-sensitive SiPM design but does not supply this
numerical PDE. The current
production electron card contains detector geometry and photon production,
and the new SiPM response card should own PDE.

The initial response model should consider a conditional PE draw for each
incident-photon cell, `N_PE | N_incident ~ Binomial(N_incident, PDE)`. This
preserves integer PE counts and the binomial detection variance while using
the already fluctuated photon input. Review this interface with the user
before adding additional sensor effects. Sample at the count level per
SiPM/bucket; the photon totals make individual-photon loops impractical.

## Cards, code and output graph

The agreed stage-2 card locations are:

```text
config/specs/response/sipm/sipm_v1.toml
config/specs/response/sipm_electronics/readout_v1.toml
config/specs/response/ion_csa/csa_v1.toml
config/specs/runs/<final_detector_response_name>.toml
```

The new SiPM agent owns the first **two** cards and the SiPM/electronics
response implementation. The ion CSA agent owns the third. A final-run card
will select one detector configuration and the three independent response
cards when both branches are ready. Current `config/specs/` contains the
detector cards, resolved run cards and the ion CSA pilot card. The SiPM
response cards and their parsers still need to be created.

`sipm_v1.toml` should state the PDE and its source/status, photon-to-PE
sampling policy and RNG seed, gas applicability, and each additional sensor
effect adopted after review (gain and its variation, dark counts, crosstalk,
afterpulsing, finite microcells/recovery, if included). The electronics card
should state the readout transfer/shaping, sampling/bucket convention, noise
model and its numerical source, dynamic range, and RNG policy. Use explicit
units in keys. Identify which quantities are measured, taken from a paper,
or user design inputs. Do not assign numerical values to effects whose
model and parameters have not been reviewed with the user.

Keep physical calculations and response types in `src/ion_signal`; put
configuration resolution and run orchestration in scripts. Record the exact
parsed cards, upstream photon-file content/identity, code hash and seeds in
each output. Changing PDE or SiPM parameters should reuse photon histories;
changing only the following electronics should reuse the SiPM-stage
realisation if it is cached. Propose whether to save an intermediate PE HDF5
and its cache identity before choosing the final output schema.

The planned combined electron-response destination is, for each sample and
gas model:

```text
results/simulations/production_pilot/<sample>/<realisation_id>/
  response/electrons/<photon_output_id>/<sipm_response_id>/electron_response.h5
```

For a combined response, `<sipm_response_id>` must depend on both the SiPM
and following-electronics cards, their relevant upstream identities, RNG
settings and model/code versions. If PE is cached as a separate stage, give
that cache an identity based only on its actual SiPM inputs so changing the
electronics can reuse it.

The response must retain event, gas, SiPM and bucket identity; incident and
detected count ledgers; sensor response; electronics noise/response and
provenance. Write versioned HDF5 atomically and mark it complete only after
validation. The final graph under `runs/<final_graph_id>/manifest.toml`
references the gas-specific response files and the ion response. Separate
gas results remain separate files. Numerical summaries belong in `results/`;
Python presents figures/tables from saved HDF5 and writes a self-contained
report in `latex_generated/` that compiles to PDF.

## Review questions and checks

The production cards do not yet specify a SiPM device, gain distribution,
dark-count rate, optical crosstalk, afterpulsing, cell count/recovery or
following-electronics noise/transfer function. Establish the first response
model and its parameter provenance with the user. Also decide the readout
window needed for dark counts and noise in zero-photon buckets, and whether
an intermediate PE cache is useful. Keep the accepted 0.25 µs incident-photon
bucket and parent-arrival timing approximation as the detector input.

Minimum tests: conditional binomial moments and PDE limits; deterministic
replay independent of event order/workers; per-event and aggregate
`N_PE ≤ N_incident` before sensor multiplication; correct zero-input and
zero-photon-bucket behaviour; signed bucket alignment; gas and pair-checksum
identity; HDF5 schema/round trip; and, for every added sensor/electronics
effect, a numerical or physical validation of its response, noise and
scaling. Run a small pilot before the four 8,000-event responses. Follow
`AGENTS.md`: Julia calculates and saves HDF5, Python presents from that
HDF5, and generated TeX/PDF gathers the parameters, assumptions and checks.
Run the relevant tests plus `make test physics`; compile LaTeX and inspect
the build log, without a visual PDF inspection unless requested.

## Verified local input checks

On 2026-10-30, all four photon files had `complete = true` and 8,000 event
groups. Every event passed the electron/photon ledger, positive sparse-cell
counts and equal sparse-column lengths. For each sample, the two gas files
had identical event-ID sets and primary-pair checksums. Their summed
incident-photon counts matched the 250 ns index totals above. Every event
also matched its earlier 2 µs counterpart in primary pairs, attachment
survivors, drift lengths, grid reference/phase and pair checksum.
