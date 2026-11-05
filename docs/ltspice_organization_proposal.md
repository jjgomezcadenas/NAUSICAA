# LTspice installation, circuits and study outputs: organization proposal

Read-only audit, 5 November 2026. This document proposes changes; no application,
run directory, model, card or manuscript was moved or deleted. The running
production job and its `notes/physics/work/julia-depot` remain untouched.

## Recommendation

Separate the simulator installation from the repository's portable scientific
inputs. Within the repository, distinguish supplied Rogelio references, our
reproductions, our design studies, and frozen collaborator releases. Keep the
existing `scripts/electronics/` task structure; it already makes most of these
distinctions. Migrate paths incrementally with compatibility links and retain
original recorded hashes and manifests.

The principal space saving is removal of an unused extracted installer, not
removal of simulation results. The current CSA production uses the native
**LTspice 17.2.4** executable unpacked beneath the manuscript's scratch folder.
That executable must first acquire a stable installation location.

## What is installed and what actually runs

| Item | Observed location / identity | Role |
|---|---|---|
| Native simulator | `notes/physics/work/ltspice-review/legacy/LTspice.pkg/Payload/Applications/LTspice.app` | Active default for CSA production and pad-history viewer. App version 17.2.4; bundle about 87 MB. |
| Native auxiliary library | `legacy/LTspice.pkg/Payload/Library/Application Support/LTspice/` | About 41 MB of installed support files. Preserve with the application until required model/symbol lookup has been checked. |
| Native executable | `…/LTspice.app/Contents/MacOS/LTspice` | SHA-256 `4036f0d4a9a94471481d498cf4e06531c74f95bf4d29dd5010e1e7888c74b209`. This same hash is recorded in the v2 acquisitions and the active `csa_v3_response_calibration/analysis.h5` calibration. |
| Newer extracted package | `notes/physics/work/ltspice-review/package/Package.pkg/Payload/Applications/LTspice.app` | About 1.3 GB; Info.plist version 26.0.2.1. Existing `LTSpice/README.md` records that this Wine runner stalled even on a resistor test and was not used for reported calculations. No active executable caller points here. |
| Downloaded installers | `LTspice_26.pkg` / `LTspice_legacy.pkg` in the same scratch directory | About 369 MB / 96 MB. Installation media, not run inputs. |
| External installation | `/Applications/` and `~/Applications/` checked | No LTspice application found in either standard location. The audit environment has no `LTSPICE` variable set. |

The executable paths and versions above were inspected from files; neither
simulator was launched. `LTSpice/` at repository root is **scientific material**,
not the application installation.

### Proposed installation convention

Use a versioned user installation, for example
`~/Applications/NAUSICAA/LTspice-17.2.4.app`, outside Git. Keep the matching native
support library and one original installer in a documented, external software
archive. A normal versioned installation under `/Applications/` is also suitable;
the essential requirement is a pinned version and recorded binary identity.

Introduce one resolver with this precedence:

1. Explicit `--simulator PATH`.
2. `LTSPICE` environment variable.
3. An ignored machine-local `config/local/tools.toml` entry.
4. A clearly reported missing-tool error with the setup instructions.

Provide a tracked `config/tools.example.toml`. Machine paths belong there, not
in detector cards. A run records the resolved executable path, version and
SHA-256. Preserve the current calibration requirement that the simulator hash
matches its saved calibration. Moving an identical binary should not change
that identity; changing the binary requires the existing calibration checks.

Do not silently select a newly installed latest version. Do not copy only the
executable while discarding its app bundle or support files.

## Active dependency map

| Entry point / source | Current mechanism | Proposed change |
|---|---|---|
| `src/pad_history_viewer/engine.py` | Tries extracted native executable first, then `/Applications/LTspice.app`; explicit CLI override supported. | Use shared resolver; preserve explicit override. |
| `scripts/production/run_ion_csa_pilot.py` | Hardcoded extracted-native default. Checks simulator hash against saved calibration. | Use shared resolver; retain hash checks and all simulation settings. |
| `scripts/production/benchmark_ion_csa_pilot.py` | Same hardcoded default. Scratch under `notes/physics/work/`. | Shared resolver and repository `work/` destination. |
| `src/ion_signal/ion_csa_pipeline.jl` | Forwards optional simulator argument to Python acquisition. | Retain forwarding and include the resolved runtime identity in provenance. |
| `scripts/viewers/view_pad_histories.py`, history/compression commands | Pass simulator to the viewer engine/acquisition. | Inherit the same resolver instead of another fallback. |
| `scripts/electronics/reproduction/{initial,spice2,rogelio_v2,rogelio_v3}/` | Acquisition commands require `--simulator`. | Keep their CLI; permit shared resolution where useful. |
| Device, bias, feedback and history runners | Explicit `--simulator` passed through workflow subprocesses and Julia tuning. | Preserve forwarding; workflow entry point resolves once. |
| Make targets `ltspice-study`, `shaping-study`, `ion-shaping-study`, `spice2-reproduction`, `spice2-enc` | Pass `$(LTSPICE)` into acquisition commands. | Document/configure common resolver; an unset variable should not become an empty path. |
| `src/ion_signal/spice_acquisition.py` | Runs native batch `executable -b deck`, converts simulator logs to readable `.run.txt`, reads `.raw`. | Retain as shared acquisition implementation. Resolver can live alongside it. |
| `src/ion_signal/spice_storage.py` and `spice.py` | Preserve immutable deck text, save waveform groups and parse circuit/raw data. | Keep in `src/ion_signal/`; no imported command-script helpers are needed for this role. |

Only three active Python files contain the extracted-app fallback: the viewer
engine, CSA pilot runner and CSA benchmark runner. Other acquisition commands
mostly take explicit executable arguments.

The Julia depot is a separate runtime dependency. `work/julia-depot`,
`work/julia_depot` and `notes/physics/work/julia-depot` coexist. Several CSA and
history launchers include the last one in `JULIA_DEPOT_PATH`; it is approximately
198 MB and is in use by current production. Consolidate depots in a later,
separate runtime maintenance step, after running jobs finish.

## Scientific material: ownership and purpose

| Existing material | Classification | Treatment |
|---|---|---|
| Original files at `LTSpice/` root | Rogelio's initial schematic, model, supplied plots and database; README identifies their original commit `7766b7d`. | Preserve unchanged as supplied reference. |
| `LTSpice/Spice2/`, including `v2/` and `v3/` | Later supplied circuit/model versions and original documents/results. | Preserve received contents, including supplied `.raw`, `.op.raw`, `.tar`, PDFs and DOCX. These are source evidence. |
| `results/ltspice/spice2_reproduction`, `rogelio_v2`, `rogelio_v3`, `noise_measurement` | Our calculations reproducing or evaluating supplied circuits. | Label as our reproduction results, not received Rogelio data. Some v2 groups deliberately use matched capacitance/temperature or modified M3; describe these separately from the supplied-condition reproduction. |
| `results/ltspice/rogelio_v2_retuned` | Our retuned design comparison using a supplied topology. | Classify under design studies, preserving its original run identity. |
| `LTSpice/variants/` | Portable circuits/model corrections derived by our transistor scans. | Keep as derived design artifacts with source and model provenance. |
| `results/ltspice/{m3_scan,m4_scan,bias_sensitivity,feedback_optimization,retention_scan,retention_scan_extension,…}` | Our design and optimization studies. | Keep reduced results, cards, source snapshots and necessary inputs. |
| `results/ltspice/csa_v3_response_calibration` and production CSA outputs | Calibration and event-response production. | Active scientific dependencies; retain. |
| `LTSpice/paper_2026-10-28`, `optimized_2026-10-29`, `v3` | Frozen, portable releases for collaborator review; approximately 81, 78 and 53 MB. | Preserve their internal layout and manifests. Repeated model/PWL files inside these packages make them independently runnable. |
| `notes/physics/tools/package_*.py` | Commands that build those releases. | Move eventually to `scripts/paper/` or `scripts/electronics/design/packaging/`, retaining old command compatibility until references are updated. |
| Notes and generated reports | Authored explanations and reproducible presentation. | Authored notes stay in their agreed source locations; final reports stay under `latex_generated/ion_readout/`. |

Several input cards refer to **decks inside previous results**, not only to
`LTSpice/variants/`. For example, `config/retention_scan.toml` uses the selected
feedback-optimization deck, model and HDF5 analysis; the pad-history cards also
consume that selection. These directories must not be treated as disposable
scratch merely because they are called `results/`.

## Proposed final tree

```text
~/Applications/NAUSICAA/LTspice-17.2.4.app    # installed runtime, outside Git
<external software archive>/ltspice/17.2.4/ # installer + identity/source record

NAUSICAA/
  LTSpice/                                 # portable circuit material only
    reference/rogelio/
      initial/
      spice2/
      v2/
      v3/
    designs/variants/                      # our derived reusable circuits
    releases/
      paper_2026-10-28/
      optimized_2026-10-29/
      csa_v3/
  config/
    local/tools.toml                       # ignored machine paths
    tools.example.toml
    studies/electronics/
      reproduction/                        # source-condition and matched studies
      devices/
      bias/
      feedback/
      readout/
      histories/
    specs/response/ion_csa/                 # production response choices retained
    specs/runs/                            # production graphs retained
  src/ion_signal/                          # acquisition, parsing and physics
  scripts/electronics/                     # existing task grouping retained
  scripts/paper/                          # manuscript/release packaging commands
  results/ltspice/
    reproduction/<study>/
    design/<study>/
    calibration/<design>/
    histories/<study>/
  work/
    spice/<study>/<run>/                    # regenerable raw solver files
    cache/
    julia-depot/                           # later runtime consolidation
    latex/
    archive/failed_runs/<case>/             # retained diagnostics, marked invalid
  releases/collaborator_archives/           # optional external/local ZIP exports
  latex_generated/ion_readout/              # published study reports
```

This is a destination map, not a request to rename every path immediately.
Keep compatibility paths for existing cards and callers. The existing production
`config/specs` organization is already appropriate and should stay intact.
Do not edit archived TOMLs or manifests to make old paths look current. Map their
recorded paths through compatibility links and document the migration separately.

## Retain / relocate / possible deletion candidates

Sizes below are approximate local disk usage, not a promise of reclaimable space.

| Candidate | Size | Evidence and proposed action |
|---|---:|---|
| Extracted version-26 `package/` | 1.3 GB | No active executable caller; README records its failed evaluation. Save installer hash/version and failure note, then this extraction is the strongest deletion candidate after approval. |
| `LTspice_26.pkg` | 369 MB | Unused installation media. Delete or move to an external software archive after recording source/hash; no scientific input depends on it. |
| Extracted native `legacy/` | 129 MB | **Retain until migration is verified.** Contains active executable and support library. Remove the extraction only after the versioned installation works with the same binary hash and existing acquisitions. |
| `LTspice_legacy.pkg` | 96 MB | Keep one copy externally for reproducibility; need not live under manuscript notes or Git. |
| Failed `rogelio_v2_incomplete_title_error/` | 289 MB | Concrete defective netlists: their first line is `M1 …`; the corrected run starts with `* Derived v2 circuit; first line is the SPICE title`. Preserve a small failure record, representative bad/corrected decks, configuration, hashes and readable logs. Consider deleting failed bulk waveforms only after that record and input-duplication checks are complete. |
| Scratch `ltspice-review/original/` | 24 KB | `.asc` and `.inc` are byte-identical to the tracked supplied originals. Its `.log` is not matched; preserve that diagnostic if consolidating. |
| Scratch Python helpers | Tens of KB | They differ from current commands/modules; they are historical versions, not demonstrated exact duplicates. Preserve them in a small audit archive or identify their Git revision before removal. |
| `work/handoff` ZIPs | 61 MB combined | Full archive has 456 files; compact archive 72. All compact entries are byte-identical to entries in the full archive. Keep releases if they were shared or represent a useful checkpoint; move exports outside notes. Do not confuse this subset relation with equivalence to today's source tree. |
| Frozen portable releases | 212 MB combined | Intentional standalone snapshots with manifests; retain. Repeated PWL/model files are not automatically waste. |
| Study `.raw`, `.log`, `.db` | Varies; multiple GB across studies | Regenerable solver products are candidates only after proving the required full-resolution waveforms/logs are retained in HDF5 or another canonical result. Supplied Rogelio files and unreduced/unique acquisition data are exceptions. Never delete these extensions globally. |

The failed v2 directory contains roughly 128 MB of `.raw`, 116.5 MB of PWL
inputs and 43.3 MB of HDF5. All 42 PWL files (116.5 MiB) and all five model
`.inc` files are byte-identical to the corresponding files in the corrected
`results/ltspice/rogelio_v2/` run. None of the 53 `.cir` files is byte-identical.
Thus the repeated failed-run inputs can be replaced by recorded references while
retaining the small defective decks and diagnosis. Its HDF5 has the same config/source/model/simulator
identities as the corrected run but a different acquisition-code hash. That is
consistent with the observed title-line fix; it is not a valid substitute for
the corrected run. Preserve the distinction in the failure record.

The audited `results/ltspice/` tree contains about 13.31 GiB of logical file
contents across 339,825 files: approximately 5,555 MiB of HDF5, 4,680 MiB of PWL
and 2,966 MiB of raw solver outputs. These counts are a local snapshot, not a
retention policy. Many per-history inputs and numerical results support production
validation; cleanup should operate by complete study, with its dependencies and
manifest, rather than by extension or directory age.

## Migration sequence and acceptance checks

1. **Finish active production.** Record the resolved runtime/depot paths and
   prevent cleanup from racing jobs. Do not move their current directories.
2. **Record installation provenance.** Hash installer, executable and relevant
   bundle/support files; retain version and download/source record.
3. **Install the pinned native runtime outside the project.** Add the shared
   resolver and local tooling card. Check existing deck execution and waveform
   agreement in a dedicated migration smoke test, then keep the old path as a
   compatibility link during transition. No such test was run in this audit.
4. **Recover the easy space.** After approval, remove unused version-26
   extraction; archive or remove its installer. Consolidate the failed v2 run
   only after saving its failure evidence and verifying copied PWL inputs.
5. **Separate source/reference/release labels first.** Add indexes and manifests
   identifying received files, our derived variants and frozen releases. Move
   one folder family at a time; verify every file hash and preserve aliases.
6. **Move study cards with their workflow.** Update only active callers/defaults.
   Validate source/model resolution and preserve the original run-card snapshots
   and identities. Run a small representative acquisition only where a changed
   execution path needs verification.
7. **Consolidate runtime scratch last.** Move notes scratch into `work/`, update
   all depot/cache/viewer defaults, and check no active code still depends on the
   old location. Split incoming solver files from saved study results when each
   workflow is maintained; do not bulk restructure all historical HDF5 runs.
8. **Verify presentation and portability.** Check catalogue links and frozen
   release manifests; compile changed reports through the shared LaTeX builder.
   Keep the original manuscript and current `manuscripts/itaca3` work untouched.

## Provenance retained for each scientific run

Keep the exact circuit/netlist, model card, PWL input (or a demonstrably complete
reproduction record), physical/study/response cards, source-code snapshot or
hashes, simulator version/hash, solver settings, readable logs, acquired or
reduced numerical results and the validation summaries. Keep source-versus-derived
identity and unsuccessful-versus-valid run status explicit. Report-path migration
must not rewrite an old run's scientific identity.
