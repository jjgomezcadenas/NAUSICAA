# LTspice runtime and study locations

The shared resolver selects the simulator from an explicit `--simulator`, the
`LTSPICE` environment variable, or `config/local/tools.toml`, in that order.
The local card is ignored by Git; `config/tools.example.toml` documents its
fields. `expected_sha256` pins the local-card executable. Existing acquisition
and calibration hash checks remain active.

```sh
PYTHONPATH=src work/event-viewer-venv/bin/python -m ion_signal.spice_runtime
make PYTHON=work/event-viewer-venv/bin/python ltspice-study
```

The first command reports the actual executable, version and SHA-256 without
launching LTspice. `make` resolves the same local configuration unless `LTSPICE`
is supplied explicitly. Commands that require `--simulator` continue to accept
an explicit executable path.

## Installation status

The currently selected, working simulator is native **17.2.4**, SHA-256
`4036f0d4a9a94471481d498cf4e06531c74f95bf4d29dd5010e1e7888c74b209`.
The selected executable is
`~/Applications/NAUSICAA/17.2.4/LTspice.app/Contents/MacOS/LTspice`.
The active calibration records the same hash. Keep the bundle named `LTspice.app`:
its resources loader requires that name.

The installer and support files are archived under
`~/Library/Application Support/NAUSICAA/ltspice/17.2.4/`.
`docs/ltspice_installation.json` records the verified copy.
A fresh calibration-deck comparison produced **60,156 bit-identical samples**
in time, V(in), V(out), and I(Is1) using the original and relocated applications;
see `docs/ltspice_runtime_validation.json`.
The original native extraction and duplicate installer have been removed after
verifying their external copies. Commands use the external application and the
canonical `work/julia-depot:$HOME/.julia` depot stack.

## Scientific files

| Purpose | Canonical location |
|---|---|
| Received Rogelio circuits/models/results | `LTSpice/reference/rogelio/{initial,spice2,v2,v3}/` |
| Our reusable circuit variants | `LTSpice/designs/variants/` |
| Frozen collaborator releases | `LTSpice/releases/{paper_2026-10-28,optimized_2026-10-29,csa_v3}/` |
| Study cards | `config/studies/electronics/{reproduction,devices,bias,feedback,readout,histories}/` |
| Production response/run cards | Existing `config/specs/response/ion_csa/` and `config/specs/runs/` |
| Saved reproduction results | `results/ltspice/reproduction/` |
| Saved design studies | `results/ltspice/design/` |
| Calibration | `results/ltspice/calibration/` |
| Saved history studies | `results/ltspice/histories/` |
| Release/manuscript commands | `scripts/paper/` |
| Generated reports | `latex_generated/ion_readout/` |

`docs/ltspice_paths.json` records each move and its verification.
Top-level `LTSpice/` and `latex_generated/` compatibility aliases have been
removed. Python and Julia read boundaries resolve historical card/manifest
paths using `src/ion_signal/project_paths.{py,jl}` and the report registry.
Moved cards and frozen release files retain their original bytes and hashes.
The pad-history resume check compares resolved paths plus content hashes and
retains its existing manifest. Nested study/result aliases remain for recorded
acquisition-deck paths.

## Cleanup performed

- Removed the unused 26.0.2.1 Wine extraction: approximately **1.218 GiB** of
  regular-file contents, with installer and executable identities recorded first.
- Archived the failed v2 acquisition under
  `work/archive/failed_runs/rogelio_v2_title_error/`, retaining its original path
  as a compatibility link.
- Replaced its 42 PWL files and five model files with links to byte-identical
  corrected-run inputs: **116.54 MiB** of duplicate contents removed.
- Retained all unique failed-run HDF5/raw data, defective decks and diagnostic
  logs. The failure manifest documents the missing SPICE title line.

These are logical file-content sizes; filesystem free-space changes may differ.
Details are in `docs/ltspice_cleanup.json` and
`docs/ltspice_failed_run_archive.json`. The version-26 installer was subsequently removed as unused; its hash and
failed-evaluation provenance remain recorded. Historical ZIP exports are
retained in `work/archive/notes_physics/handoff/`. The native installer and
support files remain in the external versioned archive.

## Validation

`make test physics PYTHON=work/event-viewer-venv/bin/python` passed 70 tests and
regenerated the physics benchmark without changing its tracked outputs. The simulator/storage tests and the separate recorded-path regression pass;
Julia pixel-area validation passes 23 checks. The script-layout check validates all 143 prior command moves,
shared acquisition/storage behavior, 130 Python commands and 14 CLI imports from a
separate working directory. Frozen inputs and moved cards are checksum-checked
against the relocation manifest.

The historical script-layout test was updated to recognize the existing,
configurable acquisition timeout: the default remains the original 120 seconds.
The acquisition implementation itself is unchanged.

## Completed root cleanup

`LTSpice/` now contains `README.md`, `reference/`, `designs/`, `releases/`,
`notes/` and `papers/`. Its 13 former top-level aliases have been removed.
`latex_generated/` contains eight subject directories, `index.html` and
`catalogue.json`; its 81 former top-level aliases have been removed.
The original root technical notes live in `LTSpice/notes/`; their PDFs are in
`latex_generated/ion_readout/rogelio_notes/` and build products in `work/latex/`.

Validation resolves 130 recorded paths, verifies 487 original file hashes and
13 production-card references, and confirms the original manuscript and
bibliography are unchanged. Presentation-only builds pass for ELM, the SiPM
model, MC selection, CSA design, Nygren and both Rogelio notes. CSA design keeps
its existing layout warning. No production data were recomputed by this cleanup.


The final tracked-move audit checked **64,840 files (3,791,193,056 bytes)**:
no file is missing and all numerical/scientific result files retain their Git
blob identity. The recorded differences are two rebuilt PDFs, one Markdown
link update and six relocated command sources; see
`docs/ltspice_tracked_moves_validation.json`.

`make note` was also attempted. Its presentation prerequisites run successfully,
but the untouched original `notes/physics/nausica.tex` fails at line 66 with an
unclosed math expression around the double-beta lifetime sentence. This source
was left unchanged. The independent generated-report builds listed above pass.

## Final scratch consolidation

No active Julia/LTspice process used the old project locations. The two Python
viewers and `work/event-viewer-venv/` were preserved. The global `~/.julia/`
installation, packages and caches were untouched.

| Scratch component | Before | After / destination |
|---|---:|---|
| `notes/physics/work/` | 885.0 MiB | 96.1 MiB of unique material in `work/archive/notes_physics/`; old directory removed |
| Canonical Julia depot | 236.1 MiB | Retained at `work/julia-depot/` |
| Failed v2 run | 171.6 MiB | 43.7 MiB at `work/archive/failed_runs/rogelio_v2_title_error/` |

Removed the native application/support extraction and duplicate native
installer only after verifying every external copy against its recorded hash.
The unused version-26 installer was also removed; its hash and failed-evaluation
record are preserved. The notes-side Julia depot contained only compiled caches,
so it was removed after validating the canonical depot's imports and HDF5 I/O.
All active command/test defaults now use `work/`; historical scratch paths are
resolved through the migration map. Unique draft text, logs, received snapshots
and historical ZIP exports are archived without content changes.

For the failed run, all **53 large raw waveforms** have every numerical array
and complete header exactly preserved in `raw.h5`. Those duplicate raw files
were removed after comparison; the **43.3 MiB HDF5**, all defective decks,
configuration/logs and **50 tiny companion raw outputs (57 KiB)** remain.
No unique failed-run numerical information was discarded. The full hash and
array comparison is `docs/ltspice_failed_raw_audit.json`.

This final pass reclaimed **965,912,767 bytes (921.2 MiB)**, including 4.3 MiB
of bit-identical migration-test outputs replaced by references to their retained
copies. Sizes describe logical file contents, not APFS free-space accounting.
The action inventory and before/after measurements are in
`docs/ltspice_scratch_cleanup.json`.

Post-cleanup Julia tests pass 18 response checks and 33 correction checks;
the canonical depot imports HDF5, Distributions, Graphs and DataFrames and
performs a successful HDF5 round trip. Python simulator/path/viewer tests pass
13 checks. The installed simulator retains its calibration-pinned hash. Both
viewers remain running; production HDF5 data have not been changed.

## Git contents after reorganization

Generated HDF5 results, current-waveform PWL files and compressed waveform CSVs
under `results/ltspice/` remain available locally and are excluded from Git.
The user approved removing their previous tracked locations on 2026-11-05.
Source decks, configuration cards, code, provenance, reports and papers remain
tracked, together with the portable inputs in `LTSpice/releases/`. Simulator
logs and build auxiliaries remain local. This changes tracking only; it does
not remove the numerical files or rewrite earlier Git history.
