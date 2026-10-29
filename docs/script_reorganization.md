# Scripts directory migration

The migration on `reorg` groups commands by scientific task using the approved
inventory. It relocates 143 files: commands, four review-content modules, two
included batch/cache helpers, and the readout HTML template. The resulting
`scripts/` tree contains 138 Python and Julia commands, including the two
input-comparison commands that were already grouped.

Use [scripts/README.md](../scripts/README.md) for workflow entry points and
[script_paths.json](script_paths.json) for every old/new path and the original
SHA-256 at checkpoint `a12a91aa`. Historical results and frozen collaborator
packages retain their recorded provenance. New runs record their current source
paths and hashes.

## Shared code and callers

`src/ion_signal/spice_acquisition.py` owns native LTspice acquisition, readable
log capture and raw parsing. Its `run` function was extracted unchanged from
the original study. `spice_storage.py` holds the unchanged v2 text-preservation,
hash and HDF5 waveform-storage functions. The viewer imports acquisition from
this shared module. Study-specific feedback helpers remain together in their
workflow directory.

The historical Julia pixel batch and diffusion-cache functions are in
`src/event_analysis/batch/`. Their including command supplies the same globals
and EventAnalysis interface as before. The current direct collection and both
carrier-population models keep their existing algorithms and HDF5 schemas.

Make targets, subprocess dispatch, Julia includes, Python imports, active
documentation and packaging paths use the new locations. Python commands resolve
the repository and shared modules when launched from another working directory.
Review content is under `docs/event_analysis_review/content/`; the analytical
readout HTML template is `assets/readout/readout.html`.

Presentation writes figures and LaTeX source. Explicit `pdf-grid-readout`,
`pdf-electron-diffusion` and `pdf-mc-selection` Make targets compile those reports.

## Validation

- `make test physics`: 52 tests pass; the physics benchmark runs through its moved command.
- Event-analysis tests: 2,151 checks pass, including the relocated diffusion cache.
- Complete CSA history pipeline: 48 checks pass using the saved matched-baseline acquisitions.
- `make test-script-layout`: all 143 destinations, extracted helper syntax trees,
  Python parsing, active path references, and 14 CLI imports from another directory.
- Fresh LTspice runs through the original and shared runner produce identical waveform arrays.
- Direct collection smoke run: all six sample/carrier workers complete; ion accounting,
  grid phase, bucket width and PE-yield validation pass.
- Restarted CSA viewer: 33 direct-input, refinement and live API checks pass at its existing URL.
- `make event-review`: HTML guide regenerates from the relocated content modules.

Julia syntax checks cover all 54 command files. The migrated source was checked
without compiling LaTeX. Scratch acquisitions, smoke readouts and comparison
files are under `notes/physics/work/`.

The optional fresh acquisition check is reproducible with:

```sh
PYTHONPATH=src python tests/script_layout/check_acquisition.py --simulator /path/to/LTspice
```
