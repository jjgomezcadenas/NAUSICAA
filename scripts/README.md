# Reproduction commands

Run commands from the repository root with Python 3.11+ and the pinned scientific
dependencies, or Julia with `--project=.`. Python launchers also resolve imports
when started by absolute path from another directory. Existing Make targets keep
their names. Julia calculates new analyses and writes HDF5; Python presents saved
results.

| Directory | Purpose | Main entry points |
|---|---|---|
| `events/` | NEXUS input, hit analysis and MC selection | `extract_hits.jl`, `analyze_mc_selection.jl` |
| `transport/` | Diffusion coefficients, cloud studies and concurrency benchmarks | `calculate_electron_diffusion.jl`, `propagate_clouds_parallel.jl` |
| `collection/` | Direct carrier sampling and physical grid/time assignment | `sample_collect_parallel.jl`, `run_pixel_charge_direct.jl` |
| `production/` | Shared realizations, detector/electronics runs and matched reports | `run_elm_matched_production.py`, `run_matched_detector.jl`, `run_ion_csa.jl`, `run_matched_sipm.jl` |
| `electronics/reproduction/` | Initial, Spice2 and Rogelio v2/v3 circuit studies | Each version's acquisition, Julia reduction and presentation commands |
| `electronics/devices/` | M3/M4 geometry and bias tuning | `run_m3_scan.py`, `run_m4_scan.py` |
| `electronics/bias/` | Bias precision study | `run_bias_sensitivity_pipeline.py` |
| `electronics/readout/` | Frames, discharge, averaging and noise covariance | `run_bucket_pipeline.py`, `run_discharge_pipeline.py` |
| `electronics/feedback/` | Feedback-device optimisation | `run_feedback_pipeline.py` |
| `electronics/histories/` | Full histories, compression, input comparisons and diagnostic sequences | `run_pad_history_pipeline.py`, `run_pad_compression.py` |
| `electronics/design/` | ADC and pixel-area estimates | `analyze_adc_estimate.jl`, `build_pixel_area.jl` |
| `benchmarks/` | Analytical signal and circuit benchmarks | `make physics readout realistic-readout` |
| `viewers/` | Interactive application launchers | `view_events.py`, `serve_cloud_viewer.py`, `view_pad_histories.py` |
| `paper/` | Manuscript figures/macros, audits and collaborator packages | `build_spice_paper.py`, `audit_paper_macros.py` |
| `validation/` | Petit references and event-analysis checkpoints | `validate_event_step*.jl`, `generate_petit_*_reference.jl` |
| `documentation/` | HTML review-guide builder | `make event-review` |

## Current matched detector and electronics workflow

The completed production is 8,000 bb0nu events with shared placement in
5–145 cm, 4 mm ion pitch / 40 ms frames, and Xe–15% He at 4 mm EL pitch /
250 ns buckets. Its selected cards and manifests are documented in
[matched production](../docs/elm_matched_production.md).

```sh
make elm-matched-production PYTHON=work/event-viewer-venv/bin/python
make elm-matched-analysis PYTHON=work/event-viewer-venv/bin/python
make elm-matched-report PYTHON=work/event-viewer-venv/bin/python
make response-viewer
```

These commands respectively run the full workflow, reduce saved responses and
present them, regenerate presentation alone, and serve the production viewer
on port 8769. The current run is already complete; its manifest selects the
existing files. The nonlinear CSA approximation is new development; the saved
production uses the SPICE/linear response.

## Earlier direct ion and electron readout studies

```sh
make grid-readout
make present-grid-readout
```

The first command samples both carrier populations directly onto their physical
grids, validates the saved sparse histories and generates the presentation. The
second regenerates presentation from saved HDF5. `make pixel-charge-1mm`,
`make pixel-charge-1mm-3mm` and `make pixel-charge-1mm-4mm` run the ion-only
direct-collection studies. To present any saved pixel study, run
`python scripts/collection/present_pixel_charge.py --results-dir RESULTS --latex-dir REPORT`.

To inspect the saved ion and PE histories event by event, run
`work/event-viewer-venv/bin/python scripts/viewers/serve_cloud_viewer.py` and
open `http://127.0.0.1:8767/`. The viewer selects one sample and one electron
gas; event numbers are positions in the file, starting at 1.

## Original ion CSA pilot publication workflow

```sh
julia --project=. scripts/production/run_ion_csa.jl \
  config/specs/runs/ion_csa_v3_pilot_all.toml
```

The Julia entry point validates and publishes the six corrected 100-event
pilot files at 2, 3 and 4 mm pitch. It saves compact input-charged pixel rows,
complete bucket histories and signed reconstructed counts. The run manifest
locates the final matrices and covariance. Branches without a saved corrected
study execute the existing acquisition and reconstruction pipeline.
`config/specs/runs/ion_csa_v3_production.toml` selects the full samples for an
explicit production run. See [interfaces, schema and controls](../docs/ion_csa_response.md).

`make ion-csa-test` runs the response, reconstruction and production tests.
`make csa-design-note` regenerates the three-pilot technical note's numerical
summary, tables and figures. Compilation is a separate step.

## CSA response to stored pad histories

```sh
python scripts/electronics/histories/run_pad_history_pipeline.py \
  --config config/studies/electronics/histories/pad_history_readout_direct.toml --simulator /path/to/LTspice
python scripts/electronics/histories/run_pad_compression.py \
  --config config/studies/electronics/histories/pad_history_viewer_direct.toml --workers 4
python scripts/viewers/view_pad_histories.py
```

The complete history pipeline screens all pads, acquires selected nonlinear
transients, checks charge and resolution, and writes its report. The compression
command evaluates the screened candidate/control histories. To regenerate their
presentations:

```sh
python scripts/electronics/histories/present_pad_histories.py
python scripts/electronics/histories/present_pad_compression.py \
  --results results/ltspice/histories/pad_compression_direct/8b0d8b859f56
```

The viewer's calculation worker remains in `src/pad_history_viewer`; its reusable
LTspice runner is `src/ion_signal/spice_acquisition.py`. Existing source hashes
and cached numerical outputs remain attached to their original studies.

## Longer CSA retention

```sh
python scripts/electronics/feedback/run_retention_comparison.py \
  --simulator /path/to/LTspice
```

This runs the preserved feedback-length studies in
`config/studies/electronics/feedback/retention_scan.toml` and
`config/studies/electronics/feedback/retention_scan_extension.toml`, screens
all direct-collection pad histories and simulates the candidate/control pads.
Julia saves calibration, noise, validation and nonlinear statistics as HDF5;
Python acquires LTspice results and presents the saved numbers.

To regenerate the combined report source and figures:

```sh
julia --project=. scripts/electronics/feedback/combine_retention_scans.jl \
  config/studies/electronics/feedback/retention_comparison.toml
python scripts/electronics/feedback/present_retention_scan.py \
  --config config/studies/electronics/feedback/retention_comparison.toml
```

Results are in `results/ltspice/design/retention_{scan,scan_extension,comparison}`;
the combined report is `latex_generated/ion_readout/retention_comparison/report.tex`.
[Study configuration and validation](../notes/physics/retention_study.md).

## Documents and path history

Presentation commands write report source and figures. Compile reports explicitly
with `make pdf-grid-readout`, `make pdf-electron-diffusion`, or
`make pdf-mc-selection`. `make note` builds the earlier physics draft;
`make -C manuscripts/itaca3` builds the current EPJC manuscript.
Review-guide content lives in `docs/event_analysis_review/content/`; the readout
HTML template is `assets/readout/readout.html`. Included historical pixel batch
and cache functions live in `src/event_analysis/batch/`.

[Full inventory](../docs/scripts_inventory.md) ·
[Old-to-new paths and original hashes](../docs/script_paths.json)
