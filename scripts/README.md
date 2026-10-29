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

## Current ion and electron readout

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

## CSA response to stored pad histories

```sh
python scripts/electronics/histories/run_pad_history_pipeline.py \
  --config config/pad_history_readout_direct.toml --simulator /path/to/LTspice
python scripts/electronics/histories/run_pad_compression.py \
  --config config/pad_history_viewer_direct.toml --workers 4
python scripts/viewers/view_pad_histories.py
```

The complete history pipeline screens all pads, acquires selected nonlinear
transients, checks charge and resolution, and writes its report. The compression
command evaluates the screened candidate/control histories. To regenerate their
presentations:

```sh
python scripts/electronics/histories/present_pad_histories.py
python scripts/electronics/histories/present_pad_compression.py \
  --results results/ltspice/pad_compression_direct/8b0d8b859f56
```

The viewer's calculation worker remains in `src/pad_history_viewer`; its reusable
LTspice runner is `src/ion_signal/spice_acquisition.py`. Existing source hashes
and cached numerical outputs remain attached to their original studies.

## Documents and path history

Presentation commands write report source and figures. Compile reports explicitly
with `make pdf-grid-readout`, `make pdf-electron-diffusion`, or
`make pdf-mc-selection`. Compile the manuscript with `make note` when requested.
Review-guide content lives in `docs/event_analysis_review/content/`; the readout
HTML template is `assets/readout/readout.html`. Included historical pixel batch
and cache functions live in `src/event_analysis/batch/`.

[Full inventory](../docs/scripts_inventory.md) ·
[Old-to-new paths and original hashes](../docs/script_paths.json)
