# Framed detector response

`run_framed_detector_response.jl` applies the MC fast trigger to bb0nu, then
calculates ion and SiPM-photon vectors for accepted events. The card is
`config/readout/bb0nu_detector_response_3mm_mainframe20.toml`. Numerical controls and selected
XY/XZ projections go to `results/readout/<run_label>/detector_response/`.
The current run uses a 20 mm frame margin on each side of the primary electrons and their delta-electron descendants.
Full event vectors remain in memory. TopMetal uses 3 mm pitch and 20 ms buckets
(3 mm in Z). SiPMs retain 4 mm pitch and 250 ns response buckets; event-display
projections sum 20 electron buckets into 3 mm in Z.

- `make test-detector-response`: frame, collection and fluctuation tests.
- `make detector-response-bb0nu PYTHON=work/event-viewer-venv/bin/python`: calculation and plots.
- `make detector-response-single-electrons PYTHON=work/event-viewer-venv/bin/python`: same detector settings, single-electron sample.
- `make detector-response-present PYTHON=work/event-viewer-venv/bin/python`: plots from saved controls.

Completed outputs are preserved; use a new run label for another calculation.
The main index is `notes/html/index.html`.

# Detector-production grid scan

The 8,000-event scan uses the 1 mm Geant4 hit tables for `bb0nu` and
`single_electrons`. The three resolved run cards are:

| Ion cells (mm) | Ion frame (ms) | Run card |
|---:|---:|---|
| 2 × 2 × 2 | 20 | `config/specs/generated/4725adea70b9/run.toml` |
| 3 × 3 × 3 | 30 | `config/specs/generated/df1414afcdc2/run.toml` |
| 4 × 4 × 4 | 40 | `config/specs/generated/91f0e47f63fc/run.toml` |

The electron card is identical in all three. At 3 mm the joint stage writes
the common pair realization, the ion history, and pure-Xe and Xe–15% He
incident-photon histories for each sample. The 2 and 4 mm stages replay the
source hits with the same pair, placement and ion-diffusion seeds; they verify
each event's pair checksum against the common realization and write ion
histories only.

Run the full pipeline from the repository root:

```sh
JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia" \
  PYTHON=python3.13 python3.13 scripts/production/run_grid_scan.py
```

To regenerate the scan index and report from the completed joint run and ion
files without repeating transport:

```sh
JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia" \
  PYTHON=python3.13 python3.13 scripts/production/run_grid_scan.py \
  --full-summary results/simulations/production_pilot/pilot_summary_238ddc5a75e0fd9b.h5
```

The numerical index is
`results/simulations/production_pilot/grid_scan_summary_238ddc5a75e0fd9b.h5`.
Python reads it to generate
`latex_generated/production_grid_scan/238ddc5a75e0fd9b/grid_scan.tex` and
its PDF. The HDF5 index contains the resolved branch-file paths, counts,
cards and provenance.

## Ion and photon pulse comparison

The saved sparse histories also feed an event-level pulse comparison. For
each event, Julia computes the mean over all occupied cells and the maximum
cell count. Ion cells are pixel–frame pairs at the three grid settings;
photon cells are SiPM–0.25 µs buckets for pure Xe and Xe–15% He. The ion
histories come from the grid-scan index and the photon histories from the
250 ns electron index. No count cut or
time integration is applied. Regenerate the numerical HDF5 reduction, browser
page, figures and PDF in one command:

```sh
JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia" \
  python3.13 scripts/production/run_readout_pulses.py
```

The outputs are `results/readout_pulses_8000/readout_pulses.h5`,
`results/readout_pulses_8000/readout_pulses.html`, and
`latex_generated/readout_pulses_8000/readout_pulses.tex`/`.pdf`.

## Electron photon histories at 250 ns

The electron-only replay uses the same 8,000-event primary-pair realisations
and detector settings as the grid scan, with `bucket_width_us = 0.25`
(0.25 mm at 1 mm/µs). It writes new pure-Xe and Xe–15% He photon files for
both samples and leaves the ion histories intact. Run from the repository
root:

```sh
JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia" \
  julia --project=. scripts/production/run_electron_250ns.jl \
  config/specs/runs/electron_250ns.toml
python3.13 scripts/production/present_electron_250ns.py \
  results/simulations/production_pilot/electron_250ns_summary_464425eac32a2c34.h5
```

The summary stores the four gas/sample paths, counts, resolved cards and
provenance. The generated report is
`latex_generated/electron_250ns/464425eac32a2c34/electron_250ns.pdf`.
The current SiPM input paths are listed in
`docs/sipm_electronics_handoff.md`. An optional final integer argument to
`run_electron_250ns.jl` limits events for a pilot and writes to `work/`.

## ELM baseline: full bb0nu Xe–He sample

`make elm-pitch4-production` runs the 8,000-event Xe–15% He baseline:
4 mm pitch, 2.5 mm holes/guides, 250 ns buckets and 12.5 ns SiPM recovery.
The photon runner selects gas branches through `selected_gases` in its run card.
Four independent Julia workers process disjoint 100-event response shards;
completed shards are validated and reused when the workflow resumes.

The controlling cards are `config/specs/runs/elm_photons_pitch4mm_ep2_full.toml`
and `config/specs/runs/elm_pitch4mm_bb0nu_full.toml`. The resulting
`results/elm_sipm_pitch4mm_production/files.toml` resolves the photon index,
response manifest and immutable histories. Numerical reduction streams quiet
noise statistics and saves event means, maxima, charge losses, clipping flags
and attachment-corrected energy widths. The report also gives the efficiency
and energy width after rejecting clipped events.

Regenerate presentation with `make elm-pitch4-production-report`. The browser
report is `results/elm_sipm_pitch4mm_production/elm_sipm_pitch4mm_production.html`;
the self-contained TeX/PDF is under
`latex_generated/electron_readout/elm_sipm_pitch4mm_production/`.
Use `JULIA="$HOME/.juliaup/bin/julia" PYTHON=work/event-viewer-venv/bin/python`
as Make overrides when those are the configured local environments.

### Framed electronics

`make electronics-bb0nu PYTHON=work/event-viewer-venv/bin/python` validates the
existing CSA calculations at 100 μs, runs the 100 accepted bb0nu events selected
by `config/readout/bb0nu_electronics_pilot_100.toml`, and presents their controls.
The input detector card and completed detector controls fix the selected events,
frames, placements and diffusion seeds. Ion and photon vectors are regenerated
in memory. Every pad and SiPM in the frame receives electronics response.

The public interfaces are in `src/electronics/ElectronicsResponse.jl` and
`notes/html/event_analysis/sim_electronics_response.html`. Existing circuit,
nonlinear calibration, SiPM and ASIC cards remain the numerical sources.
The 100-event card selects `electronics_100events_efficient/`. Outputs there are
`controls.h5`, `plot_controls.h5`, `index.html` and
`plots/` within the detector run directory. Use `make electronics-present` to
regenerate plots. Full event vectors remain in memory.

`make test-electronics` runs numerical interface tests. Peak currents above the calibration maximum are recorded. A voltage outside
the nonlinear tables stops the run. Production never launches SPICE; calibration
extension is documented in `notes/html/pending_work.html`.

CSA calibration is loaded once per run. Event-length tables, noise factors and
working arrays are reused. The linear response uses the exact exponential
update; the nonlinear model is unchanged. The report includes complete ion
processing time, one-time preparation, allocations and retained preparation.
The earlier 100-event comparison remains in `electronics_100events/`.

### Framed reconstruction and ZS controls

`make reconstruction-bb0nu PYTHON=work/event-viewer-venv/bin/python` runs the
reconstruction tests, validates the PE inverse, and processes the same 100
accepted bb0nu events through detector response, electronics and reconstruction.
The card is `config/readout/bb0nu_reconstruction_100.toml`.

`src/reconstruction/` provides separate ion and electron reconstruction functions.
Ions use measured ADC vectors, baseline subtraction, covariance-weighted linear
inversion and the direct nonlinear integral when selected by measured voltage. PE estimates
use the stationary mean correction and are rebinned by 20. ASIC clipping remains
marked. All estimates remain signed.

Outputs are `reconstruction/controls.h5`, `plot_controls.h5`, `index.html` and
`plots/` under the detector run directory. The HTML event selector shows signal
and noise distributions and threshold scans. Truth labels are used only for
these controls; no zero suppression is applied. Full event vectors stay in memory.
`make reconstruction-present` regenerates numerical summaries and presentation
without repeating detector or electronics simulation.

The direct inverse is checked with `scripts/validation/validate_direct_ion_inverse.jl`.
Its saved comparison is under `reconstruction/direct_inverse/`; regenerate that
HTML with `scripts/production/present_direct_ion_inverse.py`. The earlier
100-event ZS report used the iterative reference. Nonlinear noise-bias studies
remain deferred.

Framed SiPM runs now select fully recovered initial microcells through
`elm_pitch4mm_3015ps_tau12p5ns_recovered.toml`. Historical cards retain stationary
dark initialization. `scripts/validation/benchmark_sipm_initialisation.jl`
compares both settings on the same ten photon events, timing the complete SiPM
and ASIC calculation. Results are under `electronics_sipm_initialisation/`.
