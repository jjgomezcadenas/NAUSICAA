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
