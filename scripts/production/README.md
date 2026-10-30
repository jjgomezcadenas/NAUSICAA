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
JULIA_DEPOT_PATH="$PWD/work/julia_depot:$HOME/.julia" \
  PYTHON=python3.13 python3.13 scripts/production/run_grid_scan.py
```

To regenerate the scan index and report from the completed joint run and ion
files without repeating transport:

```sh
JULIA_DEPOT_PATH="$PWD/work/julia_depot:$HOME/.julia" \
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
photon cells are SiPM–2 µs buckets for pure Xe and Xe–15% He. No count cut or
time integration is applied. Regenerate the numerical HDF5 reduction, browser
page, figures and PDF in one command:

```sh
JULIA_DEPOT_PATH="$PWD/work/julia_depot:$HOME/.julia" \
  python3.13 scripts/production/run_readout_pulses.py
```

The outputs are `results/readout_pulses_8000/readout_pulses.h5`,
`results/readout_pulses_8000/readout_pulses.html`, and
`latex_generated/readout_pulses_8000/readout_pulses.tex`/`.pdf`.
