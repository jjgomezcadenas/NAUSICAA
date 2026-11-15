# New detector configuration: instructions for the other agent

Date: 2026-11-15. Branch: `new-config` (from `ITACA`). Commit `fecaf9e9f`.
Read `AGENTS.md` (section "Detector configuration") and
`notes/html/detector_configuration.html` first.

## What changed

- The detector is defined once, by the eight cards under `config/detector/`
  (gas, geometry, drift, voltages, el, sipm, ions, energy) and the index
  `config/detector/itaca.toml`. Design values only; derived quantities
  (potentials, yields, drift times, masses) are computed, never typed.
- `config/itaca.toml` is gone. Its content, including the coordinate
  convention you added (z = 0 at the pad plane, Geant4 offset), is in
  `config/detector/geometry.toml`; the offset is now derived from the length.
- A run is a card under `config/runs/` (schema `nausicaa.run.v1`). It names
  the detector set, changes detector values only through `[overrides]`
  (dotted keys such as `"gas.pressure_bar" = 15.0`), and holds the run
  parameters, the response-model cards (`[response]`), the present-chain
  acquisition (`[electronics]`) and the stage sections (`[hit_production]`,
  `[ions]`, `[electrons]`, `[track_analysis]`, `[hit_output]`, `[controls]`,
  `[zs_controls]`).
- The framed production scripts (`run_framed_detector_response.jl`,
  `run_hit_production.jl`, `run_framed_electronics.jl`,
  `run_framed_reconstruction.jl`) take run cards only. The cards under
  `config/readout/` no longer drive anything; do not edit them.
- Every stage writes `resolved.toml` beside its results (detector values
  after overrides, the overrides, the run card, SHA-256 of every input). A
  different configuration in an existing folder is an error.

## Your productions

Run cards that reproduce them exactly (nine overrides each: 15 bar, 293.15 K,
150 cm, 300 V/cm, 0.6 mm/us, 30 ms, diffusion 1.6 / 0.75 in card units,
15 cm/s):

| Card | Production |
|---|---|
| `config/runs/bb0nu_vh1_cf20f_full_6156.toml` | bb0nu, all events, VH1 hits of 6156 accepted events |
| `config/runs/bb0nu_vh1_cf20f_1000.toml` | bb0nu, VH1 hits of 1000 events, 3 keV ion cut |
| `config/runs/single_electrons_vh1_cf20f_full_5481.toml` | single electrons at the bb0nu energy |

Their outputs exist under `results/readout/`, so running these cards now
stops at "Output exists": that is the guard working. The pilot card
`config/runs/bb0nu_vh1_legacy_conditions_pilot100.toml` reproduced the
detector response and the hits of the bb0nu production bit for bit.

## How to run

```sh
make vh1-production RUN=config/runs/<card>.toml      # detector response, then VH1 hits
make test-itaca                                       # resolver, cards, settings (213 tests)
make itaca-detector                                   # default detector page and configuration page
julia --project=. scripts/documentation/calculate_itaca_conditions.jl config/runs/<card>.toml   # a run's derived values
```

To start a new run: copy a run card, change `run_label`, the event counts and
`output_subdir`, and put every changed detector value in `[overrides]`. A
misspelt override key or a physical key outside `[overrides]` is refused.

## Moving to the new default conditions

The default detector is 16 bar, 290 K, 360 cm diameter, 130 cm drift,
350 V/cm, Magboltz transport, 100 ms lifetime, 16.5 cm/s ions, 15 kV over
5 mm. A run at the default is a run card with an empty `[overrides]`. Do not
start it on your own: the user decides when production moves, and the hit
output of a new run label goes under a new `results/readout/<run_label>/`.

## Rules that still apply

- No production run, cluster job or long computation without the user's OK.
- No destructive step and no push without the user's OK; commit only on
  "commit".
- Do not edit `src/electronics`, `src/reconstruction`, `src/vh1`, the cards
  under `config/readout/` and `config/specs/`, or the detector cards under
  `config/detector/` (the user decides every design value; propose, then wait).
- Your uncommitted NEXUS files (`simulation/nexus/slurm/`) and the Magboltz
  folder are untouched in the tree; commit them yourself when the user says so.
- The Geant4 macro still carries 15 bar and 290 K; generating it from the
  detector set is planned, not done.
