# VH1 readout with 15 ms frames: handoff

Date: 2026-11-17. Branch: `main`. For the agent that maintains `src/vh1`
(module, cards and SPICE study; see `docs/handoffs/vh1_production.md`).

## Why

The ion-only detector iItaca (pure xenon, 16 bar, 3 mm pads) is studied in two
configurations. Both keep a cubic 3 mm lattice: one frame spans 3.0 mm.

| Configuration | Drift | Field | Ion drift velocity | Frame |
|---|---|---|---|---|
| II_L1500_D350 | 1.5 m | 350 V/cm | 15.0 cm/s | 20 ms |
| II_L1000_D500 | 1.0 m | 500 V/cm | 20.0 cm/s | 15 ms |

The velocities are user design inputs (2026-11-17), close to the calculation
with the Neves mobility in pure xenon (13.9 and 19.8 cm/s). The 20 ms frame is
the validated reference. VH1 has not been run with 15 ms frames.

## Question

Does VH1 with 15 ms frames give a reconstructed ion signal as good as with
20 ms frames, for the same charge per cell?

## What to do

1. New card beside `config/specs/response/ion_csa/vh1_device_changes_cf20f.toml`,
   same circuit, 20 fF capacitors and noise spectrum, with `frame_ms = 15.0`,
   `conversions_per_frame = 2`, `conversion_window_ms = 7.5`. Do not modify the
   existing cards.
2. State whether the circuit can swap its capacitors every 15 ms with the
   present SPICE results (gain, output window, noise spectrum), or whether the
   SPICE study must be repeated. If it must, stop and report before running it.
3. New run card from `config/runs/bb0nu_itaca_pilot_10x100_ion15.toml` (same
   NEXUS files, same seeds, 3 mm pads) with overrides: ion drift velocity
   20.0 cm/s, ion bucket 15 ms, and the new VH1 card. The reference is the
   existing pilot `bb0nu_itaca_pilot_10x100_ion15` (15.0 cm/s, 20 ms); rerun it
   only if its results are not on disk.
4. The electron arm is not needed. Skip it if the pipeline allows; otherwise
   run it unchanged and ignore it.
5. Tests: the response matrix and covariance for 15 ms frames (conditioning,
   predicted ion noise), charge conservation, and the existing VH1 tests
   unchanged.

The existing gas (Xe/He 15 %) and 1.3 m drift of those files are irrelevant to
this check: the ion spread does not depend on the velocity, and the Z cell is
3.0 mm in both runs, so the ions per pad-bucket are statistically the same.

## What to report

For 15 ms against 20 ms, same events:

- Noise after reconstruction, ions per bucket (reference: about 33 with 20 fF).
- Signal-to-noise per pad-bucket, mean and distribution (reference: 6.8 at 3 mm pitch).
- Saturated pads and converter codes out of range (reference: none in 100 events).
- Residual RMS of the reconstruction.
- Event energy from the ions: mean, FWHM, and hits above the 5 keV cut per event.
- Condition number of the response matrix. The 1.67 ms ion pulse is 22 % of a
  7.5 ms window against 17 % of a 10 ms window.
- Leakage or low-frequency noise integrated per frame, if the spectrum makes a
  difference between 15 and 20 ms visible.

A results page beside the numerical outputs under `results/readout/`, linked
from the VH1 pilots index, with the comparison table and the distributions.

## Constraints

- No change of physical design parameters in `config/detector/`; everything
  goes through run-card overrides and a new response card.
- No Slurm job, commit or push without the user's explicit OK. The 100-event
  pilot runs on the Mac.
- Scratch files in `work/`.
- Do not start SPICE automatically.

## Not in scope

The iItaca detector cards, the ion-only run mode and the emulated energy are
written separately. A 30 ms frame is no longer needed.
