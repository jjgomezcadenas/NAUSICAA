# Longer-retention CSA study

The study extends the frozen 180 nm feedback design to discharge times near
20–30 s. The user confirmed 20 and 40 ms frames. Both use centred 5 ms voltage
averages, 1.25 pF external input capacitance and a 3 fF feedback capacitor.

## Circuit and physical input

`config/retention_scan.toml` varies the feedback length from 400 to 2000 µm;
`config/retention_scan_extension.toml` adds 3000 and 4000 µm. Width is 0.36 µm
and gate voltage is zero. The M4 gate is adjusted to restore the reference
quiescent output, with four decimal places in volts. The remaining amplifier
geometry and biases, transistor card, 1.8 V supply and 20 °C temperature are
preserved. The large feedback-channel lengths are study parameters; the study
does not establish a physical layout for them.

The calibration uses the on-axis square-pad benchmark: a 0.5 mm square pad in
an infinite grounded plane, homogeneous dielectric, unit pad weighting
potential and zero at the other boundaries and infinity. Ions approach through
the finite 5 mm gap at 3 m/s. Uniform collection times fill each frame. The
integral per ion is e[1 − phi_w(g)], without renormalization or a separate
collection impulse.

## Calculation and files

- `src/ion_signal/retention_scan.jl` generates currents, tunes the rounded bias,
  reduces matched transients, integrates stationary noise, screens histories
  and gathers complete nonlinear statistics.
- `src/ion_signal/retention_validation.jl` measures zero-input numerical drift
  in volts and recovered-ion units.
- `src/ion_signal/retention_compression.jl` and
  `src/pad_history_viewer/retention_worker.jl` provide a study-specific numerical
  baseline bound. The original viewer keeps its existing kernel and caches.
- `src/ion_signal/retention_comparison.jl` combines the two preserved HDF5 studies.
- `scripts/electronics/feedback/run_retention_comparison.py` reproduces both
  studies and the combined report. The Python acquisition scripts preserve
  deck, model, configuration, simulator and source hashes.
- `scripts/electronics/feedback/present_retention_scan.py` reads HDF5 to produce
  figures, tables and the self-contained LaTeX report. It does not repeat the
  physical calculation.

Each scan saves `inputs.h5`, `operating_points.h5`, `analysis.h5`,
`baseline_validation.h5` and `presentation.h5` under its `results/ltspice`
directory. `source_snapshot/<sha256>/` retains calculation and acquisition
versions. Variant transients, PWL inputs, full screens and individual case
caches remain local and are excluded from git. Numerical HDF5 files stay local.
Decks, source snapshots, aggregate CSV tables and report sources are tracked.

## Numerical checks

Independent zero-input transients determine the baseline. Six-second tails
from 10 and 100 ions determine the small-signal exponential discharge time.
The fits are checked with two solver resolutions, three tail intervals and
two amplitudes. Fast transients check charge conservation and gain through
2000 ions. One-sided noise is integrated from 0.0001 Hz to 1 MHz through the
5 ms average and decay-corrected frame difference. Two spectral densities,
two quadrature orders and the saved 0.001 and 0.1 Hz lower-cut results check
convergence. Independent transistor variances reproduce the total variance.

Two variants have roughly 20 nV zero-current drift in the LTspice transient.
Independent controls bound its effect on frame charge. The study worker permits
100 nV only when this is below 0.005 ion; measured frame errors must be below
0.001 ion and drift below 0.01% of voltage noise RMS. The recorded numerical
values are retained. The default viewer's 1 nV check is unchanged.

## Stored histories and compression

Inputs are the corrected direct-collection histories:

| Geometry / frame | Directory |
|---|---|
| 2 mm / 20 ms | `results/pixel_charge_1mm_direct/` |
| 4 mm / 40 ms | `results/pixel_charge_1mm_4mm_direct/` |

Both contain 8000 double-beta and 8000 electron events. Every pad is screened
with the calibrated linear response, preserving zero buckets and total charge.
Histories whose predicted output falls below 0.35 V are simulated nonlinearly,
with 12 fixed-seed controls in each higher-voltage stratum. Refinement covers a
reproducible spread of cases. Stationary noise sets the deficit criterion;
no random noise is added to these mean transients.

The report counts an event once if a tested positive-count bucket loses at
least 25%, 50% or 75% of its nominal response and at least one frame ENC.
These are confirmed local-compression rates. Assessing event selection and
blob reconstruction requires the later track analysis. Each history includes
its full event interval and an additional 3 s tail.

An initial baseline check rejected two variants before simulation. Their
failure logs are preserved alongside the successful validated runs. One
isolated 800 µm case stalled at a numerical timestep near 1e−16 s. It was
rerun with Gear integration and checked with a finer grid. The source circuit,
input current and strict tolerances were preserved; the solver method is
recorded in the saved case. Summary
selection accepts only complete statistics, so partial and rejected attempts
cannot supply the report tables.

## Reproduction

From the repository root, using the pinned Python scientific environment:

```sh
python scripts/electronics/feedback/run_retention_comparison.py \
  --simulator /path/to/LTspice
```

Presentation only:

```sh
julia --project=. scripts/electronics/feedback/combine_retention_scans.jl \
  config/retention_comparison.toml
python scripts/electronics/feedback/present_retention_scan.py \
  --config config/retention_comparison.toml
```

The generated report is `latex_generated/retention_comparison/report.tex`.
`calibration.csv`, `compression.csv` and `noise_budget.csv` in the combined
results directory provide small, portable tables generated from HDF5.
`execution_sources.toml` records the source versions loaded by the long-running
batches; source snapshots also preserve file versions reported at write time.
Compilation is a separate user-requested step. Tests include
`tests/test_retention_scan.jl`, `tests/test_retention_worker.jl`, the existing
viewer checks and `make test physics`.

## Working design v3 — 2026-10-30

The user selected nominal 20 s discharge as the provisional working CSA design.
`config/csa_design_v3.toml` freezes `length_1200um`; the fitted discharge and
readout factors remain in its saved calibration. The portable handoff for
Vicente is `LTSpice/v3/`, with a matching copy in
`/Users/jjgomezcadenas/Dropbox-Personal/2026/NAUSICAA/v3/`.

Reproduce the package with `python notes/physics/tools/package_csa_v3.py`.
It includes electrically checked editable schematics, localized netlists,
model cards, PWL currents, all device dimensions and biases, operating point,
calibration, noise budget, compression statistics, report source and figures.
The manifest checksums identify every packaged source. PWL files remain local
and in the Dropbox handoff; the package's `inputs/` is excluded from git.
