# Matched ion CSA pitch scan

The 3 mm and 2 mm cards reuse the current 8,000-event bb0nu primary realization,
with ion drift uniform in 5–145 cm and complementary electron drift in the
150 cm chamber. Electron readout is reused. Ion collection replays the same
primary pairs and diffusion seeds at the requested pitch and bucket width.

The common realization is
`results/simulations/production_pilot/bb0nu/4bd2a52bf1b95a0f4ce0348587da828bad322378e1cea713cfe8c9b746ef59f3/realisation.h5`.
Each ion collection verifies event IDs, initial pair counts, pair checksums and
both drift distances against that file and the existing matched 4 mm ion file.
All events are retained. The total ion count equals the initial pair count.

## Commands and outputs

Use `scripts/production/run_matched_ion_pitch.jl 3` (or `2`) to generate the
ion collection and detector index. Then use
`scripts/production/run_ion_csa_timed.jl` with
`config/specs/runs/ion_csa_nonlinear_matched_pitch3_full.toml` (or `pitch2`).
Set `JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia"` and run Julia with
`--startup-file=no --project=.` as for the other production scripts.

These run cards select production mode, the validated nonlinear model below
0.63 V, and SPICE fallback for histories outside its calibration domain.
They preserve signed response buckets and every pixel with positive input
charge. They request no SPICE controls or refinement histories.

Numerical outputs and timing ledgers live in
`results/ion_csa_nonlinear/matched_pitch3/` and `matched_pitch2/`.
The timing driver asserts that the actual filtered response files are absent
before running, records circuit/model/kernel/calibration hashes, and separates
covariance setup and final presentation from event work. With SPICE fallbacks,
event work includes the fallback launcher's orchestration and internal report.

`work/event-viewer-venv/bin/python scripts/production/present_matched_ion_csa.py 3`
(or `2`) publishes HTML, TeX and PDF from saved numerical outputs under
`latex_generated/ion_readout/ion_csa_nonlinear_pitch3/` (or `pitch2`).

## 3 mm run and publication recovery

The 3 mm collection and fresh CSA numerical stages completed on 2026-11-06.
Their machine-readable ledgers are `collection.toml` and the execution
`timing.toml` (also preserved as `numerical_timing.toml` beside the original
preflight); generated reports take numerical values from those files.
The original publication step encountered a Julia `sendfile` error after
copying the complete corrected HDF5. The temporary copy and corrected parent
had identical size and SHA256. Numerical arrays and response tests completed
before this error.

`scripts/production/recover_ion_csa_publication.jl` verifies the original code,
circuit inputs and byte identity, applies the existing production metadata,
validates the output, atomically renames it, and registers the saved study using
`config/specs/runs/ion_csa_nonlinear_matched_pitch3_register.toml`.
Recovery and registration time are reported separately from fresh numerical
work. The complete original wall time was not recovered; the report does not
substitute a cached rerun for that measurement. The recovery script changes no
library source or physics calculation.

The 2 mm run is prepared but held pending the requested repository-audit fixes
and the subsequent instruction to launch it.
