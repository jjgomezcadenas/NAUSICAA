# Bucket-duration comparison

Checkpoint before calculations: `1d45a60` on branch `LTSpice`.

`config/bucket_study.toml` specifies the bucket durations (seconds), averaging
width, nominal discharge correction, numerical steps and final circuit.
`10ms/` and `20ms/` are independent output directories. All historical studies
remain in their original paths. To change physical parameters for a later run,
choose a **new `output` path** in the configuration. Input generation refuses to
overwrite files; acquisition checks source hashes before reusing saved records.
New ion counts can be appended to a matching acquisition.

## Result

With the final M3 ×4 / M4 ×4 circuit, external input capacitance 1.25 pF and
300 µs centred averages:

| Bucket | ENC (e RMS) | Minimum ions for S/N > 5 |
| --- | ---: | ---: |
| 10 ms | 94.61 | 474 |
| 20 ms | 112.14 | 562 |

At 100 ions/ms, direct finite-charge transients give S/N = 10.56 and 17.79.
The nominal 0.3 s discharge correction leaves residual contributions in the
following empty bucket of 3.07% and 5.75% of the originating packet's signal.
Two consecutive 1000-ion packets give second-bucket excesses of 2.82% and 5.50%
relative to isolated packets. These are deterministic history effects; the ENC
comes from stationary noise about the quiescent operating point.

Authoritative numerical results and verification are in each `analysis.h5`.
The generated report is `latex_generated/bucket_scan/report.tex`; it contains
the circuit, boundary conditions, charge convention, estimator, numerical
convergence and threshold brackets. The paper has not been updated.

## Reproduction

From the repository root, with the pinned Python environment and Julia project:

```sh
export PYTHONPATH="$PWD/src"
export JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia"
export MPLCONFIGDIR="$PWD/notes/physics/work/matplotlib-cache"
python scripts/run_bucket_pipeline.py --julia /path/to/julia --simulator /path/to/LTspice
```

The pipeline creates input files on a fresh output path and resumes matching
saved acquisitions. The current presentation compares the configured 10/20 ms
pair. The input, acquisition and reduction functions take the bucket duration
from the configuration.

Individual steps:

```sh
# Once, for a fresh output directory:
julia --project=. scripts/build_bucket_inputs.jl
python scripts/run_bucket_study.py --simulator /path/to/LTspice
julia --project=. scripts/analyze_bucket_study.jl
# Check the integer threshold predicted by the reduction:
python scripts/run_bucket_study.py --simulator /path/to/LTspice --bucket-ms 10 --counts 473 474 475
python scripts/run_bucket_study.py --simulator /path/to/LTspice --bucket-ms 20 --counts 560 561 562
julia --project=. scripts/analyze_bucket_study.jl
python scripts/present_bucket_study.py
julia --project=. tests/test_bucket_study.jl
make test physics
```

The full pipeline chooses the threshold counts from the initial reduction.
`present_bucket_study.py` reads the reduction without rerunning physics.
LaTeX compilation is a separate, user-requested step.

## Data and validation

- `input.h5`: Julia-calculated mean current per ion, finite-gap charge and
  physical/configuration/source hashes. Bin integration conserves charge.
- `configuration.toml`: the exact configuration used for acquisition.
- `decks/`: PWL currents, native LTspice circuits, model and readable run logs.
- `raw.h5`: compressed baseline, isolated and consecutive-packet transients,
  both timestep refinements and original noise spectra at 200/400 points per
  decade. The circuit and model are compared with the saved noise run before
  reusing its spectra.
- `analysis.h5`: calibrated gains, filtered noise, ENC, finite-charge S/N,
  integer-threshold verification, history contribution and convergence checks.

Both input-bin width and maximum transient timestep are halved from 5 to
2.5 µs. The new 10 ms PWL agrees with the historical input to floating-point
precision, and its calibrated gain agrees within 0.001%. Tests independently
integrate the weighting field, check current scaling and finite-path charge,
verify exponential cancellation and compare white-noise filtering with the
analytic time-domain limit. The threshold brackets have numerical errors
below one tenth of their distances from S/N = 5.
