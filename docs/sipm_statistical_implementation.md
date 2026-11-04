# Statistical SiPM response

The implementation is in `src/ion_signal/sipm_statistical.jl`, included by
`SiPMResponse`. Its derivation is in `docs/sipm_electronics_model.tex`, under
“Compact recovery-state and charge moments” and “Conditional-count statistical
implementation”. The generated model PDF is
`latex_generated/sipm_response/sipm_electronics_model.pdf`.

## Interface and state

```julia
StatisticalOptions(substeps=128, sparse_trigger_limit=16)
simulate_channel_statistical(rng, incident_photons, bucket_width_s, sensor;
                            options=StatisticalOptions())
```

One call processes one SiPM's consecutive acquisition buckets. Input photons
are integer counts. `Sensor` specifies PDE, effective microcell count, gain,
recovery time in seconds and total dark rate in Hz. The returned arrays contain
incident photons, detected PE, dark triggers, optical charge, dark charge and
total charge in coulombs. The input is preserved. The production adapter saves
the total charge as PE equivalents in the existing pixel/bucket format.

The compact state consists of the sum of microcell recovery fractions and the
sum of their squares. Detected counts are binomial; dark counts are Poisson.
At most 16 triggers in a bucket use individual arrival times and aggregate-state
charge draws. Denser buckets partition their counts into 128 intervals and use
finite-cell occupancy moments with midpoint recovery. One sampled charge
updates both output and remaining recovery. The second state moment preserves
its conditional mean loss. Bounded beta sampling approximates higher moments;
conditional first and second total-charge moments close in this state.

Initialization samples the stationary dark recovery once per microcell. The
subsequent state has two numbers per SiPM. Optical/dark charge partitioning
uses a bounded approximation to the squared-charge sum; production saves their
sum. Numerical diagnostics count state repairs, charge-variance limiting and
optical/dark-label variance limiting separately.

The existing assumptions are uniform microcell illumination, constant PDE,
exponential recharge, uniform arrivals within each input bucket and reset after
an avalanche. Crosstalk, afterpulsing and gain spread retain their card settings.
The TOML card `config/specs/response/sipm/sipm_statistical_v1.toml` selects this
model. Cards without `response_model` retain the microscopic implementation.

## Validation and reproducible commands

Use the repository's Julia environment and a Python environment with h5py and
Matplotlib. The Makefile accepts `JULIA`, `PYTHON` and `EVENT_JULIA_DEPOT` overrides.

```sh
make sipm-statistical-validation
make sipm-statistical-pilot
```

The first command runs unit and integration tests, Julia numerical comparisons,
and Python presentation. The second processes **two bb0nu events for each gas**
using the existing E/P = 3 photon inputs, then runs the existing electronics
stage. It saves the usual independent sensor and ASIC PE-per-bucket files and a
manifest under `results/sipm_statistical_pilot/`. Full production is a separate
run-card choice.

To regenerate presentation from saved numerical comparisons:

```sh
make sipm-statistical-report
```

The numerical benchmark is `results/sipm_statistical_validation/validation.h5`.
The browser report is `results/sipm_statistical_validation/validation.html`;
LaTeX, figures and PDF are under `latex_generated/sipm_statistical_validation/`.
The benchmark stores source snapshots, hashes and photon-source identities.

Synthetic tests compare fixed-count pulses, bright transients, sparse pulses
and dark-only histories at 32, 64 and 128 intervals per bucket. They measure
means, variances, adjacent-bucket covariance and total charge. Real-event checks
repeat the same saved photon event 256 times for each gas at E/P = 2.5 and 3.
An independent follow-up repeats the E/P = 2.5 pure-Xe event with 512 new seeds
and records detected-count variance alongside charge variance. Two additional
10,000-replica channel tests include binomial PDE fluctuations explicitly.
Those repetitions measure response fluctuations conditional on the photon
input; they do not measure the full energy resolution of the event ensemble.
Event timings use the median elapsed time per replica. Synthetic timings use
the ensemble elapsed time divided by its replica count.

Pilot comparison budgets are 0.05% for means and 5% for variances, or four Monte
Carlo standard errors if larger. These are numerical validation choices. The
report shows sampling errors so a passing finite-sample check is distinguishable
from a precise bound on model error.

The independent-seed check also uses the known binomial detected-count second
moment as a control variate. With fixed incident photons H, detected count N,
PDE p and event charge Q, its estimator is
`sample_variance(Q) - 0.5 * (mean((N - H*p)^2) - H*p*(1-p))`.
The subtracted term has expectation zero. The fixed coefficient 0.5 is shared
by both models; it reduces Monte Carlo scatter in their variance comparison.
The report retains raw variances and quotes uncertainty from the per-replica
estimator contributions.

## 100-event SiPM and ASIC production pilot

`make sipm-production-pilot` runs 100 bb0nu events per gas at E/P = 2.5 and 3
with the statistical sensor card and the existing ASIC card. It records SiPM
and ASIC computation/write times separately, then verifies all saved event and
cell identities, charge ledgers, clipping bounds and quiet-cell noise. Timings,
8,000-event estimates and output manifests are saved under
`results/sipm_statistical_production_pilot/`.

The driver requires fresh outputs for timing; cached histories remain available
through their manifests. `make sipm-production-pilot-report` regenerates the
[HTML report](../results/sipm_statistical_production_pilot/sipm_production_pilot.html)
and [PDF](../latex_generated/sipm_statistical_production_pilot/sipm_production_pilot.pdf)
from saved numerical results. Full production is a separate step.

The E/P = 2 extension uses `elcc_250ns_ep2.toml` with the user-specified
J = 0.030. `make electron-ep2-pilot` replays 100 photon events for both gases,
archives their files beside the common realisation, and applies the statistical
SiPM and existing 1 mA ASIC cards. It then compares them with the saved 2.5 and
3.0 pilots. The [clipping comparison](../results/sipm_statistical_production_pilot/clipping_comparison.html)
distinguishes lost charge, clipped photon-positive cells and events containing
clipping. Its before/after response spreads use expected-dark subtraction,
known-attachment correction and deposited-energy normalization, retaining
SiPM recovery nonlinearity. `make sipm-clipping-pilot-report` regenerates the
presentation from saved reductions.
