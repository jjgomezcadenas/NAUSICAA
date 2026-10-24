# Point-sample noise calculation

Reproduce from the repository root: `make noise-study`.

Inputs: `config/realistic_csa.toml`; signal: `arrival_cases_point_samples.json`.
The JSON records all input values, assumptions, hashes and numerical checks.

Frame interval: 0.01 s. Sampling aperture: 0 s. Noise band: 0.1–1e+06 Hz.
The sources are study assumptions; feedback uses the passive-resistor thermal spectrum.

| Source | Charge noise (e RMS) | Fraction of variance |
| --- | ---: | ---: |
| white_voltage | 4.148113 | 0.681130 |
| flicker_voltage | 1.886098 | 0.140818 |
| feedback | 2.120850 | 0.178053 |
| additional_current | 0.000000 | 0.000000 |
| total | 5.026153 | 1.000000 |

| Arrivals | Mean estimated charge per ion (e) | S/N per ion | N boundary for S/N = 5 | Minimum integer N for S/N > 5 |
| --- | ---: | ---: | ---: | ---: |
| uniform_per_ion | 0.374487 | 0.074508 | 67.107 | 68 |
| simultaneous_per_ion | 0.378216 | 0.075250 | 66.445 | 67 |

Checks: covariance quadrature tolerance 1e-7 versus 1e-9; independent positive spectral integral with Gauss–Legendre orders 8 and 16; variance addition. All relative discrepancies must be below 1e-6.
The ion thresholds scale each fixed arrival pattern. A different arrival history requires its own signal calculation.
