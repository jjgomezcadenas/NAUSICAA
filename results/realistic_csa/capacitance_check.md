# Topmetal-S capacitance check

Reproduce: `make capacitance-check`.

Use the published 5 pF electrode capacitance as Cin, excluding Cf, with no additional transistor capacitance. All other circuit and noise inputs remain fixed.
Source: https://arxiv.org/pdf/2010.09226, Section 2.3 (1 mm electrode). This is a benchmark substitution, not a geometrical extraction for NAUSICA.

Point readings separated by 0.01 s; assumed noise band 0.1–1e+06 Hz; passive-resistor feedback noise.

| Cin (pF) | White voltage noise (e RMS) | Flicker (e RMS) | Feedback (e RMS) | Total (e RMS) |
| ---: | ---: | ---: | ---: | ---: |
| 0.25 | 1.458 | 0.278 | 5.177 | 5.385 |
| 5 | 4.148 | 1.886 | 2.121 | 5.026 |

| Cin (pF) | Arrivals | Recovered charge / input (%) | S/N for 1000 ions | Ions for S/N = 5 (rounded up) |
| ---: | --- | ---: | ---: | ---: |
| 0.25 | simultaneous_1000 | 92.32 | 171.15 | 30 |
| 0.25 | spread_1000 | 92.42 | 169.83 | 30 |
| 5 | simultaneous_1000 | 37.88 | 75.25 | 67 |
| 5 | spread_1000 | 37.84 | 74.51 | 68 |

Checks: halve the current time step; compare covariance and independent direct noise integration; verify addition of independent source variances. Detailed values and input hashes are in capacitance_check.json.
