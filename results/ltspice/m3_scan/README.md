# M3 optimization with provisional M4 x4

From the repository root, using the pinned Python environment and Julia project:

```sh
PYTHONPATH=src python scripts/run_m3_scan.py --simulator /path/to/LTspice
PYTHONPATH=src julia --project=. scripts/tune_m3_bias.jl /path/to/python /path/to/LTspice
julia --project=. scripts/analyze_m3_scan.jl
python scripts/present_m3_scan.py
```

The additional threshold check for the provisional M3 x4 candidate is reproduced by:

```sh
PYTHONPATH=src python scripts/run_m3_scan.py --simulator /path/to/LTspice --scales 4 --gate 0.5384774789214134 --counts 473 474 475
julia --project=. scripts/analyze_m3_scan.jl
python scripts/present_m3_scan.py
```

Acquisition resumes existing runs. Move raw.h5 before a completely fresh run.
M4 is held at W=2.88 um, L=1.44 um, gate=0.556803417206 V. M3 W and L
are scaled together. Its gate restores the reference output DC voltage within
10 uV. All other dimensions, biases and noise coefficients remain unchanged.

Conditions: 1.25 pF external load plus intrinsic MOS capacitances; Cf=3 fF;
20 C; 1.8 V supply; stationary NOIMOD=2 noise over 0.1 Hz–1 MHz. Uniform
finite-gap arrivals over 10 ms; centered 300 us averages and nominal 0.3 s
discharge correction. Injection calibration: 100 electrons, checked against 10 electrons; 1 us plateau and
10 ns edges, with a 300 us rise/300 us flat-top trapezoid.

| M3 scale | W (um) | L (um) | Gate (V) | Frame ENC (e RMS) | Injection ENC (e RMS) | Gain (uV/ion) | Rise (us) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 50 | 0.5 | 0.6 | 137.88 | 89.49 | 26.54 | 5.22 |
| 2 | 100 | 1 | 0.5618814 | 103.55 | 62.76 | 26.82 | 5.16 |
| 4 | 200 | 2 | 0.5384775 | 94.61 | 55.40 | 27.02 | 5.12 |
| 8 | 400 | 4 | 0.5265441 | 92.40 | 53.46 | 27.11 | 5.11 |

Detailed results and convergence checks are in analysis.h5.
The one-electron injection tests showed numerical baseline variations; the
final injection calibration uses 100 electrons and agrees with 10 electrons
within 0.1 percent. Early timestep and injection_precise tolerance trials are
retained in raw.h5/decks for diagnosis and excluded from the final injection
ENC. scripts/diagnose_m3_injection.jl records these tests separately. The
strictest baseline trials did not converge within the acquisition timeout. Raw spectra,
transients, operating points, logs and source hashes are in raw.h5.
Standalone report and figure: latex_generated/m3_scan/.
Portable circuits: LTSpice/variants/m4_x4_m3_xN/.
The report source is generated without compiling or changing the manuscript.
