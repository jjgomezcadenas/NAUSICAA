# M4 geometry and bias scan

Reproduce from the repository root with the pinned Python environment and Julia project:

```sh
PYTHONPATH=src python scripts/run_m4_scan.py --simulator /path/to/LTspice
PYTHONPATH=src julia --project=. scripts/tune_m4_bias.jl /path/to/python /path/to/LTspice
julia --project=. scripts/analyze_m4_scan.jl
python scripts/present_m4_scan.py
```

The acquisition script resumes existing runs. Delete or move this study's raw.h5
before a completely fresh run. The original circuit and noise coefficients are
preserved. W and L of M4 are scaled together. The second scan adjusts only its
gate bias to recover the original quiescent output voltage within 10 uV.

Conditions: 1.25 pF external load plus intrinsic MOS capacitances; Cf=3 fF;
20 C; stationary NOIMOD=2 noise over 0.1 Hz–1 MHz. Uniform finite-gap ion
arrivals over 10 ms; centered 300 us averages with nominal 0.3 s discharge
correction. Injection is one elementary charge with a 1 us plateau and 10 ns
edges, filtered with a 300 us rise/300 us flat-top trapezoid. Schematic-level
simulation; parasitic layout and process corners remain to be evaluated.

| W (um) | L (um) | Gate (V) | Frame ENC (e RMS) | Gain (uV/ion) | Injection ENC (e RMS) |
|---:|---:|---:|---:|---:|---:|
| 0.72 | 0.36 | 0.620000 | 295.9 | 15.47 | 203.6 |
| 1.44 | 0.72 | 0.620000 | 133.7 | 0.27 | 85.7 |
| 1.44 | 0.72 | 0.585809 | 172.0 | 21.57 | 114.8 |
| 2.88 | 1.44 | 0.620000 | 126.9 | 0.20 | 80.3 |
| 2.88 | 1.44 | 0.556803 | 137.9 | 26.54 | 89.5 |
| 5.76 | 2.88 | 0.620000 | 126.5 | 0.19 | 79.9 |
| 5.76 | 2.88 | 0.541323 | 129.7 | 32.65 | 83.3 |
| 11.52 | 5.76 | 0.620000 | 126.4 | 0.18 | 79.9 |

Detailed results, including per-transistor noise, operating point, power,
10–90% rise time, finite-count linearity and convergence, are saved in analysis.h5.
Report source and figures: latex_generated/m4_scan/. The report has not been
compiled and the manuscript has not been changed.
