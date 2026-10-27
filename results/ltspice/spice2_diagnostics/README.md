# Revised-card capacitance and injection diagnostics

From the repository root, with the pinned Python environment and Julia project:

```sh
PYTHONPATH=src python scripts/run_spice2_diagnostics.py --simulator /path/to/LTspice
julia --project=. scripts/analyze_spice2_diagnostics.jl
python scripts/present_spice2_diagnostics.py
```

Raw data and provenance: raw.h5. Numerical reduction: analysis.h5.
Generated standalone report: latex_generated/spice2_diagnostics/report.tex.
The report source is generated without compiling or changing the manuscript.

External capacitance plus intrinsic MOS capacitances; unchanged circuit and
revised NOIMOD=2 card; 20 C; 0.1 Hz–1 MHz stationary noise. Frame: uniform
finite-gap arrivals over 10 ms, 300 us centered averages, 0.3 s nominal
discharge correction. Injection: one elementary charge over 1 us plus 10 ns
edges; trapezoid rise and flat top both 300 us in the main comparison.

| External capacitance (pF) | Frame ENC (electrons RMS) |
|---:|---:|
| 0.75 | 232.3 |
| 1.25 | 295.9 |
| 3 | 521.7 |

At 1.25 pF, injection with the frame estimator: 292.1 e RMS; with the 300 us trapezoid: 203.6 e RMS.
