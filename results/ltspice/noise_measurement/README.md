# Noise of the voltage and of the frame difference

The v3 source snapshot includes Rogelio's new noise measurement. The circuit is
unchanged from the earlier transient study. The separate directory preserves
the prior input-charge study. The noise analysis uses 0.1 Hz–1 MHz and 200/400
points per decade. LTspice `INTEG V(onoise)` in a noise analysis gives the RMS
square-integral; this is checked against independent Julia quadrature.

Julia compares the raw voltage RMS, a 300 µs rectangular average, and the
discharge-corrected difference between readings 20 ms apart. The common decay
time and frame gain come from the saved rounded-bias paper circuit. Its 20 ms
uniform-ion frame calibration and the v3 single-pulse peak are different
measurements. No new v3 gain calibration is claimed here.

From the repository root, with the pinned Python environment and Julia project:

```sh
PYTHONPATH=src python scripts/run_noise_measurement.py --simulator /path/to/LTspice
julia --project=. scripts/analyze_noise_measurement.jl
python scripts/present_noise_measurement.py
julia --project=. tests/test_noise_measurement.jl
```

The first command acquires spectra. The next two constitute the analysis and
presentation pipeline; presentation alone reads analysis.h5. The standalone
report is `latex_generated/noise_measurement/report.tex`. Compilation remains
under the user's control.
