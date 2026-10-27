# Bias-voltage sensitivity

This study varies M3 and M4 gate voltages around the current paper circuit,
and M4 around the bias-retuned Rogelio v2 circuit. The transistor geometries,
other biases, input capacitance, temperature and readout remain fixed.
Explicit rounded settings change both tuned paper voltages together.

The generated report is `latex_generated/bias_sensitivity/report.tex`.
The test of practical usefulness includes gain and amplitude dependence as
well as ENC. Every point is calibrated with its own tail time and 100-ion gain;
the 1000-ion signal checks that calibration at a higher amplitude.

## Reproduce

```sh
export PYTHONPATH=src
export JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia"
python scripts/run_bias_sensitivity_pipeline.py \
  --julia /Users/jjgomezcadenas/.juliaup/bin/julia \
  --simulator notes/physics/work/ltspice-review/legacy/LTspice.pkg/Payload/Applications/LTspice.app/Contents/MacOS/LTspice
make test physics PYTHON=python
```

`config/bias_sensitivity.toml` supplies voltages and offsets. Julia writes
`points.h5`; Python acquires operating points, AC responses, spectra and
transients into one HDF5 file per point. Julia reduces these into `analysis.h5`.
`python scripts/present_bias_sensitivity.py` regenerates the figure and
standalone LaTeX report from saved results. No compilation is performed.

Refined cases are listed in `scripts/run_bias_sensitivity_pipeline.py`:
the rounded settings, both nominal circuits, paper M3 offsets of ±0.1 and
±0.5 mV, and M4 offsets of ±1 mV in both circuits. These cases halve the
maximum timestep, tighten current/voltage/relative tolerances, double spectral
density and check two quadrature orders. Charge is checked at every point.

Raw native simulator files are ignored; HDF5 preserves their numerical
content. All decks, PWL inputs and successful run logs are tracked. The
original circuit studies and manuscript are unchanged.
