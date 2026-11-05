# Rogelio v2 assessment

Separate assessment of the circuit supplied on 27 October 2026 in
`LTSpice/Spice2/v2/180nm_professional/`. Original sources and the paper's
existing results are preserved. The generated report is
`latex_generated/rogelio_v2/report.tex`.

## Reproduce

From the repository root, using the scientific Python environment:

```sh
export PYTHONPATH=src
export JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia"
python scripts/run_rogelio_v2_pipeline.py \
  --julia /Users/jjgomezcadenas/.juliaup/bin/julia \
  --simulator notes/physics/work/ltspice-review/legacy/LTspice.pkg/Payload/Applications/LTspice.app/Contents/MacOS/LTspice
julia --project=. tests/test_rogelio_v2.jl
make test physics PYTHON=python
```

Presentation alone: `python scripts/present_rogelio_v2.py`.
The pipeline generates the LaTeX report without compiling it.

## Contents and interpretation

- `raw.h5`: native supplied and reproduced noise spectra; corrected native,
  matched-condition and enlarged-M3 operating points, spectra and transients.
- `analysis.h5`: Julia finite-path charge checks, calibrated decay, filtered
  noise, gain, ENC, nonlinear signal counts and adjacent threshold tests.
- `diagnostics.h5` and `diagnostic_analysis.h5`: voltage-ratio AC runs at
  unchanged bias, and M4 conductance from the operating-point log. The ratio
  includes the physical feedback loading; it is not an unloaded output gain.
- `solver.h5`: short zero-input tests at relaxed/default solver tolerances.
- `decks/`: exact simulation decks, input waveforms, model copies and successful
  run logs. HDF5 retains the numerical outputs; native `.raw` files are ignored.
- `../rogelio_v2_retuned/`: separate adjustment of only the v2 M4 gate voltage
  to restore the paper's output baseline, followed by the identical measurement.

The supplied model's `NOIBE` spelling is retained for the reproduction only.
Matched comparisons use the paper's corrected model (`NOIB`, closed model
parentheses). The saved Julia input waveform is reused with hashes, so both
circuits see identical uniform-arrival current and finite-path charge. The
native case retains the supplied capacitance and temperature; the matched
cases use the paper's capacitance and temperature. Both noise and signal are
processed using each circuit's calibrated tail and the same averaging and
frame-difference filter.

The inherited very tight voltage tolerance stalled the native v2 transient
in LTspice 17.2.4; paired looser/tighter settings are recorded explicitly in
the successful decks. The analysis checks their agreement as well as time
step, frequency density, quadrature, charge conservation and device-noise
variance closure. No claims of physical oscillation are made from the failed
solver attempts.

The generator checks input hashes and refuses to overwrite differing decks.
Physical calculations are in `src/ion_signal/rogelio_v2.jl`; the gate search
is controlled by Julia in `scripts/retune_rogelio_v2.jl`. Python acquires
LTspice data and presents saved results.
