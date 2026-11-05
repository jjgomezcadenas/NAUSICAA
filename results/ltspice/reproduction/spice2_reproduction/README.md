# Rogelio noise reproduction — 27 October 2026

Run the complete pipeline from the repository root:

```sh
make spice2-reproduction PYTHON=/path/to/python LTSPICE=/path/to/LTspice
```

`make spice2-plots` regenerates the presentation from saved HDF5 without
rerunning LTspice or the Julia analysis. Python acquires simulator output;
Julia validates and reduces it; Python renders the saved results.

Original inputs are in `LTSpice/Spice2`. The derived model corrects only
NMOS NOIBE to NOIB and adds the two missing closing parentheses. Original
inputs are preserved byte for byte. The circuit is unchanged: 250 fF external
capacitance, 3 fF feedback, original biases and transistor sizes, default 27 C.
The main sweep is 1 Hz–1 GHz, 1000 points per decade. A second sweep starts at
0.1 Hz to reproduce the screenshot. The previous library is run as a control.

`raw.h5` contains the simulator spectra, operating points, logs and input hashes.
`analysis.h5` contains selected frequencies, total and flicker contributions,
variance closure and consistency checks. The derived decks and text logs are
retained. Native binary output is redundant with the HDF5 acquisition.

The six-panel figure reproduces the visible scales and shapes in Rogelio's
screenshot. No numerical export was supplied for an exact cross-simulator
comparison. Both sweeps agree; the revised library leaves the operating point
and gain unchanged. The total and component noise variances close to floating
point precision, and all seven transistor flicker terms are nonzero.
All 52 Python tests and the physics checks pass.

The self-contained report source and figures are in
`latex_generated/spice2_reproduction/`. No LaTeX compilation was run, following
the user's standing instruction. No NAUSICA ENC result is updated by this
original-configuration reproduction.
