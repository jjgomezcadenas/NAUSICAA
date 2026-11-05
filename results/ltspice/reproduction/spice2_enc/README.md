# NAUSICA ENC with the revised Spice2 model

Reproduce from the repository root:

```sh
make spice2-enc PYTHON=/path/to/python LTSPICE=/path/to/LTspice
```

Python runs LTspice and saves its output to raw.h5. Julia computes baseline-
subtracted voltage averages, filtered noise, convergence checks and the next
integer counts to simulate. The acquisition loop continues until adjacent
integer counts bracket S/N=5. `make spice2-enc-plots` rebuilds the figures and
standalone LaTeX report from analysis.h5 alone.

Conditions: external pad 1.25 pF plus intrinsic MOS capacitances; Cf=3 fF;
20 C; uniform mean arrivals over 10 ms; centered 300 us averages; correction
tau=0.3 s; stationary one-sided output noise integrated over 0.1 Hz–1 MHz.
The saved PWL is the finite 5 mm approach to a 0.5 mm square pad in a grounded
plane, with no charge renormalization. Only the model noise parameters differ
from the previous 1.25 pF study. The original received files are preserved.
The corrected copy fixes NOIBE to NOIB and closes the two model parameter lists.

The gain remains 15.4694 uV/ion. The filtered RMS noise is 4.57672 mV;
ENC=295.856 electrons RMS, of which flicker contributes 289.944 and thermal
sources 58.8433 in quadrature. The verified threshold is 1485 ions:
1484 gives S/N=4.999850 and 1485 gives 5.003208. Both timestep resolutions
place the adjacent counts on the same respective sides of five. The numerical
refinement errors are at least ten times smaller than the distance to the
threshold. These rounded summaries are generated numerically in analysis.h5;
report tables and figures are produced from that file.

Charge conservation, transient resolution, spectral density, quadrature,
noise-component closure and the old/new gain agreement are checked in the
analysis. All 52 Python tests and the physics checks pass. The noise remains
an operating-point calculation using the coefficients supplied by Rogelio;
no ADC noise or independent process calibration is included.

Report source and figure: latex_generated/spice2_enc/. No LaTeX compilation
or manuscript update was performed. Previous results are retained separately.
