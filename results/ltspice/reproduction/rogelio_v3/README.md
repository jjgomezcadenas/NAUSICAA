# Rogelio v3: input charge and output amplitude

The supplied schematic and card remain unchanged in `LTSpice/Spice2/v3/`.
This study reproduces the original 40-second transient, repeats it with
waveform compression disabled, then resolves the exponential input pulse
with paired maximum timesteps and tight tolerances. The same v3 circuit is
also driven by the paper's saved 1000-ion uniform-arrival waveform.

The generated report is `latex_generated/rogelio_v3/report.tex`.

```sh
export PYTHONPATH=src
export JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia"
python scripts/run_rogelio_v3.py \
  --simulator notes/physics/work/ltspice-review/legacy/LTspice.pkg/Payload/Applications/LTspice.app/Contents/MacOS/LTspice
julia --project=. scripts/analyze_rogelio_v3.jl
julia --project=. tests/test_rogelio_v3.jl
python scripts/present_rogelio_v3.py
```

Python acquires LTspice data into `raw.h5`. Julia integrates the actual saved
input current, compares it with the behavioral `idt` accumulator, fits the
source decay, measures voltage drops and records convergence in `analysis.h5`.
Python presents the saved results. Native simulator raw files are excluded
from git; numerical content is retained in HDF5.

The original compressed current trace gives unreliable numerical quadrature
and `.meas` charge results for the tiny pulse. Disabling compression alone
restores agreement with the behavioral accumulator; the refined runs establish
the final charge and amplitude. In `analysis.h5`, `resolved_charge_integration`
identifies the refined cases. The original compressed quadrature is retained
as diagnostic evidence, not as the physical input-charge result.

The schematic's resistance label uses `V(out)/Id(M4)`. The study records
that DC ratio separately from `1/gds`, read from the transistor operating-point
log. The report is generated without compilation or manuscript edits.
