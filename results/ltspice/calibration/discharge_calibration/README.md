# CSA discharge calibrated from LTspice

The selected circuit, 20 ms bucket and 300 µs centred averages are unchanged.
The correction time is obtained by fitting the current-free output tail of a
100-ion packet between 21 and 120 ms. The 10-ion result, separate fit windows,
a late tail outside the fit window, and halved numerical steps validate it.

- Effective circuit decay: 2.454158 s; alpha = 0.9918836822 per 20 ms.
- Recalibrated gain: 27.1192 µV/ion; filtered noise: 3.05636 mV RMS.
- ENC: 112.701 e RMS; finite-charge threshold: 564 ions for S/N > 5.
- Previous-packet residual: 0.0820% of its own-bucket signal (previously 5.747%).
- For the tested eight-bucket sequence, the largest deterministic error falls
  from 303.92 to 12.23 ions. The same circuit's stationary noise is evaluated
  separately; no random noise is added to this sequence.

The previous 0.3 s results in `results/ltspice/bucket_scan/10ms/` and `20ms/`
remain unchanged. The new input, acquisition and reduction are in `20ms/`
under this study. Each run checks source hashes before resuming. A changed
physical configuration requires a new output path; input generation refuses
to overwrite an existing file.

The report and proposed manuscript macros are generated in
`latex_generated/discharge_calibration/`. They are not loaded into the paper.
The baseline is the matched zero-input circuit simulation. The input current
uses the existing central-square-pad finite-gap model and preserves its charge.
The fitted time is the decay of the full circuit response, including finite
gain and the feedback transistor.

## Reproduce

From the repository root:

```sh
export PYTHONPATH="$PWD/src"
export JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia"
export MPLCONFIGDIR="$PWD/notes/physics/work/matplotlib-cache"
python scripts/run_discharge_pipeline.py --julia /path/to/julia --simulator /path/to/LTspice
julia --project=. tests/test_discharge_calibration.jl
make test physics
```

Configuration: `config/discharge_calibration.toml`.

Individual stages:

```sh
# Input generation is for a fresh output directory.
julia --project=. scripts/build_discharge_inputs.jl
python scripts/run_discharge_calibration.py --simulator /path/to/LTspice
julia --project=. scripts/analyze_discharge_calibration.jl
python scripts/run_discharge_calibration.py --simulator /path/to/LTspice --counts 560 561 562 563 564 565
julia --project=. scripts/analyze_discharge_calibration.jl
python scripts/present_discharge_calibration.py
```

Julia performs input-current generation, the decay fit, noise integration and
charge reconstruction, saving HDF5. Python handles LTspice acquisition and
presentation. Source hashes record circuit, model, waveforms, configuration,
simulator, noise source and code. Native simulator scratch is excluded from
git; compressed HDF5, PWL sources, decks and readable logs preserve the results.
The paper and report have not been compiled.

## Manuscript assets and macro audit

The active paper now uses the calibrated values. The final-circuit averaging
scan at 180, 200, 250, 300 and 350 µs is reduced into `20ms/paper_readout.h5`.
It also contains the current and baseline-subtracted voltage, averaged samples,
threshold bracket and filtered noise-source fractions for presentation.

```sh
julia --project=. scripts/build_averaging_scan.jl
make spice-paper PYTHON=python
python scripts/audit_paper_macros.py
```

The audit checks all 41 loaded numerical macros against saved calculations,
configuration and netlist values; `20ms/macro_audit.json` records each source.
The report for the updated scan is
`latex_generated/calibrated_spice_averaging/report.tex`. The previous
`results/ltspice/averaging_scan.h5` and its report remain unchanged.
