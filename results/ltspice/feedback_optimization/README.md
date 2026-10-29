# Feedback and averaging optimization

The study scans the feedback NMOS gate voltage, width and length while retaining
the paper's remaining device dimensions, 1.25 pF external capacitance, 3 fF
feedback capacitor, 1.8 V supply and 20 C temperature. The physical input is
copied from the preserved 20 and 40 ms finite-gap studies without renormalization.
Mean populations of 120 and 400 ions are user design inputs.

The first scan holds amplifier biases fixed. The retuned scan adjusts M4's gate
to restore the reference DC output. The combined scan tests longer channels with
the feedback gate grounded. The final candidate is rerun with a four-decimal
M4 bias. Each full circuit is tested with two transient resolutions, two spectral
densities and multiple packet populations. Averaging windows are 0.3, 0.5, 1, 2
and 5 ms. Gain and noise are recalculated for every averaging window.

Reproduce from the repository root:

```sh
PYTHONPATH=src python scripts/run_feedback_pipeline.py --julia /path/to/julia --simulator /path/to/LTspice
```

The pipeline resumes matching acquisitions. A fixed-bias point whose transient
did not converge is recorded and excluded; its retuned counterpart is evaluated.
Unexpected failures without recorded acquisition metadata stop the pipeline.

- `config/feedback_optimization.toml`: primary scan and validity criteria.
- `config/feedback_combined.toml`: combined length/bias settings.
- `inputs.h5`: preserved input currents and their source hashes.
- `<variant>.h5`: complete compressed numerical SPICE output, source hashes and logs.
- `retuned_points*.h5`: DC operating-point searches (Julia determines bias settings).
- `analysis.h5`: Julia signal, noise, retention, linearity and validation results.
- `decks/inputs/`: shared current waveforms, time in seconds and current in amperes.
- `decks/<variant>/`: exact netlists, local model card and readable run logs.
- `decks/tuning/` and `tuning/`: saved operating-point trials.

Julia reduction and checks: `src/ion_signal/feedback_optimization.jl`.
Report source and scientific figure: `latex_generated/feedback_optimization/`.
Presentation alone: `python scripts/present_feedback_optimization.py`.
The report is generated directly from the saved numerical analysis. The main
paper and Vicente's manuscript package remain their previous circuit snapshot.

The selected candidate is `length_400um_gate0_rounded`: Mf W=0.36 um,
L=400 um, gate=0 V; M4 gate=0.6164 V. Its netlists reference shared PWL files
through `../inputs/`. The other physical and circuit assumptions, the scan
criteria and the longer-retention tradeoff are stated in the generated report.
