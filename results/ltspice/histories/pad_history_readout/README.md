# Complete stored pad histories through the optimized CSA

This study reads all uncut pad histories in the 2 mm / 20 ms and 4 mm / 40 ms
samples. Julia expands missing sparse cells to zero counts, propagates the
small-signal response through every frame, and selects demanding histories.
Twelve selected histories are also run directly in nonlinear LTspice, with
matched zero-input references, two resolutions and a 20 s recovery tail.

The source files contain 8000 double-beta events and 8000 electron events per
configuration. The two configurations reuse the same diffused events.
`summary.h5` holds the ensemble counts, selected-case results and ADC assessment.
`spice_analysis.h5` holds the reduced waveforms and frame estimates.

## Physical and readout assumptions

- The circuit is `length_400um_gate0_rounded` from the feedback optimization.
  Mf W=0.36 um, L=400 um, gate ground; M4 gate=0.6164 V; external input
  capacitance 1.25 pF; Cf=3 fF. Other devices retain the frozen circuit.
- Each stored bucket supplies N uniformly distributed collection times. The
  current is the uniform-arrival mean of the finite-gap single-ion current.
- The current benchmark uses a 0.5 mm square pad in a grounded plane, a 5 mm
  final gap and constant ion speed 3 m/s. The unit-pad/zero-plane weighting
  potential determines the induced charge e[1-phi_w(g)]. No renormalization
  or second collection impulse is applied. This remains a planar benchmark.
- Every event begins at the quiescent baseline. All internal empty buckets and
  all high-charge deposits remain. Event-to-event pileup is not modelled.
- Boundary readings are 5 ms centred averages. The corrected frame signal is
  D_n = -[Delta Vbar_n - exp(-frame/tau) Delta Vbar_(n-1)]. The fixed gain and
  tau are calibrated at low charge. Nonlinear histories may violate that gain.
- A 0.2 V minimum-output guard selects demanding cases. It is a diagnostic
  choice, not a calibrated saturation threshold or an exact saturation rate.
- Nominal noise references are stationary spectra at the quiescent operating
  point. This study evaluates noiseless nonlinear transients; it does not
  calculate noise at the compressed operating points.
- ADC quantization assumes independent uniform errors on averaged samples.
  ADC input noise, buffer noise and oversampling credit are not included.

## Files and interfaces

`config/pad_history_readout.toml` supplies source paths, bucket times in seconds,
pitches in metres, averaging duration, diagnostic guard and solver resolutions.
The frozen source HDF5 files remain in their event-analysis directories.

`src/ion_signal/pad_history_readout.jl`:

- `kernel(root, config, frame_s)` reads the calibrated 100-ion SPICE response
  and returns per-ion output and averaged-response coefficients at the first
  three boundaries, the discharge factor, baseline and frame gain.
- `response(counts, kernel)` propagates a complete vector of bucket counts.
  It returns raw boundary voltages, centred averages and corrected frame
  signals, including two post-input frames. An exponential tail carries the
  accumulated circuit state after the explicit early response.
- `screen(root)` groups sparse cells by event and pad, fills all missing time
  bins with zero, checks the complete collected-ion ledger and writes one row
  per pad. It also saves complete histories for the lowest linear output,
  largest bucket and closest approach to the diagnostic guard.
- `spice_inputs(root)` converts those selected histories into charge-preserving
  PWL currents at two resolutions, retaining original event times and count
  vectors. Leading zero time is shifted away for simulation and recorded.

`src/ion_signal/pad_history_spice.jl`:

- `analyze(root)` subtracts matched baselines, integrates centred averaging
  windows, calculates corrected frame signals and their linear predictions,
  checks current charge and time refinement, and measures output extrema,
  guard recovery and signed empty-frame residuals.
- `summarize(root)` saves ensemble fractions, cumulative minimum-output curves,
  selected SPICE statistics and ADC quantization/noise comparisons to HDF5.

Python acquires LTspice results and presents the saved Julia output.
`run_pad_history_spice.py` preserves the exact circuit/model/PWL/netlist hashes
and reuses matching acquired waveforms. `present_pad_histories.py` writes the
self-contained report source and plots. No LaTeX compilation is performed.

## Reproduction

From the repository root, with the Julia environment and Python scientific
requirements installed:

```sh
python scripts/run_pad_history_pipeline.py --simulator /path/to/LTspice
```

The pipeline checks stored source, generator and configuration hashes before
reusing results. Inputs changed since this study require a separate output
directory. The individual stages are:

```sh
julia --project=. scripts/screen_pad_histories.jl
julia --project=. scripts/build_pad_history_spice_inputs.jl
PYTHONPATH=src python scripts/run_pad_history_spice.py --simulator /path/to/LTspice
julia --project=. scripts/analyze_pad_history_spice.jl
python scripts/present_pad_histories.py
```

`spice_inputs.h5` is preserved on reruns: use the pipeline to reuse it. Numerical
outputs, source hashes and validation results are saved in HDF5. The report is
`latex_generated/pad_history_readout/report.tex`. The source MC files and native
LTspice binary scratch are outside the tracked study outputs.

## Findings and limits

Only a small fraction of pads crosses the linear screening guard; several
percent of 4 mm events contain at least one such pad. The twelve demanding
transistor simulations show real compression above the ground rail. A compressed
signal can therefore remain within the ADC input range. The maximum observed
positive post-collection residual is below one quiescent ENC; some histories
have a larger negative residual in the first evaluated empty bucket.

These checks do not establish full gain recovery for a new small packet during
a compressed recovery tail. A compression flag needs calibration against the
transistor response before it is used in track reconstruction. Flagged large
deposits can then be retained as high-charge locations, as requested by the
user. Reconstruction performance with those flags remains a separate study.
