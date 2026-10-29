# Scripts inventory and organization

The directory migration is implemented on branch `reorg`. The commands below
link to their current locations. The complete old-to-new map and original source
hashes are in [script_paths.json](script_paths.json); see the
[workflow guide](../scripts/README.md) for current entry points.


The current `scripts/` tree contains **138 Python and Julia commands**, grouped
by scientific task. The original inventory listed 139 top-level scripts. The
migration also includes the three direct-collection commands and two previously
grouped input-comparison commands; four content modules and two included helpers
now live alongside their documentation or reusable source code.

The active ion/PE workflow uses `collection/sample_collect_stream.jl` and
`collection/sample_collect_parallel.jl`. Ion-only pixel studies use
`collection/run_pixel_charge_direct.jl`. Their saved direct histories are the
inputs to the current electronics studies.

## Directory organization

Group by scientific task and keep each study's input generation, acquisition,
Julia analysis and Python presentation together. Use one documented complete
workflow command per study, with separate commands to regenerate presentation
from saved results.

```text
scripts/
  events/                         # readers and event reconstruction
  transport/                      # carrier diffusion and stream benchmarks
  collection/                     # clouds to pads or SiPMs and time buckets
  electronics/
    reproduction/
      initial/
      spice2/
      rogelio_v2/
      rogelio_v3/
    devices/                      # M3 and M4 geometry
    bias/                         # bias precision and sensitivity
    readout/                      # frames, discharge, averaging and noise
    feedback/                     # feedback-device optimization
    histories/                    # complete histories and compression
    design/                       # ADC and area budgets
  benchmarks/                     # analytical signal and circuit references
  viewers/                        # application launchers
  paper/                          # manuscript assembly and handoff
  validation/                     # Petit references and checkpoint reports
  documentation/                  # review-guide builder
```

Only launch/orchestration scripts should remain in this tree. Reusable numerical
code belongs in `src/event_analysis` or `src/ion_signal`; application code belongs
in its existing `src/*_viewer` package. Move imported documentation content next
to the review-guide sources rather than presenting it as executable commands.
The readout HTML template belongs with the benchmark's presentation assets.

## Which commands to start with

| Task | Present starting command | What follows |
|---|---|---|
| Import raw NEXUS data | `extract_hits.jl` | Reduced hit HDF5 files. |
| MC selection | `analyze_mc_selection.jl` | `present_mc_selection.py`. |
| Direct ion/PE readout | `collection/sample_collect_parallel.jl` | `collection/summarize_grid_readout.jl`, then `collection/present_grid_readout.py`. |
| Direct ion-only histories | `collection/run_pixel_charge_direct.jl` | `collection/present_pixel_charge.py`. |
| Six carrier clouds | `propagate_clouds_parallel.jl` | `present_cloud_propagation.py`. |
| Cloud collection | `propagate_grids_parallel.jl` | `summarize_grid_readout.jl`, then `present_grid_readout.py`. |
| Earlier pixel-charge pipeline | `analyze_mc_selection.jl --pixel-charge` | Includes `run_pixel_charge.jl` and `reduce_pixel_charge.jl`; then `present_pixel_charge.py`. |
| Original event viewer | `view_events.py` | `src/event_viewer`. |
| Paired cloud viewer | `serve_cloud_viewer.py` | `src/cloud_viewer`. |
| CSA history viewer | `view_pad_histories.py` | `src/pad_history_viewer`; cached/on-demand LTspice. |
| Feedback optimization | `run_feedback_pipeline.py` | Acquisition, tuning, Julia reduction and report. |
| Complete pad histories | `run_pad_history_pipeline.py` | Screening, selected SPICE cases, reduction and report. |
| Event compression study | `run_pad_compression.py` | Candidate/control SPICE runs and event statistics. |
| Current paper transistor results | `build_averaging_scan.jl` | `plot_averaging_scan.py`, then `build_spice_paper.py`. |

The current collection commands use direct carrier sampling. Earlier cloud
and two-stage commands reproduce their historical studies. The directory
migration preserves those calculation methods and saved outputs.

## Naming and shared code

Use consistent stage names for future commands: `run` for a complete workflow,
`inputs` for current/input generation, `acquire` for LTspice, `analyze` for Julia
reduction, and `present` for figures/tables/report source. Existing filenames can
be retained during the first move. In particular, `report_pad_compression.jl`
performs numerical analysis; `build_averaging_scan.jl` also performs analysis;
`build_spice_paper.py` only presents saved results.

`run_feedback_retuned.py` and `run_feedback_combined.py` differ in their selected
configuration. `run_feedback_candidate.py` additionally uses the rounded-case
output suffix. Replace that duplication with one acquisition command accepting
the configuration and case label.

`build_ltspice_study.py` supplies a runner imported by many later studies. Extract
that runner into a shared module before moving the scripts. Likewise, the
`run_pixel_charge.jl` and `pixel_diffusion_cache.jl` helper implementations belong
with event-analysis code, leaving thin command entry points. The older
`analyze_ltspice_study.py` and some early `build_*` scripts retain numerical Python
calculations; moving them does not by itself convert them to Julia.
`present_grid_readout.py` also calculates several means from saved event arrays;
those reductions should move to its Julia summary stage.

Presentation commands write figures and LaTeX source. Explicit `make pdf-*`
targets compile the electron-diffusion, grid-readout and MC-selection reports.

## Safe migration sequence

1. Keep this inventory as the mapping from present files to destinations.
2. Extract shared acquisition and cache helpers without changing their behaviour.
3. Move one workflow at a time, preserving filenames and user-facing Make targets.
4. Update repository-root discovery, sibling imports, Julia includes, subprocess
   calls, templates, tests, documentation and packaging paths together.
5. Record the path mapping alongside provenance; preserve historical source
   hashes and frozen collaborator packages. Changed scripts get new provenance
   rather than rewritten historical records.
6. Check each workflow's lightweight entry points and relevant integration tests.
   Reuse saved results for presentation checks; perform simulations when a moved
   acquisition path needs validation.

Archive only code whose consumers and reproduction role have been checked.
Older analytical generators are still used by build targets, and old circuit
scripts sometimes supply helpers to the current pipeline. Keep their historical
calculations reproducible.

## Complete file inventory

Roles describe the current implementation, including mixed-purpose scripts.
Each group gives its proposed destination.

### Event input and reconstruction

Directory: `scripts/events/`.

Event-file ingestion, the earlier ion-track analysis, and the current MC selection workflow. The older analysis is retained for its recorded results.

| Current file | Role | Purpose |
|---|---|---|
| [extract_hits.jl](../scripts/events/extract_hits.jl) | Input | Convert NEXUS hits to reduced HDF5 while retaining event and hit identity. |
| [analyze_hits.jl](../scripts/events/analyze_hits.jl) | Analysis | Run the earlier ion-track analysis through src/ion_signal/event_analysis.jl. |
| [plot_hit_analysis.py](../scripts/events/plot_hit_analysis.py) | Presentation | Render that earlier ion-track analysis and its LaTeX report. |
| [analyze_mc_selection.jl](../scripts/events/analyze_mc_selection.jl) | Analysis and dispatcher | Run MC selection; dispatch --pixel-charge to collection and reduction helpers. |
| [present_mc_selection.py](../scripts/events/present_mc_selection.py) | Presentation | Render saved MC selection results as figures, HTML and LaTeX. |

### Carrier transport

Directory: `scripts/transport/`.

Ion and electron diffusion, six carrier/sample streams, and their process benchmarks. The existing coarse diffusion representation needs the resolution fix identified in the time-axis audit.

| Current file | Role | Purpose |
|---|---|---|
| [calculate_electron_diffusion.jl](../scripts/transport/calculate_electron_diffusion.jl) | Analysis | Evaluate electron diffusion coefficients and widths for pure Xe and Xe–He. |
| [present_electron_diffusion.py](../scripts/transport/present_electron_diffusion.py) | Presentation | Render the saved electron diffusion calculation and compile its report. |
| [propagate_cloud_stream.jl](../scripts/transport/propagate_cloud_stream.jl) | Worker CLI | Propagate one sample and carrier stream into its own HDF5 file. |
| [propagate_clouds_parallel.jl](../scripts/transport/propagate_clouds_parallel.jl) | Pipeline | Launch six independent diffusion streams with bounded concurrency. |
| [present_cloud_propagation.py](../scripts/transport/present_cloud_propagation.py) | Presentation | Render the saved six-stream propagation summary. |
| [benchmark_three_clouds.py](../scripts/transport/benchmark_three_clouds.py) | Benchmark | Measure six-stream runtime, storage and memory for different process counts. |
| [present_three_cloud_benchmark.py](../scripts/transport/present_three_cloud_benchmark.py) | Presentation | Render the saved concurrency benchmark as a LaTeX report. |

### Carrier collection

Directory: `scripts/collection/`.

Conversion of saved clouds to detector grids, plus the earlier pixel-charge workflow. Keep the two HDF5 schemas distinct in commands and documentation.

| Current file | Role | Purpose |
|---|---|---|
| [collect_cloud_stream.jl](../scripts/collection/collect_cloud_stream.jl) | Worker CLI | Map one saved cloud stream to ion-pad or electron-SiPM time histories. |
| [propagate_grids_parallel.jl](../scripts/collection/propagate_grids_parallel.jl) | Pipeline | Launch six independent cloud-to-grid collection streams. |
| [summarize_grid_readout.jl](../scripts/collection/summarize_grid_readout.jl) | Analysis | Reduce completed grid histories to event and ensemble summaries. |
| [present_grid_readout.py](../scripts/collection/present_grid_readout.py) | Presentation and summary statistics | Render ion/PE readout summaries; currently also calculates several means and compiles the report. |
| [run_pixel_charge.jl](../src/event_analysis/batch/run_pixel_charge.jl) | Included helper | Define the earlier full-cloud transport and pixel/frame batch called by analyze_mc_selection.jl. |
| [pixel_diffusion_cache.jl](../src/event_analysis/batch/pixel_diffusion_cache.jl) | Included helper | Read and write the diffusion cache shared by pixel-grid configurations. |
| [reduce_pixel_charge.jl](../scripts/collection/reduce_pixel_charge.jl) | Analysis and included helper | Rebuild pulse distributions from saved pixel counts without repeating transport. |
| [present_pixel_charge.py](../scripts/collection/present_pixel_charge.py) | Presentation | Render saved inclusive, mean, maximum and minimum pixel-pulse distributions. |

### Interactive viewers

Directory: `scripts/viewers/`.

User-facing launch commands. The application implementations already belong in src/event_viewer, src/cloud_viewer and src/pad_history_viewer.

| Current file | Role | Purpose |
|---|---|---|
| [view_events.py](../scripts/viewers/view_events.py) | Launcher | Open the original MC event, reconstruction, transport and pad viewer. |
| [serve_cloud_viewer.py](../scripts/viewers/serve_cloud_viewer.py) | Launcher | Open the paired six-cloud browser viewer. |
| [view_pad_histories.py](../scripts/viewers/view_pad_histories.py) | Launcher | Open the linked counts, CSA-voltage and recovered-charge viewer. |

### Analytical benchmarks

Directory: `scripts/benchmarks/`.

Analytical models and reference figures. Several remain in the Makefile and manuscript build; age alone is not a reason to archive them.

| Current file | Role | Purpose |
|---|---|---|
| [build_physics.py](../scripts/benchmarks/build_physics.py) | Calculation and presentation | Generate detector transport and pad-weighting benchmark numbers and figures. |
| [build_readout.py](../scripts/benchmarks/build_readout.py) | Calculation and presentation | Generate ideal CSA/current benchmarks, an HTML demonstration and SPICE currents. |
| [build_realistic_csa.py](../scripts/benchmarks/build_realistic_csa.py) | Calculation and presentation | Generate the finite-gain linear CSA, noise and SPICE comparison benchmarks. |
| [build_gain_figure.py](../scripts/benchmarks/build_gain_figure.py) | Presentation | Draw the normalized open-loop amplifier-gain illustration. |
| [build_arrival_cases.py](../scripts/benchmarks/build_arrival_cases.py) | Calculation and presentation | Generate the uniform-arrival mean response and simultaneous-arrival cross-check. |
| [build_noise_study.py](../scripts/benchmarks/build_noise_study.py) | Calculation and presentation | Generate the illustrative linear-model frame-noise estimates. |
| [build_capacitance_check.py](../scripts/benchmarks/build_capacitance_check.py) | Historical diagnostic | Compare the former illustrative capacitance with the Topmetal-S reference load. |

### Initial transistor studies

Directory: `scripts/electronics/reproduction/initial/`.

First LTspice reproduction, shaping and capacitance studies. Some functions are imported by later scripts, so shared helpers must be extracted before any archival move.

| Current file | Role | Purpose |
|---|---|---|
| [build_ltspice_study.py](../scripts/electronics/reproduction/initial/build_ltspice_study.py) | Acquisition and shared helper | Acquire the initial circuit transients and noise; also supplies the reused LTspice runner. |
| [analyze_ltspice_study.py](../scripts/electronics/reproduction/initial/analyze_ltspice_study.py) | Analysis and presentation | Reduce the initial spectra and frame-estimator results. |
| [audit_ltspice_noise.py](../scripts/electronics/reproduction/initial/audit_ltspice_noise.py) | Diagnostic | Audit the original noise setup and a synthetic flicker-noise control. |
| [build_shaping_study.py](../scripts/electronics/reproduction/initial/build_shaping_study.py) | Acquisition and analysis | Run short-injection shaping and early ion-frame averaging studies. |
| [build_capacitance_study.py](../scripts/electronics/reproduction/initial/build_capacitance_study.py) | Acquisition and analysis | Compare external capacitance assumptions using the transistor circuit. |

### Updated Spice2 model

Directory: `scripts/electronics/reproduction/spice2/`.

Reproduction of the revised noise card, calibrated ENC, and capacitance/injection diagnostics.

| Current file | Role | Purpose |
|---|---|---|
| [reproduce_spice2.py](../scripts/electronics/reproduction/spice2/reproduce_spice2.py) | Acquisition | Acquire the updated circuit noise spectra with documented syntax corrections. |
| [analyze_spice2.jl](../scripts/electronics/reproduction/spice2/analyze_spice2.jl) | Analysis | Validate component-noise closure and reduce the updated spectra. |
| [plot_spice2.py](../scripts/electronics/reproduction/spice2/plot_spice2.py) | Presentation | Render the updated spectra and saved Julia reductions. |
| [run_spice2_enc.py](../scripts/electronics/reproduction/spice2/run_spice2_enc.py) | Pipeline and acquisition | Run the 1.25 pF ion-frame calibration and ENC study. |
| [analyze_spice2_enc.jl](../scripts/electronics/reproduction/spice2/analyze_spice2_enc.jl) | Analysis | Calculate calibrated frame gain, noise and the ion threshold. |
| [plot_spice2_enc.py](../scripts/electronics/reproduction/spice2/plot_spice2_enc.py) | Presentation | Render the calibrated revised-model ion-frame results. |
| [run_spice2_diagnostics.py](../scripts/electronics/reproduction/spice2/run_spice2_diagnostics.py) | Acquisition | Acquire capacitance and short-injection comparison cases. |
| [analyze_spice2_diagnostics.jl](../scripts/electronics/reproduction/spice2/analyze_spice2_diagnostics.jl) | Analysis | Reduce the revised-model capacitance and injection diagnostics. |
| [present_spice2_diagnostics.py](../scripts/electronics/reproduction/spice2/present_spice2_diagnostics.py) | Presentation | Render the saved diagnostic comparisons. |

### Rogelio circuit version two

Directory: `scripts/electronics/reproduction/rogelio_v2/`.

A preserved circuit-version assessment, with acquisition, bias retuning and numerical diagnostics.

| Current file | Role | Purpose |
|---|---|---|
| [run_rogelio_v2_pipeline.py](../scripts/electronics/reproduction/rogelio_v2/run_rogelio_v2_pipeline.py) | Pipeline | Reproduce v2, validate the threshold, retune and generate its report. |
| [run_rogelio_v2.py](../scripts/electronics/reproduction/rogelio_v2/run_rogelio_v2.py) | Acquisition | Acquire supplied-v2 and matched ion/noise calculations. |
| [analyze_rogelio_v2.jl](../scripts/electronics/reproduction/rogelio_v2/analyze_rogelio_v2.jl) | Analysis | Reduce the v2 or retuned-v2 results and optional diagnostics. |
| [retune_rogelio_v2.jl](../scripts/electronics/reproduction/rogelio_v2/retune_rogelio_v2.jl) | Tuning | Choose a v2 bias using saved operating-point calculations. |
| [run_v2_gate_point.py](../scripts/electronics/reproduction/rogelio_v2/run_v2_gate_point.py) | Acquisition helper | Acquire a single unchanged-geometry v2 gate-bias point. |
| [diagnose_rogelio_v2.py](../scripts/electronics/reproduction/rogelio_v2/diagnose_rogelio_v2.py) | Diagnostic acquisition | Record M4 operating parameters and voltage gain at fixed bias. |
| [check_v2_solver.py](../scripts/electronics/reproduction/rogelio_v2/check_v2_solver.py) | Numerical diagnostic | Check zero-input transient behaviour of the high-output-bias v2 circuit. |
| [compact_v2_raw.py](../scripts/electronics/reproduction/rogelio_v2/compact_v2_raw.py) | Storage utility | Losslessly compress v2 HDF5 acquisitions. |
| [present_rogelio_v2.py](../scripts/electronics/reproduction/rogelio_v2/present_rogelio_v2.py) | Presentation | Render the supplied and retuned v2 assessment. |

### Rogelio circuit version three

Directory: `scripts/electronics/reproduction/rogelio_v3/`.

Comparison of the supplied v3 input pulse and measured noise convention.

| Current file | Role | Purpose |
|---|---|---|
| [run_rogelio_v3.py](../scripts/electronics/reproduction/rogelio_v3/run_rogelio_v3.py) | Acquisition | Reproduce the supplied EXP pulse and recorded ion input. |
| [analyze_rogelio_v3.jl](../scripts/electronics/reproduction/rogelio_v3/analyze_rogelio_v3.jl) | Analysis | Reduce v3 input-charge and voltage-amplitude comparisons. |
| [present_rogelio_v3.py](../scripts/electronics/reproduction/rogelio_v3/present_rogelio_v3.py) | Presentation | Render the v3 comparison report. |
| [run_noise_measurement.py](../scripts/electronics/reproduction/rogelio_v3/run_noise_measurement.py) | Acquisition | Snapshot the updated v3 circuit and reproduce its noise measurement. |
| [analyze_noise_measurement.jl](../scripts/electronics/reproduction/rogelio_v3/analyze_noise_measurement.jl) | Analysis | Calculate the noise integrals and measurement comparisons. |
| [present_noise_measurement.py](../scripts/electronics/reproduction/rogelio_v3/present_noise_measurement.py) | Presentation | Render the saved noise-measurement analysis. |

### Transistor geometry

Directory: `scripts/electronics/devices/`.

M3 and M4 geometry sweeps and their bias tuning. Each study should expose one complete run command.

| Current file | Role | Purpose |
|---|---|---|
| [run_m3_scan.py](../scripts/electronics/devices/run_m3_scan.py) | Acquisition | Acquire M3 geometry variants at fixed W/L and unchanged bias. |
| [analyze_m3_scan.jl](../scripts/electronics/devices/analyze_m3_scan.jl) | Analysis | Reduce the M3 sweep using the shared geometry-study implementation. |
| [tune_m3_bias.jl](../scripts/electronics/devices/tune_m3_bias.jl) | Tuning | Retune the M3-study bias against its reference operating point. |
| [present_m3_scan.py](../scripts/electronics/devices/present_m3_scan.py) | Presentation | Render M3 results and export portable LTspice variants. |
| [diagnose_m3_injection.jl](../scripts/electronics/devices/diagnose_m3_injection.jl) | Numerical diagnostic | Check small-injection calibration and solver-tolerance sensitivity. |
| [run_m4_scan.py](../scripts/electronics/devices/run_m4_scan.py) | Acquisition | Acquire M4 geometry variants at fixed W/L and unchanged bias. |
| [analyze_m4_scan.jl](../scripts/electronics/devices/analyze_m4_scan.jl) | Analysis | Reduce the M4 geometry and operating-point sweep. |
| [tune_m4_bias.jl](../scripts/electronics/devices/tune_m4_bias.jl) | Tuning | Retune M4 bias for the selected geometry variants. |
| [present_m4_scan.py](../scripts/electronics/devices/present_m4_scan.py) | Presentation | Render the M4 geometry/bias comparison. |

### Bias precision

Directory: `scripts/electronics/bias/`.

Bias rounding and sensitivity studies for the paper and v2 circuits.

| Current file | Role | Purpose |
|---|---|---|
| [run_bias_sensitivity_pipeline.py](../scripts/electronics/bias/run_bias_sensitivity_pipeline.py) | Pipeline | Build bias points, acquire transients and noise, reduce and render. |
| [build_bias_sensitivity.jl](../scripts/electronics/bias/build_bias_sensitivity.jl) | Input generation | Write the configured bias offsets and decimal-rounding cases. |
| [run_bias_sensitivity.py](../scripts/electronics/bias/run_bias_sensitivity.py) | Acquisition | Acquire the configured gate-bias variations. |
| [analyze_bias_sensitivity.jl](../scripts/electronics/bias/analyze_bias_sensitivity.jl) | Analysis | Reduce bias-dependent gain, noise and operating-point results. |
| [present_bias_sensitivity.py](../scripts/electronics/bias/present_bias_sensitivity.py) | Presentation | Render the saved sensitivity results. |

### Frame duration and discharge

Directory: `scripts/electronics/readout/`.

Bucket-size comparisons, discharge calibration, averaging, neighbouring-frame noise and the pitch exercise.

| Current file | Role | Purpose |
|---|---|---|
| [run_bucket_pipeline.py](../scripts/electronics/readout/run_bucket_pipeline.py) | Pipeline | Generate inputs, acquire frame-size cases, validate and render results. |
| [build_bucket_inputs.jl](../scripts/electronics/readout/build_bucket_inputs.jl) | Input generation | Generate finite-gap inputs for the configured bucket durations. |
| [run_bucket_study.py](../scripts/electronics/readout/run_bucket_study.py) | Acquisition | Acquire bucket-duration transients and reuse matched noise spectra. |
| [analyze_bucket_study.jl](../scripts/electronics/readout/analyze_bucket_study.jl) | Analysis | Validate acquisition provenance and reduce the bucket-size study. |
| [present_bucket_study.py](../scripts/electronics/readout/present_bucket_study.py) | Presentation | Render the saved 10/20 ms comparison. |
| [run_discharge_pipeline.py](../scripts/electronics/readout/run_discharge_pipeline.py) | Pipeline | Generate decay inputs, acquire tails, calibrate and render results. |
| [build_discharge_inputs.jl](../scripts/electronics/readout/build_discharge_inputs.jl) | Input generation | Generate extended-tail and packet-sequence currents. |
| [run_discharge_calibration.py](../scripts/electronics/readout/run_discharge_calibration.py) | Acquisition | Acquire long discharge tails and ion-packet sequences. |
| [analyze_discharge_calibration.jl](../scripts/electronics/readout/analyze_discharge_calibration.jl) | Analysis | Fit discharge and calibrate the frame-difference estimator. |
| [present_discharge_calibration.py](../scripts/electronics/readout/present_discharge_calibration.py) | Presentation | Render the saved discharge fit and estimator validation. |
| [build_averaging_scan.jl](../scripts/electronics/readout/build_averaging_scan.jl) | Analysis | Reduce the selected paper-circuit averaging scan using PaperReadout. |
| [plot_averaging_scan.py](../scripts/electronics/readout/plot_averaging_scan.py) | Presentation | Render calibrated averaging-scan results from HDF5. |
| [analyze_bucket_noise_covariance.jl](../scripts/electronics/readout/analyze_bucket_noise_covariance.jl) | Analysis | Calculate stationary-noise covariance between frame estimates. |
| [present_bucket_noise_covariance.py](../scripts/electronics/readout/present_bucket_noise_covariance.py) | Presentation | Render the saved covariance and correlation results. |
| [run_pitch_exercise.py](../scripts/electronics/readout/run_pitch_exercise.py) | Pipeline | Run and report the earlier 4 mm pitch and 40 ms bucket exercise. |
| [pitch_exercise.jl](../scripts/electronics/readout/pitch_exercise.jl) | Input analysis and checks | Expose the pitch exercise's input, analysis and validation stages. |
| [present_pitch_exercise.py](../scripts/electronics/readout/present_pitch_exercise.py) | Presentation | Render the saved pitch/bucket exercise. |

### Feedback optimization

Directory: `scripts/electronics/feedback/`.

The feedback-device and averaging optimization. Three acquisition wrappers differ mainly in configuration or output suffix.

| Current file | Role | Purpose |
|---|---|---|
| [run_feedback_pipeline.py](../scripts/electronics/feedback/run_feedback_pipeline.py) | Pipeline | Run device scans, operating-point tuning and selected-candidate reporting. |
| [run_feedback_optimization.py](../scripts/electronics/feedback/run_feedback_optimization.py) | Acquisition | Acquire the configured feedback-device sweep. |
| [tune_feedback_optimization.jl](../scripts/electronics/feedback/tune_feedback_optimization.jl) | Tuning | Retune feedback-sweep operating points to the reference baseline. |
| [run_feedback_retuned.py](../scripts/electronics/feedback/run_feedback_retuned.py) | Acquisition helper | Acquire a retuned point using feedback_optimization.toml. |
| [tune_feedback_combined.jl](../scripts/electronics/feedback/tune_feedback_combined.jl) | Tuning | Retune the combined feedback variants. |
| [run_feedback_combined.py](../scripts/electronics/feedback/run_feedback_combined.py) | Acquisition helper | Acquire a retuned point using feedback_combined.toml. |
| [run_feedback_candidate.py](../scripts/electronics/feedback/run_feedback_candidate.py) | Acquisition helper | Acquire the selected rounded-bias candidate with a distinct output suffix. |
| [analyze_feedback_optimization.jl](../scripts/electronics/feedback/analyze_feedback_optimization.jl) | Analysis and checks | Reduce and validate selected feedback-optimization cases. |
| [present_feedback_optimization.py](../scripts/electronics/feedback/present_feedback_optimization.py) | Presentation | Render the feedback and averaging comparison. |

### Complete pad histories

Directory: `scripts/electronics/histories/`.

Stored-event propagation, nonlinear compression, standard-case caches, the illustrative sequence and the timing audit. Event-level conclusions from the coarse input files require the diffusion-resolution fix.

| Current file | Role | Purpose |
|---|---|---|
| [run_pad_history_pipeline.py](../scripts/electronics/histories/run_pad_history_pipeline.py) | Pipeline | Screen all stored pad histories, simulate selected cases, reduce and render. |
| [screen_pad_histories.jl](../scripts/electronics/histories/screen_pad_histories.jl) | Analysis | Propagate every pad with the calibrated linear kernel and select demanding cases. |
| [build_pad_history_spice_inputs.jl](../scripts/electronics/histories/build_pad_history_spice_inputs.jl) | Input generation | Generate charge-preserving currents for the selected complete histories. |
| [run_pad_history_spice.py](../scripts/electronics/histories/run_pad_history_spice.py) | Acquisition | Acquire matched-baseline nonlinear transients at two resolutions. |
| [analyze_pad_history_spice.jl](../scripts/electronics/histories/analyze_pad_history_spice.jl) | Analysis | Reduce selected transients, empty-frame residuals and ADC quantization estimates. |
| [present_pad_histories.py](../scripts/electronics/histories/present_pad_histories.py) | Presentation | Render the complete-history screening and nonlinear validation report. |
| [run_pad_compression.py](../scripts/electronics/histories/run_pad_compression.py) | Pipeline and acquisition | Run the conservative nonlinear candidate/control study and generate its report. |
| [report_pad_compression.jl](../scripts/electronics/histories/report_pad_compression.jl) | Analysis and checks | Build numerical report tables and validate event compression accounting. |
| [present_pad_compression.py](../scripts/electronics/histories/present_pad_compression.py) | Presentation | Render event-loss fractions and affected-charge distributions. |
| [cache_standard_pad_histories.py](../scripts/electronics/histories/cache_standard_pad_histories.py) | Acquisition cache | Precompute the initial pad responses for the standard-event menus. |
| [run_pad_sequence_example.py](../scripts/electronics/histories/run_pad_sequence_example.py) | Pipeline and acquisition | Run the supplied ten-bucket sequence and paired controls. |
| [pad_sequence_example.jl](../scripts/electronics/histories/pad_sequence_example.jl) | Input generation and analysis | Build example currents or reduce its SPICE waveforms. |
| [present_pad_sequence_example.py](../scripts/electronics/histories/present_pad_sequence_example.py) | Presentation | Render the ten-bucket mean response and probe comparison. |
| [audit_history_time.jl](../scripts/electronics/histories/audit_history_time.jl) | Diagnostic | Audit source z, event arrival-time axes and diffusion-cell resolution. |
| [present_history_time_audit.py](../scripts/electronics/histories/present_history_time_audit.py) | Presentation | Render the timing and spatial-resolution audit. |

### ADC and circuit area

Directory: `scripts/electronics/design/`.

Practical digitisation and approximate layout calculations.

| Current file | Role | Purpose |
|---|---|---|
| [analyze_adc_estimate.jl](../scripts/electronics/design/analyze_adc_estimate.jl) | Analysis | Calculate ADC range and resolution estimates from saved study inputs. |
| [present_adc_estimate.py](../scripts/electronics/design/present_adc_estimate.py) | Presentation | Render the ADC estimate. |
| [build_pixel_area.jl](../scripts/electronics/design/build_pixel_area.jl) | Analysis | Calculate the device geometry and approximate per-pixel area budget. |
| [present_pixel_area.py](../scripts/electronics/design/present_pixel_area.py) | Presentation | Render the area budget and its layout assumptions. |

### Paper notes and collaborator packages

Directory: `scripts/paper/`.

Paper assembly, numerical-macro checks, the Ramo explanation and calculation handoff.

| Current file | Role | Purpose |
|---|---|---|
| [build_spice_paper.py](../scripts/paper/build_spice_paper.py) | Paper presentation | Generate manuscript macros, tables and figures from calibrated transistor results. |
| [audit_paper_macros.py](../scripts/paper/audit_paper_macros.py) | Audit | Check generated manuscript macro values against their numerical sources. |
| [check_ramo_note.jl](../scripts/paper/check_ramo_note.jl) | Validation | Calculate independent weighting-potential and current checks for the Ramo note. |
| [write_ramo_note.py](../scripts/paper/write_ramo_note.py) | Presentation | Generate the Spanish clarification of the Ramo calculation. |
| [package_rogelio_calculations.py](../scripts/paper/package_rogelio_calculations.py) | Packaging | Create collaborator ZIP packages from tracked calculation files. |

### Reference and checkpoint validation

Directory: `scripts/validation/`.

Petit reference fixtures and the staged event-analysis validation workflow.

| Current file | Role | Purpose |
|---|---|---|
| [generate_petit_reference.jl](../scripts/validation/generate_petit_reference.jl) | Reference generation | Regenerate graph and backbone fixtures from the original Petit code. |
| [generate_petit_diffusion_reference.jl](../scripts/validation/generate_petit_diffusion_reference.jl) | Reference generation | Regenerate pinned diffusion fixtures with an injected replay RNG. |
| [generate_petit_collection_reference.jl](../scripts/validation/generate_petit_collection_reference.jl) | Reference generation | Regenerate occupancy references from the original simple-RECO functions. |
| [validate_event_step1.jl](../scripts/validation/validate_event_step1.jl) | Validation | Validate real-file reading, arithmetic and deterministic reconstruction. |
| [report_event_step1.py](../scripts/validation/report_event_step1.py) | Presentation | Render the saved reader/reconstruction checkpoint. |
| [validate_event_step2.jl](../scripts/validation/validate_event_step2.jl) | Validation | Validate ion production and diffusion with numerical scans. |
| [report_event_step2.py](../scripts/validation/report_event_step2.py) | Presentation | Render the saved ion/diffusion checkpoint. |
| [validate_event_step3.jl](../scripts/validation/validate_event_step3.jl) | Validation | Validate the hits-to-diffusion-to-selection-to-pads chain. |
| [report_event_step3.py](../scripts/validation/report_event_step3.py) | Presentation | Render the saved collection and selection checkpoint. |
| [scan_event_analysis.jl](../scripts/validation/scan_event_analysis.jl) | Validation and analysis | Run the paired selection scans and event ensemble checkpoint. |
| [report_event_step4.py](../scripts/validation/report_event_step4.py) | Presentation | Render the saved paired scans and ensemble checkpoint. |

### Review documentation

Directory: `scripts/documentation/`.

The HTML review-guide builder and its imported content modules. Content modules are not standalone commands.

| Current file | Role | Purpose |
|---|---|---|
| [build_event_review.py](../scripts/documentation/build_event_review.py) | Documentation builder | Build the event-analysis HTML guide, Petit assessment and roadmap. |
| [event_design_content.py](../docs/event_analysis_review/content/event_design_content.py) | Imported content | Supply reviewed design narratives and interface proposals. |
| [petit_review_content.py](../docs/event_analysis_review/content/petit_review_content.py) | Imported content | Supply the source-based Petit assessment and adaptation choices. |
| [event_viewer_content.py](../docs/event_analysis_review/content/event_viewer_content.py) | Imported content | Supply the event-viewer implementation and interface documentation. |
| [event_review_voxels.py](../docs/event_analysis_review/content/event_review_voxels.py) | Historical content | Preserve the superseded early voxel-interface proposal. |

### Presentation asset

| Current file | Purpose | Proposed location |
|---|---|---|
| [templates/readout.html](../assets/readout/readout.html) | HTML template used by the analytical readout demonstration. | With the readout presentation assets, outside executable scripts. |
