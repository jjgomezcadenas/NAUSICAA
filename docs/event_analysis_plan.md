# Event analysis: implementation roadmap

The canonical design and decision record are in the
[HTML documentation](event_analysis_review/index.html). Read Purpose and physics,
Event flow, Reuse and development, Interface drafts, Configuration and outputs,
and Validation and roadmap. The Petit source assessment and the full decision
record remain embedded as appendices.

## Current status

Step 1 is implemented: the v3 hit reader preserves MC time and actual primary
IDs; the local Petit adaptation supplies Float32 binning, graphs, backbone and
truth paths. Both reduced samples are regenerated and verified. The HTML links
the code docstrings, regression evidence and generated checkpoint report.
Step 2 is implemented: ion production, thermal widths, explicit RNG sampling and
Petit Monte Carlo diffusion with full charge accounting. The numerical checkpoint
covers all events for ion production, 40 real-event diffusion cases, Gaussian
statistics and histogram resolution/padding scans. Results are in
results/event_analysis/step2_validation.h5 and its generated PDF report.
Step 3 is implemented: exactly-one-component selection, integer ion thresholds,
full retained-charge pad histories, explicit axes and selectable relative-time
order. The integration checkpoint saves 100 validation events, decisions and
accepted histories. The subset acceptance is recorded without tuning cuts.
Step 4 is implemented: central-pad ensemble alignment, exact integer numerators,
accepted/processed denominators, HDF5 round trips and twelve paired settings on
100 validation events. Selection is stable across 100, 200 and 400 histogram
bins; pad histories retain measurable resolution sensitivity. Review the working
alignment and required pad-response accuracy before the full-sample run.

Interactive viewer: the user has approved the two-tab MC energy / MC tracks
scope, Petit defaults, repository layout and explicit Apply controls. Viewer
step 1 interfaces were accepted as implementation guidance. Step 2 is implemented:
the local browser displays energy maps and occupied 3D voxels, with navigation,
camera controls and session restore. Launch with `make view-events`; chapter 7
contains usage and validation details. Browser review is accepted, including
event 2 with independent axis scales and independent projection zoom.
Step 3 is implemented: MC graph inspection uses the existing Petit adaptation,
with a component table, switchable layers, backbones and endpoints. The graph
and applied TOML settings are saved in HDF5 and reused across display rebins.
Next: review MC tracks, then interactive graph parameter controls (step 4).

The central flow is: read hits → produce ions → diffuse → retain the full
diffusion histogram → reconstruct and select exactly one connected track →
fill pad/time histories from the full histogram of accepted events. The
TOML-controlled graph threshold affects reconstruction. Pad histories retain
the accepted event's charge below that threshold.

Reconstruction uses b = 3 sigma_T and maximum connection distance = 6 sigma_T
throughout the drift range, with both factors in TOML. This replaces the earlier
fixed-bin proposal and removes Petit's 50 cm switch from the baseline. Physical
pad pitch remains independently fixed at 2 mm.

F12 is resolved: retain explicit Float32 edge vectors in computational interfaces
and the initial saved-axis schema. Compact axis storage can be considered later
without changing the computational interface.

## Development sequence

1. **Review contracts.** Review types, function signatures, ownership, configuration and saved schemas in the HTML. Settle the listed open contracts and approve the implementation scope.
2. **Input and local extraction · step 1 complete.** Extended reduced HDF5 with MC time and actual primary IDs. Copied the required Petit functions into src/event_analysis with licence and provenance, retaining their algorithms. Updated the root Julia dependencies.
3. **Deterministic reconstruction · step 1 complete.** Adapted voxelisation, graphs, backbone and truth reference to Float32. Added geometry, threshold and indexing tests. Compared copied stages with the pinned Petit reference.
4. **Ion production and diffusion · step 2 complete.** Implemented explicit RNG boundaries, common W-value, per-hit rounding, thermal widths and histogram accounting. Validated statistics, reproducibility and charge counts, with generated figures and a PDF checkpoint.
5. **Selection and pad histories · step 3 complete.** Connected exact ion-threshold reconstruction to single-track selection and the full retained diffusion table to pad histories. Validated grid phase, extended coverage, relative-time axes and exact charge conservation.
6. **Ensembles and paired scans · step 4 complete.** Save exact ensemble numerators, both denominators, aligned histories and provenance in HDF5. Compare numerical and selection settings on paired events. Python generates PNGs, tables and a compiled LaTeX report.
7. **Production review and full samples.** Review ensemble alignment and the measured histogram sensitivity, select the required pad-response accuracy, then run both full samples and measure runtime/memory.

## Documentation and review

The HTML includes the physical assumptions, data flow, Petit reuse map, complete
draft contracts, configuration example, HDF5 proposal, examples, test strategy
and open decisions. Review one type or function at a time. Extensive docstrings
are part of each implementation increment. Step-one numerical validation is saved in results/event_analysis/step1_validation.h5
and presented in latex_generated/event_analysis_step1/report.pdf.

Rebuild with `make event-review PYTHON=/path/to/python3.11-or-newer`.
The builder reads `docs/event_analysis_review/content/event_design_content.py` for the design and this
roadmap, `docs/event_analysis_review/content/petit_review_content.py` for the source assessment/decision
record, and `docs/event_analysis_review/template.html` for presentation.
The default Petit checkout is `../Petit`; the builder also accepts `--petit-root`.
Refresh the existing HTML tab after rebuilding.

The [previous plan](event_analysis_review/archive/event_analysis_plan_before_consolidation.md)
and [previous review](event_analysis_review/archive/review_before_consolidation.html)
are historical snapshots. Their earlier API sketches and workflow are superseded
by the current HTML design.
