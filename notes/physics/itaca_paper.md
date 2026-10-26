# NAUSICA for ITACA: paper plan and work status

Working plan agreed with J. J. Gómez-Cadenas. Updated 26 October 2026.

Manuscript: [nausica.tex](nausica.tex). This document records the intended paper,
the work already available, and what remains. It does not replace the manuscript
or authorise wholesale rewriting of the reviewed text.

## Current manuscript — 26 October 2026

The author approved replacing the finite-gain analytical electronics section
with the transistor SPICE study. This update supersedes the earlier status and
implementation descriptions below. The event/track section remains the next task.

- Abstract and Sections 1–2 retain their reviewed wording.
- Section 3 keeps the concept and focusing text. The frame duration is now
  explicitly 10 ms, the 300 us averaging window is introduced, and 1/Cf is
  identified as the ideal gain. Averaging needs repeated reads in the acquisition
  window or an averaging stage before the ADC; hardware scheduling remains to
  be selected in the complete design.
- Section 4 keeps the weighting potential, single-ion current and its figure,
  plus the ideal CSA explanation. The illustrative voltage plot and discharge-time
  table are removed from the active text. The sum of ion currents supplies SPICE.
- Section 5 is now **CSA response and noise**: circuit, uniform arrivals, frame
  readout, and noise/sensitivity. It uses 1.25 pF external pad capacitance plus
  intrinsic transistor capacitances, Cf=3 fF, 10 ms frames and 300 us averages.
  The 0.3 s correction is an estimator parameter, not an imposed circuit decay.
  Saved SPICE results give ENC=40.47 ions RMS and a verified threshold of
  203 ions for S/N > 5. Pulse injection stays in the code, outside the paper.
- Appendix A remains. The active finite-gain appendix has been replaced by
  `appendix_spice_readout.tex`, covering transistor dimensions, averaged-sample
  covariance and noise integration. The original `appendix_realistic_csa.tex`
  and original models remain available as references.
- New circuit, ion-waveform and noise figures, device table and value macros
  are produced by `scripts/build_spice_paper.py` / `make spice-paper` directly
  from the saved 1.25 pF dataset. The script checks the reproduced gain and
  noise, the adjacent-count threshold bracket and schematic connectivity/biases.
  Source hashes are in `generated/spice_paper_provenance.json`.
- The reviewed wording before this reorganisation is preserved in
  `archive/2026-10-26_pre_spice/`. New/reworked passages await author review and
  are not marked green. The old green list below records historical approval.

The 52 numerical tests and `make physics` pass. The asset generator reproduces
the saved gain/noise and threshold. Active LaTeX inputs, graphics, custom macros,
references and environment nesting were checked without compiling the manuscript.

### Remaining before final noise claims

The last manuscript plot has been replaced by ENC versus averaging duration
at fixed 1.25 pF, using the saved SPICE simulations. The new analysis follows
the Julia/HDF5/Python workflow (`make spice-averaging`). ENC falls from
41.80 to 40.18 electrons RMS between 180 and 350 us; at 300 us it reproduces
40.47 electrons RMS. Both signal calibration and filtered noise are recomputed
at each width. The raw/filtered noise spectrum remains a repository diagnostic.
This new paragraph and figure caption await author review.

The supplied-model execution reports zero transistor 1/f contributions. The
original 250 fF noise audit reproduces the scale of Rogelio's PDF; the diagnostic
nonzero-KF copy verifies that the reader detects flicker terms. Compare the
simulator/model settings and spectra with Rogelio. No arbitrary extra noise has
been added. The manuscript states the model condition beside the ENC result.

### Next section

Develop the ion signal from simulated event tracks, then reconstruct noisy
tracks with the selected response. Use the event reader/analysis conventions in
AGENTS.md (Julia and HDF5, preserving particle-physics units until the SI
calculation interface). This turn does not implement or modify the event work.
The response range, frame-history correction and digitisation scheme must be
checked with those inputs. The 1.25 pF pad load remains a design estimate.

## Earlier plan and review history

The following material records the original staged plan and earlier approvals.
For current circuit values, manuscript structure and completed work, use the
26 October status above.

## Purpose and scope

Present the conceptual design of NAUSICA and demonstrate, with a simplified
quantitative model, that its readout can reconstruct ion tracks above electronic
noise while preserving their characteristic topology.

The complete calculation is:

**Track energy deposits → ion drift and diffusion → arrivals at each pad →
induced current → CSA voltage → noisy digitised samples → reconstructed track.**

The paper should derive the important design choices from this calculation.
Signal formation is an essential part of the study, but the final product is a
complete conceptual readout design supported by reconstructed tracks.

The abstract states the intended successful outcome of the completed paper.
That outcome still needs to be established by the event-level calculation.
The initial performance study will use signal events throughout the fiducial
volume. Repeating the background selections with a full electronics response
would be valuable, but may be a separate paper; it is not a prerequisite for
the present study of track preservation.

NAUSICA may be used as part of MARS or to instrument the full cathode. The ITACA
paper presents MARS; the full-cathode alternative is introduced here. Do not
attribute both alternatives to the ITACA paper or develop the MARS mechanics in
this paper. A separate study may develop the alternative implementation.

## Approved green text — updated 24 October 2026

The author identifies the following prose as green and approved. This records
his assessment, not an independent Pangram test.

| Manuscript passage | Status |
| --- | --- |
| Abstract | Green; preserve the current wording |
| Section 1, Introduction | Green; preserve the current wording |
| Section 2, The ITACA detector, including Principle of operation | Green; reviewed and approved |
| Section 3 opening paragraphs, before Charge collection and focusing | Green; author-revised and approved |
| Section 3, Charge collection and focusing | Green; approved through the end of this subsection |
| Section 3, Pixel electronics and readout | Green; author-revised and approved on 24 October; preserve current wording |
| Parameters of the design and Section 4 | Not covered by this green approval; continue review |
| Section 5 opening paragraph, before Finite gain and response time | Green; author wording with approved typo corrections and removal of repeated gain/bandwidth phrase; retain final SPICE sentence |
| Section 5, Finite gain and response time | Green apart from the newly authorized Topmetal-S capacitance paragraph; full active subsection in nausica.tex approved, including the gain figure, current-to-voltage response and Laplace-transform derivation of the integration limit; preserve current wording |
| Section 5, Response to uniform arrivals | Replaces the previously green two-case example at the author's explicit request; fixed 5 pF, free N and normalized figure/table. Revised text awaits renewed green approval |
| Section 5, Noise in the recovered charge | Green with authorized N notation and strict S/N threshold update; full active subsection approved, including the correction defining $g_n$ as the mean estimated charge per ion; preserve current wording |
| Remaining Section 5 reference draft and appendices | Draft; continue author review |

The segment approved on 23 October starts at `\subsection{Principle of operation}`
and ends immediately before `\subsection{Pixel electronics and readout}`.
Checkpoint commit `3cf0db3` preserves the manuscript and comments before cleanup.
Comment-only lines in this interval were subsequently removed; the active text
was preserved exactly. Comments outside the interval remain for later review.

Handle green prose with care. Do not insert longer sentences, reformulate it,
or change its meaning without the author's approval. As clarified on 24 October,
all manuscript edits, including typos and small word-level corrections, require
explicit instruction. A request to produce a figure does not authorise manuscript
edits or insertion of that figure. A green designation does not
mean pending calculations, including the Appendix C collection study, have
already been completed.

## Agreed working assumptions

- The positive carrier is molecular Xe_2^+, with charge +e.
- The detector context is xenon at 15 bar, a bulk field of 200 V/cm and a
  maximum drift length of 1.5 m. The nominal ion drift speed is 0.1 m/s.
- With the pads at ground, the upper IFG electrode is at +8 kV and the lower
  electrode at +3 kV. There are 5 kV across the 10 mm IFG and 3 kV across the
  5 mm gap to the pads: nominal fields of 5 and 6 kV/cm respectively.
- Describe the focusing geometry and voltages, and assume collection efficiency
  approximately unity for the initial track study. State this once as a model
  assumption. A detailed focusing-field solution is not required before the
  first track calculation, and has not been completed.
- The existing single-ion calculation uses a square pad in an infinite grounded
  plane, with a central trajectory and constant velocity. This is a stated
  analytical approximation, not a solution of the complete IFG geometry.
- Physical fields determine trajectories; weighting fields determine induced
  signals. Keep the distinction and the relevant boundary conditions explicit.
- Preserve integrated signal charge, q times the change in weighting potential.
  Do not add a second collection impulse or renormalise a finite-path integral.
- Use sigma for spatial diffusion and resolution. The reviewed detector section
  describes an average spread of order 1 mm. Retain the actual drift-distance
  dependence in the calculation. The energy ROI is ±1 sigma; the quoted NEXT
  energy-resolution measurement remains expressed in FWHM.
- Separate total bulk drift time, event arrival span, arrival duration on one
  pad, single-ion pulse duration, and feedback discharge time.

## Proposed paper structure

### Abstract

**Status:** Reviewed and accepted by the author. Preserve the current text.
It already describes the full track-to-readout study and its intended result.
Do not reopen the CSA wording without a request from the author.

### 1. Introduction

**Purpose:** Explain the physics motivation, the limitations of electron-track
imaging, the ITACA proposal, and the specific contribution of NAUSICA.

**Completed:** Reviewed with the author. ITACA is named before NAUSICA; the
second ITACA paper is cited for the earlier topology study. Spatial resolution
uses sigma. The author wishes to retain the strong statement about preserving
the information needed for topology discrimination.

**Remaining:** Complete the final section roadmap when the new structure is in
place. Avoid another general rewrite.

### 2. The ITACA detector

**Purpose:** Supply the detector conditions and operating sequence needed to
understand the ion readout.

**Completed:** Reviewed sentence by sentence. Pressure dependence, scintillation
and t0, diffusion notation, and the role of the three event coordinates have
been clarified. The text explicitly allows NAUSICA within MARS or across the
full cathode. The sentence assuming detector movement and a particular gate
window is commented out. The parameter table is retained.

**Remaining:** Keep this section stable while developing the quantitative
sections. Revisit only inconsistencies revealed by the calculations.

### 3. The NAUSICA concept

**Purpose:** Introduce the device before discussing numerical specifications.

Start with the need to measure ion charge without gas amplification. Explain
the TopMetal choice, with the existing references, in a short paragraph. Then
describe the collection geometry and the pad–CSA–sampling–ADC readout chain.
Give the focusing voltages and the near-unit collection assumption.

Identify the quantities the study must determine: pad dimensions, pitch,
instrumented area, feedback capacitance and discharge time, sampling interval,
noise, and digitisation. Distinguish nominal geometries used for calculations
from parameters ultimately selected by the study.

**Available:** The present section contains the basic concept, references,
geometry and architecture. It is largely inherited from ITACA and must be
reorganised. The existing TopMetal image is only a starting point.

**Drafted on 23 October:** Section 3 now follows the agreed concept, collection,
electronics and design-parameter sequence. Existing wording has been reused;
the original section remains archived as comments at the end of the source,
after `\end{document}`.
The opening and focusing text have since been reviewed and approved as green;
the remaining subsections await review. The schematic was subsequently replaced; see the figure update below.

**Remaining:** Continue author review from Pixel electronics and readout. Keep its distinction between pad width and pitch, nominal geometry
and calculated fields, and equivalent circuit and transistor implementation.

Move the track extent and occupancy figures out of this section. Do not open
with the current list of numerical requirements or announce inherited
electronics values as final specifications.

**Figure update, 23 October:** A new original vector schematic now replaces
the inherited image. It shows the IFG geometry and nominal voltages beside
an equivalent CSA with feedback capacitance, effective MOS resistance,
sampling and ADC. Editable source, vector PDF and PNG preview are in
`figures/`; the previous figure block is preserved in checkpoint `3cf0db3`.
The replacement is included in the approved segment. The author-revised
opening paragraphs were left byte-for-byte unchanged. Treat paragraphs
the author identifies as green as protected prose: obtain approval before
introducing longer sentences or substantial wording changes.

### 4. Signal formation and response of one pad

**Purpose:** Follow one ion from its approach to the pad to the CSA output.

Keep the essential weighting-potential result, induced-current expression,
integrated charge and pulse timescale in the main text. State the approximation
and move the detailed derivation to an appendix. Introduce the analytical CSA
response here, including finite feedback discharge, and show current and voltage
together. This defines the response to be summed over arrivals in a real event.

**Available:** The single-ion physics model and the linear CSA model are
implemented, with physical and numerical tests. The existing paper contains
the weighting-potential derivation and current results. An offline analytical
readout interface and example current files are available.

**Drafted on 23 October:** Section 4 and Appendix A are now in the manuscript,
awaiting author review. The section gives the induced current, finite-gap charge,
pulse duration and linear CSA response. The current figure has full-gap and near-pad panels, with distance and
time-to-collection axes and the final 90% of full-path charge shaded. A separate
figure shows CSA voltage, and a generated table compares illustrative discharge times.
Appendix A contains the square-pad derivation, charge-duration definitions,
CSA solution, charge recovery between samples and numerical checks. The model
is the existing grounded-plane benchmark, not an IFG field calculation.

The readout generator now creates `generated/readout_values.tex`,
`generated/single_ion_current.pdf`, `generated/single_ion_response.pdf` and
`generated/single_ion_response_table.tex`. `make note` regenerates these as well
as the physics outputs. All 22 existing tests passed. The three finite-feedback
examples were also checked against independent continuous quadrature; their
largest relative peak-voltage difference was 0.046%. No PDF visual inspection
was performed.

**Approved corrections and figure integration, 24 October:** Retained the author's
wording and full-path footnote. Fixed the listed grammar and LaTeX errors,
kept the earlier gap definition and removed the duplicate, and retained the
approximate 200 microseconds without the duplicated precise value. Inserted the
two-panel current figure before the CSA subsection, added its text reference,
and corrected the separate voltage figure caption and table reference. The
claim of validation in Appendix C is commented out with a TODO until that study
exists. These specific edits were explicitly authorised after review; they do
not authorise further manuscript edits.

The previous manual-like analytical sections are preserved as comments after
`\end{document}` for reuse in the later track study and Appendix C. Sections
1–3 and the abstract were left unchanged by this drafting step. Appendix C's
quantitative collection study remains to be done.

**Remaining:** Review Section 4 and Appendix A with the author. Keep example
capacitances and discharge times separate from the eventual design choices.

### 5. Beyond the ideal CSA model

**Editing location:** The author is incorporating the reviewed wording directly
in `nausica.tex`. `section_beyond_ideal.tex` is the unchanged reference draft;
its input is currently commented out. Apply approved corrections to the active
text in `nausica.tex`, not to the reference draft. The active section now contains the opening paragraph, finite-gain discussion,
frame-charge derivation, the uniform-arrival study with free N and the noise subsection. The noise subsection explains voltage fluctuations, propagation
through the estimator, sample covariance and ion-count sensitivity. Its equation
labels now supply the references used by Appendix B.

**Purpose:** Calculate the response and sampled-charge noise for specified ion
arrivals at one pad, extending the ideal CSA to finite gain, bandwidth and input
capacitance. Keep signal and noise transfer functions consistent.

**Drafted, 24 October:** [section_beyond_ideal.tex](section_beyond_ideal.tex) is
retained unchanged as the reference draft; reviewed text is incorporated directly
in the active manuscript. [appendix_realistic_csa.tex](appendix_realistic_csa.tex)
is Appendix B; the pending focusing appendix is now C. The section defines
controlled N-ion arrivals, a two-sample frame-charge estimator, its bias and its
noise covariance. The implementation and independent checks are complete. The uniform mean response per ion and numerical noise threshold are included
in the active manuscript. Simultaneous mid-frame arrivals appear only as a
brief cross-check. Unaffected reviewed prose is preserved.

**Gain figure:** `scripts/build_gain_figure.py` generates the normalized
open-loop magnitude plot from `CSA.amplifier_gain`, with logarithmic f/fa and
linear |A|/A0 axes and the 1/sqrt(2) point marked. `make note` regenerates
`generated/amplifier_gain.pdf`; no particular amplifier performance is implied
by this normalized curve.

**Current uniform-arrival study:** The author requested a single main case,
N ions arriving uniformly over a 10 ms frame, with N left free until the noise
calculation. Input capacitance is fixed at 5 pF from the start, citing Topmetal-S.
`uniform_ion_waveform` in `src/ion_signal/realistic_csa.py` calculates the mean
current per ion from a continuous uniform arrival density and the finite-gap
single-ion current. `scripts/build_arrival_cases.py` produces a two-panel
current/voltage figure, a charge table normalized by N, macros and waveform CSVs.
`arrival_cases_point_samples.json` records the uniform response and the
simultaneous mid-frame cross-check. `make noise-study` regenerates both the
signals and noise, without compiling LaTeX. The earlier discrete 1000-ion
benchmarks remain diagnostic cases in `make realistic-readout`, not the basis
of the paper's main example. The archived `section_beyond_ideal.tex` is unchanged.

**Implementation:** A delegated agent has rewritten the useful circuit assumptions
from Rogelio's HTML in Python, preserving the original as provenance. Shared
specification: [realistic_csa_spec.md](readout/realistic_csa_spec.md). Parameters
are illustrative until supported by a device model. The rewrite removes the
arbitrary resistance cap, unsupported fixed stability number and inconsistent
noise formulas rather than treating the HTML output as validation.

**Available outputs:** Input currents, finite-gain and ideal voltage waveforms,
true and recovered frame charge, noise RMS and adjacent-sample covariance,
template-calibrated S/N and conditional ion-count sensitivity. The initial
implementation passed 37 tests, including 15 circuit/noise tests covering physical
charge, circuit limits, numerical convergence and independent noise integrals.
The parent independently checked the configured noise against a direct positive
spectral integral. LaTeX compiles without warnings; no PDF visual check was made.
Do not infer a track threshold or foundry performance from these benchmarks.

**Noise subsection — green, 24 October:** Reviewed and approved by the author,
including the definition of $g_n$ as the mean estimated charge per ion. Added directly to `nausica.tex`, following the
approved two-case study. It keeps the frame-boundary readings used there and
starts from the fluctuation in each voltage reading before deriving the charge
variance. No new noise numbers were inserted: the existing noise benchmarks use
a 10 microsecond aperture, so the point-sample numerical study must be presented
with its stated source spectra and frequency band. The reference draft stays
unchanged. No LaTeX compilation was run.

**Current point-sample noise calculation:** Reproduce with `make noise-study`.
The script `scripts/build_noise_study.py` reads the uniform mean signal per ion
and computes the noise with the same fixed 5 pF circuit and frame-boundary
readings (no aperture). Noise is 5.02615 e RMS in the nominal charge estimator;
white voltage, flicker and feedback contributions are 4.14811, 1.88610 and
2.12085 e RMS. Amplifier voltage noise accounts for 82.19% of the variance.
The uniform mean estimated charge per ion is 0.374487 e, giving S/N=0.074508 N.
The strict condition S/N > 5 requires N > 67.107, hence at least 68 ions.
The simultaneous mid-frame cross-check requires 67 ions. These are mean-signal
thresholds for the assumed spectra, not detection probabilities for individual
random arrival histories. Finite-gain recovery is about 38%; point-sample ideal
recovery is 99.99% after rounding. All quoted numbers and tables are generated.
The two-case noise sensitivity table is no longer included; the main result is
an equation and ion-count threshold in the prose. All 41 tests and physics checks
pass, including independent uniform-arrival integration, conserved finite-path
charge, scaling, timestep refinement and independent noise integration. Maximum
noise quadrature comparison discrepancy is 1.3e-11 relative. No LaTeX compilation.

**Historical capacitance comparison (superseded in the active paper):** At the author's request, use the
published 5 pF value directly before considering any area scaling. The source is
Mei et al., arXiv:2010.09226, Section 2.3: approximately 5 pF to ground for a
1 mm diameter electrode. In this benchmark 5 pF replaces total Cin, excluding
Cf, without an additional unknown transistor capacitance. Other circuit and
noise assumptions remain fixed. `make capacitance-check` runs the tracked
`scripts/build_capacitance_check.py`; JSON and Markdown results are saved as
`results/realistic_csa/capacitance_check.*`. Point samples, 10 ms frames, the
configured 0.1 Hz–1 MHz noise band, passive-resistor feedback noise and the
central grounded-plane pad trajectory are retained. Recovered fractions are
37.88% / 37.84% (simultaneous / uniform); S/N for 1000 ions is 75.25 / 74.51;
rounded template thresholds at S/N = 5 are 67 / 68 ions. Noise in the nominal
charge estimator is 5.026 e RMS, with white voltage / flicker / feedback terms
4.148 / 1.886 / 2.121 e RMS. Amplifier voltage noise now dominates. The similar
numerical charge noise does not mean unchanged sensitivity: the signal is
smaller because of finite gain. Signal timestep refinement changes results by
less than 1.6e-6 relative; independent noise integration agrees within 8e-14.
The optional `make capacitance-check` remains a historical diagnostic comparing
0.25 and 5 pF with 1000-ion packets. Its generated values are no longer loaded
by the paper, and it is not part of `make note`. The active paper instead uses
5 pF in `config/realistic_csa.toml` throughout. No LaTeX compilation was made.

The notebook containing the source papers is NEXT & NEXEL & ITACA:
https://notebook.google.com/notebook/8beb347c-cdfd-40c6-a101-fb9a575faa82 .
Its source PDFs are also local in `/Users/jjgomezcadenas/Dropbox-Personal/2026/ITACA/`.
Liang et al., JINST 19 (2024) C04004, Section 2.2 reports a corrected guard-ring
injection capacitance of 255 +/- 10 fF (not total detector capacitance) and ENC
120–140 e; earlier Topmetal-S noise reports were underestimated by calibration.
Do not equate that measured ENC to this calculation's estimator noise.

**Active-paper provenance audit, 24 October:** Reviewed source and generators;
no compilation or visual PDF inspection. The revised uniform-response plot was
checked separately as a PNG; the manuscript PDF was not inspected. At the time of the audit the active paper had six figures
(seven image files) and four tables. The numerical noise subsection now adds only the noise-input table, generated
by `scripts/build_noise_study.py`; the sensitivity comparison table was removed.
The active paper therefore has six figures and five tables. The text explicitly
states that these inputs are not derived from the actual NAUSICA circuit. Commented archival material after the
document is excluded.

| Active item | Tracked source / generator | Status |
| --- | --- | --- |
| Figure 1: two ITACA diffusion panels | `img/SOURCES.txt` and tracked PNGs | Copied assets; recorded hashes verified. No Python, Julia or notebook generator found in the local ITACA manuscript directory. Generator/data still needed for local reproduction. |
| Figure 2: IFG and pixel schematic | `figures/nausica_concept.tex`, `figures/build_nausica_concept.sh` | Reproducible TikZ drawing; numerical labels are literal, not imported from detector config. |
| Figures 3–4: single-ion current and ideal CSA voltage | `scripts/build_readout.py` | Generated from tracked physical model and config. |
| Figure 5: amplifier gain | `scripts/build_gain_figure.py` | Generated from tracked amplifier model. |
| Figure 6: uniform mean response per ion | `scripts/build_arrival_cases.py` | Generated with point samples and recorded signal provenance. |
| Table 1: ITACA detector parameters | Inline in `nausica.tex` | Inherited design values and derived quantities mixed; no table generator. Needs a sourced config and generated derived quantities. |
| Table 2: single-ion CSA response | `scripts/build_readout.py` | Generated table. |
| Table 3: arrival/circuit parameters | `scripts/build_arrival_cases.py` | Table layout inline; all numerical entries use generated macros. |
| Table 4: input and recovered charge | `scripts/build_arrival_cases.py` | Generated table. |
| Section 4 / Appendix A numerical results | `scripts/build_physics.py`, `scripts/build_readout.py` | Generated macros, except rounded 200 microseconds typed in prose (computed duration 193.51 microseconds). |
| Section 5 recovery fractions | `scripts/build_arrival_cases.py` | Revised uniform example uses generated macros for the ideal and finite-gain recovery fractions. |
| Numerical noise threshold, revised after the audit | `scripts/build_noise_study.py` | JSON/CSV/Markdown and LaTeX values/table generated; build dependency included. |

Other manually entered derived results include the 0.45 mm flux-based minimum
hole diameter; detector volume, anode voltages and drift times; the mapping of
1–2 mm to 10–20 ms. The former 1000-ion frame estimate has been removed
from the main calculation. Some
underlying drift calculations already exist, but these prose values do not use
generated macros. Published numbers (NEXT resolution, decay energies, ITACA
rejection) require citations rather than a new calculation. Design inputs
require recorded sources. Abstract claims about simulated and reconstructed
tracks remain goals; the event-level calculation is still outstanding.
No green prose or figure content was changed during this audit.

**Specs to SPICE is a distinct task:** [specs_to_spice.md](readout/specs_to_spice.md)
defines the numerical package and comparison protocol for Rogelio. It includes
charge-checked PWL currents, circuit connections/parameters, sampling, transient
convergence and noise integration. A behavioral netlist tests the analytical
model; Rogelio's transistor calculation tests its physical adequacy. No message
or package is sent without an explicit user instruction.

### 6. The ion signal from a track

**Purpose:** Establish the spatial and temporal input the electronics must handle.

Describe the Geant4 event sample, energy deposits, ion production, transport,
diffusion and assignment to pads. Calculate containment, charge per pad, charge
per sampling interval, arrival duration, and repeated arrivals at a pad. Use
the single-ion response to obtain pad currents and voltages.

**Available:** ITACA figures showing track extent, occupancy and per-frame
signals have been copied locally. Bulk transport and diffusion calculations
exist. Specified one-pad arrival distributions can already be convolved with
the single-ion current and CSA response.

**Remaining:** Obtain the underlying event-level numerical data or a reproducible
simulation producing them. Imported figures are not a substitute for energy
deposits and positions. Implement event transport and pad assignment, including
the stated collection assumption and appropriate fluctuations. Recalculate the
distributions for each pitch and sampling interval being studied. Preserve
orientations, scattering, nonuniform deposition and endpoint blobs.

The current illustrative 260-ion packet is not a simulated event or proof of a
charge-per-frame distribution. A straight 10 cm track along drift may be used
as a diagnostic extreme, not as the event ensemble.

### 7. The NAUSICA design

**Purpose:** Select and assemble the conceptual design using Sections 4--6.

Relate pitch and sampling to image quality; capacitance to gain and range;
feedback discharge to charge retention and temporal reconstruction; and noise
and digitisation to recoverable track signals. Treat these choices together.
A sampled CSA voltage includes the history of arriving charge and is not simply
the charge delivered during the latest frame.

Text moved from Section 4 for use here when choosing feedback time and sampling:

> The single-ion current lasts much less than a readout frame. The feedback must
> retain its contribution until the voltage is sampled. The discharge time must
> therefore be considered together with the sampling interval and the arrival
> distribution on each pad.

The Section 4 superposition paragraph now introduces multiple arrivals at a pad,
separately from retention between samples. Section 7 is not yet in the active
manuscript, so the moved text is retained here for its draft.

Describe the proposed circuit implementation and compare the analytical model
with Rogelio's circuit results when available. Finish with a justified parameter
table, multiplexing and readout architecture, and a credible tiling concept.

**Available:** Candidate values and architecture are in the current text.
Analytical scans over feedback time and sampling phase already exist. Rogelio's
supplied calculator provides an initial circuit concept.

**Remaining:** Establish the noise model, including temporal correlations,
and the usable voltage range. Compare realistic input currents in the
analytical and SPICE models. Derive the parameter choices rather than retaining
the copied specifications by default. Clarify what is supported by circuit
calculations and what remains a design target.

The existing 16 × 16 cm² area can describe a track-sized readout unit; it must
not silently fix the architecture to MARS. Full cathode instrumentation would
use a different number and arrangement of units. The inherited assertion of
dead-area-free tiling needs an actual layout justification.

### 8. Track reconstruction and performance

**Purpose:** Demonstrate that the selected design preserves the ion track.

Explain baseline treatment, recovery of charge from sampled CSA voltages,
thresholding and three-dimensional reconstruction. Include electronic noise
and digitisation consistently with the selected readout model.

**Available:** The physical and analytical electronics components provide the
starting point. There is no completed noisy event reconstruction yet.

**Remaining:** Implement and validate the event-level chain. Present original,
transported and reconstructed tracks, with representative pad waveforms.
Quantify charge recovery, false hits, continuity, spatial resolution and
endpoint-blob preservation across drift distances and orientations. Study
sensitivity to noise and the most important design parameters, including
collection losses if needed to establish operating margin.

Use a circuit-supported analytical response for event production where possible;
running SPICE for every pad of every event is not required. A full comparison of
signal and background selections can be developed later.

### 9. Conclusions

**Purpose:** State the resulting conceptual design and demonstrated imaging
performance, followed by the specific work needed for implementation.

**Remaining:** Write after the calculations and design choices are established.
Separate achieved results from further electrostatic, circuit and experimental
validation. Avoid repeating the development plan in the published conclusions.

## Where Rogelio's work enters

| Stage | Role | Present status |
| --- | --- | --- |
| Section 4 | Ideal current-to-voltage response | Implemented and in the manuscript |
| Section 5 | Finite gain, bandwidth and sampled-charge noise | Draft and tested Python rewrite complete; numerical examples await author review |
| Section 7 | Transistor-level design and SPICE response to realistic pad currents | Circuit calculations/results still need to be supplied or connected and assessed |
| Section 8 | Achievable response and noise parameters for event reconstruction | Pending circuit comparison and noise model |

The exchange is: our track model supplies realistic input currents; Rogelio's
circuit work establishes the achievable response; the track reconstruction
tests whether that response preserves the image.

The supplied HTML contains analytical approximations, not a SPICE engine or a
foundry-validated circuit simulation. It is preserved unchanged as
[rogelio_reference.html](readout/rogelio_reference.html). The adapted
[analytical model](readout/analytical_model.html) replaces its calculations with
the repository's tested physical-current and linear-CSA models. This distinction
does not imply that Rogelio has no other preliminary calculations; those must
be examined separately when available.

## Existing implementation and artifacts

Paths below are relative to this plan.

| Item | Location and status |
| --- | --- |
| Finite-gain CSA and noise | [realistic_csa.py](../../src/ion_signal/realistic_csa.py), [configuration](../../config/realistic_csa.toml), [benchmark outputs](../../results/realistic_csa/benchmarks.json) |
| Specs to SPICE | [handoff recipe](readout/specs_to_spice.md), [PWL inputs and behavioral netlists](../../results/realistic_csa/spice/) |
| Single-ion physics and transport | [physics.py](../../src/ion_signal/physics.py) |
| Arrival benchmarks and linear CSA | [readout.py](../../src/ion_signal/readout.py) |
| Numerical inputs | [detector.toml](../../config/detector.toml), [readout.toml](../../config/readout.toml) |
| Reproducible generators | [build_physics.py](../../scripts/build_physics.py), [build_readout.py](../../scripts/build_readout.py) |
| Physical/numerical tests | [test_physics.py](../../tests/test_physics.py), [test_readout.py](../../tests/test_readout.py) |
| Generated model results | [physics_benchmarks.json](../../results/physics_benchmarks.json), [readout benchmarks](../../results/readout/benchmarks.json) |
| SPICE input examples | [PWL current files](../../results/readout/spice/) |
| Analytical interface and assumptions | [readout README](readout/README.md), [offline interface](readout/analytical_model.html) |
| Imported ITACA figures | [img](img/), with [provenance](img/SOURCES.txt) |
| Bibliography and notation | [references.bib](references.bib), [itaca_defs.tex](itaca_defs.tex) |

The existing CSA is an ideal amplifier with linear feedback, obeying
`dV/dt + V/tau = -I/Cf`, with `tau = Rf*Cf`. It does not currently predict
transistor noise, finite gain or bandwidth, saturation, power or stability.
The arrival benchmarks include a single ion, a uniform packet, a diffused
compact packet and two separated packets. They are not complete tracks.

Tests cover boundary conditions, independent quadrature, charge conservation,
scaling, numerical convergence and analytical CSA limits. The physical and
readout tests passed in the earlier implementation work; the latest manuscript
build after the Section 2 review also passed without warnings. No tests are
being newly claimed by writing this plan.

## Candidate parameters that still require justification

| Parameter | Existing candidate or example | What must establish the choice |
| --- | --- | --- |
| Pad size and pitch | 0.5 mm square pads, 2 mm pitch | Collection model, capacitance and track sampling |
| Readout unit area | 16 × 16 cm² | Track containment and chosen deployment |
| Feedback capacitance | 5 fF in inherited text; 3 fF in the adapted calculator | Gain, range, noise and circuit feasibility |
| Feedback discharge time | Approximately 1 s in inherited text; a broad scan in the model | Pad arrival distributions, reconstruction and circuit behaviour |
| Sampling interval | 10 ms | Longitudinal image quality and noise after charge recovery |
| Noise | 50-electron RMS ENC target in inherited text | Defined noise bandwidth/correlations, circuit results and track recovery |
| Digitisation | 12 bits over approximately 1 V | Actual output range, noise, clipping and quantisation |

In particular, a one-second feedback time is not the full bulk drift time, and
a 10-ion cut used in an occupancy plot is not an established electronics
threshold. These inherited statements must be resolved in the quantitative
design section rather than repeated as specifications.

## Disposition of the current later sections

**Focusing update:** The main text now introduces the uniform-field flux
estimate, which motivates 0.5 mm holes at 2 mm pitch with the existing
5 kV drop across the 10 mm IFG. Higher voltages are not the preferred
design direction. Appendix C is reserved for the quantitative focusing
and collection study that will test the near-unit-efficiency assumption.
That calculation remains pending; it is distinct from the flux estimate.
Keep weighting-potential derivations in Appendix A and place any separate
numerical-validation appendix after Appendix C.

| Current material | Destination |
| --- | --- |
| Scope and definitions | Essential definitions at first use in Sections 3–6 |
| Detector definition and voltages | Concise geometry in Section 3; detailed inputs where needed |
| Focusing entrance argument | Brief estimate in Section 3; quantitative collection study in Appendix C |
| Transport and diffusion | Section 6, referring back to Section 2 instead of rederiving it twice |
| Weighting potential and single-ion signal | Essential results in Section 4 |
| Central-pad derivation and shape comparison | Appendix A, with selected results retained in Section 4 |
| Pulse-duration definitions and numerical checks | Duration derivations in Appendix A; separate validation appendix if needed |
| Implications for the next physics stage | Replace with the actual event calculation in Section 6 |
| Validation, reproduction and remaining work | Scientific validation after Appendix C if needed; commands and implementation instructions in repository documentation |

Keep model assumptions and boundary conditions beside the results they qualify.
Do not move all limitations to an appendix. Remove repeated warnings, historical
design discussions, configuration hashes and development instructions from the
paper narrative. Retain useful technical material in the repository.

## Sources and outstanding inputs

- Local ITACA source: `/Users/jjgomezcadenas/Papers/itaca/`, particularly
  `itacaRevisited.tex`, `srcx/`, the bibliography in `pool/`, and the figure
  directories. Reuse selected assets locally so the manuscript stays portable.
- Main reference: `itaca2026`; first ITACA paper: `Gomez-Cadenas:2025ect`.
  Other relevant keys include `Nygren_2018`, `mei2020topmetalcmosdirectcharge`,
  `An2016TopmetalII` and `You2021`.
- Event data: locate the numerical deposits/positions underlying the ITACA
  study or establish how to regenerate them. The imported plots alone do not
  provide this input.
- Circuit inputs: obtain Rogelio's applicable netlists, device assumptions,
  simulation conditions and output/noise results. Existing PWL files provide
  initial physical input currents for comparison; event-derived examples follow
  from Section 6.
- The title currently uses NAUSICAA while the expanded acronym in the body uses
  NAUSICA. Settle the final spelling with the author rather than silently
  changing the reviewed text.

## Working method and next step

Section 5 and Appendix B are drafted, and the parallel Python implementation
has passed its tests and independent review. The next task is author review of
the section and illustrative numerical cases. The specs-to-SPICE handoff is
generated but has not been executed in SPICE; the transistor comparison remains
a separate task. Author review of Appendix A also
remains pending; the event-data study follows as Section 6. The author
has shortened the electronics discussion and closed Section 3 with its
parameter summary. The discarded roadmap paragraph has been removed and
the commented original Section 3 moved to the end of the source. Preserve
the reviewed wording; do not reorganise the entire manuscript in one pass.

Preserve the author's voice. The author reports the abstract and reviewed
Sections 1–2 as acceptable to Pangram; this plan does not independently assess
that score. Do not edit manuscript text without explicit instruction, including
minor corrections. Review requests are read-only. Producing a figure does not
authorise inserting it or changing prose. For substantive edits, show the
original and replacement one at a time,
wait for approval, apply the agreed change, and then proceed. Comment out
proposed removals for review before deleting material when appropriate.

**Writing rule explicitly requested by the author:** Use plain wording in the
first draft, rather than waiting for the author to ask for a simpler version.
Prefer familiar words and one idea per sentence. Explain the physical meaning
and connections between equations, but do not repeat their mathematics in
elaborate prose. Technical precision does not justify unnecessarily complicated
language. For example, prefer "following the same trajectory" to "following the
same approach", and "The integral adds the voltage signals produced by charge
delivered at different times" to "Each charge increment contributes a delayed
copy of this response". Preserve the author's approved wording. This rule does
not authorise unsolicited manuscript edits.

Keep physical calculations in `src/ion_signal`, use SI internally, and generate
paper values and figures reproducibly. All calculations must be documented in
Python scripts tracked by git, not left in scratchpads (author instruction,
24 October). Run the relevant physical tests when
changing calculations. Compile LaTeX only when the author explicitly requests it; do not compile
automatically after manuscript edits. When compilation is requested, inspect the
build log. Do not visually inspect the PDF unless explicitly requested. Keep scratch
work under `work/`. Do not commit or push without a new user request.

Update this plan as sections and calculations are completed, preserving the
distinction between inherited ITACA results, new analytical benchmarks, circuit
results, and demonstrated event-level performance.

### 2026-10-24: shaping checkpoint and capacitance comparison

- Saved the complete 180–350 us shaping study in commit `58f2fad` on LTSpice.
- Ran the authorized 1.25 pF external-load estimate at 300 us, retaining device
  capacitances and the existing 0.5 mm pad current. See
  `results/ltspice/capacitance/pad_1.25pF/README.md` and its reproducing script.
- Area scaling is a design assumption; a 0.5 mm pad capacitance is not yet
  established by an electrostatic/layout calculation. The model has zero
  flicker noise. All 52 numerical tests and the physics checks pass.
- No paper prose changed or LaTeX compilation performed.
