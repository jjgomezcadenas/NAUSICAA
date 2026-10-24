# NAUSICA for ITACA: paper plan and work status

Working plan agreed with J. J. Gómez-Cadenas. Updated 24 October 2026.

Manuscript: [nausica.tex](nausica.tex). This document records the intended paper,
the work already available, and what remains. It does not replace the manuscript
or authorise wholesale rewriting of the reviewed text.

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
| Section 5, Finite gain and response time | Green; full active subsection in nausica.tex approved, including the gain figure, current-to-voltage response and Laplace-transform derivation of the integration limit; preserve current wording |
| Section 5, Response to simultaneous and uniform arrivals | Green; full active subsection approved, including the parameter table, current/voltage figure, recovered-charge table and concluding recovery fractions; preserve current wording |
| Section 5, Noise in the recovered charge | Green; full active subsection approved, including the correction defining $g_n$ as the mean estimated charge per ion; preserve current wording |
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
frame-charge derivation, the two arrival-case study and the approved noise subsection. The noise subsection explains voltage fluctuations, propagation
through the estimator, sample covariance and ion-count sensitivity. Its equation
labels now supply the references used by Appendix B.

**Purpose:** Calculate the response and sampled-charge noise for specified ion
arrivals at one pad, extending the ideal CSA to finite gain, bandwidth and input
capacitance. Keep signal and noise transfer functions consistent.

**Drafted, 24 October:** [section_beyond_ideal.tex](section_beyond_ideal.tex) is
retained unchanged as the reference draft; reviewed text is incorporated directly
in the active manuscript. [appendix_realistic_csa.tex](appendix_realistic_csa.tex)
is Appendix B; the pending focusing appendix is now C. The section defines
controlled k-ion arrivals, a two-sample frame-charge estimator, its bias and its
noise covariance. The implementation and independent checks are complete. The point-sample two-case results are now included in the active manuscript.
The noise derivation is now in the paper and approved as green. Numerical noise
results and the other benchmark cases remain outside it pending author review.
Existing reviewed manuscript prose is unchanged.

**Gain figure:** `scripts/build_gain_figure.py` generates the normalized
open-loop magnitude plot from `CSA.amplifier_gain`, with logarithmic f/fa and
linear |A|/A0 axes and the 1/sqrt(2) point marked. `make note` regenerates
`generated/amplifier_gain.pdf`; no particular amplifier performance is implied
by this normalized curve.

**Two-case study inserted after author approval of the draft:** The active paper
now presents 1000 simultaneous mid-frame arrivals and 1000 evenly spaced arrivals
over 10 ms. `scripts/build_arrival_cases.py` generates a four-panel current/voltage
figure, a charge-recovery table and LaTeX values from the stored physical
waveforms. This comparison uses point voltage samples at the frame boundaries,
with no noise; it does not reuse the 10 microsecond aperture averages from the
noise benchmarks. Point-sample results and provenance are saved separately in
`results/realistic_csa/arrival_cases_point_samples.json`. The input-charge
integrals are independently checked against the waveform CSVs. `make arrival-cases`
regenerates these outputs, and `make note` includes this dependency. No LaTeX
compilation was run for this insertion, following the author's instruction.

**Implementation:** A delegated agent has rewritten the useful circuit assumptions
from Rogelio's HTML in Python, preserving the original as provenance. Shared
specification: [realistic_csa_spec.md](readout/realistic_csa_spec.md). Parameters
are illustrative until supported by a device model. The rewrite removes the
arbitrary resistance cap, unsupported fixed stability number and inconsistent
noise formulas rather than treating the HTML output as validation.

**Available outputs:** Input currents, finite-gain and ideal voltage waveforms,
true and recovered frame charge, noise RMS and adjacent-sample covariance,
template-calibrated S/N and conditional ion-count sensitivity. All 37 tests pass, including 15 new circuit/noise tests covering physical
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
paper values and figures reproducibly. Run the relevant physical tests when
changing calculations. Compile LaTeX only when the author explicitly requests it; do not compile
automatically after manuscript edits. When compilation is requested, inspect the
build log. Do not visually inspect the PDF unless explicitly requested. Keep scratch
work under `work/`. Do not commit or push without a new user request.

Update this plan as sections and calculations are completed, preserving the
distinction between inherited ITACA results, new analytical benchmarks, circuit
results, and demonstrated event-level performance.
