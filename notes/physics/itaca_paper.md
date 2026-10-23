# NAUSICA for ITACA: paper plan and work status

Working plan agreed with J. J. Gómez-Cadenas. Updated 23 October 2026.

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

## Approved green text — 23 October 2026

The author identifies the following prose as green and approved. This records
his assessment, not an independent Pangram test.

| Manuscript passage | Status |
| --- | --- |
| Abstract | Green; preserve the current wording |
| Section 1, Introduction | Green; preserve the current wording |
| Section 2, The ITACA detector, including Principle of operation | Green; reviewed and approved |
| Section 3 opening paragraphs, before Charge collection and focusing | Green; author-revised and approved |
| Section 3, Charge collection and focusing | Green; approved through the end of this subsection |
| Pixel electronics and readout, and subsequent text | Not covered by this approval; continue review |

The latest approved segment starts at `\subsection{Principle of operation}`
and ends immediately before `\subsection{Pixel electronics and readout}`.
Checkpoint commit `3cf0db3` preserves the manuscript and comments before cleanup.
Comment-only lines in this interval were subsequently removed; the active text
was preserved exactly. Comments outside the interval remain for later review.

Handle green prose with care. Do not insert longer sentences, reformulate it,
or change its meaning without the author's approval. Pure grammar and small
word-level corrections may be treated as minor. A green designation does not
mean pending calculations, including the Appendix B collection study, have
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
the original section remains archived as comments after the active draft.
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

**Remaining:** Turn this material into a compact paper section, add the CSA
equations and relevant output plots to the manuscript, and select a small set
of instructive examples. Keep example capacitances and discharge times clearly
separate from the eventual design choices.

### 5. The ion signal from a track

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

### 6. The NAUSICA design

**Purpose:** Select and assemble the conceptual design using Sections 4 and 5.

Relate pitch and sampling to image quality; capacitance to gain and range;
feedback discharge to charge retention and temporal reconstruction; and noise
and digitisation to recoverable track signals. Treat these choices together.
A sampled CSA voltage includes the history of arriving charge and is not simply
the charge delivered during the latest frame.

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

### 7. Track reconstruction and performance

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

### 8. Conclusions

**Purpose:** State the resulting conceptual design and demonstrated imaging
performance, followed by the specific work needed for implementation.

**Remaining:** Write after the calculations and design choices are established.
Separate achieved results from further electrostatic, circuit and experimental
validation. Avoid repeating the development plan in the published conclusions.

## Where Rogelio's work enters

| Stage | Role | Present status |
| --- | --- | --- |
| Section 4 | Analytical current-to-voltage response, informed by the supplied CSA calculator | Adapted model implemented; manuscript integration pending |
| Section 6 | Transistor-level design and SPICE response to realistic pad currents | Circuit calculations/results still need to be supplied or connected and assessed |
| Section 7 | Achievable response and noise parameters for event reconstruction | Pending circuit comparison and noise model |

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
design direction. Appendix B is reserved for the quantitative focusing
and collection study that will test the near-unit-efficiency assumption.
That calculation remains pending; it is distinct from the flux estimate.
Keep weighting-potential derivations in Appendix A and place any separate
numerical-validation appendix after Appendix B.

| Current material | Destination |
| --- | --- |
| Scope and definitions | Essential definitions at first use in Sections 3–5 |
| Detector definition and voltages | Concise geometry in Section 3; detailed inputs where needed |
| Focusing entrance argument | Brief estimate in Section 3; quantitative collection study in Appendix B |
| Transport and diffusion | Section 5, referring back to Section 2 instead of rederiving it twice |
| Weighting potential and single-ion signal | Essential results in Section 4 |
| Central-pad derivation and shape comparison | Appendix A, with selected results retained in Section 4 |
| Pulse-duration definitions and numerical checks | Duration derivations in Appendix A; separate validation appendix if needed |
| Implications for the next physics stage | Replace with the actual event calculation in Section 5 |
| Validation, reproduction and remaining work | Scientific validation after Appendix B if needed; commands and implementation instructions in repository documentation |

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
  from Section 5.
- The title currently uses NAUSICAA while the expanded acronym in the body uses
  NAUSICA. Settle the final spelling with the author rather than silently
  changing the reviewed text.

## Working method and next step

Work through the sections in order. The next task is author review of
Pixel electronics and readout. The preceding green text is approved and should
remain stable. Do not reorganise the entire manuscript in one pass.

Preserve the author's voice. The author reports the abstract and reviewed
Sections 1–2 as acceptable to Pangram; this plan does not independently assess
that score. Pure grammar and small word-level corrections can be treated as
minor. For substantive edits, show the original and replacement one at a time,
wait for approval, apply the agreed change, and then proceed. Comment out
proposed removals for review before deleting material when appropriate.

Keep physical calculations in `src/ion_signal`, use SI internally, and generate
paper values and figures reproducibly. Run the relevant physical tests when
changing calculations. Compile LaTeX and inspect the build log after manuscript
edits; do not visually inspect the PDF unless explicitly requested. Keep scratch
work under `work/`. Do not commit or push without a new user request.

Update this plan as sections and calculations are completed, preserving the
distinction between inherited ITACA results, new analytical benchmarks, circuit
results, and demonstrated event-level performance.
