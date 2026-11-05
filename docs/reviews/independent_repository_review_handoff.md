# Independent repository review handoff

Prepared 2026-11-05 for an agent working directly with the user.

## Assignment and authority

Review the organization, documentation and code of NAUSICAA/ITACA as an
independent reviewer. The user wants a detailed, evidence-based assessment
while another agent develops a fast nonlinear CSA response. Report directly
to the user. You are not a sub-agent of the implementation agent.

The first pass is read-only for source code, scientific configurations,
manuscripts and production outputs. Write your review artifacts and any small
diagnostic programs/results in the locations below. Propose corrections;
implementation of those corrections is a separate decision. Ask the user when
a physical assumption or consequential choice is genuinely unclear. Historical
handoffs contain production authorizations for earlier tasks; this assignment
is a review, not an instruction to repeat those productions.

## Baseline and concurrent work

- Repository: `/Users/jjgomezcadenas/Projects/NAUSICAA`.
- Remote: `https://github.com/jjgomezcadenas/NAUSICAA.git`.
- Reviewed baseline: commit `251118d41` on `main`, pushed to `origin/main`.
  Record the full SHA using `git rev-parse 251118d41` in your report.
- Active implementation branch: `csa_nonlinear`, created from that baseline.
- Review preparation on this branch updates `AGENTS.md`, adds this handoff,
  and marks two historical handoffs as historical. Read these preparation
  changes as current guidance; distinguish them from the committed baseline.
- At handoff preparation, nonlinear CSA implementation has not started.

Do not switch branches, reset files, stage another agent's changes or commit
in the shared checkout. Record `git status --short` at the beginning. For a
stable view, inspect baseline files using `git show 251118d41:path/to/file`.
For test runs requiring writes, use an isolated checkout/worktree or a scratch
copy under `work/repository_review/`, pinning the baseline. Record the exact
revision and local modifications tested. A separate checkout does not include
ignored numerical data, the Python environment, local tool cards or Julia
caches: reference the existing resources read-only, and explicitly configure
test output locations. Avoid copying many gigabytes merely to inspect them.

The implementation agent will modify CSA library code, scripts, tests, cards
and model documentation on `csa_nonlinear`. Your initial findings should stay
anchored to the baseline. A second review of the nonlinear change can follow
when the user identifies a completed implementation commit.

## Deliverables

Write:

- `docs/reviews/repository_review_2026-11-05.md` — complete findings.
- `docs/reviews/repository_review_2026-11-05.html` — readable, self-contained
  browser version with navigation by topic and severity, and links to evidence.
- `work/repository_review/` — diagnostic programs, command logs and temporary
  results. Keep large data and build auxiliaries out of Git.

Each finding should have a stable ID, severity, category, precise file/function
or baseline line reference, evidence, practical consequence, proposed remedy
and validation method. Distinguish confirmed defects, supported concerns,
documentation gaps and optional improvements. State uncertainty explicitly.
Include the commands you actually ran, their outcomes, and coverage limitations.
Record important parts that were inspected and found consistent as well.
Conclude with a short prioritized action list for discussion, not an automatic
rewrite plan. Group findings that have one underlying cause.

This is a code/documentation review: Markdown plus HTML is sufficient. Any new
scientific numerical study should follow the project's Julia → HDF5 → Python
presentation → generated LaTeX/PDF workflow.

## Reviewed entry material

Read these in order. This list was checked against the current checkout,
selected cards, manifest and implementation on 2026-11-05. It replaces the
earlier broad reading list of plans and handoffs.

1. `AGENTS.md` and this handoff: conventions, current scope and review contract.
2. `docs/elm_matched_production.md`: completed matched workflow, configuration,
   output selection and numerical/presentation commands. Its listed cards and
   `results/elm_matched_production/manifest.toml` identify the actual run.
3. `docs/ltspice_setup.md`: current installation, canonical paths, completed
   cleanup, validation records and Git policy. Migration counts and test
   outcomes are dated records, not instructions to redo the migration.
4. `docs/sipm_statistical_implementation.md`: current statistical interfaces,
   illuminated-cell wrapper and validation methodology. Its older pilots are
   labelled as separate configurations; use the matched cards for production.

`README.md` and `scripts/README.md` now distinguish the current matched workflow
from earlier benchmark/pilot commands. They are navigation aids. Inspect
`Makefile` and `Project.toml` as implementation artifacts when reviewing
workflows and dependencies; a dependency declaration is not evidence that its
entire workflow has been tested.

The earlier detector pipeline plan, CSA/SiPM handoffs, event-analysis roadmap,
diffusion audit, migration proposal, script inventory and decisions chronology
are removed from required onboarding. They remain repository review subjects,
with their historical content preserved. The assignment does not rely on them
to establish current production settings, completion or authorization.

For equations, follow the source-specific pointers in the sections below to
`docs/elm.tex`, `notes/physics/csa_design.tex` and
`docs/sipm_electronics_model.tex`. Review the formulas alongside their active
implementation and selected cards; this handoff's documentation reconciliation
does not claim a new independent physics validation of those notes.

Preparation checks underlying this list:

| Material | Evidence checked |
| --- | --- |
| Current production | Combined completion manifest, matched realization/detector/response cards, controller and Make targets |
| LTspice setup | Runtime resolver and local configuration, canonical filesystem locations, recorded migration and cleanup audits |
| Statistical SiPM interfaces | Constructors and functions in `sipm_statistical.jl`, `Sensor` in `sipm_response.jl`, selected illuminated-cell card |
| README hit-table contract | `event_io.jl`: v3 schema, Float32 MC time/energy, stored primary IDs, reader signatures |
| Script guide | Current matched controller and viewer targets; distinction from the earlier direct-readout/pilot commands |

This preparation corrected stale current-status claims and interface examples.
It did not rerun the physical production or replace the independent review.

## Physical and computational context

One MC ionization track produces shared primary electron/positive-ion pairs.
Primary counts fluctuate with the agreed Fano prescription and are shared
between the two arms; they are not independently resampled. Preserve the
existing treatment of hits with mean charge below one pair when tracing the
sampler. Trace its implementation and tests rather than substituting the old
Poisson(1) suggestion from discussion history.

One longitudinal placement gives complementary ion/electron drift lengths.
Diffusion differs by carrier and gas. Sample carriers directly into physical
transverse cells and longitudinal buckets. The production files contain sparse
histories and provenance, not all individual carrier coordinates. Histogram
centres are for display, not a replacement for physical sampled positions.

Detector-stage ion histories and electron photon histories are independent
products linked to a common realization. The response stage applies the CSA
to ions, and PDE/SiPM/electronics to photons. Changing one arm's geometry or
electronics should invalidate the appropriate descendants, while preserving
the unaffected arm. Two gas prescriptions exist; a particular run can select
only Xe–He. Preserve model/gas identity and photon-band metadata.

NEXUS hits retain Float32 x/y/z in mm, energy in MeV and MC time in ns. Signal
formation/electronics convert explicitly to SI where needed. Inspect actual
numerical precision and documented error checks; generic preferences for
Float64 or SI in event tables are not findings by themselves.

Track reconstruction stages are A: MC; B: diffused; C: discretized readout.
Check their actual implementation and documentation separately. A/B selection
does not gate later propagation. Do not reintroduce earlier abandoned topology
or barycentre prescriptions without tracing the final reviewed design.

The central-square-pad weighting model is an adopted approximation. Weighting
potential and applied electric field are separate boundary-value problems.
Input charge obeys integral(I dt) = q Delta(phi_w); no extra collection impulse
or silent renormalization belongs in that accounting.

Electron attachment acts only on electrons. Event-energy comparisons correct
the mean attenuation using the known drift; residual attachment fluctuations
remain. Distinguish RMS/mean from FWHM and clipping selections from intrinsic
resolution. For pulse means/maxima, keep the agreed common positive-input
channel–bucket cells across stages and retain signed reconstructed values.

Physics formulas should cite primary literature or documented local
derivations; Petit supplies implementation lineage. User-supplied values and
idealizations should be labelled as design inputs. Scientific review should
identify concrete inconsistencies and quantify relevance where possible.

## Current completed production

The latest coherent sample is 8,000 bb0nu events, Xe–15% He only:

| Quantity | Current selection |
| --- | --- |
| Source | 1 mm G4 steps, `events/0nubb_hitdf.h5` |
| Chamber | 150 cm; shared ion drift uniform in 5–145 cm |
| Electron drift | 150 cm minus ion drift |
| Ion cells | 4 mm pitch; 40 ms buckets |
| Electron cells | 4 mm pitch; 250 ns buckets; 1 mm/µs |
| EL | 15 bar; E/P = 2 kV/cm/bar; 5 mm depth; 2.5 mm holes |
| Guides | 2.5 mm diameter; 0.5 mm radial xenon gap; nominal 100 mm length |
| Optics | Analytical hole acceptance; TPB 0.65; forward factor 0.5; documented meridional capture; ideal configured propagation/coupling |
| Lifetime | 15 ms for electron attachment |
| SiPM | S14160-3015PS reference; 3 × 3 mm² active area; 15 µm cells; PDE 0.32 |
| Cell accounting | 39,984 physical; 21,808 illuminated by the guide |
| SiPM recovery | Provisional 12.5 ns statistical model |
| CSA | Existing v3 calibration; 0.63 V guard; 0.5 ENC model tolerance |

The user accepts a compact matrix implementation at 4 mm pitch. The commercial
package envelope is not a new blocker for this design study. Optical idealizations
and provisional recovery assumptions remain documented assumptions to inspect.

Start from these cards:

- `config/specs/realisations/elm_matched_full.toml`
- `config/specs/detector/configurations/elm_matched_full.toml`
- `config/specs/detector/electrons/elm_4mm_250ns_ep2_shared.toml`
- `config/specs/detector/ions/pitch4mm_bucket40ms.toml`
- `config/specs/runs/elm_matched_detector_full.toml`
- `config/specs/runs/elm_matched_csa_full.toml`
- `config/specs/runs/elm_matched_sipm_full.toml`
- `config/specs/runs/elm_matched_sipm_response_full.toml`
- `config/specs/runs/elm_matched_response_viewer.toml`

Read `results/elm_matched_production/manifest.toml` and the linked stage
manifests for actual file paths and identities. Important local artifacts:

- `results/elm_matched_production/detector_index.h5`
- `results/elm_matched_production/ion_csa/359613d487dfeee8/manifest.toml`
- `results/elm_matched_production/response/77ed411cb39be020/manifest.toml`
- `results/elm_matched_production/elm_matched_production.html`
- `latex_generated/electron_readout/elm_matched_production/elm_matched_production.pdf`

The full run is complete. No new full production is needed for the review.
Earlier productions retain their own geometry, lifetime, drift range, J values
and gas choices. Do not combine histories based only on matching event IDs.

## Directory organization to inspect

| Area | Role |
| --- | --- |
| `src/event_analysis/` | Julia event, sampling, transport and collection calculations |
| `src/ion_signal/` | Signal formation, CSA, SiPM and electronics; shared simulator/path helpers |
| `src/cloud_viewer/` | Production viewer server, HDF5 adapters and browser UI |
| `scripts/` | Task-grouped workflows; Julia reduction and Python presentation |
| `config/specs/` | Realization, detector, response and run cards |
| `config/studies/electronics/` | Relocated circuit-study cards |
| `results/` | Numerical results, manifests, presentation figures/HTML |
| `latex_generated/` | Generated reports grouped by scientific subject; catalogue/index |
| `notes/physics/` | Earlier manuscript, physics notes and shared bibliography |
| `manuscripts/itaca3/` | Self-contained EPJC manuscript and copied dependencies |
| `papers/` | Tracked reference literature, subject to size checks |
| `LTSpice/reference/` | Received circuits, models and reference material |
| `LTSpice/designs/` | Reusable circuit variants |
| `LTSpice/releases/` | Frozen portable collaborator packages and their inputs |
| `LTSpice/notes/`, `LTSpice/papers/` | Circuit notes and literature |
| `work/` | Ignored caches, scratch, archives and test outputs |

The recent reorganization removed top-level alias forests in `LTSpice/` and
`latex_generated/`. Some nested result aliases support immutable historical
paths. Inspect `src/ion_signal/project_paths.py`, `project_paths.jl`,
`docs/ltspice_paths.json`, `docs/path_migrations.toml` and the report registry.
Assess whether active callers resolve paths correctly, whether fresh checkout
behavior is documented and whether retained aliases are necessary.

About 3.6 GB of formerly tracked generated HDF5/current-waveform data were
removed from Git tracking with explicit user approval. They remain local.
Logs and build auxiliaries are ignored. Frozen release inputs stay tracked.
Do not mistake these intentional exclusions for data loss. Git history was not
rewritten, so old blobs remain in historical commits. Missing data in a fresh
clone should lead to clear prerequisites/reproduction instructions, not silent
substitution of another run.

The installed simulator is external to the repository:
`~/Applications/NAUSICAA/17.2.4/LTspice.app/Contents/MacOS/LTspice`.
Support/installers live in
`~/Library/Application Support/NAUSICAA/ltspice/17.2.4/`.
`config/local/tools.toml` is local/ignored; `config/tools.example.toml` is the
tracked template. Verify the selected executable without moving/reinstalling it.

## Code review priorities

1. Trace public workflow commands through cards, libraries and outputs. Check
   default paths, stage dependencies, resume behavior and completeness markers.
2. Audit shared primary realization and complementary placement, per-arm seed
   streams, identity/hash inputs, units and sparse-cell accounting.
3. Check attachment, Fano/EL/optical sampling and SiPM statistical response
   against the documented formulas and explicit approximations.
4. Trace CSA selection, acquisition, waveform reduction, noise covariance and
   final reconstruction. Check units, state/history dependence and provenance.
5. Review tests for independent numerical/physical checks and missing important
   failure modes. Distinguish synthetic tests from validation against saved data.
6. Inspect library/script boundaries, duplicated logic and public interfaces.
   Prioritize maintainability changes with an identifiable practical benefit.
7. Check report regeneration, generated values, stale output references and
   document status. Protect author-written manuscript content.
8. Verify viewer event ordinal vs MC ID, arm identity matching, geometry axes,
   charge sums, signed values, stages, PNG labels and configuration selection.

The latest viewer maps longitudinal buckets to source z as:

`z_center = reference_z_mm + drift_direction * (bucket + 0.5) * bucket_width_mm`.

Electron direction is -1; ion direction is +1. It reverses electron image rows
into increasing z and uses a common initial extent with 25% extra width, equal
spatial scale and independent zoom/rotation. Arrival times remain in 3D hover.
Relevant regression: `tests/test_production_viewer.py`.

## CSA performance context and planned nonlinear work

The current branch will add an optional fast nonlinear response while preserving
the existing SPICE reference and threshold. Proposed card mode names are
`hybrid_spice` and `hybrid_nonlinear`; they are proposed interfaces, not modes
already implemented in the review baseline.

Current selection uses the minimum voltage predicted by the linear model for
the whole pad history. BELOW 0.63 V triggers SPICE, together with controls.
Larger signals drive the output downward. The 0.5 ENC criterion applies to
response accuracy; it is distinct from charge quadrature and compression flags.

Latest saved run evidence:

- 697,133 pad histories; 29,920 direct transistor histories including controls.
- Initial acquisition phase approximately 1,106 s; complete CSA chain about
  1,734 s. These are different scopes, not separately measured retry costs.
- Six histories triggered a hardcoded 1e-7 relative input-charge quadrature
  check. Their maximum discrepancy was about 0.033%. The retry changed the
  maximum timestep from 250 µs to 2.5 µs and restarted Julia reduction.
- Previously reduced histories were cached. It did not rerun entire events.
- Twelve additional histories were independently refined for solver validation.
- Existing final corrected fast-path validation was about 0.434 ENC; the raw
  control summary gives about 0.438 ENC. Check which stage a number describes.

Read the saved evidence rather than repeating these numbers uncritically:
`results/elm_matched_production/ion_csa/359613d487dfeee8/execution/pitch4/filtered/5bcfd782ef02bcad/`.
This includes `plan.h5`, `calibration.h5`, `jobs.toml`,
`acquisition_summary.json`, `fast_path_validation.toml`, and
`spice_cases/<key>/{input.toml,current.pwl,raw.h5,response.h5}`.

The implementation agent's scope is a compact stateful nonlinear model fitted
to saved histories, with event-level calibration/validation separation, the
existing finite-gap current prescription, voltage averaging and reconstruction,
unchanged noise treatment, explicit calibration identity and SPICE fallback
outside its validated domain. The intended first application is matched 4 mm
bb0nu, followed by checks at 20/30/40 ms. A memoryless per-bucket amplitude
correction is not the agreed model. Assess the eventual implementation against
held-out histories and final corrected outputs, not fit residuals alone.

No need to decide a new error tolerance or commission another full run during
this repository review. Identify any implementation issue that could affect
calibration and report it to the user promptly.

## Validation procedure and environment

First inspect commands and targets. Some report targets regenerate tracked
artifacts, and production targets can run for many minutes or write gigabytes.
Run broad writing/build checks in the isolated review workspace. Keep command
logs and record environment/version information. Do not stop the running viewers.

Available local Python environment: `work/event-viewer-venv/bin/python`.
Julia is normally `julia` or `~/.juliaup/bin/julia`, with project dependencies
and `JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia"` from the main checkout.
For an isolated checkout, use the absolute existing depot/environment paths.
The user allows no PDF visual inspection unless explicitly requested; compile
and inspect logs for LaTeX checks.

Useful commands from a suitably configured repository root:

```sh
git status --short
git rev-parse HEAD
make test physics PYTHON=work/event-viewer-venv/bin/python
make test-script-layout VIEWER_PYTHON=work/event-viewer-venv/bin/python
make ion-csa-test
make sipm-test
make elm-matched-test
work/event-viewer-venv/bin/python -m unittest discover -s tests -p test_production_viewer.py
PYTHONPATH=src work/event-viewer-venv/bin/python -m ion_signal.spice_runtime
```

The baseline reorganization recorded passing Python physics/path/simulator
checks and Julia response/correction checks. The viewer coordinate regression
passed five Python tests and was checked visually on event 7. These are prior
results, not a claim that you reran them or that they cover everything.

`make note` builds the earlier `notes/physics/nausica.tex`, not the new EPJC
manuscript. A pre-existing unmatched-math failure around line 66 was reported
in that older draft during migration; verify current status and leave author
text untouched in this review. `manuscripts/itaca3/` has its own Makefile and
was compiled successfully using the supplied EPJC template. Do not change its
format based on the older manuscript's appearance.

The live production viewer is at `http://127.0.0.1:8769/`; an older viewer may
also be running on 8768 and a smoke viewer on 8774. Inspect the production card
before comparing displays. Avoid launching a second server on an occupied port.

## Completion criteria

The review is complete when the report states its tested baseline and coverage,
provides actionable evidence for each finding, separates active/historical work,
and identifies a manageable priority order for the user. Do not claim a full
scientific audit of files or workflows you did not inspect. No source rewrite,
data deletion, history rewrite, publication or production rerun is part of this
handoff. The user will decide which follow-up changes to authorize.
