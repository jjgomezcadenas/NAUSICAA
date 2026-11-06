# NAUSICAA/ITACA independent repository review, 2026-11-05

Independent review of the organization, documentation and code of the NAUSICAA/ITACA
repository, performed on 2026-11-05 for the user under
`docs/reviews/independent_repository_review_handoff.md`. The review is anchored to the baseline
commit `251118d414e1ea28ac39f6b83e38e0b85c6b36d1` (`251118d41`, "Complete matched ELM production
and reorganize simulation artifacts", on `main`). The checkout under review was branch
`csa_nonlinear` at `d0b6476457f0e013f470617f9dec8a84341b3719`, which differs from the baseline
only in guidance documents (`AGENTS.md`, `README.md`, `scripts/README.md`, `docs/*.md` and the
handoff). No tracked implementation file differs between HEAD and the baseline, so every source
reference below applies to both. During the review the concurrent CSA agent added untracked files
(`src/ion_signal/ion_csa_nonlinear.jl`, `scripts/production/calibrate_ion_csa_nonlinear.jl`,
`tests/test_ion_csa_nonlinear.jl`, `config/specs/response/ion_csa/nonlinear_calibration_v1.toml`,
`docs/ion_csa_nonlinear_model.html`, `results/ion_csa_nonlinear/`) and, towards the end of the review, modified
tracked CSA files (`Makefile`, `src/ion_signal/ion_csa_{pipeline,response,summary}.jl`,
`scripts/production/{ion_csa_pilot.jl,reduce_matched_analysis.jl,run_ion_csa_pilot.py}`,
`tests/test_ion_csa_pipeline.jl`). Those changes were not read or reviewed; every line reference in this report is
to the baseline (`git show 251118d41:<path>`), which equals the working tree as it stood when each file was read.

The review is read-only for source, cards, manuscripts and production outputs. Tests that write
were executed in an isolated `git worktree` pinned to the baseline under
`work/repository_review/baseline_worktree/`, with the existing Julia depot, Python environment and
NEXUS hit tables referenced read-only. Scratch programs, command logs and the five sub-review notes
are under `work/repository_review/`. Nothing was committed, staged, deleted or published. Five
parallel read-only sub-reviews (realisation/sampling, CSA chain, SiPM chain, viewer, organization)
fed this report; their claims about code were spot-checked against source, and the two with the
largest consequences (S-01 source extents, W-01 deck hashing) were independently confirmed.

**Severity scale.** *critical*: wrong physics or wrong numbers in a current product. *major*:
incorrect or misleading behaviour, or a gap that can plausibly produce wrong results or block
reproduction. *minor*: robustness, clarity or maintainability with a concrete practical benefit.
*info*: observation or optional improvement. **Categories**: confirmed defect, supported concern,
documentation gap, optional improvement. Each finding states its uncertainty. No critical finding
was identified; no finding changes the numbers of the completed 8,000-event matched production.

## Contents

1. Baseline, environment and method
2. Summary of the assessment
3. Findings
4. Inspected and found consistent
5. Commands run and outcomes
6. Coverage limitations
7. Prioritized action list

## Baseline, environment and method

| Item | Value |
|---|---|
| Baseline reviewed | `251118d414e1ea28ac39f6b83e38e0b85c6b36d1` (`main`, pushed) |
| Checkout HEAD | `d0b6476457f0e013f470617f9dec8a84341b3719` on `csa_nonlinear` |
| `git status --short` at start | clean (0 lines), 2026-11-05T09:13Z |
| HEAD vs baseline | 8 files, documentation only (`git diff --stat 251118d41 HEAD`) |
| Isolated test checkout | `git worktree add --detach work/repository_review/baseline_worktree 251118d41` |
| Python | `work/event-viewer-venv/bin/python` 3.13.13; h5py 3.16.0, numpy 2.4.4, matplotlib 3.10.9 |
| Julia | 1.12.6 (`~/.juliaup/bin/julia`); `Manifest.toml` julia_version 1.12.6; depot `work/julia-depot:~/.julia` |
| TeX | TeX Live 2026, latexmk 4.88 |
| Simulator | native LTspice 17.2.4 at `~/Applications/NAUSICAA/17.2.4/…`, SHA `4036f0d4…` (from `config/local/tools.toml`, not launched) |
| Running viewers | ports 8768 (legacy grid), 8769 (production card), 8774 (smoke); none started or stopped |

Method: read the handoff entry material in order; trace the matched workflow controller
(`run_elm_matched_production.py` → `run_matched_detector.jl` → concurrent `run_matched_sipm.jl` and
`run_ion_csa.jl` → `assemble_matched_response.jl` → `reduce_matched_analysis.jl` →
`present_matched_analysis.py`) through cards, libraries and saved outputs; run the test suites and report
builds in the isolated worktree; delegate five parallel read-only sub-reviews with explicit questions and
cross-check their claims against source and saved data.

## Summary of the assessment

The completed 8,000-event matched production is coherent and traceable at the level the handoff asks
about: one shared primary-pair ledger feeds both arms, the complementary placement is exact to
7.6e-6 cm, identities are verified across six histories per event, charge is conserved exactly
(897,306,359 ions sampled equals pairs), the EL yield, optical acceptance, attachment survival and
SiPM moments agree with the documented formulas and with independent recomputation, the CSA input
current preserves `∫I dt = qΔφ_w`, and the signed three-frame reconstruction and noise covariance
match their derivation. Python (70) and Julia (≈2,300 checks) tests pass once two fresh-checkout
prerequisites are supplied, the EPJC manuscript and the matched report compile.

The most consequential items for discussion are:

1. **S-01** – NEXUS deposits up to 2.5 m apart are propagated as if inside the 1.5 m chamber, and
   the ASIC acquisition window (global photon span ±4 buckets over all illuminated channels) makes the
   electron energy width after electronics a function of bookkeeping rather than physics. This needs a
   user decision on fiducialisation and on the energy estimator definition.
2. **W-01 / W-02** – the CSA identity checks do not cover the deck, model and kernel waveform that are
   actually simulated, nor the simulator build. Harmless today, but directly relevant while the
   nonlinear model is fitted to saved histories: report to the implementation agent now.
3. **V-01 / V-02** – the viewer shows photon-cell-restricted electron totals without saying so, and
   renames the gas identity `xe_he15` → `xe_he` without cross-checking HDF5 metadata.
4. **R-03** – the superseded Oliveira Figure 6 J anchors are still enforced by the card validator,
   the interactive form and two tests; the production avoided them via `user_specified`.
5. **T-01** – fresh-checkout test failures (`work/`, `events/`) and data-dependent tests without
   prerequisite messages.
6. **W-09 / O-04** – absolute machine paths in every production manifest, and documented commands
   that die inside HDF5 in a fresh clone instead of naming the missing prerequisite.
7. **O-02 / O-03** – Git tracking policy: 516 MB of tracked HDF5 (with inconsistent sibling ignore
   rules) and about 25,000 generated SPICE decks tracked while their inputs and outputs are ignored.

Documentation is extensive but layered: the three entry documents are current, while several older
documents (`docs/csa_4mm_production.md`, `docs/sipm_electronics_model.tex` implementation paragraphs,
`scripts/production/README.md`, `config/electron_transport.toml`) describe superseded states without
being marked historical.

## Findings

Finding counts by severity and category:

| Severity | confirmed defect | supported concern | documentation gap | optional improvement | Total |
|---|---|---|---|---|---|
| critical | 0 | 0 | 0 | 0 | 0 |
| major | 3 | 5 | 3 | 0 | 11 |
| minor | 1 | 5 | 10 | 9 | 25 |
| info | 0 | 0 | 5 | 8 | 13 |
| **Total** | 4 | 10 | 18 | 17 | 49 |

Index:

| ID | Sev. | Category | Topic | Title |
|---|---|---|---|---|
| W-01 | major | supported concern | Ion CSA chain | Simulated CSA deck, model and kernel waveform are not covered by the calibration identity checks |
| W-02 | major | documentation gap | Ion CSA chain | Simulator identity and acquisition summaries are not recorded in tracked production outputs |
| W-03 | minor | confirmed defect | Ion CSA chain | Hard-coded CSA calibration prehistory duplicates the card value |
| W-04 | minor | supported concern | Ion CSA chain | CSA study resume is fragile: late `jobs.toml`, newest-by-mtime directory selection and `latest.toml` pointer |
| W-05 | minor | optional improvement | Ion CSA chain | Charge-quadrature retry restarts the whole Julia reduction once per failing history |
| W-06 | minor | optional improvement | Ion CSA chain | Numerical knobs of the CSA pipeline are hard-coded rather than card-controlled |
| W-07 | minor | documentation gap | Tests and validation | Python CSA checkers are pinned to pilot sizes; no independent validation record for the matched run |
| W-08 | minor | documentation gap | Documentation and reports | CSA production documentation describes the earlier 4 mm run, not the matched run |
| W-09 | major | confirmed defect | Paths, Git tracking and organization | Absolute machine paths are embedded in tracked manifests, plans and job files and consumed verbatim |
| W-10 | minor | supported concern | Workflow and provenance | Reduction and presentation read cards from the working tree instead of the saved run provenance |
| W-11 | info | optional improvement | Ion CSA chain | Code-hash coverage lists omit two included sources |
| W-12 | info | optional improvement | Ion CSA chain | The fast CSA model is level-independent linear superposition; direct SPICE carries nonlinear state across buckets |
| R-01 | minor | supported concern | Realisation, transport and sampling | Common-realisation code hash is whole-file and asymmetric between arms |
| R-02 | minor | documentation gap | Realisation, transport and sampling | Fano sampler is not Bernoulli for all sub-unity means, contrary to its docstring and the formula note |
| R-03 | major | documentation gap | Realisation, transport and sampling | Superseded Oliveira Figure 6 J anchors remain in the card validator, the interactive form and tests |
| R-04 | minor | optional improvement | Workflow and provenance | Matched validation and manifest hard-code 5/145/150 cm instead of reading the realisation card |
| R-05 | info | documentation gap | Documentation and reports | Grid-reference conventions are only partly documented |
| R-06 | info | documentation gap | Documentation and reports | AGENTS.md ELM design values differ from the production card |
| R-07 | info | optional improvement | Realisation, transport and sampling | Stream identity embeds the reference card's placement values |
| S-01 | major | supported concern | SiPM and electronics chain | Source hits outside the 150 cm chamber enter the electron chain, and the ASIC acquisition window scales with the event's z extent |
| S-02 | minor | documentation gap | SiPM and electronics chain | Energy-width estimators and dark subtraction are undocumented; RMS/mean is not FWHM |
| S-03 | minor | documentation gap | Documentation and reports | The SiPM derivation note describes the microscopic model and old cards, not the statistical production |
| S-04 | minor | optional improvement | Tests and validation | Statistical SiPM tests lack dense-path, illuminated-wrapper and regime-limit checks |
| S-05 | minor | optional improvement | Workflow and provenance | SiPM production manifest lacks top-level model, band and source identities |
| S-06 | info | documentation gap | SiPM and electronics chain | Attachment survival uses one reference drift per event for all hits |
| V-01 | major | documentation gap | Production viewer | Electron-panel totals are restricted to photon-occupied cells without saying so |
| V-02 | major | supported concern | Production viewer | Gas identity is renamed `xe_he15` → `xe_he` by the assembler and viewer and never cross-checked against HDF5 metadata |
| V-03 | minor | supported concern | Production viewer | Viewer start-up fails with a raw traceback when referenced data are absent |
| V-04 | minor | optional improvement | Production viewer | `CloudSource` is dead code with a different indexing convention |
| V-05 | minor | optional improvement | Tests and validation | Viewer regression covers synthetic fixtures only and is outside `make test` |
| V-06 | minor | documentation gap | Production viewer | PNG export names use the file position, not the MC ID, and re-render at a different aspect |
| V-07 | info | optional improvement | Paths, Git tracking and organization | Duplicated helpers across the three viewers and hard-coded labels |
| V-08 | info | optional improvement | Production viewer | Minor robustness notes in the viewer server |
| T-01 | major | confirmed defect | Tests and validation | Julia test suites fail in a fresh checkout because `work/` and `events/` are assumed to exist |
| T-02 | minor | optional improvement | Tests and validation | `make test` runs only Python; there is no aggregate Julia test target, and one Julia test writes into `results/` |
| D-01 | minor | documentation gap | Documentation and reports | `scripts/production/README.md` is stale and contradicts the Makefile |
| D-02 | info | documentation gap | Documentation and reports | `make note` fails on an unclosed math expression in the older physics note (author text) |
| D-03 | info | optional improvement | Documentation and reports | Status prints and card flags use misleading names |
| D-04 | info | optional improvement | Paths, Git tracking and organization | Presentation script imports helpers from another script |
| O-02 | major | supported concern | Paths, Git tracking and organization | Large numerical binaries are tracked, with an inconsistent sibling policy |
| O-03 | major | supported concern | Paths, Git tracking and organization | About 25,000 generated SPICE deck files are tracked while their inputs and outputs are ignored |
| O-04 | major | confirmed defect | Paths, Git tracking and organization | Documented commands fail inside HDF5 in a fresh clone instead of naming the missing prerequisite |
| O-05 | minor | documentation gap | Paths, Git tracking and organization | Nested result and config aliases are load-bearing but documented as deck-only; the Julia path map is derived |
| O-06 | minor | documentation gap | Documentation and reports | Stale counts and statements in navigation documents |
| O-07 | minor | documentation gap | Documentation and reports | Historical documents without a status marker |
| O-08 | minor | supported concern | Documentation and reports | Manuscript generated values are a one-time copy with no regeneration or check rule |
| O-09 | minor | optional improvement | Paths, Git tracking and organization | Library code lives under `scripts/production/` with duplicated helpers |
| O-10 | info | optional improvement | Paths, Git tracking and organization | Makefile resolves the simulator on every expansion and errors loudly when unconfigured |
| O-11 | info | documentation gap | Paths, Git tracking and organization | Environment reproducibility note |

### Topic: Workflow and provenance

#### W-10 · Reduction and presentation read cards from the working tree instead of the saved run provenance

| Field | Value |
|---|---|
| Severity | minor |
| Category | supported concern |
| Topic | Workflow and provenance |
| Location | `scripts/production/reduce_matched_analysis.jl:57-58,33`; `scripts/production/present_matched_analysis.py:95` |

**Evidence.** `reduce_matched_analysis.jl` embeds `read(joinpath(ROOT,"config/specs/response/ion_csa/csa_0p5_enc_4mm.toml"))`
and the reconstruction card from disk, although the CSA manifest already stores `run_card_toml` and the
response file stores the resolved cards. `present_matched_analysis.py` appends
`config/specs/realisations/elm_matched_full.toml` read from disk into the "Resolved cards" block. The
reduction also requires `run_timings.toml`, which only the Python controller writes, and hard-codes the
`tau_12p5/xe_he15` group name.

**Consequence.** If a card is edited after the run, the report's “Resolved cards” no longer describe the run that produced the numbers, while claiming to. Running the Julia stages by hand (without the controller) makes `make elm-matched-analysis` fail on the missing timings file.

**Proposed remedy.** Take card texts from the saved manifests/HDF5 metadata; treat `run_timings.toml` as optional; derive the gas/tau group from the manifest.

**Validation.** Edit a card copy in a scratch checkout and regenerate: the report block must still show the original text.

**Uncertainty.** None.


#### R-04 · Matched validation and manifest hard-code 5/145/150 cm instead of reading the realisation card

| Field | Value |
|---|---|
| Severity | minor |
| Category | optional improvement |
| Topic | Workflow and provenance |
| Location | `scripts/production/run_matched_detector.jl:45`; `scripts/production/assemble_matched_response.jl:42,69` |

**Evidence.** `@assert 5<=li<=145 && 5<=le<=145 && isapprox(li+le,150;atol=1e-4)` and `"chamber_length_cm"=>150.,"minimum_drift_cm"=>5.,"maximum_drift_cm"=>145.` are literals; the realisation card holds the same values.

**Consequence.** Manifest provenance is transcribed rather than derived; a card with other bounds fails or is mislabelled.

**Proposed remedy.** Read `min_ion_drift_cm`, `max_ion_drift_cm`, `chamber_length_cm` from `spec.realisation.data`.

**Validation.** Smoke run with a shifted card.

**Uncertainty.** None.


#### S-05 · SiPM production manifest lacks top-level model, band and source identities

| Field | Value |
|---|---|
| Severity | minor |
| Category | optional improvement |
| Topic | Workflow and provenance |
| Location | `scripts/production/sipm_production.jl:74-87,302-308`; `scripts/production/sipm_pipeline.jl:135-138,333-337` |

**Evidence.** `nausicaa.sipm_production_manifest.v1` records run id, card, summary and per-shard sample/gas/paths/SHA-256 but not sensor/electronics model ids, `photon_band` or the source photon SHA at top level (the pilot manifest does record model ids). The band is checked at run time but only stored inside the embedded card text.

**Consequence.** Model/band identity must be recovered from embedded card text or shard metadata.

**Proposed remedy.** Add `sensor_stage_id`, `asic_stage_id`, `photon_band`, `gas` and `source_sha256` to the manifest root.

**Validation.** Smoke run; inspect the manifest.

**Uncertainty.** None.


### Topic: Realisation, transport and sampling

#### R-01 · Common-realisation code hash is whole-file and asymmetric between arms

| Field | Value |
|---|---|
| Severity | minor |
| Category | supported concern |
| Topic | Realisation, transport and sampling |
| Location | `scripts/production/run_common_realisation.jl:18-25`; `scripts/production/run_ion_readout.jl:31-41` |

**Evidence.** `run_common` hashes all of `src/event_analysis/production_readout.jl` and asserts `metadata/code_hash`
on reuse, so an edit to an unrelated function (for example `write_photon_history`) makes reuse of a valid
`realisation.h5` fail with a bare `@assert`, while the path (card-only identity) is unchanged. The ion arm's
`check_realisation` verifies schema, source, card, stream and seed but not the code hash.

**Consequence.** Operators must delete a valid ledger by hand after unrelated edits; conversely the ion arm accepts a ledger produced by any code version.

**Proposed remedy.** Hash only the ledger-determining functions (`sample_fano_pairs`, `sample_integer_moments`, `sample_drift_cm`, `event_seed`) or include the code hash in the directory name; make both arms apply the same check.

**Validation.** Touch a comment in `production_readout.jl` and rerun `run_common` against the existing ledger (expected assertion failure today).

**Uncertainty.** Inferred from code; not exercised to avoid writing into results/.


#### R-02 · Fano sampler is not Bernoulli for all sub-unity means, contrary to its docstring and the formula note

| Field | Value |
|---|---|
| Severity | minor |
| Category | documentation gap |
| Topic | Realisation, transport and sampling |
| Location | `src/event_analysis/production_readout.jl:145-186`; `docs/production_formulas.tex:59` |

**Evidence.** The docstring says "for means below one, the integer draw is Bernoulli(mean)" and the note states
`N_i ~ Bernoulli(μ_i)` for `0≤μ_i<1`. The implementation always calls `sample_integer_moments`, which uses
the adjacent-integer law only when `F·μ ≤ f(1−f)` (line 153), i.e. for `μ ≤ 1−F = 0.85`; for `0.85<μ<1`
it mixes with `Binomial(2, μ/2)` (lines 154-159), so such a hit can yield two pairs (verified in source).

**Consequence.** Mean and variance targets are met (variance F·μ rather than μ(1−μ)); numerically negligible for 1 mm Geant4 steps (hits of 18.6–21.9 eV). Contract mismatch only.

**Proposed remedy.** Either special-case `mean<1 → Bernoulli` as documented, or correct the docstring and `production_formulas.tex`.

**Validation.** Add a test at μ=0.9 asserting the intended support.

**Uncertainty.** None on the code; the policy choice is the user's.


#### R-03 · Superseded Oliveira Figure 6 J anchors remain in the card validator, the interactive form and tests

| Field | Value |
|---|---|
| Severity | major |
| Category | documentation gap |
| Topic | Realisation, transport and sampling |
| Location | `scripts/production/specs.jl:93-100`; `scripts/production/run_pilot.py:24-27` (`FIGURE6_J_READINGS`), `:121-133`, `:156-161`; `tests/test_production_form.py:26-27`; `tests/event_analysis/test_production_specs.jl:11`; `docs/detector_response_pipeline_plan.md:157` |

**Evidence.** `j_anchors=Dict(1.5=>0.1,2.0=>0.05,3.0=>0.01)` is enforced whenever a card declares
`j_status="oliveira_fig6_user_reading"`; the interactive pilot form defaults to J=0.01 at E/P=3 with that
status and prints "E/P 1.5 → 0.1; 2 → 0.05; 3 → 0.01" (seen in the unit-test output). `AGENTS.md` and
`docs/production_formulas.tex:150-153` record the corrected 2026-11-04 readings 0.030/0.015/0.025 at
2.0/2.5/3.0 (minimum ≈0.01 near 2.8). The matched card bypasses the check with `j_status="user_specified"`
and `j_variance_ratio=0.030`, so the production used the intended value.

**Consequence.** A future card declaring the Figure 6 status with the corrected 0.030 at E/P=2 is rejected, while the superseded 0.05 is accepted and labelled as a Figure 6 reading; the tests pin the superseded values.

**Proposed remedy.** Add a dated anchor set (keep the historical values for historical cards) in one shared place used by both `specs.jl` and `run_pilot.py`; update the form text, the plan document and the two tests.

**Validation.** Card with `oliveira_fig6_user_reading`, E/P=2, J=0.030 must load; J=0.05 must be rejected or labelled historical.

**Uncertainty.** None.


#### R-07 · Stream identity embeds the reference card's placement values

| Field | Value |
|---|---|
| Severity | info |
| Category | optional improvement |
| Topic | Realisation, transport and sampling |
| Location | `scripts/production/specs.jl:49-50` |

**Evidence.** `stream_values` is the reference card minus `samples`/`max_events_per_sample`, so it still contains `min_ion_drift_cm=10` and `chamber_length_cm` of `config/specs/generated/df1414afcdc2/realisation.toml`. Anchoring works (verified: old and new ledgers share stream identity `7b0a08c6…`).

**Consequence.** The primary stream identity is not conceptually placement-free; a new card without `primary_reference` and a different drift range gets a new pair stream. Documented as intentional compatibility (`specs.jl:40`).

**Proposed remedy.** None required; mention in the formula note.

**Validation.** n/a

**Uncertainty.** None.


### Topic: Ion CSA chain

#### W-01 · Simulated CSA deck, model and kernel waveform are not covered by the calibration identity checks

| Field | Value |
|---|---|
| Severity | major |
| Category | supported concern |
| Topic | Ion CSA chain |
| Location | `scripts/production/run_ion_csa_pilot.py:43-46`; `src/ion_signal/ion_csa_pipeline.jl:23-51` (`check_calibration_resources`); `config/csa_v3_response_calibration.toml:4-5`; `src/ion_signal/pad_history_readout.jl:11-14`; `src/ion_signal/ion_csa_response.jl:108,384-388` |
| Related | W-02, W-11 |

**Evidence.** The direct-SPICE deck body and transistor model used in production are read from
`results/ltspice/csa_v3_response_calibration/decks/length_1200um/30ms_single_N100_r2.cir` and
`.../corrected.inc` (`run_ion_csa_pilot.py:43-46`). The Julia identity check hashes a different pair
of files, the calibration card's `source_deck` and `model`, which point to
`results/ltspice/feedback_optimization/decks/length_400um_gate0_rounded/{40ms_single_N400_r2.cir,corrected.inc}`
(`ion_csa_pipeline.jl:27-28,41-45`; `csa_v3_response_calibration.toml:4-5`). The linear-kernel waveform
file `variants/length_1200um.h5` read through `cfg["circuit_results"]` is not hashed either. The
`deck_sha256`/`model_sha256` attributes exist only in `raw.h5`/`response.h5` under the gitignored
`spice_cases/`, and `finish` copies only three datasets, so the final `ion_response.h5` carries no
deck identity. Today both model copies are content-identical (verified with `shasum`), and
`git status` shows the tracked deck files unchanged.

**Consequence.** Editing the length_1200um deck or model, which is plausible while the nonlinear
model is being fitted to saved histories, would make direct SPICE and the linear kernel describe
different circuits while every identity check stays green. The published response could not be
traced to the simulated circuit body. This is the one implementation issue found that could affect
the nonlinear calibration work in progress, and should be communicated to the implementation agent.

**Proposed remedy.** Hash `decks/<variant>/*.cir`, `decks/<variant>/corrected.inc` and `variants/<variant>.h5`
in `check_calibration_resources` (or fold them into `inputs_sha256`); compare the `raw.h5`
`model_sha256` with the calibration `model_sha256` in `reduce_cases`; write deck/model/simulator SHAs
into the `finish` metadata. Keep the existing checks.

**Validation.** Extend the `tests/test_ion_csa_pipeline.jl` fixture (lines 116-122) with a copy of
the decks directory and assert an `ArgumentError` after mutating the deck; `h5dump -a` of the response
metadata after a smoke run should show the new SHAs.

**Uncertainty.** Low on the code path (confirmed in two places). No present numerical impact.


#### W-02 · Simulator identity and acquisition summaries are not recorded in tracked production outputs

| Field | Value |
|---|---|
| Severity | major |
| Category | documentation gap |
| Topic | Ion CSA chain |
| Location | `src/ion_signal/spice_runtime.py:60-61`; `scripts/production/run_ion_csa_pilot.py:28-33`; `src/ion_signal/ion_csa_pipeline.jl:206-239,352` |
| Related | W-01, W-09 |

**Evidence.** The LTspice executable is pinned by SHA-256 only in the Python pilot (`spice_runtime.py:60-61`
against `config/local/tools.toml`, and `run_ion_csa_pilot.py:28-33` against `analysis.h5`
`simulator_sha256`, which match: `4036f0d4…`). `IonCSAProduction.settings` and the final response
metadata (54 keys) have no `simulator_*` entry; `acquisition_summary.json` and
`fast_path_validation.toml` hold the simulator identity and the 0.4375/0.4336 ENC validation numbers
but sit under `execution/`, which `run` gitignores (`ion_csa_pipeline.jl:352`).

**Consequence.** A reader of the tracked manifest or the final HDF5 cannot tell which simulator
build produced the direct histories, and the validation evidence quoted in the handoff lives only on
the production machine.

**Proposed remedy.** Store simulator version/SHA in `plan.h5` attributes and in the `finish` metadata; copy
`acquisition_summary.json` and `fast_path_validation.toml` into the tracked run directory next to
`manifest.toml`.

**Validation.** Smoke run; inspect metadata keys and the tracked directory listing.

**Uncertainty.** None.


#### W-03 · Hard-coded CSA calibration prehistory duplicates the card value

| Field | Value |
|---|---|
| Severity | minor |
| Category | confirmed defect |
| Topic | Ion CSA chain |
| Location | `src/ion_signal/ion_csa_response.jl:108,117`; `config/csa_v3_response_calibration.toml:12`; `src/ion_signal/pad_history_readout.jl:14,19` |

**Evidence.** `calibration()` passes `"prehistory_s"=>.01` to `P.kernel` and `B.pwl(.03,…,.01,.15,2.5e-6)`,
duplicating `prehistory_s = 0.01` in the calibration card. The kernel alignment `t.-=c["prehistory_s"]`
would use 0.01 even if the card changed; the only guard is the loose `@assert abs(interp(t,v,3frame)/(raw[3]*alpha)-1)<1e-4`.

**Consequence.** Latent coupling: a changed calibration prehistory forces regeneration (card SHA) but the kernel taps could be misaligned without a check tripping.

**Proposed remedy.** Read `prehistory_s` from the calibration card in `calibration()`.

**Validation.** Kernel unit test with a prehistory different from 0.01.

**Uncertainty.** Low; the defect is latent (no card change has occurred).


#### W-04 · CSA study resume is fragile: late `jobs.toml`, newest-by-mtime directory selection and `latest.toml` pointer

| Field | Value |
|---|---|
| Severity | minor |
| Category | supported concern |
| Topic | Ion CSA chain |
| Location | `src/ion_signal/ion_csa_response.jl:166,217,224,278,305`; `scripts/production/run_ion_csa_pilot.py:41,100`; `scripts/production/run_ion_csa_correction.py:10`; `scripts/production/assemble_matched_response.jl:21` |

**Evidence.** `prepare` returns early when `plan.h5` exists (line 166) but writes `jobs.toml` after moving the
plan (217, 224); an interruption between the two leaves a study the Python pilot cannot resume
(`FileNotFoundError` at `run_ion_csa_pilot.py:100`). The pilot and the correction step pick the newest
`plan.h5`/`summary.h5` by mtime rather than by the identity-derived study id. `reduce_cases` skips an
existing `response.h5` without comparing its stored `raw_sha256` with the current `raw.h5` (278, 305).
`assemble_matched_response.jl:21` reads `ion_csa/latest.toml`, a last-writer pointer, instead of the
run identity derived from `elm_matched_csa_full.toml`.

**Consequence.** Safe in the matched flow (one branch per directory, `simulate_branch` requires exactly
one study) but not in shared results directories or after a second CSA run in the same directory, where
the assembly could pair the electron arm with a different CSA run that passes the event-identity checks.

**Proposed remedy.** Regenerate `jobs.toml` when absent; select studies by identity; compare `raw_sha256` before
skipping; make the assembly resolve the CSA manifest from the run card identity (or verify that
`latest.toml` `run_identity` matches the card-derived identity).

**Validation.** Fixture test: two studies in one directory; interrupted prepare; mutated raw.h5.

**Uncertainty.** Low; behaviour inferred from code, not exercised.


#### W-05 · Charge-quadrature retry restarts the whole Julia reduction once per failing history

| Field | Value |
|---|---|
| Severity | minor |
| Category | optional improvement |
| Topic | Ion CSA chain |
| Location | `results/elm_matched_production/csa.log:609-907`; `src/ion_signal/ion_csa_response.jl:285-291`; `scripts/production/run_ion_csa_pilot.py` (retry loop) |
| Related | W-06 |

**Evidence.** The production log shows six cycles of `Julia stage reduce` → `ERROR: LoadError: Charge
quadrature failed: … error 2.5e-5` → `Refining current quadrature: <case> step 2.5e-06` → restart. Each
cycle re-ran the reduction stage for one failing case (errors 1.3e-7 … 3.3e-4 against the 1e-7 bound).
Previously reduced histories were cached, so no event was recomputed, but the stage restarted six times.

**Consequence.** Wall-clock and log noise; a run with many marginal cases would restart many times. The stage boundary is also confusing when reading the log (each `ERROR` is expected behaviour).

**Proposed remedy.** Collect all quadrature failures in one reduction pass, refine them together, then re-run once; log the retry as a warning rather than a top-level `ERROR`.

**Validation.** Replay the six saved cases with the batched retry; identical `response.h5` content.

**Uncertainty.** None on the observed behaviour.


#### W-06 · Numerical knobs of the CSA pipeline are hard-coded rather than card-controlled

| Field | Value |
|---|---|
| Severity | minor |
| Category | optional improvement |
| Topic | Ion CSA chain |
| Location | `src/ion_signal/ion_csa_response.jl:167,256-257,285-291,404,445`; `scripts/production/run_ion_csa_pilot.py:71-74,90-93,99,116` |

**Evidence.** PWL charge check 1e-8, SPICE charge check 1e-7, refinement tolerance 0.01 ENC, storage
threshold 0.001 ENC, `maxframes=64`, 8 iterations, 12 refinements, 120/240 s timeouts with gear fallback,
retry step `min(oldstep/4, current_step_s/2^k)` with at most 3 attempts. All are covered by the code
hash, so identity is intact.

**Consequence.** Changing a tolerance requires a code edit and changes the code hash of every study; the values are not visible in cards or reports.

**Proposed remedy.** Move them to the response card with their current values as defaults, and print them in the report's resolved-card block.

**Validation.** Smoke run with the defaults reproduces the current identity-independent numbers.

**Uncertainty.** None.


#### W-11 · Code-hash coverage lists omit two included sources

| Field | Value |
|---|---|
| Severity | info |
| Category | optional improvement |
| Topic | Ion CSA chain |
| Location | `src/ion_signal/ion_csa_response.jl:74`; `src/ion_signal/ion_csa_pipeline.jl:220-230` |

**Evidence.** The hashed file lists omit `src/event_analysis/EventAnalysis.jl` (included at `ion_csa_response.jl:9`) and `src/ion_signal/spice_runtime.py` (resolves the simulator).

**Consequence.** Edits to those files do not invalidate cached studies although they can change behaviour.

**Proposed remedy.** Add both files to the lists (accepting one-time identity change) or document why they are excluded.

**Validation.** Touch a comment in `EventAnalysis.jl`; the production code hash must change.

**Uncertainty.** None.


#### W-12 · The fast CSA model is level-independent linear superposition; direct SPICE carries nonlinear state across buckets

| Field | Value |
|---|---|
| Severity | info |
| Category | optional improvement |
| Topic | Ion CSA chain |
| Location | `src/ion_signal/pad_history_readout.jl:9-40`; `src/ion_signal/ion_csa_response.jl:195-215,375`; `src/ion_signal/ion_csa_correction.jl:144-147`; `src/ion_signal/ion_csa_pipeline.jl:126` |

**Evidence.** The kernel is the N=100 minus N=0 waveform divided by −100 with 5 ms averages at 0, T, 2T and an
exponential tail `α=exp(−T/τ)`; `response` propagates `oldr=α·oldr+raw[3]·q`, so baseline recovery is
modelled but gain and τ are level-independent. Direct SPICE simulates the whole pad history from the
quiescent point (`V(out)[0]==0.8863611214728842`, equal to the stored baseline); each (event, pad) history
is independent with no inter-pad coupling; noise is stationary and level-independent. Example
`spice_cases/bb0nu_p4_0_-10_4`: linear minimum 0.546 V, direct `maximum_error_ENC` 0.957, i.e. about 1 ENC
compression at 0.55 V. `model_used` codes are `{0,1}`; the model-version strings are compared by equality at
`ion_csa_response.jl:65-66`; the three-frame inverse is applied to every mean with small-signal taps.

**Consequence.** Context for the nonlinear work rather than a defect: an alternate fast mode plugs in at `P.response`
(line 195), needs a new `model_used` code, a new `model_version`, and must be validated post-inverse as
`validate_controls` already does for controls. A memoryless per-bucket correction would not reproduce the
state carried by direct SPICE.

**Proposed remedy.** None for the baseline; the second review should check the nonlinear implementation against these insertion points.

**Validation.** n/a

**Uncertainty.** None.


### Topic: SiPM and electronics chain

#### S-01 · Source hits outside the 150 cm chamber enter the electron chain, and the ASIC acquisition window scales with the event's z extent

| Field | Value |
|---|---|
| Severity | major |
| Category | supported concern |
| Topic | SiPM and electronics chain |
| Location | `events/0nubb_hitdf.h5` (label `NEXT_PRECDR_XENON`); `src/event_analysis/production_readout.jl:245`; `scripts/production/sipm_pipeline.jl:110-123`; `scripts/production/reduce_elm_sipm_pilot.jl:110-131` |

**Evidence.** The NEXUS sample spans x, y, z ∈ [−2000, 2000] mm. Of 8,995 source events, 15 have a z extent above
1,500 mm (65 above 1,000 mm; median 42.5 mm, 99th percentile 895 mm), and 117 have more than 1% of their
energy farther than 750 mm from the median z (`work/repository_review/logs/06_source_extent_check.log`).
Carriers are sampled "into physical grid cells without clipping"; the SiPM acquisition window is the global
first…last photon bucket ±4 over all illuminated channels. Saved ASIC shards show bucket spans 71…10,287
(median 223), 14 events span more than 6,008 buckets (1.5 m at 0.25 mm) and 309 events have reference drift
plus span above 150 cm. Event 6544 (electron drift 143.9 cm): 61% of its 1.52 M photons sit 2.3 m from the
reference; 374 channels × 9,447 buckets = 3.53 M acquired cells, expected dark 618,306 PE exceeding the
486,859 detected PE, summed noise ≈6,657 PE. Per-event acquired cells range 3,597 … 3.53 M (median 40,879).

**Consequence.** (a) The ASIC-stage energy observables are sums over the whole acquisition rectangle, so their widths
(0.294%/0.567% all events, 0.292% retained, vs 0.202% after the SiPM) are dominated by noise and dark integrated
over a window whose size is a bookkeeping policy, not detector physics. (b) Deposits that would lie outside any
1.5 m chamber are read out as if inside; the deposited-energy normalisation uses all hits, so the ratio hides it.
The ion arm is affected in the same way for z extent (no clipping) but has no window-noise term.

**Proposed remedy.** Ask the user whether hits outside the drift volume should be removed or flagged at the detector stage
(fiducial selection), and whether an ASIC energy estimator restricted to a defined signal region (photon-positive
cells plus a fixed margin) should accompany the full-window sum. State the all-acquired-cells definition in the
report.

**Validation.** Recompute `readout_corrected` on photon-positive cells and compare widths; count events whose hits exceed the chamber.

**Uncertainty.** The numerical evidence is firm; whether out-of-chamber hits are intentional (NEXUS volume larger than the TPC) is a design question for the user.


#### S-02 · Energy-width estimators and dark subtraction are undocumented; RMS/mean is not FWHM

| Field | Value |
|---|---|
| Severity | minor |
| Category | documentation gap |
| Topic | SiPM and electronics chain |
| Location | `scripts/production/reduce_elm_sipm_pilot.jl:10,110,159-176`; `scripts/production/present_matched_analysis.py:93`; `docs/sipm_electronics_model.tex` §Dark subtraction |
| Related | S-01 |

**Evidence.** Widths are `100·std/mean` (population RMS/mean). The retained subset is `clipped_acquired_cells==0`
including dark/noise-only channels. `dark_expected = acquired_cells·R_d·Δt` depends on the acquisition size
(S-01). Corrections are applied once (`detected/survival/energy`, `(after−dark_expected)/survival/energy`), but
the four estimator formulas and the RMS-vs-FWHM factor (2.355 for a Gaussian) appear in neither the .tex nor the .md.

**Consequence.** The `AGENTS.md` target “0.5% at 2.5 MeV” is FWHM-like; comparing it with RMS/mean numbers without the factor misleads.

**Proposed remedy.** Add the formulas and the RMS/FWHM statement to `docs/sipm_electronics_model.tex` and the matched report definitions paragraph.

**Validation.** n/a

**Uncertainty.** None.


#### S-06 · Attachment survival uses one reference drift per event for all hits

| Field | Value |
|---|---|
| Severity | info |
| Category | documentation gap |
| Topic | SiPM and electronics chain |
| Location | `scripts/production/run_electron_250ns.jl:193-198`; `scripts/production/reduce_elm_sipm_pilot.jl:109` |
| Related | S-01 |

**Evidence.** Survival `exp(−t/τ)` with `t` from the event reference drift is applied per hit (binomial), and the reducer inverts exactly that law; units checked (cm→mm→µs→ms). Hits displaced in z (S-01) share the reference survival.

**Consequence.** Consistent and unbiased in expectation; an upstream simplification worth stating with the fiducial decision.

**Proposed remedy.** One sentence in `elm.tex` near line 494.

**Validation.** n/a

**Uncertainty.** None.


### Topic: Production viewer

#### V-01 · Electron-panel totals are restricted to photon-occupied cells without saying so

| Field | Value |
|---|---|
| Severity | major |
| Category | documentation gap |
| Topic | Production viewer |
| Location | `src/cloud_viewer/production.py:142-157`; `src/cloud_viewer/web/app.js:20`; `src/cloud_viewer/web/index.html:42` |

**Evidence.** `px,py,pb` are read from the photon file and `values=group[field][()][positions]` samples the sensor/ASIC arrays
only at photon cells. Live event 1: viewer `readout` sum 492,196.6 PE over 14,402 cells, while the ASIC file's
`readout_pe_total` is 504,741.6 over 89,250 acquired cells (sensor 491,582.8 vs 504,589.8). Meanwhile
`clipped_buckets=count_nonzero(a['saturated'])` counts all acquired cells. Titles and the card meta show the
number with no qualifier.

**Consequence.** A reader comparing the viewer with the production summary sees a ≈2.5% deficit (dark and noise in quiet cells) with no explanation. The selection itself follows the `AGENTS.md` common-cells rule.

**Proposed remedy.** Add “photon-positive cells only” to the electron title/meta and footnote; optionally show `readout_pe_total`; state that clipping counts cover all acquired cells.

**Validation.** `work/event-viewer-venv/bin/python -I work/repository_review/viewer/compare_live.py` reproduces both numbers.

**Uncertainty.** Low.


#### V-02 · Gas identity is renamed `xe_he15` → `xe_he` by the assembler and viewer and never cross-checked against HDF5 metadata

| Field | Value |
|---|---|
| Severity | major |
| Category | supported concern |
| Topic | Production viewer |
| Location | `scripts/production/assemble_matched_response.jl:72` (`"default_gas"=>"xe_he"`); `src/cloud_viewer/production.py:38-40,155`; `src/cloud_viewer/web/app.js:7,105,136`; `config/specs/runs/elm_matched_response_viewer.toml:7` |

**Evidence.** `gas = 'xe_he' if shard['gas'] == 'xe_he15' else shard['gas']` is the only place the manifest key meets the card
key; the label "Xe + 15% He" is hard-coded for anything not `pure_xe`. Every data file carries the identity
(`sensor/metadata/gas_model = b'xe_he15'`, `events/N/gas_model`, photon `metadata/branch`), but the viewer reads
none of them, and the per-event identity checks (`pair_checksum`, drifts) are gas-independent. Today the mapping
is correct (verified live).

**Consequence.** A mis-written manifest gas field or a sensor shard pointing at `photons_pure_xe.h5` of the same realisation would display and export (`…_xe_he_…png`) as Xe–He without error. The UI/PNG key differs from the production key, against the `AGENTS.md` rule to preserve gas identity in outputs.

**Proposed remedy.** Use the manifest key (`xe_he15`) end to end, or keep the alias but assert `gas_model`/`branch` metadata equal `shard['gas']` in `event_payload` and compare `sensor['metadata/source_sha256']` with the production manifest.

**Validation.** Fixture test with `gas_model=b'pure_xe'` in the sensor file and `gas="xe_he15"` in the manifest expecting `ValueError`.

**Uncertainty.** Low on the code path; the trigger needs a mis-written manifest.


#### V-03 · Viewer start-up fails with a raw traceback when referenced data are absent

| Field | Value |
|---|---|
| Severity | minor |
| Category | supported concern |
| Topic | Production viewer |
| Location | `src/cloud_viewer/production.py:17,26,31`; `src/cloud_viewer/server.py:197`; `README.md` (paired production viewer section) |
| Related | W-09 |

**Evidence.** The card is repo-relative, but the SiPM manifest and `sensor['metadata/source_path']` are absolute, and `results/simulations/` is gitignored; `h5py.File` at `__init__` raises `FileNotFoundError` before the server binds, with no handler in `serve()`.

**Consequence.** `make response-viewer` is documented as the entry point without stating that it only works where the production data exist.

**Proposed remedy.** Check existence of all referenced files in `ProductionSource.__init__` and raise one `ValueError` listing them; state the data prerequisite in the README.

**Validation.** Launch with a card pointing at a non-existent `ion_response.h5`.

**Uncertainty.** None.


#### V-04 · `CloudSource` is dead code with a different indexing convention

| Field | Value |
|---|---|
| Severity | minor |
| Category | optional improvement |
| Topic | Production viewer |
| Location | `src/cloud_viewer/server.py:20-114,156,194-210`; `tests/event_viewer/check_cloud_viewer.py:8,12` |

**Evidence.** No instantiation of `CloudSource` exists; only `CloudSource.project` is reused; its `event_payload` is 0-based while `GridSource`/`ProductionSource` are 1-based. `check_cloud_viewer.py` is a manual script against `results/grid_readout_direct_8000`.

**Consequence.** Two indexing conventions in one file; nothing checks `ProductionSource` on real data.

**Proposed remedy.** Remove `CloudSource` (move `project` to a module function) or document it as the retired view; add a real-data smoke test for `ProductionSource` guarded by file presence.

**Validation.** n/a

**Uncertainty.** Low.


#### V-06 · PNG export names use the file position, not the MC ID, and re-render at a different aspect

| Field | Value |
|---|---|
| Severity | minor |
| Category | documentation gap |
| Topic | Production viewer |
| Location | `src/cloud_viewer/web/app.js:160-161,173,179` |

**Evidence.** `String(current).padStart(5,'0')` yields `…_xe_he_event00001_response_readout.png` for MC ID 0; the MC ID appears only in the drawn meta string. `Plotly.toImage(div,{width:565,height:550})` renders at a different aspect than the on-screen plot, so with `scaleanchor` the exported range differs from what the user framed.

**Consequence.** Exported figures are hard to match to MC events and may not show the framed region.

**Proposed remedy.** Include `mc<id>` in the filename and headline; export at the on-screen aspect or document the re-framing.

**Validation.** Export event 1 and compare with the screen.

**Uncertainty.** Low (reasoned from code, not executed in a browser).


#### V-08 · Minor robustness notes in the viewer server

| Field | Value |
|---|---|
| Severity | info |
| Category | optional improvement |
| Topic | Production viewer |
| Location | `src/cloud_viewer/server.py:266,199`; `src/cloud_viewer/production.py:16,19,51,118` |

**Evidence.** The handler catches only `ValueError, KeyError, OSError`, so a `StopIteration` from an unknown operating-point id yields an HTTP 500; the positional default `operating_point="legacy"` collides with the card `id` semantics; `sample` is checked against shards but not against `ion_response.h5 /metadata/sample`.

**Consequence.** Unclear errors on mis-configured cards.

**Proposed remedy.** Broaden the handler, default to `None`, add the sample check.

**Validation.** n/a

**Uncertainty.** None.


### Topic: Tests and validation

#### W-07 · Python CSA checkers are pinned to pilot sizes; no independent validation record for the matched run

| Field | Value |
|---|---|
| Severity | minor |
| Category | documentation gap |
| Topic | Tests and validation |
| Location | `tests/production/check_ion_csa_response.py:15`; `tests/production/check_ion_csa_correction.py:36`; `tests/production/check_ion_csa_pipeline.py` |

**Evidence.** `len(ids)==100` and `events==200` are asserted, so the checkers cannot run on the 8,000-event
matched run; `find results/elm_matched_production -name independent_validation.json` returns nothing.
The Julia `validate_response` did run inside `run` (`ion_csa_pipeline.jl:361`) and replays every sparse cell.
None of the `tests/production/check_*.py` or `tests/event_analysis/check_*.py` scripts are wired into a Make target.

**Consequence.** The independent Python cross-check of the published matrices exists but was not applied to the current production, and its result is not kept.

**Proposed remedy.** Parametrize the expected counts from the manifest; run `check_ion_csa_pipeline.py` on the matched manifest and keep the JSON next to it; add a Make target that runs the check scripts.

**Validation.** Run the checker on `results/elm_matched_production/ion_csa/359613d487dfeee8/manifest.toml`.

**Uncertainty.** None.


#### S-04 · Statistical SiPM tests lack dense-path, illuminated-wrapper and regime-limit checks

| Field | Value |
|---|---|
| Severity | minor |
| Category | optional improvement |
| Topic | Tests and validation |
| Location | `tests/test_sipm_statistical.jl`; `tests/test_elm_sipm_pilot.jl:18-48`; `work/repository_review/sipm/check_illuminated.jl` |

**Evidence.** Tested: occupancy moments by enumeration, cancellation, `recover_totals`, beta moments, two-trigger quadrature,
binomial PDE, charge conservation, reproducibility. Not in `tests/`: dense path (>16 triggers) vs microscopic
reference (only in the benchmark script), `simulate_channel_illuminated` with finite τ against a two-sensor
reference, label-partition moments, bucket-width change, saturation. The sub-review's scratch script covers these
and agrees with the microscopic model (e.g. 13,286.2±0.8 vs 13,286.6±0.8 PE for the densest production bucket;
−0.12% midpoint bias only at 0.4 triggers per cell per substep, outside the production regime).

**Consequence.** Regressions in the dense or two-population paths would not be caught by `make sipm-test`.

**Proposed remedy.** Promote checks A, C, D, E, F of the scratch script to tests with Monte Carlo tolerances; record the regime limit in the .md.

**Validation.** `make sipm-test` runs them.

**Uncertainty.** None.


#### V-05 · Viewer regression covers synthetic fixtures only and is outside `make test`

| Field | Value |
|---|---|
| Severity | minor |
| Category | optional improvement |
| Topic | Tests and validation |
| Location | `tests/test_production_viewer.py:14-116`; `Makefile` (`response-overview-test`) |

**Evidence.** Five tests pass (0.06 s, tempfile only): edge mapping for ±1 direction, reversed electron rows, relative times,
one-based position → MC IDs, 1.25× ranges, unknown-gas rejection, per-stage sums, titles, clipping fields, drift
mismatch detection. Not covered: the photon-cell restriction (V-01), gas mapping for `pure_xe`/two gases,
`chamber_length_cm` sum check, `pair_checksum` mismatch, the HTTP layer, `info()`, relative-path resolution,
missing files, all of `app.js`.

**Consequence.** Display-level regressions (labels, PNG names, gas alias) are untested.

**Proposed remedy.** Add the listed cases; include the file in `make test` discovery (it is already a `test_*.py` module, so `unittest discover` runs it; the Makefile target only duplicates it).

**Validation.** `make test`

**Uncertainty.** None.


#### T-01 · Julia test suites fail in a fresh checkout because `work/` and `events/` are assumed to exist

| Field | Value |
|---|---|
| Severity | major |
| Category | confirmed defect |
| Topic | Tests and validation |
| Location | `tests/test_ion_csa_response.jl:5`; `tests/test_ion_csa_correction.jl:5`; `tests/test_ion_csa_pipeline.jl:79`; `tests/test_matched_realisation.jl`; `tests/event_analysis/test_production_specs.jl:7`; `scripts/production/specs.jl` (`realpath` of sample paths) |

**Evidence.** In the baseline worktree (fresh checkout, no ignored data) the three CSA suites errored with
`IOError: mktempdir("<root>/work"): no such file or directory` (14, 24 and 0 passes respectively); after
`mkdir work` all pass (`work/repository_review/logs/02_julia_tests_worktree.log`, `05_julia_tests_rerun.log`).
`test_matched_realisation.jl` and `tests/event_analysis/runtests.jl` errored with
`IOError: realpath(".../events/electron_0nubb_energy_hitdf.h5")` because the `events/` directory holds only
ignored files; with the hit tables linked read-only they pass (8/8 and 2052/2052).
`tests/test_elm_streaming_reduction.jl` needs `results/elm_sipm_pitch4mm_pilot/photons_ep2.h5` (ignored) and fails
with an HDF5 "file may not exist" message. Python `unittest discover` passed 70 tests in the same worktree.

**Consequence.** A fresh clone cannot run `make ion-csa-test`, `make elm-matched-test` or `make test-event-core` without undocumented manual steps, and the failures are raw IOErrors rather than prerequisite messages.

**Proposed remedy.** `mkpath(joinpath(ROOT,"work"))` in the tests (or a shared test helper); check sample-file existence in `load_specs` and raise a message naming `make extract-hits`; guard data-dependent tests with a clear skip/prerequisite message; document which tests need local data.

**Validation.** Re-run the suites in a fresh worktree without `work/` and `events/`.

**Uncertainty.** None (reproduced).


#### T-02 · `make test` runs only Python; there is no aggregate Julia test target, and one Julia test writes into `results/`

| Field | Value |
|---|---|
| Severity | minor |
| Category | optional improvement |
| Topic | Tests and validation |
| Location | `Makefile` (`test`, `ion-csa-test`, `sipm-test`, `elm-matched-test`, `test-event-core`, `response-overview-test`); `tests/test_sipm_response.jl:186-188` |

**Evidence.** `test:` runs `unittest discover -s tests`; Julia suites are spread over seven separate targets. `test_sipm_response.jl` writes `results/sipm_response_validation/validation.h5` and a `.gitignore` into `results/` during a unit test (observed in the worktree).

**Consequence.** No single command validates the repository; a unit test leaves artefacts in the results tree.

**Proposed remedy.** Add a `test-julia` aggregate target and a `check` that runs both; write the validation file under `work/`.

**Validation.** `make test test-julia` in a fresh worktree.

**Uncertainty.** None.


### Topic: Paths, Git tracking and organization

#### W-09 · Absolute machine paths are embedded in tracked manifests, plans and job files and consumed verbatim

| Field | Value |
|---|---|
| Severity | major |
| Category | confirmed defect |
| Topic | Paths, Git tracking and organization |
| Location | `results/elm_matched_production/manifest.toml`; `results/elm_matched_production/ion_csa/359613d487dfeee8/manifest.toml` (`calibration_file`, `responses[].path`); `results/elm_matched_production/response/77ed411cb39be020/manifest.toml` (shard `sensor`/`asic`); `results/ion_csa_*/*/jobs.toml`; `scripts/production/assemble_matched_response.jl:62-69`; `src/ion_signal/ion_csa_pipeline.jl:290`; `src/ion_signal/ion_csa_response.jl:278,439`; `scripts/production/reduce_matched_analysis.jl:26,32`; `src/ion_signal/project_paths.py:resolve_recorded_path` |
| Related | V-03, O-04 |

**Evidence.** 133 tracked files contain `/Users/jjgomezcadenas/...`; 80,040 occurrences sit in 70 TOML files, mostly
`folder = "/Users/.../spice_cases/..."` lines in `results/ion_csa_*/*/jobs.toml` (one file has 29,339 lines, 5.5 MB).
The combined manifest stores absolute paths for every history, the SiPM manifest for every shard, the CSA
manifest for `calibration_file`; 15 HTML provenance pages and `latex_generated/electron_readout/sipm_response/sources.json`
also contain them. None occur in `src/`, `scripts/`, `tests/`, `Makefile` or `config/`. Consumers use the paths
verbatim: `folder=job["folder"]` → `isfile(rawpath) || error("Missing acquisition …")`, `h5open(graph["ion_response"])`,
`check_selection` opening `ledger.source_path`. `resolve_recorded_path` re-roots an absolute path only when it lies
under the current root, so a clone elsewhere gets the foreign path back unchanged. The viewer card is the exception
(relative, via `relpath`).

**Consequence.** `make elm-matched-analysis`, any CSA re-reduction and the viewer fail on a relocated checkout or another user's machine with an HDF5 open error rather than a prerequisite message, even when the data are copied; provenance files leak the user's home layout.

**Proposed remedy.** Write repository-relative paths at the writers (the `spice_runtime.source_path()` convention exists) and resolve with `resolve(root,path)` at read time; for existing manifests add a read-boundary re-rooting (strip any `…/NAUSICAA/` prefix or store `recorded_root`) without rewriting recorded bytes.

**Validation.** Run `reduce_matched_analysis.jl` and `make response-viewer` from a copy of the repository at a different path; `git ls-files | xargs grep -l /Users/` must not grow after the change.

**Uncertainty.** Low for the cited consumers; not every consumer of `sensor`/`sipm_manifest` keys was enumerated.


#### V-07 · Duplicated helpers across the three viewers and hard-coded labels

| Field | Value |
|---|---|
| Severity | info |
| Category | optional improvement |
| Topic | Paths, Git tracking and organization |
| Location | `src/event_viewer/server.py:23-67` vs `src/pad_history_viewer/engine.py:28-70`; `src/cloud_viewer/server.py:38-50` vs `production.py:95-98`; `scripts/viewers/serve_cloud_viewer.py:7,14`; `web/app.js:7-8,105,136`; `web/index.html:16-23` |

**Evidence.** `write_toml` and `Worker` (≈40 lines) are duplicated; projection by `np.add.at` appears twice; ROOT is computed three times and `sys.path` inserted twice; sample/gas labels are hard-coded although `/api/info` supplies them. `src/event_viewer` and `src/pad_history_viewer` exist beside `src/cloud_viewer` but are not in the handoff's directory table.

**Consequence.** Maintenance cost; label drift between viewers.

**Proposed remedy.** A small shared `viewers/common.py`; let `/api/info` carry display labels; list the two other viewers in the organization documentation.

**Validation.** n/a

**Uncertainty.** None.


#### D-04 · Presentation script imports helpers from another script

| Field | Value |
|---|---|
| Severity | info |
| Category | optional improvement |
| Topic | Paths, Git tracking and organization |
| Location | `scripts/production/present_matched_analysis.py:11-12` |

**Evidence.** `from present_elm_sipm_pilot import table,tex_table,esc,fmt` relies on the script directory being on `sys.path`; `build_report` is imported after a `sys.path.insert` of `scripts/documentation`.

**Consequence.** Script-to-script coupling; moving or renaming a presenter breaks another.

**Proposed remedy.** Move the shared HTML/TeX table helpers into a small module under `scripts/documentation/` or `src/`.

**Validation.** n/a

**Uncertainty.** None.


#### O-02 · Large numerical binaries are tracked, with an inconsistent sibling policy

| Field | Value |
|---|---|
| Severity | major |
| Category | supported concern |
| Topic | Paths, Git tracking and organization |
| Location | `git ls-files results` (28,428 files, 657 MB logical); `results/three_cloud_sample_500/.gitignore`; `.gitignore` |
| Related | O-03 |

**Evidence.** 63 tracked `.h5` files total 516 MB: `results/three_cloud_sample_500/electrons_xe_he.h5` 48.7 MB,
`results/event_analysis/step4_scans.h5` 47.2 MB, `three_cloud_sample_500/bb0nu_ion.h5` 39.1 MB, `electrons_ion.h5`
30.0 MB, two `results/ion_tracks/*/analysis.h5` of ≈20 MB, and 24 `.h5` under `results/three_cloud_benchmark/`
(257 MB, the same physics repeated for 1/2/4/6 workers as a timing benchmark). In `three_cloud_sample_500/` the
ignore rules exclude `bb0nu_pure_xe`, `bb0nu_xe_he` and `electrons_pure_xe` but track their three siblings;
`event_analysis/step3_validation.h5` is ignored while `step4_scans.h5` is tracked. Four 6.3 MB `ion_csa.html`
pages with embedded figures are tracked. `papers/` totals 86.8 MB with a 13.4 MB maximum, within the 20 MB rule.

**Consequence.** Clone and fetch cost of about 0.6 GB for results alone, and the ignore rules do not express the actual policy. The 2026-11-05 untracking covered `results/ltspice` only.

**Proposed remedy.** User decision per directory; the most defensible candidates are the `three_cloud_benchmark` files and the four `three_cloud_sample_500` files. Untrack with `git rm --cached` plus ignore rules; do not rewrite history.

**Validation.** `git ls-files -z results | xargs -0 ls -l | awk '$5>5e6'` after the change.

**Uncertainty.** Whether some files are intentionally retained report inputs (`step4_scans.h5` feeds a report) was not traced per consumer.


#### O-03 · About 25,000 generated SPICE deck files are tracked while their inputs and outputs are ignored

| Field | Value |
|---|---|
| Severity | major |
| Category | supported concern |
| Topic | Paths, Git tracking and organization |
| Location | `results/ltspice/histories/**` (10,536 `.cir`, 8,683 `.inc`, 8,502 `.toml`, 57.5 MB); `.gitignore` comment "source decks are tracked" |
| Related | O-02 |

**Evidence.** The per-pad-history decks are generated by the pipeline from tracked cards; the PWL current files they drive (`/results/ltspice/**/*.pwl`) and the `.h5` outputs are ignored, so the tracked decks are neither complete inputs nor runnable in a fresh clone.

**Consequence.** Slow `git status`/grep, noisy diffs and documentation path checks; the tracking policy is ambiguous.

**Proposed remedy.** Decide whether decks are provenance (keep a hash manifest and ignore the decks) or inputs (then track the PWL files too); state the choice in `.gitignore`.

**Validation.** `git ls-files results/ltspice/histories | awk -F. '{print $NF}' | sort | uniq -c`.

**Uncertainty.** Not verified that every `.cir` references an ignored `.pwl`.


#### O-04 · Documented commands fail inside HDF5 in a fresh clone instead of naming the missing prerequisite

| Field | Value |
|---|---|
| Severity | major |
| Category | confirmed defect |
| Topic | Paths, Git tracking and organization |
| Location | `scripts/production/run_elm_matched_production.py`; `scripts/production/present_matched_analysis.py:25`; `scripts/production/reduce_matched_analysis.jl:26,32`; `tests/test_rogelio_v2.jl:4,24,32`; `tests/test_rogelio_v3.jl`; `tests/test_adc_estimate.jl`; `tests/test_bucket_study.jl:49` |
| Related | T-01, W-09, V-03 |

**Evidence.** `make elm-matched-analysis` and `make elm-matched-report` open ignored `results/elm_matched_production/*.h5`
and the absolute `ion_response` path with no existence check; four Julia data tests open ignored files through
`results/ltspice/<alias>` symlinks and are not invoked by any Make target. The handoff asks that missing data lead
to clear prerequisite instructions.

**Consequence.** In a clone the documented commands die with `unable to open file` and no hint of which target regenerates the data.

**Proposed remedy.** A small shared "require results" helper (Python and Julia) that names the regenerating target from `docs/generated_reports.json`; data tests print a clear skip message when inputs are absent.

**Validation.** Clone to a temporary directory and run `make elm-matched-report` and `make sipm-test`.

**Uncertainty.** Code paths confirmed; not executed in a clone without data (the worktree test in T-01 covers the Julia suites).


#### O-05 · Nested result and config aliases are load-bearing but documented as deck-only; the Julia path map is derived

| Field | Value |
|---|---|
| Severity | minor |
| Category | documentation gap |
| Topic | Paths, Git tracking and organization |
| Location | `src/ion_signal/project_paths.py:44`; `src/ion_signal/project_paths.jl`; `docs/path_migrations.toml`; `scripts/documentation/catalogue_reports.py:193`; `src/ion_signal/{noise_measurement.jl:14,23,bucket_noise_covariance.jl:22,discharge_calibration.jl:8,bias_sensitivity.jl:8,pitch_exercise.jl:8,16,29}` |

**Evidence.** The Python and Julia maps are identical (217 entries). `resolve_recorded_path` rewrites only strings beginning
with `LTSpice/` or `latex_generated/`; card keys such as `output = "results/ltspice/retention_scan"` and the
hard-coded `results/ltspice/...` paths in six `src` files rely on 37 tracked symlinks under `results/ltspice/` and
23 under `config/` (184 alias references in 53 TOML, 368 in 20 JSON manifests). `docs/ltspice_setup.md` says aliases
remain only for recorded acquisition-deck paths. `docs/path_migrations.toml` is written only by `make report-catalogue`.

**Consequence.** Removing an alias breaks active cards and library code; editing a JSON registry without running the catalogue desynchronises Julia from Python.

**Proposed remedy.** Extend the prefix filter to `results/` and `config/`, route the `src` hard-codes through the resolver, then retire the aliases; add a test asserting the TOML equals the JSON-derived map.

**Validation.** `find results LTSpice config -type l` empty after migration with tests green.

**Uncertainty.** Low.


#### O-09 · Library code lives under `scripts/production/` with duplicated helpers

| Field | Value |
|---|---|
| Severity | minor |
| Category | optional improvement |
| Topic | Paths, Git tracking and organization |
| Location | `scripts/production/sipm_production.jl` (included by 11 scripts); `sipm_pipeline.jl` (by 3); `specs.jl`; `src/ion_signal/ion_csa_pipeline.jl:11-14`; `src/ion_signal/pad_compression.jl:86` vs `scripts/production/sipm_pipeline.jl:165`; four `run_elm_*` scripts (`archive_photons`) |
| Related | D-04, V-07 |

**Evidence.** `str/schema/require/resolve` helpers are defined in four places; `atomic_h5` exists in two variants with different semantics (the script version also sets `metadata/complete=true`); `archive_photons` is defined four times in two variants. `AGENTS.md` places calculations in `src/`.

**Consequence.** Behavioural drift between copies; the SiPM pipeline cannot be tested as a library without including a script.

**Proposed remedy.** Move the SiPM pipeline modules and `specs.jl` to `src/ion_signal/` (or `src/production/`), keep one `atomic_h5` with an explicit `complete` flag.

**Validation.** Existing SiPM tests pass after the move; code hashes change once.

**Uncertainty.** None.


#### O-10 · Makefile resolves the simulator on every expansion and errors loudly when unconfigured

| Field | Value |
|---|---|
| Severity | info |
| Category | optional improvement |
| Topic | Paths, Git tracking and organization |
| Location | `Makefile:4`; `src/ion_signal/spice_runtime.py:50-53` |

**Evidence.** `LTSPICE ?= $(shell … spice_runtime --path)` is recursively expanded, so the resolver runs only when a simulator target expands it (verified read-only: it reads `Info.plist` and hashes the file). When nothing is configured it raises `FileNotFoundError` pointing to `docs/ltspice_setup.md`.

**Consequence.** Fine for non-simulator targets; on a fresh clone a simulator target prints a Python traceback rather than a one-line hint.

**Proposed remedy.** Catch the error in `--path` mode and print the hint to stderr with a non-zero exit.

**Validation.** `make -n ltspice-study` without `config/local/tools.toml`.

**Uncertainty.** Not tested without the user's local card.


#### O-11 · Environment reproducibility note

| Field | Value |
|---|---|
| Severity | info |
| Category | documentation gap |
| Topic | Paths, Git tracking and organization |
| Location | `work/event-viewer-venv`; `requirements.txt`; `config/tools.example.toml` |

**Evidence.** The venv's `pip freeze` lists conda-built file URLs, so it is not recreatable from the freeze; `requirements.txt` pins numpy 2.4.4, scipy 1.17.1, matplotlib 3.10.9, h5py 3.16.0, plotly 6.9.0 and is consistent with the installed versions. `python -m ion_signal.spice_runtime` reports LTspice 17.2.4 with the pinned SHA.

**Consequence.** None today; a new machine must rebuild from `requirements.txt`, which is adequate.

**Proposed remedy.** State in the README that `requirements.txt`, not the venv freeze, is the reproducible specification.

**Validation.** n/a

**Uncertainty.** None.


### Topic: Documentation and reports

#### W-08 · CSA production documentation describes the earlier 4 mm run, not the matched run

| Field | Value |
|---|---|
| Severity | minor |
| Category | documentation gap |
| Topic | Documentation and reports |
| Location | `docs/csa_4mm_production.md`; `docs/elm_matched_production.md:16` |

**Evidence.** `csa_4mm_production.md` quotes 29,339 histories, 715,442 pixels and 0.39523/0.39136 ENC from the earlier
run. The matched run has 697,133 pixel histories, 29,650 histories below 0.63 V plus 510 controls (29,920
direct), 6 charge retries, 12 refinements. The two validation numbers in the handoff describe different
stages: 0.437537906541256 ENC is the pre-inverse scalar error of the filtered study
(`fast_path_validation.toml` and `acquisition_summary.json`, `assess` in `ion_csa_response.jl:319-339`),
0.4336132729993806 ENC is the post-inverse control error of the corrected study
(`validate_controls` in `ion_csa_correction.jl`). `elm_matched_production.md` states only the qualitative policy.

**Consequence.** Readers can attribute the earlier run's statistics to the matched production, or compare the two ENC numbers as if they measured the same stage.

**Proposed remedy.** Generate these counts and both ENC numbers from the saved summaries into the matched report/HTML (they are already in `matched_analysis.h5/validation`), and label the stage of each.

**Validation.** Regenerate the report; numbers appear from HDF5, not prose.

**Uncertainty.** None.


#### R-05 · Grid-reference conventions are only partly documented

| Field | Value |
|---|---|
| Severity | info |
| Category | documentation gap |
| Topic | Documentation and reports |
| Location | `scripts/production/run_detector_pilot.jl:10-17`; `scripts/production/run_ion_readout.jl:99-100`; `scripts/production/run_electron_250ns.jl:192-193`; `docs/production_formulas.tex:185-188`; `AGENTS.md` |

**Evidence.** Production uses the undiffused source-energy-weighted x/y centre, the track minimum z with direction +1 for
ions and the track maximum z with direction −1 for electrons, with an independent random transverse phase per
arm (saved per event: event 0 shares x/y (1203.98, 307.48), z −1465.4 ion / −1396.2 electron). The formula note
documents the ion convention only; `AGENTS.md` describes the interactive viewer's convention and says batch
studies keep their recorded one.

**Consequence.** Readers must infer the electron reference and the per-arm phase from code.

**Proposed remedy.** One paragraph in `production_formulas.tex` or `elm.tex` stating both references and the per-arm phase.

**Validation.** n/a

**Uncertainty.** None.


#### R-06 · AGENTS.md ELM design values differ from the production card

| Field | Value |
|---|---|
| Severity | info |
| Category | documentation gap |
| Topic | Documentation and reports |
| Location | `AGENTS.md` (electron EL grid: 5 mm holes, 10 mm pitch, E/P=3, 2 mm bucket); `config/specs/detector/electrons/elm_4mm_250ns_ep2_shared.toml:6-12`; `config/electron_transport.toml:19-26`; `docs/elm.tex:601-603` |

**Evidence.** The production card uses 2.5 mm holes, 4 mm pitch, E/P=2, 0.25 µs buckets and a square 3 mm SiPM; `elm.tex` records the post-pilot selection; `electron_transport.toml` still holds the earlier values.

**Consequence.** New readers of `AGENTS.md` start from superseded ELM geometry.

**Proposed remedy.** Point `AGENTS.md` to `elm.tex`/the matched cards as the current design and mark the earlier values as historical.

**Validation.** n/a

**Uncertainty.** None.


#### S-03 · The SiPM derivation note describes the microscopic model and old cards, not the statistical production

| Field | Value |
|---|---|
| Severity | minor |
| Category | documentation gap |
| Topic | Documentation and reports |
| Location | `docs/sipm_electronics_model.tex:293-296,361-363,506,598-600`; `scripts/production/sipm_production.jl:140-142,167-169`; `src/ion_signal/sipm_statistical.jl:141-154` |

**Evidence.** The note says the production follows individual avalanches and the compact model is proposed, quotes 0.5 PE
noise and says production saves the sum Q with the `sipm_v1`/`readout_v1` cards. The matched run uses the
statistical two-population model (`simulate_channel_illuminated`, not described in the note), 3.5417 PE noise, the
ELM card, and saves `optical_charge_pe_equiv` and `dark_charge_pe_equiv` separately; the "SiPM optical charge"
energy stage uses the optical component, which depends on the approximate label partition.

**Consequence.** The authoritative derivation document contradicts what is saved and analysed.

**Proposed remedy.** Update the quoted paragraphs; state that the optical/dark split is an approximation used as a diagnostic.

**Validation.** n/a

**Uncertainty.** None.


#### D-01 · `scripts/production/README.md` is stale and contradicts the Makefile

| Field | Value |
|---|---|
| Severity | minor |
| Category | documentation gap |
| Topic | Documentation and reports |
| Location | `scripts/production/README.md`; `Makefile:5` (`EVENT_JULIA_DEPOT`) |

**Evidence.** The README uses `JULIA_DEPOT_PATH="$PWD/work/julia_depot:…"` (underscore) and `python3.13`, and cites
`latex_generated/production_grid_scan/238ddc5a75e0fd9b/` and `latex_generated/electron_250ns/…`, which no longer
exist after the subject reorganization; the Makefile uses `work/julia-depot` (hyphen). Both depots exist locally
(`work/julia_depot` 228 MB, `work/julia-depot` 239 MB).

**Consequence.** Following the README creates or reuses a second depot and points to removed report paths.

**Proposed remedy.** Rewrite the README against the current Make targets (or replace it with a pointer to `scripts/README.md` and `docs/elm_matched_production.md`); remove the unused depot after confirming nothing references it.

**Validation.** Path-existence check of the README's backticked paths.

**Uncertainty.** None.


#### D-02 · `make note` fails on an unclosed math expression in the older physics note (author text)

| Field | Value |
|---|---|
| Severity | info |
| Category | documentation gap |
| Topic | Documentation and reports |
| Location | `notes/physics/nausica.tex:66`; `notes/physics/itaca_defs.tex:12` |

**Evidence.** Reproduced in the worktree: `! Missing $ inserted. l.66 ...on of $\bbonu\ lifetimes in excess of $10^`. `\bbonu` is defined with `\ensuremath`, so the surrounding `$…$` opens a math group that the `$10^{28}$` closes early. The `manuscripts/itaca3` manuscript and the matched report both compile (16 pages; `work/repository_review/logs/03_nausica_latexmk.log`, `04_report_builds.log`).

**Consequence.** The `make note` target fails; the active manuscript is unaffected.

**Proposed remedy.** Author decision: drop the `$` around `\bbonu` or rewrite the sentence. Not changed in this review.

**Validation.** `make note` after the edit.

**Uncertainty.** None.


#### D-03 · Status prints and card flags use misleading names

| Field | Value |
|---|---|
| Severity | info |
| Category | optional improvement |
| Topic | Documentation and reports |
| Location | `scripts/production/sipm_production.jl:254`; `scripts/production/specs.jl:67` (`pilot_only`); `config/specs/runs/elm_matched_detector_full.toml` (`pilot_only = true`); `config/specs/runs/elm_matched_sipm_response_full.toml` (`event_ids = [0]`) |

**Evidence.** The SiPM progress line prints `bb0nu/xe_he15: 99 / 7999` (last event ID of the shard over the last ID of the run, not a count). The 8,000-event production cards must carry `pilot_only = true` ("Pilot runner requires pilot_only=true") and the response template card lists `event_ids = [0]`, which the production run overrides.

**Consequence.** Readers of logs and cards misjudge progress and scope.

**Proposed remedy.** Print counts; rename or drop the `pilot_only` requirement for full cards; comment the template role of `event_ids`.

**Validation.** n/a

**Uncertainty.** None.


#### O-06 · Stale counts and statements in navigation documents

| Field | Value |
|---|---|
| Severity | minor |
| Category | documentation gap |
| Topic | Documentation and reports |
| Location | `docs/generated_reports.md` (layout table); `docs/scripts_inventory.md:9`; `docs/decisions.md:13,21`; `scripts/README.md:44-45` vs `docs/elm_matched_production.md:16` |

**Evidence.** The layout table says 52 ion_readout and 9 electron_readout collections (82) while the registry has 53 and 13 (87); the inventory says 138 commands while 233 `.py/.jl` files exist under `scripts/`; `decisions.md` cites Petit paths without prefix; the two READMEs' phrases "the saved production uses the SPICE/linear response" and "selected nonlinear histories ... receive transistor simulation" read as contradictory although both are true. All 29 documented Make targets exist.

**Consequence.** Minor misdirection for new readers.

**Proposed remedy.** Generate the counts in the catalogue step; prefix Petit paths; one reconciling sentence about SPICE for selected histories versus the planned fast nonlinear model.

**Validation.** n/a

**Uncertainty.** None.


#### O-07 · Historical documents without a status marker

| Field | Value |
|---|---|
| Severity | minor |
| Category | documentation gap |
| Topic | Documentation and reports |
| Location | `docs/detector_response_pipeline_plan.md`; `docs/response_production_spec.md`; `docs/script_reorganization.md`; `docs/ion_csa_handoff.md`; `docs/sipm_electronics_handoff.md` |

**Evidence.** These describe plans, a 3 mm/30 ms first case, a migration in present tense, and 2026-10-30 authorizations, without a status line. `docs/event_analysis_plan.md`, `docs/ltspice_organization_proposal.md`, `docs/ion_csa_response.md` and `docs/sipm_production_handoff.md` are properly marked.

**Consequence.** A reader may take superseded settings or authorizations as current.

**Proposed remedy.** One-line "Status: historical, superseded by X (date)" header on each.

**Validation.** n/a

**Uncertainty.** None.


#### O-08 · Manuscript generated values are a one-time copy with no regeneration or check rule

| Field | Value |
|---|---|
| Severity | minor |
| Category | supported concern |
| Topic | Documentation and reports |
| Location | `manuscripts/itaca3/itaca3.tex:13-20,384,544,582`; `manuscripts/itaca3/provenance.json`; `manuscripts/itaca3/Makefile` |

**Evidence.** The manuscript `\input`s `generated/*.tex` copied from `notes/physics/generated/`; no script writes `manuscripts/itaca3/generated`; 10 of 11 files are identical today and `spice_devices_table.tex` already differs (a stale `LTSpice/variants/` comment). Author-written content is correctly separate from generated includes.

**Consequence.** Manuscript macros can silently diverge from regenerated note values, against the "LaTeX includes generated values rather than copies" convention.

**Proposed remedy.** Add a `sync`/`check` rule to the manuscript Makefile that copies from `notes/physics/generated/` with hash comparison, or `\input` through a path variable.

**Validation.** `make -C manuscripts/itaca3 check` after regenerating the notes.

**Uncertainty.** None.


## Inspected and found consistent

Items inspected and found consistent with their documentation and with each other (file references
in the sub-review notes under `work/repository_review/*/agent_report_*.md`):

- **Workflow controller.** `run_elm_matched_production.py` runs the detector stage, then the SiPM and
  CSA arms concurrently with per-stage logs, then the assembler; `run_timings.toml` records wall times
  and the runner SHA. `run_matched_detector.jl` cross-checks event IDs, both drifts, `initial_pairs` and
  `pair_checksum` across the common ledger, ion and photon files before marking the index complete;
  `assemble_matched_response.jl` repeats the check across six histories (plus sensor and ASIC shards),
  asserts `gas=="xe_he15"` per shard, and writes the manifest and viewer card atomically. The combined
  manifest, `detector_files.toml`, `sipm_files.toml` and both stage manifests report `complete = true`
  with 8,000 events and `gases = ["xe_he15"]`.
- **Primary pairs.** Per-hit mean `1e6·E/W` with W = 21.9 eV, F = 0.15; both arms replay the same seed
  and error on mismatch; old (`c25084b9…`, drift ≥10 cm) and new (`4bd2a52b…`, 5–145 cm) ledgers share
  stream identity `7b0a08c6…` with identical 8,000 pair counts and checksums; the new ion drift equals the
  old minus 5 cm to 7.6e-6 (same uniform variate rescaled), confirming the primary-reference anchoring.
- **Placement.** `le = 150 − li` in Float32; li ∈ [5.0012, 144.971), KS distance 0.0079 against uniform.
- **Seeds.** Distinct labels `pairs`, `placement`, `ion_phase`, `electron_phase`, `ion_diffusion`,
  `attachment`, `xe_he15_diffusion_el`; arm code hashes include only that arm's runner so an arm-only edit
  invalidates only its descendants.
- **Diffusion.** Electron σ = d·√(L_cm/P_bar) with the card's d (xe_he15: 1.6/0.75 mm·√bar/√cm) equal to
  `config/electron_transport.toml`; ion σ = 10·√(2kT L/E) mm; scaling tests exist.
- **Direct sampling.** Floor on Float32 edges; electrons mirrored about the reference; Int64 checked sums;
  `IonPadHistory` enforces Σcount = initial pairs; buckets 0.1 m/s × 40 ms = 4 mm and 1 mm/µs × 0.25 µs =
  0.25 mm as stored; viewer edge reconstruction inverts the writer exactly (live event 1 checked).
- **Attachment.** Binomial per hit with `exp(−t/τ)`, electrons only; survivors/pairs 0.9518 vs analytic 0.9507.
- **EL and optics.** Yield 139·(2.0−0.80)·15·0.5 = 1251 photons/electron equals the saved ratio; optical
  acceptance 0.0118746 equals file metadata and `elm_values.tex`; transmission = 1 idealisations labelled
  as baseline/design inputs in `elm.tex` and the card.
- **Units.** NEXUS tables Float32 mm/MeV/ns with unit attributes checked; explicit conversions at SI
  interfaces; no silent Float32/Float64 mixing found in accumulators.
- **CSA sign and guard.** Positive ions pull the output below the 0.886 V baseline; candidates are
  histories with linear minimum below `initial_spice_guard_V = 0.63` (card line 25) plus controls; the
  0.5 ENC tolerance is used only in `assess` and `validate_controls`.
- **CSA input current and PWL.** `[F(t)−F(t−T)]/T` with `F(u)=e[φ(clamp(−vu,0,g))−φ(g)]`, ending at
  collection; no second impulse; no renormalisation; 5 µs grid on the ramps; `.tran` max step 250 µs
  with 125 µs refinement and 2.5 µs retry; the saved charge error of the inspected case is 4.4e-16.
- **Acquisition.** Batch LTspice with 120 s timeout and gear/240 s fallback recorded as
  `solver_method`; raw parser rejects stepped/FastAccess files; cache key per (sample, pitch, event, pad);
  `--reuse-cases` requires equal PWL, model, simulator and step SHAs; atomic writes; deterministic seeds.
- **Reduction and noise.** Exact piecewise-linear 5 ms averages; gain 4.1479e-5 V/ion; taps
  [0.01734, 1, 0.01504] removed exactly by the inverse; covariance from the one-sided ASD with |F|² and
  variance pinned to the full band; Cholesky per history length; `M⁻¹CM⁻ᵀ` with symmetry/positivity
  asserts; lag-1 correlation −0.273; control variance 3410 vs 3415 e²; signed estimates and no zero
  suppression enforced by the card.
- **SiPM statistics.** Occupancy moments, `recover_totals`, Samuelson bounds, stationary dark
  initialisation Beta(1, R_dτ/M), sparse/dense paths, illuminated/unilluminated split 21,808/18,176 of
  39,984 with dark rate apportioned by count, independent substreams; τ in s, dark in Hz, q₀ = 5.7678e-14 C,
  noise 3.5417 PE, clip 4334.38 PE; ASIC clips after noise and keeps negatives (minimum saved −18.4 PE).
  Independent two-sensor microscopic references agree (dense ELM bucket, state continuity, bucket-width
  change, optical/dark partition). Shard 0001 ledgers consistent between sensor and ASIC; PDE fraction
  0.320002.
- **Pulse statistics.** Photon-positive cells used at every electron stage and `n>0` cells at every ion
  stage; signed values; one mean and one maximum per event, as `AGENTS.md` prescribes.
- **Viewer.** One-based position → sorted MC IDs with an init-time assertion that the ion ID set equals
  every shard ID set; per-event identity checks between arms; drift sum vs chamber length; equal spatial
  scale (`scaleanchor`, `aspectmode data`); arrival time in hover data; live totals equal independent
  recomputation (ion response 110,617.13 with 41 negative cells; input 112,216 = `initial_pairs`).
- **Builds.** `manuscripts/itaca3` compiles (16 pages) with its own Makefile; the tracked matched report TeX
  compiles through `build_report.py`; LTspice setup documentation matches the directory layout
  (`LTSpice/{README.md,reference,designs,releases,notes,papers}`, eight subject directories under
  `latex_generated/`).

## Commands run and outcomes

All commands were run from the main checkout (read-only) or from the baseline worktree (writes). Logs
are under `work/repository_review/logs/`; sub-review scratch under `work/repository_review/{realisation,csa,sipm,viewer,organization}/`.

| # | Command (abridged) | Where | Outcome |
|---|---|---|---|
| 0 | `git status --short; git rev-parse HEAD 251118d41; git diff --stat 251118d41 HEAD` | main | clean; HEAD `d0b6476…`; baseline `251118d4…`; 8 doc files differ (`00_git_status_start.log`) |
| 1 | `git worktree add --detach work/repository_review/baseline_worktree 251118d41` | main | worktree at baseline, clean |
| 2 | `python -m unittest discover -s tests -v` (PYTHONPATH=src) | worktree | 70 tests OK, 0.43 s (`01_python_unittest_worktree.log`); wrote `results/sipm_response_validation/validation.h5` |
| 3 | 12 Julia suites (`test_ion_csa_*`, `test_sipm_*`, `test_matched_realisation`, `event_analysis/runtests`, `test_event_io`, `test_elm_streaming_reduction`, `test_response_overview`) | worktree, no `work/`, no `events/` | 7 pass; 3 CSA suites error on `mktempdir("…/work")`; matched realisation and event-analysis error on missing `events/*.h5`; streaming reduction errors on missing ignored data (`02_julia_tests_worktree.log`) |
| 4 | same suites after `mkdir work events` and read-only symlinks to `events/*.h5` | worktree | all pass except `test_elm_streaming_reduction.jl` (needs `results/elm_sipm_pitch4mm_pilot/photons_ep2.h5`) (`05_julia_tests_rerun.log`) |
| 5 | `latexmk -pdf … nausica.tex` | worktree | fails at line 66, `Missing $ inserted` (`03_nausica_latexmk.log`) |
| 6 | `build_report.py latex_generated/electron_readout/elm_matched_production/elm_matched_production.tex`; `make -C manuscripts/itaca3` | worktree | both compile; itaca3 16 pages (`04_report_builds.log`) |
| 7 | `make -n elm-matched-production elm-matched-analysis elm-matched-report elm-matched-test response-viewer ion-csa-test sipm-test test-script-layout` | worktree | all targets resolve to the documented scripts |
| 8 | h5py scan of `events/0nubb_hitdf.h5` z extents | main | 8,995 events; 15 with z span > 1,500 mm; volume ±2,000 mm (`06_source_extent_check.log`) |
| 9 | Sub-review scripts: `realisation/{inspect_h5,check_reference}.py`, `realisation/run_tests.sh`; `sipm/check_illuminated.jl`; `viewer/{inspect_h5b,compare_live}.py`; `python -m unittest tests.test_production_viewer` | main (read-only) / scratch | see sub-review notes; matched-realisation 8/8, grid-readout 27/27, production-readout 51/51, sipm-statistical 41/41, viewer 5/5 |
| 10 | `curl http://127.0.0.1:8769/api/info` and `/api/event` for events 1 and 7 | main | payloads equal independent recomputation; invalid queries return HTTP 400 |
| 11 | `lsof` on 8767/8768/8769/8774 | main | three viewers running; none touched |

Not run: any `make` target that writes `results/` or `latex_generated/` in the main checkout; LTspice;
any production or pilot script; `make note` (its prerequisites regenerate tracked outputs; the LaTeX step
was run directly instead).

## Coverage limitations

- The review inspected the matched 4 mm Xe–He production chain in depth. Earlier productions (pixel-charge
  studies, grid scans, operating-point pilots, CSA 2/3 mm runs, retention and device studies) were
  inspected only where their code is shared with the current chain or where documentation claims about
  them were checked. No claim is made about their numbers.
- The `src/event_analysis/petit/` reconstruction code (backbone, graphs, truth; stages A/B/C) was not
  audited beyond confirming that its tests pass (1,644 parity checks). The handoff's A/B/C item remains
  open for a later pass.
- The manuscript text in `manuscripts/itaca3/` and `notes/physics/nausica.tex` was compiled but not reviewed
  for scientific content; author text was left untouched.
- Sub-review conclusions about saved data (shard spans, live viewer totals, ledger comparisons) were
  produced by scratch scripts kept under `work/repository_review/`; the two with the largest consequences
  (S-01 extents, W-01 deck hashing) were independently re-derived here.
- Tests were run in a worktree with the main checkout's Julia depot and Python environment; a truly clean
  machine (no depot, no venv) was not tested.
- The untracked nonlinear-CSA files that appeared during the review were not read; the handoff reserves them
  for a second review once the user identifies a completed commit.

## Prioritized action list

Prioritized for discussion, not an automatic rewrite plan:

1. **Decide the fiducial policy (S-01).** Whether to drop or flag source hits outside the drift volume, and
   whether to add a signal-region ASIC energy estimator beside the full-window sum. This changes reported
   electron energy widths and affects both arms' z extents; nothing else in this list depends on it.
2. **Close the CSA identity gap now (W-01, W-02, W-03).** Hash the simulated deck/model/kernel file and record
   simulator identity in tracked outputs before the nonlinear calibration is fitted. Communicate W-01 to the
   implementation agent immediately.
3. **Fix the gas-identity path and viewer labelling (V-02, V-01, V-06).** Use `xe_he15` end to end or assert
   HDF5 metadata; label photon-cell restriction; put MC IDs in PNG names.
4. **Retire the superseded J anchors (R-03)** in `specs.jl`, `run_pilot.py`, both tests and the plan document.
5. **Make a fresh checkout testable and self-explaining (T-01, T-02, O-04, V-03, D-01).** Create `work/` in
   tests, give prerequisite messages for missing `events/` and ignored results (one shared helper naming the
   regenerating target), add an aggregate Julia test target, stop the unit test writing into `results/`,
   rewrite `scripts/production/README.md`.
6. **Portability of manifests (W-09, O-05).** Relative paths through `ProjectPaths` at the writers; extend the
   resolver's prefix filter so the retained aliases can be retired.
7. **Git tracking policy (O-02, O-03).** User decision per directory: the timing-benchmark HDF5 files and the
   generated SPICE decks are the main candidates; untrack without rewriting history.
8. **Documentation reconciliation (W-08, S-02, S-03, R-02, R-05, R-06, O-06, O-07, O-08).** Generate the matched
   CSA statistics and both ENC numbers from saved results; add the estimator formulas and RMS/FWHM statement;
   update the SiPM note's implementation paragraphs; fix or implement the Bernoulli contract; document both grid
   references; point `AGENTS.md` at the current ELM design; mark the five historical documents; add a manuscript
   `sync`/`check` rule for generated values.
9. **Robustness and organization refinements (W-04, W-05, W-06, W-11, R-01, R-04, R-07, S-04, S-05, S-06, V-04,
   V-05, V-07, V-08, D-03, D-04, O-09, O-10, O-11).** Lower priority; each has an identified practical benefit and
   a validation method in its finding.
