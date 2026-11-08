# NAUSICAA working conventions

This repository develops reproducible ion and electron readout models for the
ITACA TPC, including detector transport, optical collection and electronics.
Establish the physics before implementing electronics.

- Xe2+ means the molecular ion Xe_2^+, with charge +e, not atomic Xe^{2+}.
- Extend event readout to both carrier populations: a single-electron or
  double-electron ionisation track produces a trail of electrons and positive
  ions. Work to 2026-10-28 simulated the positive ions; electron drift to the
  anode is now part of the simulation scope. Trace the implemented electron
  diffusion laws and pure-Xe/Xe–He coefficient choices before implementation; preserve
  carrier and gas-model identity in configurations, caches and outputs.
  Electron readout operates at 15 bar; the Xe–He prescription represents 15% He
  (user design inputs). The production cards under config/specs/ select the
  operating conditions; config/electron_transport.toml records the earlier
  transport study. Cite diffusion formulas to papers or local derivations.
- Sample one event location along the 1.5 m TPC: ion and electron drift
  lengths are complementary, L and L_max-L. Sample each diffused carrier
  directly into its physical grid cell and longitudinal bucket; save sparse
  readout histories, seeds and provenance. Geometry, chamber length, lifetime,
  fields and time buckets are configurable. The current matched ELM production
  is specified below; earlier 10 mm electron pitch and E/P = 3 cards remain
  historical configurations. Store entry/exit electrode potentials in TOML.
  Route electrons geometrically to holes. Compute direct VUV collection
  acceptance analytically for the assumed electron path and sample that
  fraction. Save photons reaching each SiPM per time bucket. Apply SiPM PDE and
  electronics in a separate response stage. At 1 microsecond per mm the current
  0.25 mm electron bucket is 250 ns. SiPM PDE is TOML-configurable.
  These are user design inputs;
  record the adopted gas-dependent EL yield and optical boundary conditions.
- For EL variance, use the user's corrected approximate readings of Oliveira
  et al., Figure 6 (2026-11-04): J ≈ 0.030 at E/P = 2.0,
  J ≈ 0.015 at 2.5, and J ≈ 0.025 at 3.0 kV/cm/bar. The minimum
  is near E/P = 2.8, J ≈ 0.01. Use the same J prescription for both gases.
  Preserve historical J = 0.01 at E/P = 3 productions with their original cards;
  new comparison cards explicitly specify the corrected values.
- Use SI internally for transport, fields and electronics. Preserve particle-
  physics units in NEXUS event tables and hit DataFrames: positions in mm and
  deposited energy in MeV, with the original Float32 precision. Do not convert
  these tables to m or J. Convert explicitly at the interface to a calculation
  that requires SI. Show units in configuration keys, columns and figure axes.
- Implement the event reader and analysis in Julia. Keep HDF5 as the event-data
  storage format, including reduced hit tables; retain particle_id and hit_id.
- Event-analysis hit types use separate x_mm, y_mm and z_mm vectors. Coordinates,
  energy_MeV and time_ns use Float32; the same concrete EventHits type represents
  source and vertex-centred coordinates, labelled by CoordinateFrame. Identify
  hits by (event_id, particle_id, hit_id) within a sample. Use Float32 as the
  proposed numerical baseline and validate reduction errors against Float64
  references before finalising accumulators; the energy-resolution scale is
  0.5% at 2.5 MeV (12.5 keV).
- Store MC time, actual primary IDs and the common primary generation vertex in
  reduced HDF5 files. Keep hits in source coordinates through ion diffusion.
  Distinguish reconstruction A on pure MC, B after diffusion before pad
  discretisation, and C after diffusion and pad discretisation. Preserve and
  compare their components, backbones and decisions; A/B selections do not gate
  propagation to later stages. Interactive pad collection uses the exact
  retained carriers' ion-weighted x/y centre and minimum z as its grid
  reference, with one random transverse phase per event realisation. Extend
  the grid to cover all retained diffusion charge, including ions below
  reconstruction thresholds.
  Existing selected-event batch studies retain their recorded reference convention.
- Electronics track observables are ions per pixel per frame;
  time-integrated pixel charge belongs to energy analysis. Say "pixels with
  charge". Inclusive, event-mean and event-maximum pulse distributions use all
  positive pixel-frame cells. The n_min cut applies only to the minimum histogram:
  compute the true event minimum over all positive cells, then include the event
  only if that minimum > n_min. This display selection causes no charge loss.
  Current pixel-charge cards set n_min = 0, so all four distributions have no
  charge cut and the minimum histogram includes every event with pulses.
- Response mean/max comparisons select positive input channel–bucket cells and
  use the same cells at every response stage. Keep signed reconstructed values
  and apply no amplitude or signal-to-noise cut. Each event contributes one
  mean and one maximum.
- Production provenance: shared v2 primary ledgers hash the dedicated sampling,
  seed and ordered-hit-reading code; both arms verify that identity. Preserve RNG
  streams and historical ledgers when migrating; replay-verify counts, checksums
  and drift placement. CSA production verifies the actual deck/model/kernel
  against acquired calibration hashes. Report regeneration uses saved cards.
  See docs/reviews/provenance_fixes_2026-11-06.md.
- The active readout replays source hits and deterministic diffusion seeds when
  pitch or bucket width changes; it does not save individual carrier coordinates.
  Histogram cell centres are valid for display and historical comparison, not
  physical pad/time assignment. Version new HDF5 readout schemas and retain compact
  binned regression fixtures and provenance for audit. Record source contents, cards, seeds and code hashes
  so the sampled realization is reproducible.
- General analysis workflow: Julia performs physical calculations, data analysis
  and reduction, and saves numerical results and provenance to HDF5. Python
  reads these results to produce final figures, tables and LaTeX macros without
  repeating the physical calculations. Presentation can be regenerated without
  rerunning the analysis; provide a command to run the complete pipeline too.
- Always present analysis results in a generated, self-contained LaTeX report
  that compiles to PDF and gathers the relevant figures, tables, parameters,
  assumptions and validation results. Generate values from the saved results;
  do not transcribe numbers by hand. Keep numerical outputs in results/ and
  generated analysis reports in latex_generated/.
- When requirements, physical assumptions or analysis choices are unclear,
  consult the user before proceeding with the uncertain part of the work.
- Evaluate user-proposed methods on their merits. Offer a better alternative
  when one comes to mind, with its rationale; otherwise explain why the
  proposal works. Do not agree by default or invent complexity merely to
  offer an alternative.
- For event-analysis development, trace the loaded Petit implementation and
  its callers before proposing types or algorithms. Present existing behaviour,
  types/interfaces and the justified adaptation in the HTML review guide first.
  Review those prerequisites together before implementing the next interfaces.
- The revised design has 5 kV across the 10 mm IFG and 3 kV across the 5 mm
  pad gap, with 200 V/cm in the bulk. Older paper voltage schemes are superseded.
- Physical electric fields and weighting fields are different boundary-value
  problems. Fixed geometry weighting potentials do not depend on applied bias.
- State boundary conditions and approximations alongside every numerical result.
  The central square pad in a grounded plane is a benchmark, not a solved IFG.
- Preserve total signal charge: integral I dt = q Delta(phi_w). Do not append a
  second collection impulse. Do not silently renormalize finite-path integrals.
- A 10 cm straight track along drift is an extreme test, not the event ensemble.
- Keep event-analysis calculations in src/event_analysis and signal formation
  and electronics in src/ion_signal; scripts reproduce HDF5 analysis results
  and the Python presentation step. Existing Python/JSON benchmarks
  remain available; the Julia-to-HDF5-to-Python workflow applies to new analyses.
  LaTeX includes generated values rather than copies.
- Use meaningful numerical/physical tests: independent quadrature, boundary
  conditions, charge conservation, scaling and convergence.
- Run `make test physics` with a Python >=3.11 environment containing the pinned
  scientific dependencies. Run `make note` to compile the PDF with latexmk.
- Correct clear typos automatically. Keep prose edits small and show wording
  changes for review; preserve the author's voice and avoid substantial
  reformulation unless explicitly requested.
- In documentation, review pages and explanations, state what each type or
  function represents and does. Remove defensive exclusions and unsolicited
  contrasts (for example, "It does not contain the hit cloud"). Keep wording
  focused on purpose, data, units, behaviour and relevant interface contracts.
  Express input preservation and array ownership directly; document concrete
  validation errors where they are part of the contract.
- For LaTeX changes, compile and check the build log, but do not visually
  inspect the PDF unless the user explicitly requests it.
- Bound production file copies and temporary publication files. Use the shared
  `GuardedFiles` helper instead of `Base.cp` for production copies: bounded writes,
  exact size/SHA-256 verification, at most 10% growth for publication metadata,
  free-space checks with a 1 GiB reserve, atomic rename and cleanup on failure.
  For generated datasets, derive an explicit byte budget from the output schema;
  input-file size is not a valid universal bound. Investigate unexpected growth
  immediately. See docs/reviews/file_copy_failure_2026-11-08.md.
- Keep scratch files in work/.
- Storage cleanup records live in `docs/reviews/storage_cleanup_2026-11-08.md`
  and `results/storage_cleanup/2026-11-08/`. Retired corrected responses resolve
  through checksum-verified `.canonical.toml` receipts to published files;
  original and published hashes remain distinct. Keep these receipts with the
  retained outputs. Legacy binned replay fixtures are in
  `results/validation/direct_replay_fixture/`. Use `work/julia-depot` for the
  active local Julia depot.
- Cite primary sources in notes/physics/references.bib. Record user-provided
  revisions as design inputs, not as published measurements.
- Cite physics formulas to primary papers or documented local derivations,
  and point to our implementation when describing the calculation used.
  Petit establishes implementation lineage, not physical provenance.
- When the user says "commit", include the completed work in the current flow.
  Inspect staged, unstaged and untracked files; commit relevant code, tests,
  TOML, TeX, notes, documentation, papers and generated report artifacts.
  Check sizes and file types before staging. Keep large datasets and binary
  simulation outputs outside Git, and leave unrelated work untouched. Ask the
  user if the boundary of the current work is unclear. Never configure or use
  Git LFS without the user's explicit approval.
- Track and commit changes under papers/, including PDFs, with the relevant
  work. Check each file's size; ask before committing a paper larger than 20 MB
  rather than silently omitting it.
- Do not change the existing remote or publish/push without a user request.

## Current production and active development (2026-11-08)

- Read docs/elm_matched_production.md and its linked cards/manifests for the
  completed 8,000-event bb0nu production: ions at 4 mm pitch / 40 ms, Xe–15% He
  photons at 4 mm pitch / 250 ns. One shared realization places ion drift
  uniformly in 5–145 cm of the 150 cm chamber, with complementary electron
  drift. Shared primary-pair counts and checksums are preserved across arms.
  Electron attachment uses the configured lifetime (currently 15 ms); ions
  receive no electron-attachment loss. These are configurable design inputs.
- Matched 8,000-event ion productions at 3 mm / 30 ms and 2 mm / 20 ms are also
  complete, preserving the same primary pairs and 5–145 cm drift placement.
  See docs/matched_ion_csa_pitch_scan.md. The 2 mm publication was recovered
  using guarded copies after a confirmed Julia 1.12.6 sendfile-counter overflow;
  its numerical results were preserved.
- The current ELM has 2.5 mm diameter holes/guides, 5 mm EL depth, E/P = 2
  kV/cm/bar at 15 bar, and 0.5 mm radial xenon gaps around the guides. The
  optical model includes analytical hole acceptance, TPB conversion 0.65,
  forward fraction 0.5 and the documented meridional guide capture. See
  docs/elm.tex and config/specs/detector/electrons/elm_4mm_250ns_ep2_shared.toml.
  Photons saved at the SiPM are blue after TPB in this configuration; preserve
  the photon-band metadata. The SiPM uses the statistical recovery model,
  with PDE 0.32 and provisional 12.5 ns recovery in the current setup.
- The current matched run is bb0nu / Xe–He only. Earlier pure-Xe, electron
  background, 2/3 mm and other optical configurations remain separate runs.
  Historical handoff instructions describe their original tasks; consult the
  current manifest before selecting a production input or claiming completion.
- The CSA nonlinear development keeps the existing 0.63 V guard and 0.5 ENC
  validation target. New production defaults to `hybrid_nonlinear`: use the
  calibrated stateful model when the predicted minimum output falls BELOW the
  guard. Run cards select `mode = "production"` (default for new runs) or
  `mode = "validation"`. Production reuses the validated calibration and runs
  SPICE only for histories outside its domain. Validation additionally runs
  SPICE controls and solver-refinement histories. Keep the inexpensive
  provenance, charge, storage and noise consistency checks in both modes.
  `response_model = "hybrid_spice"` with validation mode selects the reference.
  Saved-study registration preserves its recorded mode. The calibration uses
  disjoint events for fitting and validation; 20/30/40 ms histories and a
  100-event 4 mm production pilot passed. See docs/ion_csa_nonlinear_model.html
  and latex_generated/ion_readout/ion_csa_nonlinear/index.html. Changes to the
  circuit or calibrated input-current prescription require revalidation.
- Production-viewer geometry uses source x/y/z in mm, inverting the saved
  electron bucket direction. Both arms start with shared ranges 25% larger
  than their combined occupied extent; zoom and rotation remain independent.
- Canonical LTspice locations, local simulator configuration and historical
  path resolution are documented in docs/ltspice_setup.md. Generated HDF5,
  PWL and compressed waveform CSV files under results/ltspice/ remain local
  and untracked. Logs and build auxiliaries are excluded. Portable release
  inputs under LTSpice/releases/ remain tracked. Do not delete local numerical
  histories when applying Git exclusions; they are calibration evidence.
- The current publication manuscript is manuscripts/itaca3/itaca3.tex (EPJC).
  The user also maintains the earlier notes/physics/nausica.tex draft. Preserve
  author edits and distinguish paper text from generated numerical reports.

## Related ITACA sources

- Keep `/Users/jjgomezcadenas/Papers/itaca/` as the local source directory for
  future reuse of manuscript material, figures and available data.
- Main manuscript: `itacaRevisited.tex`; its active sections are in `srcx/`.
- Shared notation: `pool/NextDefs.tex`; packages: `pool/NextLatexPackagesBase.tex`;
  bibliography: `pool/NextRefs.bib`. A selected, self-contained notation subset
  is imported into `notes/physics/itaca_defs.tex` and loaded by `nausica.tex`.
- Figures: `img/`, `det/`, `gxe/`, `trks/`, and `fig2/`; editable detector and
  MARS diagrams: `tikz/` and `mars/`. Check the active manuscript's references
  to identify the intended version of an asset.
- Consult this directory when figures or data are needed; verify provenance and
  whether underlying numerical data exist before treating a plot as a dataset.
  Copy selected assets into this repository when used so builds stay portable.
- Older ITACA operating values do not supersede the revised NAUSICAA design
  inputs in `config/detector.toml`.
