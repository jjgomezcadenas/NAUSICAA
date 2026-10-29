# NAUSICAA working conventions

This repository develops a reproducible physical input-current model for a
positive-ion CMOS detector. Establish the physics before implementing electronics.

- Xe2+ means the molecular ion Xe_2^+, with charge +e, not atomic Xe^{2+}.
- Extend event readout to both carrier populations: a single-electron or
  double-electron ionisation track produces a trail of electrons and positive
  ions. Work to 2026-10-28 simulated the positive ions; electron drift to the
  anode is now part of the simulation scope. Trace Petit's electron diffusion
  laws and pure-Xe/Xe–He coefficient choices before implementation; preserve
  carrier and gas-model identity in configurations, caches and outputs.
  Electron readout operates at 15 bar; the Xe–He prescription represents 15% He
  (user design inputs). Store these in config/electron_transport.toml and derive
  operating-pressure diffusion coefficients from Petit's pressure scaling.
- Sample one event location along the 1.5 m TPC: ion and electron drift
  lengths are complementary, L and L_max-L. Sample each diffused carrier
  directly into its physical grid cell and longitudinal bucket; save sparse
  readout histories, seeds and provenance. The electron EL grid has 5 mm
  diameter holes, 10 mm pitch and 5 mm depth, with E/P = 3 kV/cm/bar.
  Store entry/exit electrode potentials in TOML.
  Route electrons geometrically to holes. Compute direct VUV collection
  acceptance analytically for the assumed electron path and sample that
  fraction. Save only photoelectron histories per SiPM and time bucket;
  use 1 microsecond per mm and an initial 2 mm bucket (2 microseconds at this
  speed). SiPM PDE is TOML-configurable, initially 0.30 (user input).
  These are user design inputs;
  record the adopted gas-dependent EL yield and optical boundary conditions.
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
- The active readout replays source hits and deterministic diffusion seeds when
  pitch or bucket width changes; it does not save individual carrier coordinates.
  Histogram cell centres are valid for display and historical comparison, not
  physical pad/time assignment. Version new HDF5 readout schemas and retain old
  binned outputs for audit. Record source contents, cards, seeds and code hashes
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
- Keep scratch files in work/.
- Cite primary sources in notes/physics/references.bib. Record user-provided
  revisions as design inputs, not as published measurements.
- Keep large datasets and binary simulation files outside Git. Never configure
  or use Git LFS without the user's explicit approval.
- Do not change the existing remote or publish/push without a user request.

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
