# NAUSICAA: ion signal physics

Reproducible calculations and a LaTeX note for positive molecular Xe2+ ions
collected by a focused CMOS pad plane. Each ion carries **+e**.

## First milestone

The implemented model covers detector geometry and the revised nominal fields,
constant-mobility transport, thermal diffusion, and an analytical single-ion
weighting-potential benchmark. It does **not** yet solve the actual IFG field,
simulate realistic event tracks, or evaluate the feedback transistor.

The revised design uses 200 V/cm in the bulk, 5 kV across a 10 mm IFG, and
3 kV across the final 5 mm gap. The pads are at 0 V, IC at +3 kV and gate at
+8 kV. These values supersede the older voltage scheme in the ITACA preprint.

## Reproduce

Use Python >=3.11 (the tested environment is recorded in results/physics_benchmarks.json):

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
make PYTHON=.venv/bin/python test physics readout
make PYTHON=.venv/bin/python note
```

The last command needs a TeX installation with `latexmk`, `pdflatex`, and
`bibtex`. It builds `notes/physics/build/nausica.pdf`. The LaTeX project can also
be opened in Codex at `notes/physics`; run the Python generation step first.
If dependencies are already available, set `PYTHON` to that interpreter.

- `config/detector.toml`: SI inputs, source labels, and unresolved assumptions.
- `src/ion_signal/physics.py`: reusable equations, current convention and models.
- `scripts/build_physics.py`: figures, LaTeX tables/macros and numerical outputs.
- `tests/test_physics.py`: independent surface integration, charge accounting,
  boundary values, drift/diffusion consistency, scaling and convergence.
- `notes/physics/nausica.tex`: the physics note and derivations.
- `notes/physics/references.bib`: primary sources.
- `notes/physics/generated/`: generated figures, tables and numerical macros.
- `results/physics_benchmarks.json`: numerical results, configuration hash and
  software versions.
- `docs/decisions.md`: provenance, limitations and next milestones.

Generated benchmark outputs are checked in so the note is inspectable without
first installing Python. Regenerate them after changing any input or equation.
The note's explanatory prose describes the baseline design; review it alongside
the generated tables whenever changing that design.
Build intermediates and scratch files are ignored. No random sampling is used
in milestone one. Future event simulations must record seeds and data provenance.

## Julia event-analysis core

The reviewed event types and coordinate translation live in
`src/event_analysis/`. Load the module in Julia with:

```julia
include("src/event_analysis/EventAnalysis.jl")
using .EventAnalysis
```

`EventInfo` holds event-level information. `EventHits` stores separate Float32
position, energy and MC-time vectors in mm, MeV and ns, plus integer identifiers.
`CoordinateFrame` labels source and vertex-relative positions; `center_hits`
subtracts the true vertex and returns independently owned arrays. Each public
interface has an extensive docstring with an executable example.

Run `make test-event-core` for these tests. The next design review starts with
the offline [Petit assessment](docs/event_analysis_review/index.html): active
implementation, existing types/interfaces, and the proposed adaptation. Earlier
interface pages are retained as superseded pending review.
Run `make event-review` with Python 3.11+ to rebuild the assessment from the
sibling `../Petit` checkout. Select Python with `PYTHON=/path/to/python3` if
necessary; the builder also accepts `--petit-root /path/to/Petit`.
Source excerpts and working-tree hashes are embedded for offline reading.
Open the HTML directly in a browser and refresh after rebuilding. The
[analysis plan](docs/event_analysis_plan.md) records the review status.

## Julia event hit tables

The event reader uses Julia and HDF5. Install the pinned Julia environment with
`make julia-env`, run `make test-events`, then run `make extract-hits` to convert
the two NEXUS samples in `events/`. These targets keep new Julia packages and
caches in `work/julia-depot` and also consult the user's existing depot.

The outputs are `events/0nubb_hitdf.h5` and
`events/electron_0nubb_energy_hitdf.h5`. Each contains column datasets under
`/hits`: `event_id` (Int64), `particle_id` and `hit_id` (Int32), and `x_mm`,
`y_mm`, `z_mm`, `energy_MeV` (Float32). Coordinates and deposited energies retain
the original mm and MeV units, precision and row order. All particles and hits
are retained; no topology cuts, ion production or transport are applied.
The original units and identifier conventions follow the
[NEXUS output format](https://github.com/next-exp/nexus/wiki/Output-format).

The constant volume label is stored once in `/metadata/label`; nonconstant
labels cause an error rather than silent information loss. Deposit time is
omitted, treating the event as instantaneous for later ion transport.
`/events` stores `event_id`, `first_row` (one-based), `n_hits`, and `energy_MeV`
(Float64 sum of the original deposits), plus `vertex_x_mm`, `vertex_y_mm`,
`vertex_z_mm` (Float32) and `n_primaries`. The vertex is copied from the common
initial position of primary particles in `/MC/particles`; missing or distinct
primary vertices are errors. Schema v2 adds these event-level fields without
changing any hit values. Use `read_event_metadata(path)` to obtain this table.
Track analysis shifts coordinates to this generation vertex, before any cuts;
it does not use the unfiltered energy barycenter. Events must occupy contiguous rows;
nonconsecutive event IDs are supported. `/metadata` records the original file
name and SHA-256, software/schema versions, reader and manifest hashes, units,
and the original configuration as ordered key/value columns with repeated
keys retained. Particle ancestry remains available in the original HDF5 via
`event_id` and `particle_id`.

Each numeric hit column uses lossless shuffle + deflate level 3 compression,
with 4096-row chunks; conversion reads 65536 rows at a time. A temporary file
under `work/` is compared against every original hit and the per-event energy
sums before being renamed to the destination. Original files are not modified.
Existing outputs require an explicit `--overwrite` flag.

For a different input or interactive Julia use:

```sh
export JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia"
julia --project=. scripts/extract_hits.jl events/example.next.h5
julia --project=.
```

```julia
include("src/ion_signal/event_io.jl")
using .EventIO

df = read_hitdf("events/0nubb_hitdf.h5"; event_id=0)
all_hits = read_hitdf("events/0nubb_hitdf.h5")
original_event = read_nexus_hits("events/0nubb.next.h5"; event_id=0)
```

The reduced reader selects an event using its row index; the original reader
must scan `/MC/hits`. Missing event IDs raise `KeyError`. Event IDs are local to
each source: add a source identifier when combining DataFrames from different
files. Hit IDs are local to a particle and volume, not globally unique.

## Ion-track analysis and generated report

`make ion-tracks PYTHON=/path/to/python3` runs Julia analysis, then Python
presentation and LaTeX compilation. `make ion-report PYTHON=/path/to/python3`
regenerates only PNGs, tables, macros and the report PDF from saved results.
Python needs the pinned `requirements.txt`; Julia uses `Project.toml` and
`Manifest.toml`. Run `make test-events test-ion-analysis` for the Julia checks.

Inputs are the v2 hit tables above. `config/ion_tracks.toml` specifies the
samples, a nominal mean pair-production energy of 21.9 eV, a vertex-relative
XY cut of +/-150 mm, a 20 ms frame, and the grid margin. Pitch and drift speed
come from `config/detector.toml`. The initial analysis computes continuous ion
expectations `1e6 * energy_MeV / w_eV`, without ionization fluctuations or
transport. Each event is shifted to the primary generation vertex before cuts;
there is no rotation or barycenter recentering. Both samples use the same grids.
The XY cut is inclusive; Z has no additional cut. Each bin is left-closed and
right-open, except that the final bin includes its upper edge.

`results/ion_tracks/<sample>/analysis.h5` contains:

- `/events`: original event energy and expected ions, mean ions/hit, generation
  vertices, retained hit counts/energy/ions, and discarded energy/ions.
- `/hits`: expected ions and a selection mask aligned with the input hit rows.
- `/bins`, `/mean`: common edges, relative longitudinal time edges, the XY
  matrix and Z vector, averaged over all original saved events including zeros.
- `/diagnostics`: complete uncut X/Y histograms, both hit counts and ion weights.
- `/example`: selected hits and projections for the configured example event.
- `/summary`, `/selection`, `/validation`, `/metadata`, `/parameters`: raw
  summaries, selected totals and losses, conservation checks and provenance.

No discarded charge is restored by renormalization. Matrices are `[x,y]` in
Julia and `[y,x]` in h5py; `axis_order_on_disk` records the latter explicitly.
The mean projections describe initial selected ionization, not collection at
physical pads or joint pad/frame occupancy. Relative time uses a uniform speed
and is not an absolute arrival-time prediction.

`scripts/plot_hit_analysis.py` reads only saved analysis HDF5 files. PNG figures
go beside each result; `latex_generated/ion_tracks/report.pdf` gathers the
calculation, raw summaries, retained totals, tail diagnostics, figures,
validation and source provenance. The central Z detail and the full Z vector
are both plotted. All macros and tables in that directory are generated;
the canonical physics bibliography supplies the cited reference. Build scratch
files stay in `work/ion_tracks_latex`.

## Interpretation

The square-pad benchmark assumes a 0.5 mm square surrounded by an infinite
grounded conductor, a central normal trajectory and uniform local speed.
It is an analytical reference, not a prediction of the full perforated IFG.
The finite final-gap charge is kept as q Delta(phi_w); it is not silently
normalized to e. Full-geometry fields, measured transport and event samples
are the next steps before realistic packet distributions and electronics.

## Analytical front-end model

`make readout` generates a self-contained interactive model at
`notes/physics/readout/analytical_model.html`, numerical results in
`results/readout/benchmarks.json`, and SPICE input-current files in
`results/readout/spice/`. It combines the finite-gap square-pad current with
specified one-pad arrival distributions and an ideal CSA with linear feedback.
The original calculator supplied by Rogelio is preserved under
`notes/physics/readout/`; its transistor and noise estimates are not used.
See `notes/physics/readout/README.md` for assumptions, controls, and SPICE handoff.

## Uniform-arrival sensitivity calculation

`make noise-study` generates the paper's mean response per ion for uniform
arrivals, with fixed input capacitance 5 pF adopted from Topmetal-S. The ion count
N remains free; the noise calculation solves for the minimum N with S/N > 5.
The simultaneous mid-frame pulse is a numerical cross-check. All numerical
values, the normalized current/voltage figure, and charge/noise tables are
produced by `scripts/build_arrival_cases.py` and `scripts/build_noise_study.py`.
Their inputs are `config/detector.toml` and `config/realistic_csa.toml`; results
and provenance are in `results/realistic_csa/arrival_cases_point_samples.json`
and `noise_point_samples.json`. Signal and noise use point voltage readings.

`make realistic-readout` regenerates the separate explicit diagnostic pulses
and SPICE comparison files, using the same circuit but the acquisition aperture
specified in the configuration. `make capacitance-check` retains the historical
0.25-versus-5 pF diagnostic; it is not used by the active manuscript.
