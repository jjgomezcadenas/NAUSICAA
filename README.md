# NAUSICAA / ITACA: detector and electronics simulation

LTspice installation and scientific files: [setup and layout](docs/ltspice_setup.md).

Generated notes: [browse the catalogue](latex_generated/index.html).
Directory layout and regeneration: [generated reports](docs/generated_reports.md).

The current workflow simulates shared ion/electron production, diffusion,
ion-pad and photon-SiPM collection, CSA response and statistical SiPM/electronics
response. The completed matched run contains 8,000 bb0nu events with 4 mm ion
and EL pitch, 40 ms ion frames and 250 ns electron buckets, for Xe–15% He.
See [matched production](docs/elm_matched_production.md) for its cards,
manifests and commands. Run `make response-viewer` for the paired display.

New CSA production defaults to the validated fast nonlinear response, with
SPICE only for domain fallback. Run cards select `mode = "production"` or
`mode = "validation"`; validation repeats SPICE controls. See the [model](docs/ion_csa_nonlinear_model.html)
and [validation report](latex_generated/ion_readout/ion_csa_nonlinear/index.html).
`make ion-csa-nonlinear-pilot` runs the 100-event integration case and report;
set `mode = "validation"` and `response_model = "hybrid_spice"` for the reference mode.

The publication manuscript is `manuscripts/itaca3/itaca3.tex`; build it with
`make -C manuscripts/itaca3`. Earlier analytical studies remain available below.
Positive molecular Xe2+ denotes Xe_2^+; each ion carries **+e**.

## Original analytical benchmark

This benchmark covers detector geometry and its nominal fields,
constant-mobility transport, thermal diffusion, and an analytical single-ion
weighting potential. Event transport and transistor-level response are separate
implemented workflows described in the current production documentation.

The revised design uses 200 V/cm in the bulk, 5 kV across a 10 mm IFG, and
3 kV across the final 5 mm gap. The pads are at 0 V, IC at +3 kV and gate at
+8 kV. These values supersede the older voltage scheme in the ITACA preprint.

## Reproduce the analytical benchmark

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
- `scripts/benchmarks/build_physics.py`: figures, LaTeX tables/macros and numerical outputs.
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
Build intermediates and scratch files are ignored. This analytical benchmark
uses no random sampling; the event simulations record seeds and data provenance.

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

Run `make test-event-core` for these tests. The offline
[Petit assessment](docs/event_analysis_review/index.html) records the design
review and implementation lineage, including historical interface proposals.
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
`y_mm`, `z_mm`, `energy_MeV`, `time_ns` (Float32). Coordinates and deposited energies retain
the original mm and MeV units, precision and row order. All particles and hits
are retained; no topology cuts, ion production or transport are applied.
The original units and identifier conventions follow the
[NEXUS output format](https://github.com/next-exp/nexus/wiki/Output-format).

The constant volume label is stored once in `/metadata/label`; nonconstant
labels cause an error. Schema v3 preserves the original MC deposit time in ns.
`/events` stores `event_id`, `first_row` (one-based), `n_hits`, and `energy_MeV`
(Float32 source-order sum of the original deposits), plus `vertex_x_mm`, `vertex_y_mm`,
`vertex_z_mm` (Float32) and `n_primaries`. The vertex is copied from the common
initial position of primary particles in `/MC/particles`; missing or distinct
primary vertices are errors. Schema v3 also stores actual primary particle IDs
and their per-event row ranges. Use `read_event_metadata(path)` for the event
table and `read_primary_ids(file; event_id=...)` on an open HDF5 file for primary IDs. Production
transport uses the stored source coordinates. Events must occupy contiguous rows;
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

Original NEXUS files are stored in `~/NAUSICAA/`; reduced `_hitdf.h5` files
remain in this project's `events/` directory. Conversion of archive inputs
writes to `events/`, and the viewer reuses those reduced tables. Set
`NAUSICAA_RAW_DIRECTORY` to use a different archive location.

For a different input or interactive Julia use:

```sh
export JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia"
julia --project=. scripts/events/extract_hits.jl ~/NAUSICAA/example.next.h5
julia --project=.
```

```julia
include("src/ion_signal/event_io.jl")
using .EventIO

df = read_hitdf("events/0nubb_hitdf.h5"; event_id=0)
all_hits = read_hitdf("events/0nubb_hitdf.h5")
original_event = read_nexus_hits(expanduser("~/NAUSICAA/0nubb.next.h5"); event_id=0)
```

The reduced reader selects an event using its row index; the original reader
must scan `/MC/hits`. Missing event IDs raise `KeyError`. Event IDs are local to
each source: add a source identifier when combining DataFrames from different
files. Hit IDs are local to a particle and volume, not globally unique.

## Earlier ion-track analysis and generated report

`make ion-tracks PYTHON=/path/to/python3` runs Julia analysis, then Python
presentation and LaTeX compilation. `make ion-report PYTHON=/path/to/python3`
regenerates only PNGs, tables, macros and the report PDF from saved results.
Python needs the pinned `requirements.txt`; Julia uses `Project.toml` and
`Manifest.toml`. Run `make test-events test-ion-analysis` for the Julia checks.

This study reads reduced hit tables. `config/ion_tracks.toml` specifies the
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

`scripts/events/plot_hit_analysis.py` reads only saved analysis HDF5 files. PNG figures
go beside each result; `latex_generated/ion_readout/ion_tracks/report.pdf` gathers the
calculation, raw summaries, retained totals, tail diagnostics, figures,
validation and source provenance. The central Z detail and the full Z vector
are both plotted. All macros and tables in that directory are generated;
the canonical physics bibliography supplies the cited reference. Build scratch
files stay in `work/ion_tracks_latex`.

## Paired production viewer

The independent repository-review assignment and current development boundary
are documented in [the review handoff](docs/reviews/independent_repository_review_handoff.md).

The paired production viewer is launched with `make response-viewer` and opens
at `http://127.0.0.1:8769/`. Its default card is
`config/specs/runs/elm_matched_response_viewer.toml`, written after the matched
4 mm ion and Xe--He production completes. It selects the completed CSA response
and SiPM shards explicitly. To inspect another production, pass
`RESPONSE_VIEWER_CARD=path/to/card.toml` to `make response-viewer`.

Event numbers are one-based positions in the file; each panel also shows the MC
event ID and complementary drift distances. Gas and ion-pitch choices come from
the selected production. Electron stages are photons, detected PE, SiPM charge
and electronics response; ion stages are input charge and corrected CSA charge.
Every event is checked for matching primary-pair identity and placement across
the arms. Clipping is shown alongside the electron signal. Each projection zooms
independently, the 3D plot rotates, and **Produce PNG** exports the displayed views.
Both arms use source x, y and z in mm. The saved grid reference, longitudinal
bucket width and drift direction map each bucket back to its spatial interval;
electron buckets are reversed into increasing z. Initial ranges span the combined
occupied cells with 25% extra width (12.5% on each side), using equal spatial
scales. Relative arrival time remains available when hovering over a 3D point.

## Interpretation

The square-pad benchmark assumes a 0.5 mm square surrounded by an infinite
grounded conductor, a central normal trajectory and uniform local speed.
It is an analytical reference, not a prediction of the full perforated IFG.
The finite final-gap charge is kept as q Delta(phi_w); it is not silently
normalized to e. The production CSA uses this documented weighting approximation
with simulated event histories and the calibrated transistor response.

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
produced by `scripts/benchmarks/build_arrival_cases.py` and `scripts/benchmarks/build_noise_study.py`.
Their inputs are `config/detector.toml` and `config/realistic_csa.toml`; results
and provenance are in `results/realistic_csa/arrival_cases_point_samples.json`
and `noise_point_samples.json`. Signal and noise use point voltage readings.

`make realistic-readout` regenerates the separate explicit diagnostic pulses
and SPICE comparison files, using the same circuit but the acquisition aperture
specified in the configuration. `make capacitance-check` retains the historical
0.25-versus-5 pF diagnostic; it is not used by the active manuscript.
