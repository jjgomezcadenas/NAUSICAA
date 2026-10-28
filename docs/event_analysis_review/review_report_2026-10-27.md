# Review of the event-analysis design guide

Reviewed document: `docs/event_analysis_review/index.html` (built 2026-10-27 13:27 UTC).
Reference source: Petit checkout at `../Petit`, git HEAD `305b32e`.
Reviewer: Claude (Fable 5.1), 2026-10-27.

## 1. Scope and method

The guide describes a Julia re-implementation, inside `src/event_analysis`, of
the ITACA event pipeline found in Petit: hit reading, ion production, Monte Carlo
diffusion, voxelisation, graph construction, single-track selection, backbone
extraction, a truth comparator and per-pad, per-slice ion histories. The guide
has six design chapters and three appendices (Petit implementation, Petit
interfaces, decision record).

This review did the following.

- Read the full text of the guide (about 19,000 words).
- Recomputed every SHA-256 hash listed in the provenance section against the
  current Petit working tree and the NAUSICAA files. All 29 hashes match.
- Read the Petit functions that the guide's claims rest on: `voxelize_event`,
  `build_tracks_kdtree`, `build_tracks_knn`, `make_tracks`,
  `find_extremes_combined`, `walk_track_from_extremes`, `get_raw_path`,
  `compute_mc_path`, `diffuse_xyz_image_mc`, `transform_hits_df`,
  `sigma_t_ion_mm`, `get_sigma`, `get_voxel_size_and_distance`,
  `get_energy_threshold`, `select_events_itaca`, `is_single_track` and
  `compute_occupancy_analysis`.
- Read the implemented NAUSICAA types (`types.jl`, `coordinates.jl`), the
  reduced-HDF5 reader (`src/ion_signal/event_io.jl`), both configuration cards
  and the Makefile targets.
- Recomputed the thermal-width table and the reconstruction-binning rule.

## 2. Overall assessment

The guide is accurate about Petit. Every behavioural statement I checked
against the source is correct, including the subtle ones: the two different
longitudinal-width choices in `get_sigma` and `select_events_itaca`, the
`ldrift` unit mismatch between the `DiffusionParams` docstring and the batch
runner, the three different W-values in play (25 eV in diffusion, 22 eV in
occupancy, 21.9 eV in the NAUSICAA card), the stale assertions in Petit's test
file, the kNN defaults by entry point, and the fact that pad occupancy in Petit
only sees thresholded reconstruction bins.

The design decisions are sound and well separated from the draft contracts. The
charge ledger (expected, sampled, retained, outside, graph-retained,
below-threshold) is a genuine improvement over Petit. The pad-history input
change (full diffusion histogram, not the thresholded track) is the right call.
The "approved" and "draft" labels are applied consistently.

The guide is ready to serve as the review surface for the interface round. The
findings below are things I would settle or change before the green light.
They are ordered by how much they affect the physics result, then by how much
they affect the implementation.

## 3. Findings

### 3.1 Physics and algorithm design

**F1. The common drift distance per event is a poor approximation at the short
end of the range.** The design samples one L per event and uses one σ for all
hits, while the hits keep their z spread. A 2.5 MeV track in 15 bar xenon spans
several centimetres in z. At L = 10 cm a hit 5 cm closer to the pads than the
event reference has σ smaller by a factor √(5/10) ≈ 0.71, and a hit 5 cm
farther has σ larger by √(15/10) ≈ 1.22. The guide acknowledges the
approximation in one sentence. Because diffusion is already sampled per hit,
using L_i = L_event + (z_i − z_ref) per hit costs nothing extra and removes the
approximation. This also makes the longitudinal structure of the time histories
physically consistent with Δt = Δz/v. I recommend changing the approved width
decision to per-hit drift distances, with the common-L mode retained as the
Petit-comparison option. Whether z increases towards or away from the pads is
the same sign convention the guide already lists as open.

**F2. The 50 cm breakpoint in `get_voxel_size_and_distance` is a factor-2
discontinuity in reconstruction settings.** Recomputed with the detector card
(T = 293.15 K, E = 200 V/cm):

| L / cm | σ / mm | bin width b / mm | max edge d / mm | d / b |
|-------:|-------:|-----------------:|----------------:|------:|
| 10.0   | 0.503  | 1.51 | 3.02 | 2.0 |
| 49.9   | 1.123  | 3.37 | 6.74 | 2.0 |
| 50.1   | 1.125  | 1.69 | 2.53 | 1.5 |
| 100.0  | 1.589  | 2.38 | 3.58 | 1.5 |
| 150.0  | 1.947  | 2.92 | 4.38 | 1.5 |

Across 50 cm the bin width halves and the edge cutoff drops by 2.7×. The
connectivity rule also changes qualitatively: d = 2b reaches face, edge, corner
and second-face neighbours (32 cells), while d = 1.5b reaches only face and
edge neighbours (18 cells). This rule was tuned in Petit for electrons at 15 bar
with different widths. Inheriting it unchanged for ions will produce a step in
single-track acceptance versus L that is an artefact of the settings, not of
the physics. The guide plans to "check acceptance around the breakpoint". I
would go further: make a continuous rule (one scale factor, or a smooth
interpolation) the baseline, and keep the Petit two-branch rule as the
comparison. At minimum, the acceptance study should be presented versus L with
the breakpoint marked.

**F3. Rebinning the diffusion histogram into reconstruction bins causes
aliasing.** The diffusion histogram has 100 bins over the event extent plus
padding, so about 0.5–1.2 mm per bin for typical events. Reconstruction bins
are 1.5–3.4 mm. The ratio is not an integer, so some reconstruction cells
receive one diffusion-cell centre and neighbouring ones receive two or three.
The 10-ion threshold then acts on a modulated count. The guide's "three grids"
section does not mention this interaction. Two remedies, in order of
preference:

1. Histogram the sampled ion positions directly into the zero-anchored
   reconstruction grid of width b and, separately, directly into the pad/slice
   grid. Both grids are known before sampling (see F4). The intermediate
   histogram D, nbins, nsigma and the "outside" ion loss all disappear. Petit's
   D can be kept as an optional debug product or for parameter scans.
2. If D is kept, make the diffusion bin width an integer fraction of b and of
   the pad pitch, instead of fixing nbins = 100.

The decision record approves retaining D "as the initial reference". I am not
proposing to abandon the Petit comparison, only to make the direct histograms
the production path once parity with Petit has been shown.

**F4. The "selected track reference" for the pad grid is redundant.** A grid
anchored at the track barycentre and then shifted by a uniform offset in
[0, pitch) has exactly the same distribution as a grid anchored at the origin
and shifted by a uniform offset in [0, pitch). The reference therefore adds
nothing and creates a dependency of pad counting on reconstruction. Dropping it
means the pad phase can be sampled before diffusion, which is what F3's first
remedy needs. It also removes the "pad coverage" open contract: the pad
histogram simply covers the sampled range.

**F5. Threshold on integer ion counts, not on Float32 energy.** The diffusion
table carries an integer `ions` column, but voxelisation sums `energy` and
`make_tracks` compares energy ≥ threshold_keV with threshold = n_ions × W/1000.
In Float32, (6 + 4) × 21.9e-6 summed in two steps is not guaranteed to equal
10 × 21.9e-6, so exactly-at-threshold cells can flip. The guide lists this as
an open contract ("review any comparison tolerance"). The clean fix is for the
adapted `voxelize_event` to sum the integer `ions` column alongside energy and
for the selection branch to threshold on that integer. Exact equality is then
trivial and no tolerance is needed. Petit's energy threshold stays available for
the parity tests.

**F6. W-value wording.** The card says 21.9 eV; the decision record's example
says "For W = 22 eV this is 0.22 keV". With 21.9 eV it is 0.219 keV. Use the
card value in the example or say "about 0.22 keV".

### 3.2 Numerical and reproducibility

**F7. Exact edge-set parity between Float32 and Float64 kNN graphs will fail
for reasons that are not bugs.** On a regular lattice, a cell's 12 edge
neighbours are all at distance √2 b. With k = 10 the kNN query must pick 4 of
them, and NearestNeighbors breaks ties by whatever order the tree returns.
Float32 and Float64 bin centres round differently in the last bit, so the
chosen 4 will differ. The validation table says "compare edge sets". For kNN
on binned data, compare connected components and component energies instead,
and state in the plan that edge sets are compared only for the KDT radius
graph. The same tie issue reaches the backbone: the weighted Dijkstra and the
BFS double sweep will pick different equal-cost paths. The plan already says
"record genuine tie-dependent differences"; it should say that these are
expected, not exceptional.

**F8. `make_tracks` has a 1.0 mm default edge cutoff that silently disables kNN
connectivity on grids coarser than 1 mm.** `build_tracks_knn` defaults to
`max_distance = Inf`, but the `make_tracks` wrapper defaults to
`max_distance_mm = 1.0` for every method. The simple reconstruction script
avoids this by passing Inf explicitly. The adapted `make_tracks` draft
keeps `max_distance_mm` as a required keyword, which is correct. Record this
trap in the source map so that a parity test does not accidentally compare
against a Petit call with the default.

**F9. Diffusion histogram memory.** The dense array is `nbins³` Int64: 8 MB at
100 bins, 216 MB at 300 bins (the simple-reco setting the guide mentions as a
comparison). If D is retained, use Int32 counts and document the memory bound.
With F3's direct histograms the issue disappears.

**F10. RNG stream derivation is still undefined.** The guide proposes streams
derived from run seed, sample, event and stage, and says "specify the stable
derivation algorithm before implementation". Agreed. A concrete proposal to
close it: `Xoshiro(hash((seed, sample_id, event_id, stage)))` with the
stage as a fixed small integer, Julia's `hash` of a tuple of integers and a
string being stable within a Julia minor version. Record the Julia version
with the results, as the guide already requires.

### 3.3 Interfaces and persistence

**F11. Per-event HDF5 groups for axes and paths will not scale.**
`/samples/<id>/axes/<event_id>/` and `/paths/<event_id>/…` create several
groups per event. With tens of thousands of events this makes the file slow to
open and to iterate. The axes are uniform grids and are fully determined by
(offset_x, offset_y, pitch, first pad index, number of pads, z origin, slice
width, number of slices), which are scalars that belong in the event table.
Paths should be one columnar table per sample (event_id, path_order, x, y, z,
s, kind) rather than one group per event.

**F12. `compute_occupancy_analysis` draft signature versus F4.** The draft
takes explicit edge vectors. If the grid is defined by (offset, pitch, slice
width, z origin), pass those scalars and let the function derive the integer
pad indices with `floor`. This matches the zero-anchored convention of
`voxelize_event` and makes the "integer-pitch relabelling" test trivial.

**F13. Ownership of the walk result.** The draft proposes that endpoint rows
refer to the returned copied `path_voxels` rather than being views into
`track.voxels`. Agreed; this is the only Petit return that leaks a view. Make
the extremes plain integer indices into `path_voxels` (1 and end) rather than
rows; the rows are then derivable and the NamedTuple stays small.

**F14. `Tracks.diffusion` carries settings that do not belong to a track.**
Petit stores `DiffusionParams` inside each `Tracks` because its scripts had
nowhere else to put it. In the adapted design the per-event settings are saved
in the event record. Keeping the field is fine for parity, but the HDF5 writer
should not read it; the event record is the authority.

**F15. Reduced-HDF5 schema extension.** The current converter explicitly
omits `time` (`metadata/omitted_fields`) and stores only the `primary == 1`
vertex. The guide's increment 2 adds `time_ns` and actual primary IDs. The
double-beta sample has two primaries; the converter's `n_primaries` field
already counts them, so the IDs are available. Bump the schema string and have
`read_event` reject the old schema with a message that names the regeneration
command, as the draft already says.

### 3.4 Validation plan

**F16. Running Petit as a live reference is expensive.** Petit's
`Project.toml` depends on Plots, IJulia, Pluto, TestImages, GLM and others.
Loading `Petit.jl` for a parity test would drag all of that into CI. Generate
the reference outputs once, under Petit's environment, for the fixed fixtures
(straight, curved, two-arm, branched, separated components, exact-boundary
positions), store them as HDF5 in `tests/event_analysis/fixtures/` with the
Petit commit hash, and compare against those files. Re-generate only when the
pinned Petit commit changes.

**F17. Petit's own test file cannot serve as evidence without repair.** The
guide documents the three stale assertions (old keyword names, three-field
`Tracks`, 60,000 carriers). Agreed. Do not copy these tests; write the
adapted contract tests from the drafts and keep the Petit fixture geometry.

**F18. The energy-resolution criterion needs an operational form.** "Compare
energy errors against the 0.0125 MeV resolution scale" is the right scale, but
Float32 accumulation of 10⁵ ions × 2.19e-5 MeV has relative error of order
1e-7, far below it. State the acceptance criterion as a fraction of that
scale (for example 1e-3 × 12.5 keV) so that the test is meaningful rather than
trivially passed.

### 3.5 Documentation

**F19. Duplicated material.** The thermal-width table appears in chapter 1 and
in appendix A with the same numbers. The `get_sigma` excerpt appears twice in
appendix A. The `voxelize_event` excerpt appears in appendices A and B. The
builder could render each once and link.

**F20. Builder-derived numbers are labelled as such.** The thermal table says
it is computed by the documentation builder from the detector card. This is
good practice and consistent with the AGENTS.md rule of not transcribing
numbers by hand. Keep it.

**F21. Provenance is complete and current.** All listed hashes match, the
Petit HEAD is recorded, the MIT licence text is embedded, and the line counts
quoted (1,512 lines across the two graph files) are correct.

## 4. Verified claims

The following statements in the guide were checked against the Petit source
and are correct.

- `voxelize_event`: `floor(x / b)`, centre `(i + 0.5) b`, dictionary energy
  sums, dictionary-order rows, `electrons = round(E × ef)` with
  ef = 40,000 per MeV.
- `make_tracks`: keV → MeV conversion, methods KDT/kNN/kNN_mutual, descending
  energy sort, `≥` threshold.
- `build_tracks_knn`: `k_actual = min(k, N−1)`, undirected union of directed
  neighbour relations, `d ≤ max_distance` filter, single-vertex and
  all-isolated trivial returns.
- `find_extremes_combined`: topology first with 0.85 early exit, dense/sparse
  split at average degree 6, the edge-energy, energy, spatial, curvature and
  optional MST order, 1.1× and 1.2× length comparisons, sparse branch 0.8
  threshold.
- `walk_track_from_extremes`: extremes are `DataFrameRow` views; `path_voxels`
  is a copy.
- `compute_mc_path`: hard-coded primary IDs 1 and 2, `filter` then `sort!`
  on the copy, consecutive-bin grouping, energy-weighted group positions, arm
  reversal and concatenation.
- `diffuse_xyz_image_mc`: per-carrier `randn()` on the global RNG, extent
  plus nsigma × σ padding, `nbins` per axis, out-of-range samples dropped,
  energy = count / ef.
- `transform_hits_df`: drops time, label, particle_id, hit_id; adds integer
  carriers; no translation.
- `sigma_t_ion_mm`: √(2 (kB/e) T L / E) × 10 with L in cm and E in V/cm. The
  three tabulated values (0.503, 1.589, 1.947 mm) reproduce.
- `get_sigma("ion")`: σ_l = 0. `select_events_itaca`: σ_l = σ_t, marked
  "approximation!" in the source, with `lion = lmax − levt`.
- `get_voxel_size_and_distance`: breakpoint at 50 cm, factors (3, 2) and
  (1.5, 1.5), truth width 0.5 mm hard-coded.
- `get_energy_threshold`: 10 ions → 0.25 keV via 40 ions/keV.
- `is_single_track` in `filter_single_tracks.jl`: kNN, k = 10, 1 keV, CLI
  defaults 2 mm bins and 2 mm cutoff.
- `compute_occupancy_analysis`: grid centred on the energy-weighted
  barycentre, 2 mm z slices from the track minimum, W = 22 eV, half-open
  slices, per-slice matrices discarded, only summaries returned.
- Petit entry-point defaults: simple reco kNN with Inf cutoff; batch runner
  and single-track analysis KDT.
- Petit's `runtests.jl` includes only `test_voxels_tracks.jl`.

## 5. Recommended changes before the green light

Ordered by priority.

1. Decide on per-hit drift distance (F1). This is the one change that alters
   the approved physics model; it should be a deliberate decision.
2. Replace the two-branch binning rule with a continuous baseline, keeping the
   Petit rule for comparison (F2).
3. Decouple the pad grid from reconstruction (F4) and histogram sampled ions
   directly into reconstruction and pad grids (F3); retain D as an optional
   product for Petit parity.
4. Threshold on integer ions (F5).
5. Rewrite the graph-parity criterion as component parity for kNN (F7) and
   state the tie-breaking expectation.
6. Flatten the HDF5 layout to columnar tables with grid scalars in the event
   record (F11, F12).
7. Fix the RNG stream derivation (F10) and the Float32 acceptance criterion
   (F18) as concrete text in the contracts.
8. Store Petit reference outputs as fixtures rather than importing Petit
   (F16).
9. Editorial: W-value example (F6), duplicated excerpts (F19).

The remaining open contracts listed in chapter 5 (pad coverage, boundary
assignment, time origin and sign, ensemble axes, graph energy versus integer
count) are each addressed by one of the findings above: pad coverage by F4,
boundary assignment and ensemble axes by F12, graph energy by F5. Time origin
and sign remains a decision for the user and is tied to F1.
