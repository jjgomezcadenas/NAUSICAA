> Historical snapshot before documentation consolidation. The [current HTML design](../index.html) is authoritative for new work.

# Event tracks, ion transport and pad arrivals

## Status and working agreement

The review has restarted from Petit. Open the offline
[Petit assessment](../index.html) and review, in order:

1. Active implementation and its actual call paths.
2. Petit types, interfaces and existing test fixtures.
3. Proposed reuse, changes and compatibility validation.

Review these prerequisites before defining the next NAUSICAA types, signatures
and tests. The API sketches and implementation sequence below are historical
proposals, **superseded pending this review**. The earlier HTML is retained in
[previous_interfaces.html](../previous_interfaces.html).

Agreed during the adaptation review: express units in field/parameter names
and documentation, preserving positions in mm, deposited energy in MeV and MC
time in ns. Clarify Petit's drift-length convention explicitly: the ITACA
runner supplies cm, so an adapted name such as `drift_length_cm` records that
unit. This agreement concerns naming and documentation, preserving the values.

Input boundary approved: retain the existing Julia/HDF5 reader and reduced hit
files, add MC time and primary-particle IDs to the reduced schema, and keep
`EventHits` as the event input. Adapt its columns to Petit's DataFrame layout
(for example `x_mm` to `x`, `energy_MeV` to `energy`, `time_ns` to `time`),
preserving values and units.

Voxelisation and binning approved: adapt Petit's `voxelize_event` directly,
binning the supplied coordinates with `floor(coordinate / bin_width)` and
summing deposited energy per occupied bin. Return a DataFrame with event ID,
geometric bin centres and summed energy. Pass `bin_width_mm` as a parameter.
Use Float32 and validate numerical differences against Petit's Float64 results.
This supersedes the VoxelGrid/VoxelCloud type proposals. The path result,
pixel-count data flow and diffusion decisions are recorded as approved below.

Carrier conversion approved: retain deposited energy for binning and graph
construction. Calculate expected ion counts in the ion-production stage using
`mean_ions = energy_MeV * 1f6 / ionisation_energy_eV`, with the ionisation energy
supplied by configuration. Retain fractional expectations; any later integer
sampling is a separate explicit step.

Track representation approved: retain Petit's `Tracks` structure (voxel
DataFrame, graph, component membership and diffusion/binning settings) and
`make_tracks` returning tracks ordered by decreasing deposited energy. Use
Float32 numerical columns and document units. Preserve the index conventions:
local graph vertex `i` refers to row `i` of the track's voxel table; stored
component indices refer to the filtered parent table. Load graph settings from
TOML cards and pass them explicitly to functions.

Graph algorithms approved: retain Petit's KDT radius graph, kNN and mutual-kNN
implementations, with method and parameters selected through TOML. The baseline
is now kNN single-track selection with a finite edge-distance cutoff, following
`filter_single_tracks.jl`; this supersedes the earlier KDT baseline choice.
Checked reconstruction defaults:
`itacaScripts/itaca_simple_reco.jl` uses kNN (k = 10, default distance cutoff Inf);
`scripts/itaca_track_reco_mt.jl` and `scripts/itaca_single_track_analysis.jl`
default to KDT and allow kNN variants for RECO. The approved NAUSICAA
reconstructed-data baseline is kNN with k = 10 and a finite distance cutoff.

Backbone extraction approved: retain Petit's combined endpoint/path algorithm,
including its existing branch conditions and energy weighting, with MST as a
separate comparison. Validate whether the Float32 adaptation preserves endpoint
and path decisions on shared fixtures.

MC truth comparator approved: retain Petit's primary-hit time ordering,
consecutive-bin grouping with energy-weighted positions, and reversal/joining
of the two primary sequences for double beta. Obtain primary IDs from
`EventInfo`. Use this trajectory as a truth reference for assessing the
graph-derived backbone of the diffused event, expressing both paths in a common
coordinate frame for comparison.

Whole-event selection approved: follow the reconstructed-event branch in
`scripts/itaca_track_reco_mt.jl`, selecting kNN. Set drift/diffusion parameters,
convert energy to carriers, diffuse, voxelise, build tracks, require exactly
one component, then extract the backbone. The separate raw-MC
`itacaScripts/filter_single_tracks.jl:is_single_track` illustrates the kNN
criterion; its pre-diffusion ordering is not the adopted pipeline. Bin the
diffused distribution, apply the bin-energy threshold, build the kNN graph
with a finite maximum edge distance, and accept exactly one connected component.
Extract the backbone from that accepted component. Keep bin width, k, edge
cutoff and energy threshold in TOML. The raw-MC filter uses k = 10 and a 1 keV
bin threshold, with CLI defaults of 2 mm bins and a 2 mm edge cutoff. Those
are prefilter settings; review diffused-event settings alongside diffusion.

Checked event preparation: Petit fixes or uniformly samples a per-event drift
distance, computes common diffusion widths from it and retains the input hit
coordinates. `transform_hits_df` changes columns and adds carrier counts;
it does not translate coordinates. The approved transverse alignment keeps
event coordinates unchanged and shifts the pad grid as described below.

Checked drift bounds: the batch runner accepts configurable `lmin`/`lmax` in
cm and samples L directly (or uses fixed L when they are equal). Our requested
range is obtained with 10 and 150 cm using the same sampler; the detector card
already records a 1.5 m bulk reference. A second routine,
`select_events_itaca`, samples `levt` between `lmin + lbuff` and `lmax - lbuff`,
then uses `lion = lmax - levt` and `lele = levt`. With its default lmin = 0,
lmax = 200 and lbuff = 10 (cm as consumed by the diffusion formulas), lmax acts
as a detector length. These values are keyword arguments, not fixed geometry.
The proposed drift adaptation is therefore configuration of the existing sampler.

Checked pad occupancy: `itacaScripts/itaca_simple_reco.jl` calls
`compute_occupancy_analysis` for exactly-one-track events, after optional
diffusion, voxelisation and graph construction. It centres a pixel grid on the
energy-weighted track barycentre and slices z in 2 mm bins from the track minimum.
Each voxel contributes `energy_MeV * 1e6 / 22` ions to its pixel/slice. It returns
per-pixel totals and slice summaries; each full slice matrix is temporary.
Pixel width is taken from the graph maximum distance (voxel width if infinite).
This is geometric square-cell assignment. The archived `scriptBkup/mipa.jl`
field-line prototype is separate from the event occupancy path.

Pad histories approved: retain every pad's ion counts by longitudinal slice,
with event, pad and slice indices and axis information in HDF5. Occupied entries
can represent the full sequences; the saved axes define zero-count slices.
Calculate totals, averages and extrema from these histories. The sequences
provide the charge distribution needed for subsequent pad-signal calculations.

Pad pitch approved: configure the physical pad pitch (2 mm for this study)
independently of the graph's maximum connection distance in TOML.

Time bins approved: use `Δt = Δz / v_ion` with the configured ion drift velocity.
Save both longitudinal bin edges and corresponding time-bin edges in the output.
For example, a 2 mm slice corresponds to 20 ms at 0.1 mm/ms.

Common W-value approved: use `ionisation_energy_eV` for both diffusion carrier
conversion and pad counting. Petit's diffusion default is 40,000 carriers/MeV
(25 eV), whereas its occupancy calculation uses 22 eV; the adaptation uses the
same configured value in both stages. Validate count accounting and final-bin
inclusion.

Transverse grid offset approved: define the pixel-grid edges around the selected
track and shift them by independent uniform x and y offsets in `[0, pitch_mm)`.
The pitch is 2 mm for this study. Keep event coordinates unchanged and use the
same offset pair for every longitudinal slice of an event. Save the offsets
with the event results. An integer-pitch shift only changes pad labels.
Event preparation uses the existing drift-distance sampler; diffusion widths
are approved below; the sampling/binning procedure is the next review item.

Diffusion widths checked: `sigma_t_ion_mm` implements the one-coordinate thermal
width `sigma² = 2(kB/e)TL/E`. The batch runner's `get_sigma("ion", L)` sets
longitudinal width to zero; the loaded `select_events_itaca` instead sets
`sigma_l_ion = sigma_t_ion`. Both use one drift distance for the event. With
the detector card's T = 293.15 K and uniform bulk E = 200 V/cm, thermal widths
are 0.503, 1.589 and 1.947 mm at L = 10, 100 and 150 cm. These are bulk,
singly charged ion estimates. Under isotropy and constant v = 0.1 mm/ms, the
corresponding time widths are 5.03, 15.89 and 19.47 ms.

Diffusion widths approved: retain Petit's formula using configured T and E and
one sampled L per event; use equal x, y and z widths, following
`select_events_itaca` and the detector card's thermal isotropic approximation.
Keep zero longitudinal width as a configurable comparison. Validate the formula,
L/T/E scaling, zero-drift limit, unit conversion and longitudinal settings.
Sampling approved: retain Petit's per-ion Gaussian Monte Carlo
as the initial reference, with an explicit RNG and the approved common W-value.
Benchmark runtime and memory on representative events before selecting an
alternative. Petit rounds expected counts per hit to obtain integer carriers;
this conversion is approved below. Review histogram boundaries explicitly. Preserve fractional
expected counts as the charge reference. Validate displacement means, variances,
cross-axis covariances, seeded reproducibility and generated/histogram counts.

Integer carrier conversion approved: retain Petit's per-hit
`round(Int, mean_ions)` for Monte Carlo sampling, using the configured W-value.
Preserve fractional expectations and record the event-level difference between
summed integer counts and summed expectations. Evaluate the rounding effect on
the input samples.

Intermediate histogram approved: retain Petit's dense 3D counting
histogram over event hit bounds padded by `nsigma * sigma` on each axis, with
`nbins` bins per axis. Return occupied bin centres/counts/energy for subsequent
graph voxelisation. Configure nbins and nsigma in TOML; the batch defaults are
100 and 3, while simple reconstruction uses 300 and 3. Start comparisons at the
batch defaults, benchmark runtime and memory, and check track-selection and
pad-history convergence as bin count and padding increase. Record generated,
retained and out-of-range counts. Final settings follow those results.

Reconstruction binning approved: retain the batch helper
`get_voxel_size_and_distance` as the initial reference. For L ≤ 50 cm, use
bin width b = 3 sigma_t and maximum graph edge distance d = 2b; for L > 50 cm,
use b = 1.5 sigma_t and d = 1.5b. Apply this to the approved kNN reconstruction.
Read the breakpoint and factors from TOML and save evaluated b/d per event.
Pad pitch remains independent at 2 mm. Check acceptance on both sides of the
50 cm scale change as part of graph/binning sensitivity studies.

Graph-input threshold approved: retain the batch runner's starting
threshold of 10 ions per reconstruction bin, configured in TOML as
`[graph] min_ions_per_bin = 10`. Convert to
keV with `threshold_ions * ionisation_energy_eV / 1000`, using the approved
common W-value. Petit uses 40 ions/keV (10 ions = 0.25 keV); at W = 22 eV the
consistent threshold is 0.22 keV. Retain bins at or above the threshold before
kNN construction. Measure removed charge and single-track acceptance versus
threshold.

kNN neighbour count approved: retain k = 10, the default in both
the batch and simple reconstruction scripts, configurable in TOML. Each bin
considers up to k nearest other bins, with edges limited by the approved maximum
distance. Preserve Petit's reduction to min(k, N − 1) for small inputs and its
undirected edge rule (either endpoint can select the other). Assess acceptance
versus k alongside the other graph parameters.

Backbone path representation approved: retain Petit's ordered local
vertex indices (rows in the selected track's voxel table) and path DataFrame
with x, y, z and cumulative distance s in mm, using Float32 numerical columns.
Keep the walk result's total length and algorithm confidence score. Validate
index/row correspondence, cumulative distances, total length and single-bin/
empty-path handling. Exact signatures and ownership contracts follow in the
interface review.

Pad-count input approved: Petit fills occupancy from `track.voxels`,
using reconstruction-bin centres and charge surviving the graph threshold.
For physical pad histories, select events with the approved graph
procedure, then filling pad/slice counts from each accepted event's full
intermediate diffusion histogram, before reconstruction rebinning and the
graph threshold. Use the selected track for the grid reference and topology
diagnostics. This preserves subthreshold charge and the finer histogram
positions for pad assignment; finite diffusion-histogram losses remain in the
charge accounting. The TOML graph threshold controls reconstruction; pad
histories include the accepted event's subthreshold charge.

Reference validation proposed for review: compare each adapted stage with loaded
Petit functions using identical inputs/settings and fixed small fixtures with
known geometry (straight, curved, two-arm and separated components). Compare bin
energies, graph edges/components and backbone indices/lengths, accounting for row
order and equivalent path reversal. Use identical diffused tables for graph
comparisons. Test diffusion with seeded reproducibility and Gaussian/count
statistics, and quantify Float32 differences against Float64. Validate approved
changes with their physical expectations and charge accounting. Then evaluate
acceptance, pad histories and performance on representative events from both
HDF5 samples. Exact function contracts and tests follow this prerequisite review.

Retain `diagnose_skeleton_coverage` (energy fraction near discrete backbone
vertices) and `compute_track_distances` (main-to-satellite closest bin-pair
distances) as diagnostics. The earlier maximum original-hit distance to
continuous backbone segments is deferred. Validate against Petit's existing
single-track filter and report acceptance versus graph/binning settings.

Record each design decision in the HTML for final review. Implementation awaits
the user's explicit green light after that review.

`CoordinateFrame`, `EventInfo`, `EventHits` and `center_hits` remain implemented
in `src/event_analysis/`, with tests in `tests/event_analysis/`. Their integration
with Petit is assessed in the new guide. The corrected pipeline order is:
set drift/diffusion parameters, diffusion, voxelisation, kNN graph, exactly-one-component
selection, then backbone extraction. Diffusion changes connectivity; the
single-track condition therefore acts on the diffused event. The undiffused MC
path provides a truth reference. The approved binned output is a DataFrame.
The transverse grid offset, diffusion method and graph settings are approved;
the backbone path representation and full-diffusion pad-count input are approved.

Run `make event-review` with Python 3.11+ to regenerate the assessment from the
sibling Petit checkout. The builder accepts `--petit-root` for another path.
`scripts/petit_review_content.py` contains the assessment and
`scripts/build_event_review.py` extracts source excerpts and provenance;
`docs/event_analysis_review/template.html` provides the presentation. This is
a source assessment; execution of Petit and numerical compatibility checks are
proposed work. `make test-event-core` separately runs the existing group-one tests.

This document records the agreed scope and the proposed implementation sequence.
Before each increment, we will review signatures, input/output contracts, physical
assumptions, failure behaviour and tests together.

New analysis code will be Julia, with physical calculations expressed as pure
functions wherever possible. Functions will have extensive docstrings. Julia
will save numerical results to HDF5; Python will produce final figures, tables
and LaTeX macros from those saved results. A generated, self-contained LaTeX
report will compile to PDF. Julia supplies the numerical physics results.

## Objective

Starting from the NEXUS samples, construct a backbone from the undiffused MC
hits, select events by their compatibility with that backbone, and calculate
the ion charge arriving at each pad as a function of time. The pad sequences
will subsequently feed the existing signal and electronics model and track
reconstruction.

The first deliverable stops at MC backbones, per-hit distances and an
event-selection scan. Diffusion and pad arrivals follow after review of that
deliverable.

## 1. Event data and coordinate conventions

- Preserve positions in mm, deposited energy in MeV and the source Float32
  precision in the hit tables and analysis hit types. Float32 is the proposed
  baseline for numerical calculations and accumulations; validate accumulation
  error against Float64 references before finalising the reduction routines.
  The relevant energy-resolution scale is 0.5% at 2.5 MeV, or 12.5 keV.
- Add the original MC `time_ns` to the reduced hit tables, retaining `event_id`,
  `particle_id`, `hit_id`, `x_mm`, `y_mm`, `z_mm` and `energy_MeV`.
- Store event-level information in `EventInfo`: the event ID, generation vertex
  and primary particle IDs read from the MC particle table. Electron-background
  events have one primary electron; double-beta events have two.
- Preserve hit identity and row order in per-hit results, including backbone distances.
- Centre each event by subtracting its true generation vertex, preserving its
  orientation and relative geometry.
- Distinguish MC deposition time, used to order truth histories, from ion
  arrival time, calculated from transport.
- Version the reduced HDF5 schema and define the handling of older files when
  reviewing the reader changes. Original NEXUS files remain the source of truth.

## 2. MC graph, backbone and event selection

### Graph construction

Follow Petit's existing MC procedure initially:

1. Voxelise all undiffused event hits, summing energy in each voxel.
2. Retain graph nodes above a configurable voxel-energy threshold.
3. Connect nodes within a configurable Euclidean radius using a spatial index.
4. Find connected components and select the component with the greatest total
   energy for backbone extraction. Retain the other components for diagnostics.

Starting settings are 0.5 mm voxels, a 1.25 mm connection radius and a 1 keV
voxel-energy threshold, configured through TOML. The subsequent rejection test
uses every original hit. Components remain separate during path extraction.
Specify deterministic node ordering, component tie-breaking and path tie-breaking
when reviewing the interfaces.

### Backbone definition

Adapt Petit's existing topology/energy-aware endpoint and path method as the
baseline. Its ordered voxel centres define a piecewise-linear 3D backbone.

For a double-beta event the intended path joins the two electron stopping
regions through the common-origin region. For a single electron it follows the
principal trajectory from its start to its stopping region. Time-ordered primary
histories provide an independent diagnostic of the graph-derived backbone.

Keep the MST-based method as an explicit comparison. Record the actual path method and relevant
diagnostics. Represent inability to construct a backbone with an explicit status,
separate from rejection by the distance cut.

### Hit distances and selection

For every original hit, calculate the minimum Euclidean distance to the line
segments of the backbone, including their endpoints. Treat a zero-length
segment as a point; define empty and single-node backbone handling explicitly.

For a segment from `a` to `b` and hit position `h`:

```text
t = clamp(dot(h-a, b-a) / dot(b-a, b-a), 0, 1)
d_segment = norm(h - (a + t*(b-a)))
d_hit = minimum(d_segment over backbone segments)
D_max = maximum(d_hit over all original hits)
```

Accept an event exactly when a valid backbone exists and
`D_max <= max_hit_distance_mm`. A single hit outside this neighbourhood rejects
the whole event. Accepted events keep all their hits and energy. This supersedes
the provisional rectangular XY hit-removal cut. Construct the graph and evaluate
distances from the original event.

Start with 20 mm and scan 10, 20 and 30 mm using the same saved distances. Save
per-hit distances, the maximum distance and responsible hit, component summaries,
backbone geometry, acceptance and failure status. Examine acceptance and the
stability of distributions for both input samples. Report unconditional
acceptance and conditional averages over accepted events separately.

## 3. Placement, detector and transport

Consolidate the existing detector inputs.
The ITACA fiducial dimensions are 3.2 m diameter and 1.5 m height; the active
local study uses a periodic grid without radial placement or fiducial cuts.
Use the revised fields and local geometry from `config/detector.toml`.

For each accepted vertex-centred event:

- Sample the vertex drift distance uniformly over 100–1500 mm.
- Define local positive z away from the grid, giving hit drift distance
  `L_i = L_vertex + z_i` and the corresponding arrival-time structure.
- Sample independent transverse grid phases uniformly over `[0, pitch_mm)`.
  An integer-pitch translation only relabels pads, so larger offsets are
  redundant. The nominal pitch is 2 mm.
- Preserve event orientation.
- Use reproducible random placement, independent of processing order, selection
  threshold and thread count. Use the same placements for comparative studies.

Review the convention for hits whose translated drift distance is negative
before implementation. The 1500 mm endpoint specifies the vertex-placement range.

Use the existing thermal Gaussian diffusion model. A common width evaluated at
the vertex distance is the existing simple approximation; evaluating widths at
individual hit distances is also inexpensive. Scope the baseline choice and
longitudinal-diffusion treatment explicitly before implementing transport.

For speed, integrate Gaussian probabilities over nearby collection cells to
obtain expected ion counts. Retain a small
individual-ion Monte Carlo comparison for validation. Define numerical tail
tolerances and account explicitly for any omitted charge.

Assume ideal focusing and collection: each collection cell feeds its corresponding
hole and pad with unit collection efficiency in this approximation. Use the
existing hole/gap transport and signal model at the later response interface.

Use one explicit W value throughout ionisation and transport, initially the
existing 21.9 eV configuration value. Expected charge may be fractional;
individual-ion and ionisation fluctuations are a separate later choice.

## 4. Output and downstream interface

The pad-arrival output is a sparse table with an explicit schema, conceptually:

```text
event_id, pad_ix, pad_iy, time_bin, expected_ions
```

Save grid origin/phase, pitch, bin edges, drift placement and time-origin
conventions that fully specify the positions and times.
Each pad has an ordered time sequence. Time bin width is configurable; the
initial proposed 20 ms corresponds to 2 mm at the assumed 0.1 mm/ms ion speed.
Retain a sufficiently clear bin-integrated interface for later current and
electronics calculations. Counts represent arrivals integrated over each bin.

Store results under `results/` and generated presentation under
`latex_generated/`. Save configuration, source identities/hashes, schema and
algorithm versions, random-placement information, selection accounting and
validation summaries. The presentation step must be rerunnable without repeating
graph construction or transport. Compile LaTeX and inspect its build log; visual
PDF inspection is only performed if requested.

## 5. Julia software design

### Project location and module boundary

Event analysis belongs in `src/event_analysis/`, separate from signal formation
and electronics in `src/ion_signal/`. Use one Julia module, `EventAnalysis`, with
one public entry file and shared types. Internal files organise implementation.

Proposed layout (files are introduced as their milestones are implemented):

```text
docs/event_analysis_plan.md
src/event_analysis/
    EventAnalysis.jl       # module, imports, includes and public exports
    types.jl               # reviewed data/result/config types and validation
    coordinates.jl         # vertex subtraction and coordinate conventions
    voxelization.jl        # 3D hit binning and energy aggregation
    graphs.jl              # neighbourhood graph and connected components
    backbones.jl           # adapted Petit endpoint/path algorithms
    distances.jl           # point/segment and hit/backbone distances
    selection.jl           # event classification and cut scans
    placement.jl           # explicit sampling boundary; deterministic geometry
    diffusion.jl           # expected Gaussian charge sharing
    arrivals.jl            # pad/time distributions
    io.jl                  # NEXUS/reduced HDF5 adapters and result writers
    config.jl              # TOML parsing into validated typed settings
    pipeline.jl            # effectful orchestration of the computational core
    projections.jl         # existing energy/ion summaries and projections
tests/event_analysis/
    runtests.jl
    test_types.jl
    test_io.jl
    test_coordinates.jl
    test_voxelization.jl
    test_graphs.jl
    test_backbones.jl
    test_distances.jl
    test_selection.jl
    # placement/transport tests added with their milestone
scripts/                   # thin existing/new Julia runners; Python reporting
config/event_analysis.toml  # graph/selection/placement settings
config/detector.toml        # shared physical detector inputs
results/event_analysis/
latex_generated/event_analysis/
```

The directory reorganisation is a separate, behaviour-preserving increment.
Migrate the existing `src/ion_signal/event_io.jl` and event-analysis functions
into this module, separating their I/O and pipeline responsibilities. Update
script/test includes, Make targets and README examples together. Consolidate the
implementation in the new `EventAnalysis` module. Preserve
existing entry-point behaviour and tests during migration; introduce the new
event-level selection explicitly in the subsequent increment. Existing saved
results and cards remain readable/reproducible as the legacy analysis, not
inputs silently interpreted with the new selection policy.

Leave Python signal/electronics modules and `averaging_scan.jl` in `ion_signal`.
Update AGENTS.md's source-location convention when implementing the move, so
that it distinguishes event analysis from signal formation.

### Root Project.toml and dependency scope

Keep one root Julia environment and adapt the required Petit functions locally.
Retain the existing direct module-inclusion loading model for this increment:
scripts include `src/event_analysis/EventAnalysis.jl` once and use
`.EventAnalysis`.

The existing Project.toml has DataFrames, HDF5, Dates, SHA, TOML and Test.
Group one needs Base and these existing standard-library dependencies for tests
and their validation record, so Project.toml and Manifest.toml remain unchanged.
Review and add further dependencies at their corresponding milestone:

| Dependency | Scope |
| --- | --- |
| DataFrames, HDF5 | Retain for table/storage boundaries |
| Dates, SHA, TOML | Retain for provenance and configuration |
| Test | Retain in the environment for direct test-runner scripts |
| Graphs | Add for graph storage, components and path algorithms |
| NearestNeighbors | Add for radius-neighbour graph construction |
| LinearAlgebra | Declare directly; geometry and existing averaging analysis use it |
| SparseArrays | Declare for the adapted Petit sparse edge-weight implementation |
| Random | Declare for the later explicit placement/Monte Carlo boundary |
| SpecialFunctions | Proposed addition only at diffusion milestone for Gaussian CDF integration |

Add compatibility bounds for new external packages using versions resolved and
validated. Regenerate Manifest.toml through Julia's package manager and check the
existing Julia scripts alongside the new tests.

The current project declares Julia 1.11 compatibility, while its manifest was
generated with Julia 1.12.6. Record the version used for validation and review
version support explicitly.

Project.toml and source changes belong to the implementation increment.

### Functional core and explicit effects

- Computational functions derive their results from typed data and explicit
  parameters, preserving caller-owned inputs. I/O, plotting and progress reporting
  belong to the application layer.
- Parse TOML and validate configuration at the application boundary, then pass
  typed parameter objects into the computational core.
- Keep DataFrames at table/I/O boundaries where convenient. Review lightweight
  typed structures or arrays for geometry and graph computation, agreeing their
  contracts before implementation.
- Separate random placement generation from deterministic transport. Supply an
  explicit RNG to the sampling boundary and pass the sampled placement into
  pure functions. Document that RNG consumption is an effect.
- Keep HDF5 readers/writers, command-line orchestration and presentation separate.
- Distinguish invalid input errors from expected analysis outcomes such as no
  usable backbone, using explicit result statuses.
- Adapt only the required Petit algorithms, record their provenance and preserve
  applicable licence notices. Verify behavioural compatibility on agreed fixtures.

### Proposed concrete types for the first milestone

These declarations propose concrete element types and named result structures
for joint review.

```julia
const Point3MM = NTuple{3,Float32}
@enum CoordinateFrame source_frame vertex_frame
@enum BackboneMethod petit_combined mst_comparison
@enum BackboneStatus backbone_ok no_graph_nodes no_usable_path
@enum SelectionStatus event_accepted distance_rejected backbone_failed

struct EventInfo
    event_id::Int64
    vertex_mm::NTuple{3,Float32}
    primary_particle_ids::Vector{Int32} # one for electron background, two for double beta
end

struct EventHits
    event_id::Int64
    frame::CoordinateFrame
    x_mm::Vector{Float32}            # N; x coordinates
    y_mm::Vector{Float32}            # N; y coordinates
    z_mm::Vector{Float32}            # N; z coordinates
    energy_MeV::Vector{Float32}       # N; original values
    time_ns::Vector{Float32}          # N; original MC time
    particle_id::Vector{Int32}        # N
    hit_id::Vector{Int32}             # N; identified within its particle
end

struct VoxelGrid
    size_mm::Float32
    origin_mm::NTuple{3,Float32}
end

struct VoxelCloud
    event_id::Int64
    frame::CoordinateFrame
    grid::VoxelGrid
    x_mm::Vector{Float32}            # V; geometric cell-center coordinates
    y_mm::Vector{Float32}
    z_mm::Vector{Float32}
    energy_MeV::Vector{Float32}       # V; sums, before graph threshold
end

struct GraphSpec
    connection_radius_mm::Float32
    min_voxel_energy_keV::Float32
end

struct VoxelGraph
    graph::SimpleGraph{Int}
    node_to_voxel::Vector{Int}        # graph node → VoxelCloud index
    components::Vector{Vector{Int}}   # graph node IDs; decreasing component energy
    component_energy_MeV::Vector{Float32}
end

struct BackboneSpec
    method::BackboneMethod
    use_energy_weighting::Bool
    use_edge_energy_weighting::Bool
    dense_degree_threshold::Float32
end

struct Backbone
    event_id::Int64
    points_mm::Matrix{Float32}        # 3 × K, ordered polyline
    voxel_indices::Vector{Int}        # K, indices in VoxelCloud
    component_rank::Int              # baseline: 1, highest-energy component
    method::BackboneMethod
    path_variant::Symbol             # validated label for actual Petit branch
end

struct BackboneResult
    status::BackboneStatus
    backbone::Union{Nothing,Backbone}
end

struct HitDistances
    event_id::Int64
    distance_mm::Vector{Float32}      # N; original hit order
    particle_id::Vector{Int32}        # N; same order as EventHits
    hit_id::Vector{Int32}             # N; same order as EventHits
    farthest_index::Int              # index into those N-element vectors
end

struct EventSelection
    event_id::Int64
    status::SelectionStatus
    cut_mm::Float32
    max_distance_mm::Union{Nothing,Float32}
    farthest_hit::Union{Nothing,Tuple{Int32,Int32}} # (particle_id, hit_id)
end
```

#### Group two: voxelisation proposal

Review `VoxelGrid`, then `VoxelCloud`, then `voxelize_hits`. The proposed cloud
uses separate coordinate vectors, carries its coordinate frame and stores the
grid specification used to produce it. This extends the initial cloud sketch
with the conventions established during group-one review.

- `size_mm` is a finite positive cubic cell size. `origin_mm` is a finite
  `(x, y, z)` grid corner in the input hit frame. The analysis cell size is
  chosen in the backbone stability study, independently of sensor pitch.
- Define integer cell keys by `floor(Int, (coordinate - origin) / size)` on
  each axis. Cells include their lower boundary and exclude their upper one.
  Review numerical boundary behaviour using exact edges and adjacent Float32
  inputs, including negative coordinates and translated grids.
- Store occupied cells in ascending lexicographic `(ix, iy, iz)` order. Their
  positions are geometric cell centres, `origin + (index + 0.5) * size`,
  represented as Float32. Each occupied cell contains at least one hit, including
  cells whose deposits sum to zero. Every hit contributes once to its cell.
- Keep the original `EventHits` alongside `VoxelCloud` for the later
  hit-to-backbone distance calculation.
- Proposed energy reduction: accumulate in Float32, ordering hits within each
  cell by `(particle_id, hit_id)` for reproducibility under input permutation.
  Compare with independent Float64 sums and report discrepancies in MeV and
  relative to the 0.0125 MeV resolution scale before accepting the reduction.
- Preserve the input event ID and coordinate frame; the pipeline supplies
  vertex-centred hits. Returned arrays are independently owned. The later graph
  stage applies its energy threshold to the voxel energies.

These contracts and their tests remain proposals until jointly reviewed.

#### CoordinateFrame enum

The declaration `@enum CoordinateFrame source_frame vertex_frame` creates the
type `CoordinateFrame` and two named values:

| Value | Coordinate convention |
| --- | --- |
| `source_frame` | Original simulation coordinates, as read from the NEXUS hit table. |
| `vertex_frame` | Coordinates relative to the true generation vertex: source position minus vertex position. The vertex is at (0, 0, 0). |

The enum belongs in `src/event_analysis/types.jl`. The `EventHits` field
`frame::CoordinateFrame` holds one of these values and specifies the coordinate
system of `x_mm`, `y_mm` and `z_mm`. The position, energy, time and identifier arrays are
fields of `EventHits`. `center_hits` maps `source_frame` positions into
`vertex_frame` positions.

Contract details:

- `EventHits` stores x, y, z, energy and MC time as Float32 vectors of equal
  length. The same concrete type represents source and vertex-centred data;
  `frame` records the coordinate system. Centring subtracts the Float32 vertex
  coordinates and preserves the other fields.
- A hit is identified by `(event_id, particle_id, hit_id)` within its sample.
  `HitDistances` retains the particle and hit IDs in the original hit order;
  `EventSelection.farthest_hit` records the selected `(particle_id, hit_id)` pair.
  Ensemble I/O records sample identity.
- This key was checked in both original NEXUS samples: 2,833,001 unique hits in
  `0nubb.next.h5` and 2,732,400 in `electron_0nubb_energy.next.h5`, with zero
  duplicate keys in either file.
- Coordinate matrices may be assembled at spatial-index boundaries. Hit fields
  remain separate vectors; internal graph/backbone representations are reviewed
  with those functions.
- Constructors validate shapes, finite values, nonnegative energies and ID
  consistency. Event inputs are nonempty; graph filtering may yield zero nodes.
  Explicit result statuses handle this expected analysis outcome.
- `backbone_ok` requires a nonempty backbone; failure requires `nothing`.
  A single-node backbone is a point. Invalid coordinates raise input errors.
  Missing distances use `nothing`.
- Arrays inside a Julia `struct` are mutable. Pure public functions preserve
  input values and return independently owned arrays. Review allocation changes
  after measuring representative workloads.
- Use deterministic tie rules. Distances retain original row order; ties for
  farthest hit use lexicographic `(particle_id, hit_id)` order. Graph/path tie rules are part of the
  backbone implementation review and must be tested.

### Typed signatures and associated tests for joint review

Return types below state intended contracts. Validate inputs explicitly.

```julia
center_hits(hits::EventHits, info::EventInfo)::EventHits
voxelize_hits(hits::EventHits, grid::VoxelGrid)::VoxelCloud
build_voxel_graph(voxels::VoxelCloud, spec::GraphSpec)::VoxelGraph
extract_backbone(voxels::VoxelCloud, graph::VoxelGraph,
                 spec::BackboneSpec)::BackboneResult
point_segment_distance(h::Point3MM, a::Point3MM, b::Point3MM)::Float32
hit_backbone_distances(hits::EventHits, b::Backbone)::HitDistances
classify_event(event_id::Int64, result::BackboneResult,
               distances::Union{Nothing,HitDistances}, cut_mm::Float32)::EventSelection
```

| Function | Required tests before accepting its implementation |
| --- | --- |
| `center_hits` | Vertex subtraction, translation invariance, frame/ID validation, unchanged energy/time/IDs, input preservation |
| `voxelize_hits` | Known cell assignments including negative coordinates and boundaries; independent energy sums; every hit contributes once; canonical ordering |
| `build_voxel_graph` | Known radius edges and components, threshold equality, node/voxel correspondence, all-nodes-filtered case, deterministic energy ties |
| `extract_backbone` | Straight/curved/two-arm fixtures, branches and shortcuts, agreement with agreed Petit fixtures, explicit empty/failure/single-node outcomes |
| `point_segment_distance` | Orthogonal projection, endpoints, collinear points, zero-length segment, known analytic distances |
| `hit_backbone_distances` | Continuous-polyline distances, complete original-hit coverage, independent exhaustive reference, stable farthest-hit identity |
| `classify_event` | Whole-event rejection, equality passes, all-original-hit coverage, monotonic cut scan, inconsistent IDs/results fail validation |

`classify_event` permits `distances=nothing` only for an unsuccessful backbone
result, returning `backbone_failed` with no geometric maximum. Successful results
require matching event IDs and valid distances. The distance routine is called
only for a valid backbone. Review these invariants together with constructors.

Before each implementation, finish reviewing the exact field constraints,
examples, errors, test fixtures and tolerances. 

### Later transport data contracts

Review the transport interfaces together with the diffusion choices. The
proposed minimal contracts are:

- `EventPlacement`: `event_id::Int64`, `vertex_drift_mm::Float32`,
  `grid_phase_mm::NTuple{2,Float32}`.
- `PlacementSpec`: `drift_range_mm::Tuple{Float32,Float32}` and
  `pitch_mm::Float32`. Sampling explicitly consumes `rng::AbstractRNG`.
- `DetectorSpec`: validated Float32 fields with unit suffixes for pitch,
  temperature, field, speed and local IFG/pad geometry; parsed once from the
  shared detector configuration. Exact fields and SI conversion boundary are
  reviewed with transport.
- `PadTimeCounts`: `event_id::Int64`, `pad_ix::Vector{Int64}`,
  `pad_iy::Vector{Int64}`, `time_bin::Vector{Int64}` and
  `expected_ions::Vector{Float32}`; equal lengths, unique cell/bin keys, sorted
  by pad then time. A separate typed axis record carries bin edges and origin.
- `ChargeAccounting`: Float32 input, collected and explicitly omitted expected
  ion counts, with an agreed numerical tolerance.

The next review will define the full detector, diffusion, binning and output
types and settle whether direct integration into pad/time bins is one pure
operation.

### Docstring requirements

Each public computational function must have an extensive Julia docstring with:

1. The signature and physical/geometric purpose.
2. Every argument, including units, dimensions, conventions and defaults.
3. The returned type and the meaning, units and ordering of every output field.
4. The algorithm, equations and relevant Petit provenance where applicable.
5. Assumptions, boundary conditions and numerical approximations.
6. Invalid-input behaviour, expected failure statuses and degenerate cases.
7. Mutation, determinism and RNG behaviour.
8. A small executable example with interpretable output.

Internal helpers also need documentation proportional to their mathematical
role, especially where a convention or approximation is not evident from code.

## 6. Validation strategy

Tests will check independent geometric and physical expectations.

### Reader and coordinates

- Round-trip preservation of original time, x/y/z coordinates, energy and hit
  identities, with Float32 numerical fields.
- Uniqueness of `(event_id, particle_id, hit_id)` within each sample.
- Float32 energy sums versus Float64 references on representative events; report
  errors in keV and relative to the 12.5 keV energy-resolution scale.
- Primary metadata and ordering against the source particle/hit tables.
- Vertex subtraction, common-translation invariance and energy preservation.
- Explicit behaviour for old schemas, missing fields and invalid metadata.

### Graph, backbone and event selection

- Straight, curved and two-arm synthetic tracks with known intended paths.
- Detached deposits, connected side branches and close approaches that could
  create shortcuts; inspect backbone fidelity separately from cut acceptance.
- Energy conservation in voxelisation and correct threshold/component handling.
- Analytic point-to-segment distances, endpoint projection, duplicate points,
  single-node and missing backbones.
- Compare any accelerated distance search against exhaustive segment distances.
- One distant original hit rejects the event, even below the graph energy
  threshold; accepted events retain every hit and their original energy.
- Boundary equality at the cut and monotonic acceptance at 10, 20 and 30 mm.
- Determinism under input reordering and explicit tie cases.
- Resolution and voxel-threshold studies to separate graph instability from
  sensitivity to the physical selection cut. Check voxel-grid phase sensitivity
  and characterise the dependence on orientation.
- Representative NEXUS comparisons with Petit and with primary MC histories.

### Placement, diffusion and pad arrivals

- Reproducible placement independent of event order and threading.
- Integer-pitch translations preserve charge/time sequences up to pad relabelling.
- Zero-diffusion assignment, including a documented cell-boundary convention.
- Gaussian widths and cell probabilities against analytic cases and Monte Carlo
  with statistical tolerances; verify expected sharing at a cell boundary.
- Account for input charge, output charge and any explicitly truncated tails.
- Consistent arrival-time ordering, bin-edge handling and spatial/time conversion.
- Convergence with integration tolerances and time resolution; benchmark runtime
  on representative events before choosing an optimisation.

## 7. Implementation and review sequence

1. **Agree the data contracts and first function signatures/tests.** Start with
   reader/schema changes, centring and voxelisation. Then agree graph/backbone
   and distance interfaces before implementing them.
2. **Deliver MC topology selection.** Add time, construct backbones, save per-hit
   distances and scan 10/20/30 mm. Produce acceptance tables, maximum-distance
   distributions and selected event displays, including rejected/ambiguous cases.
3. **Review the selection results together.** Confirm backbone settings and
   the default distance cut before using the selected ensemble downstream.
4. **Agree transport interfaces and tests.** Settle placement conventions,
   diffusion options, grid mapping and pad/time output contracts.
5. **Deliver detector arrivals.** Implement and validate fast diffusion and
   collection; compare pad-charge distributions across the saved selection cuts.
6. **Scope readout and reconstruction separately.** Feed the reviewed arrival
   distributions into the existing response model, then study reconstructed
   backbones and the tighter reconstructed-track proximity requirement.

At each increment, run relevant Julia tests and required repository checks for
affected functionality, save the numerical evidence, and regenerate the report.
Review additional physics or API decisions together before implementing them.
