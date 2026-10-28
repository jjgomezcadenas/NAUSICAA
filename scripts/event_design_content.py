"""Canonical design narrative and interface drafts for the offline HTML guide.

The builder also derives the short Markdown roadmap from this module. Scientific
decisions are labelled separately from API/schema proposals and measured results.
"""
from pathlib import Path
import math
import re
import json
import tomllib
from event_viewer_content import viewer_page


MILESTONES = [
    ('Review contracts', 'Review types, function signatures, ownership, configuration and saved schemas in the HTML. Settle the listed open contracts and approve the implementation scope.'),
    ('Input and local extraction · step 1 complete', 'Extended reduced HDF5 with MC time and actual primary IDs. Copied the required Petit functions into src/event_analysis with licence and provenance, retaining their algorithms. Updated the root Julia dependencies.'),
    ('Deterministic reconstruction · step 1 complete', 'Adapted voxelisation, graphs, backbone and truth reference to Float32. Added geometry, threshold and indexing tests. Compared copied stages with the pinned Petit reference.'),
    ('Ion production and diffusion · step 2 complete', 'Implemented explicit RNG boundaries, common W-value, per-hit rounding, thermal widths and histogram accounting. Validated statistics, reproducibility and charge counts, with generated figures and a PDF checkpoint.'),
    ('Selection and pad histories · step 3 complete', 'Connected exact ion-threshold reconstruction to single-track selection and the full retained diffusion table to pad histories. Validated grid phase, extended coverage, relative-time axes and exact charge conservation.'),
    ('Ensembles and paired scans · step 4 complete', 'Save exact ensemble numerators, both denominators, aligned histories and provenance in HDF5. Compare numerical and selection settings on paired events. Python generates PNGs, tables and a compiled LaTeX report.'),
    ('Production review and full samples', 'Review ensemble alignment and the measured histogram sensitivity, select the required pad-response accuracy, then run both full samples and measure runtime/memory.'),
]

CONFIG = '''# Proposed event-analysis card; names/schema await interface review.
# Detector values are read from config/detector.toml.
[inputs]
detector_card = "config/detector.toml"
samples = ["events/0nubb_hitdf.h5", "events/electron_0nubb_energy_hitdf.h5"]

[ionization]
ionisation_energy_eV = 21.9 # Current config/ion_tracks.toml: ionization.w_eV

[drift]
min_cm = 10.0
max_cm = 150.0

[diffusion]
longitudinal = "thermal" # "zero" is the approved comparison
nbins = 100
nsigma = 3.0

[graph]
method = "kNN"
k = 10
min_ions_per_bin = 10
binning = "sigma_scaled"
bin_sigma_scale = 3.0
edge_sigma_scale = 6.0

[truth]
bin_width_mm = 0.5 # Petit reference value; explicit contract for review

[pads]
z_bin_width_mm = 2.0
random_grid_offset = true
# xy pitch comes from detector.geometry.pitch_m (currently 2 mm).

[random]
seed = 2026 # Example run seed; stream derivation is an interface-review item.

[output]
results_dir = "results/event_analysis"
report_dir = "latex_generated/event_analysis"
'''


def roadmap_markdown():
    items = '\n'.join(f'{i}. **{title}.** {body}' for i, (title, body) in enumerate(MILESTONES, 1))
    return f'''# Event analysis: implementation roadmap

The canonical design and decision record are in the
[HTML documentation](event_analysis_review/index.html). Read Purpose and physics,
Event flow, Reuse and development, Interface drafts, Configuration and outputs,
and Validation and roadmap. The Petit source assessment and the full decision
record remain embedded as appendices.

## Current status

Step 1 is implemented: the v3 hit reader preserves MC time and actual primary
IDs; the local Petit adaptation supplies Float32 binning, graphs, backbone and
truth paths. Both reduced samples are regenerated and verified. The HTML links
the code docstrings, regression evidence and generated checkpoint report.
Step 2 is implemented: ion production, thermal widths, explicit RNG sampling and
Petit Monte Carlo diffusion with full charge accounting. The numerical checkpoint
covers all events for ion production, 40 real-event diffusion cases, Gaussian
statistics and histogram resolution/padding scans. Results are in
results/event_analysis/step2_validation.h5 and its generated PDF report.
Step 3 is implemented: exactly-one-component selection, integer ion thresholds,
full retained-charge pad histories, explicit axes and selectable relative-time
order. The integration checkpoint saves 100 validation events, decisions and
accepted histories. The subset acceptance is recorded without tuning cuts.
Step 4 is implemented: central-pad ensemble alignment, exact integer numerators,
accepted/processed denominators, HDF5 round trips and twelve paired settings on
100 validation events. Selection is stable across 100, 200 and 400 histogram
bins; pad histories retain measurable resolution sensitivity. Review the working
alignment and required pad-response accuracy before the full-sample run.

Interactive viewer: the user has approved the two-tab MC energy / MC tracks
scope, Petit defaults, repository layout and explicit Apply controls. Viewer
step 1 interfaces were accepted as implementation guidance. Step 2 is implemented:
the local browser displays energy maps and occupied 3D voxels, with navigation,
camera controls and session restore. Launch with `make view-events`; chapter 7
contains usage and validation details. Browser review is accepted, including
event 2 with independent axis scales and independent projection zoom.
Step 3 is implemented: MC graph inspection uses the existing Petit adaptation,
with a component table, switchable layers, backbones and endpoints. The graph
and applied TOML settings are saved in HDF5 and reused across display rebins.
Next: review MC tracks, then interactive graph parameter controls (step 4).

The central flow is: read hits → produce ions → diffuse → retain the full
diffusion histogram → reconstruct and select exactly one connected track →
fill pad/time histories from the full histogram of accepted events. The
TOML-controlled graph threshold affects reconstruction. Pad histories retain
the accepted event's charge below that threshold.

Reconstruction uses b = 3 sigma_T and maximum connection distance = 6 sigma_T
throughout the drift range, with both factors in TOML. This replaces the earlier
fixed-bin proposal and removes Petit's 50 cm switch from the baseline. Physical
pad pitch remains independently fixed at 2 mm.

F12 is resolved: retain explicit Float32 edge vectors in computational interfaces
and the initial saved-axis schema. Compact axis storage can be considered later
without changing the computational interface.

## Development sequence

{items}

## Documentation and review

The HTML includes the physical assumptions, data flow, Petit reuse map, complete
draft contracts, configuration example, HDF5 proposal, examples, test strategy
and open decisions. Review one type or function at a time. Extensive docstrings
are part of each implementation increment. Step-one numerical validation is saved in results/event_analysis/step1_validation.h5
and presented in latex_generated/event_analysis_step1/report.pdf.

Rebuild with `make event-review PYTHON=/path/to/python3.11-or-newer`.
The builder reads `scripts/event_design_content.py` for the design and this
roadmap, `scripts/petit_review_content.py` for the source assessment/decision
record, and `docs/event_analysis_review/template.html` for presentation.
The default Petit checkout is `../Petit`; the builder also accepts `--petit-root`.
Refresh the existing HTML tab after rebuilding.

The [previous plan](event_analysis_review/archive/event_analysis_plan_before_consolidation.md)
and [previous review](event_analysis_review/archive/review_before_consolidation.html)
are historical snapshots. Their earlier API sketches and workflow are superseded
by the current HTML design.
'''


def design_pages(root: Path, section, table, code):
    detector = tomllib.loads((root / 'config/detector.toml').read_text())
    ions = tomllib.loads((root / 'config/ion_tracks.toml').read_text())
    gas, bias, geom = detector['gas'], detector['bias'], detector['geometry']
    temp, field, speed = gas['temperature_K'], bias['bulk_field_V_m'], gas['bulk_speed_m_s']
    w = ions['ionization']['w_eV']
    pitch = geom['pitch_m'] * 1000
    sigma_rows = []
    for length_cm in (10, 100, 150):
        sigma_m = math.sqrt(2 * 1.380649e-23 * temp * (length_cm / 100) / (1.602176634e-19 * field))
        sigma_rows.append((str(length_cm), f'{sigma_m * 1000:.3f}', f'{sigma_m / speed * 1000:.2f}'))

    guide = section('The scientific question', '''
<p>Given a simulated energy-deposition track in xenon, determine whether its
diffused topology forms one connected track and calculate how many positive ions
arrive at each sensor pad as a function of time. The saved pad histories supply
the charge input for the signal and electronics model in <code>src/ion_signal</code>.</p>
<p>The two input ensembles are double-beta events and single-electron backgrounds.
Topology selection and pad-charge distributions are studied together, with
selection efficiencies and charge accounting reported for each sample.</p>
<p><strong>Reading order:</strong> physics → event flow → reuse strategy → concrete
interface drafts → configuration/output → validation. The three appendices retain
the inspected Petit source, its original interfaces and the discussion’s decision record.</p>''', 'design-purpose')
    guide += section('Design status', table(['Category', 'Status'], [
        ('Approved design', 'Petit algorithms, Float32 event arithmetic, drift sampling, isotropic thermal diffusion, Monte Carlo/rounding, graph settings, grid phase and full-diffusion pad histories.'),
        ('Implemented steps 1–3', 'Foundation, v3 reader, Petit reconstruction and diffusion; integer-threshold single-track selection; full-charge sparse pad histories, grid coverage, relative time and per-event projections.'),
        ('Draft for joint review', 'Adapted signatures, ownership details, HDF5 layout, RNG stream scheme, pad-axis coverage and relative-time convention.'),
        ('Validation evidence', 'Pinned Petit graph/backbone fixtures, independent geometry and ownership tests, complete reduced-file verification and real-event software checks. See the step-one checkpoint below.'),
    ]) + '<p>Steps 1–4 are complete. Next: review the ensemble convention and measured histogram sensitivity before full-sample production.</p>', 'design-status')
    guide += section('Physical quantities and units', table(['Quantity', 'Meaning and convention'], [
        ('Xe₂⁺', 'Molecular xenon ion with charge +e.'),
        ('Hit energy and position', 'Original NEXUS energy in MeV; positions in mm; MC deposition time in ns. Preserve the source Float32 precision.'),
        ('W-value', f'Mean energy to create an electron–ion pair. The existing nominal input is {w:g} eV in config/ion_tracks.toml. The adapted key is ionisation_energy_eV; use the same value in all charge conversions.'),
        ('Ion expectation and sampling', 'mean_ions = energy_MeV × 10⁶ / W_eV; round each hit’s expectation to an integer for the approved Monte Carlo. Save both expectation and rounding difference.'),
        ('SI interface', 'Detector cards and field/signal calculations use SI. Convert explicitly when constructing event-analysis parameters: pitch to mm, field to V/cm where required by Petit, speed to mm/ms.'),
        ('Arithmetic', 'Float32 physical arrays and reductions; integer identifiers and ion counts. Quantify numerical differences against Float64. The energy-resolution reference is 0.5% at 2.5 MeV = 12.5 keV.'),
    ]), 'design-units')
    guide += section('Bulk transport model', '''
<p>Sample one drift distance L uniformly from 10 to 150 cm per event. Use that L
to calculate common diffusion widths for all its hits. Preserve the source
coordinates and axis directions. This approximates the event’s transport with
a common drift distance; its longitudinal track structure remains in the hit coordinates.</p>
<p>Review finding F1 is resolved: retain this common-distance approximation for
the study. The per-hit drift refinement was considered and declined.</p>
<div class="equation">σ² = 2 (k<sub>B</sub>T / e) L / E<br>σx = σy = σz = σ &nbsp; · &nbsp; Δt = Δz / v<sub>ion</sub></div>
<p>This is thermal isotropic diffusion in a uniform bulk field for charge +e.
Petit already contains the isotropic choice in <code>select_events_itaca</code>.
The batch runner’s zero longitudinal width remains a configurable comparison.</p>'''
        + table(['Input from detector.toml', 'Current value'], [
            ('Temperature', f'{temp:g} K — room-temperature design input'),
            ('Uniform bulk field', f'{field / 100:g} V/cm'),
            ('Constant drift speed', f'{speed:g} m/s = {speed:g} mm/ms — operating assumption'),
            ('Bulk reference length', f'{geom["bulk_reference_length_m"]:g} m'),
        ]) + table(['L / cm', 'One-coordinate σ / mm', 'Longitudinal σ / v / ms'], sigma_rows)
        + '<p class="caption">Illustrative values are calculated by the documentation builder from the current detector card and the displayed formula. They are uniform-field thermal estimates, not measured transport or an executed event study.</p>', 'design-transport')
    guide += section('Collection and time histories', f'''
<p>Use a periodic square collection grid with pitch {pitch:g} mm. The selected
track supplies a transverse reference; independently shift grid edges in x and
y by uniform offsets in [0, pitch). Use the same grid for all slices of that event.</p>
<p>Geometric cell assignment represents the approximation that ions in each
collection cell reach its corresponding hole/pad. Pitch determines collection
cells; the physical metal pad width is a separate detector dimension. The pad
histories are the input to subsequent signal formation through the IFG and pad gap.</p>
<p>Slice the full diffused distribution in 2 mm longitudinal bins. At the configured
speed, each bin spans {2 / speed:g} ms. Preserve both spatial and time axes with
the per-pad counts. The relative-time origin and drift-axis sign are explicit
interface-review items below.</p>
<p>The revised detector design specifies 5 kV across the 10 mm IFG and 3 kV
across the 5 mm pad gap. The event-stage Gaussian widths use the bulk field;
local transport and weighting-field response belong to the signal calculation.</p>''', 'design-collection')

    flow = section('One event, two consumers of the diffused charge', '''
<div class="pipeline" role="img" aria-label="Hits are diffused; a reconstruction branch selects events and an accepted-event branch counts all diffused charge in pads.">
<div class="pipeline-node">HDF5 → EventHits + EventInfo → ion counts + sampled L</div>
<div class="pipeline-arrow">↓</div><div class="pipeline-node highlight">Diffuse ions → full occupied-bin table D + charge accounting</div>
<div class="pipeline-arrow">↓ keep D available</div>
<div class="pipeline-branches"><div><b>Reconstruction</b><p>D → reconstruction bins → ion threshold → kNN → connected components</p><p>Exactly one component → accepted track → backbone and diagnostics</p></div>
<div><b>Pad histories</b><p>Accepted event → use all of D → shifted pad grid + z slices → pad counts</p><p>The graph threshold controls selection; counts below it remain in D.</p></div></div>
<div class="pipeline-arrow">↓</div><div class="pipeline-node">HDF5 results → Python figures/tables → generated LaTeX → PDF</div></div>
<p>The time-ordered primary trajectory is a parallel truth reference built from
the original hits. Compare it to the reconstructed backbone in the same source frame.</p>''', 'flow-overview')
    flow += section('Stages and data products', table(['Stage', 'Operation', 'Output'], [
        ('Read', 'Load one sample/event; retain hits and primary identity.', 'EventHits + EventInfo; source/sample provenance.'),
        ('Truth reference', 'Sort primary hits by MC time; group consecutive visits to the same bin; join reversed first arm and second arm for two primaries.', 'Truth path table and actual primary IDs.'),
        ('Production', 'Calculate E_event, mean ions per hit and total; round each hit for sampling.', 'Expected and integer count vectors; rounding difference.'),
        ('Drift', 'Sample L; evaluate widths and graph scales.', 'Per-event drift/diffusion/reconstruction settings.'),
        ('Diffusion', 'Sample three Gaussian displacements per ion; accumulate in padded 3D histogram.', 'Full occupied-bin table D; sampled/in-range/out-of-range counts.'),
        ('Reconstruction', 'Rebin D; keep bins at or above graph.min_ions_per_bin; build graph.', 'Tracks ordered by energy; component multiplicity and threshold accounting.'),
        ('Selection', 'Accept exactly one connected component.', 'Decision for every event; accepted track for backbone extraction.'),
        ('Pad counts', 'For accepted events, use D with the shifted grid and z/time bins.', 'Sparse pad/slice counts, axes and grid offsets.'),
        ('Reduction', 'Accumulate event energy/ion summaries, acceptance, pad maps and z projections in Julia.', 'HDF5 event records and ensemble reductions.'),
        ('Presentation', 'Read saved numerical reductions in Python.', 'PNGs, tables, macros and standalone LaTeX/PDF.'),
    ]), 'flow-stages')
    flow += section('Three grids with different purposes', table(['Grid', 'Definition', 'Use'], [
        ('Diffusion histogram', 'Event extent plus nsigma × σ padding; initially 100 bins on each axis.', 'Intermediate representation of sampled ion positions; converged numerically.'),
        ('Reconstruction grid', 'Petit floor(coordinate / b) convention with zero-anchored boundaries; b = 3σt and maximum graph distance = 6σt = 2b throughout the drift range.', 'Graph vertices and connectivity.'),
        ('Collection grid', '2 mm pitch in x/y with one random phase per event; 2 mm z slices.', 'Pad charge and arrival-time histories.'),
    ]) + '<p>Graph voxelisation receives the diffusion table; pad counting receives that same original diffusion table. Changes to reconstruction binning therefore affect acceptance, while the accepted event’s charge input stays available at diffusion-histogram resolution.</p>', 'flow-grids')
    flow += section('Charge accounting', code('''N_expected = sum(mean_ions_per_hit)
N_sampled = sum(integer_ions_per_hit)
rounding_difference = N_sampled - N_expected

N_sampled = N_diffusion_retained + N_diffusion_outside
N_diffusion_retained = N_graph_retained + N_graph_below_threshold

# For accepted events with complete pad-axis coverage:
sum(pad_slice_counts) = N_diffusion_retained''') + '''
<p>These equalities separate rounding, finite diffusion bounds and the graph
threshold. Each is saved with its own name. Pad charge conservation uses the
full retained diffusion count. Event rejection is recorded as an event decision.</p>
<p>Example: two occupied diffusion cells contribute 6 and 14 ions to different
reconstruction bins. At threshold 10, the graph contains the 14-ion bin. If the
event passes the one-component condition, its pad histories contain all 20 ions.
This is the approved change from Petit’s selected-track occupancy input.</p>''', 'flow-accounting')
    flow += section('Ensemble reductions', '''
<p>Save per-event energy, total expected ions, mean ions per hit and total sampled
ions. A per-event pad history retains its event-specific axes and grid offset.
For ensemble maps, the working convention labels the pad containing the selected
track’s x/y reference as (0,0). Existing cells retain their integer counts.
Arrival bin zero starts at each history’s first saved arrival edge.</p>
<p>Julia sums exact Int64 counts on shared axes covering every cell. Float32
means use two explicit denominators: accepted events, and all processed events
with rejected events contributing zero selected charge. Both appear in figures
and HDF5. With no accepted events the conditional mean is NaN; the processed
mean is zero. Review this working alignment alongside the step-four results.</p>''', 'flow-ensemble')

    reuse = section('Development philosophy', '''
<p>This is a focused adaptation of Petit. Copy the required loaded functions,
preserve their algorithms and MIT attribution, and make small reviewable changes.
Maintain a source map recording the original file, function, source hash and
reason for each adaptation. Existing graph and backbone machinery supplies the
numerical core.</p>
<p>Deterministic functions derive results from explicit inputs and preserve those
inputs. Random sampling receives an explicit RNG. HDF5 access, TOML parsing,
execution order and report generation live at the application boundary. Typed
parameters carry validated values into computation.</p>
<p>Review each type and function together: purpose, fields/signature, units,
ownership, examples and tests. Implement in small increments that produce
checkable output. Keep physical assumptions visible beside numerical results.</p>''', 'reuse-philosophy')
    reuse += section('What carries over', table(['Petit source', 'Reuse', 'Targeted changes'], [
        ('src/voxels.jl', 'voxelize_event; cell keys, centres and energy sums.', 'Float32; explicit units; keep ion count accounting consistent.'),
        ('src/track_building.jl', 'Tracks, DiffusionParams, KDT/kNN/mutual-kNN, components and energy ordering.', 'Float32 physical fields; explicit parameter units; preserve local/parent index conventions.'),
        ('src/track_extreme_finding.jl', 'TrackCoords, combined backbone algorithm, energy weighting, traversal and MST comparison.', 'Float32 allocations/annotations; document output ownership.'),
        ('src/track_path_reconstruction.jl', 'get_raw_path and endpoint comparisons.', 'Float32 paths; documented empty/single-point behaviour.'),
        ('src/mctrack.jl', 'Time ordering, consecutive-bin grouping and two-arm joining.', 'Read actual primary IDs from EventInfo.'),
        ('src/itaca_functions.jl', 'Thermal formula, Monte Carlo histogram and the low-drift graph-scale prescription.', 'Explicit RNG; configured detector/W; isotropic default; count accounting; apply b = 3σt and d = 6σt at all drift distances.'),
        ('itacaScripts/itaca_simple_reco.jl', 'Pixel/slice assignment and occupancy reductions; satellite distance diagnostic.', 'Full diffusion input; random grid offset; independent pitch; retain per-slice counts; expose reusable functions.'),
    ]) + '<p>The graph-building and extreme-finding files together contain 1,512 lines including documentation in the inspected checkout. Most of that algorithmic work is reusable. The final extraction diff will give the exact reused/changed scope.</p>', 'reuse-map')
    reuse += section('New integration work', table(['Area', 'Deliverable'], [
        ('Existing reader', 'Add MC time and actual primary IDs to reduced HDF5; construct existing EventHits/EventInfo.'),
        ('Orchestration', 'Resolve cards, choose per-event random streams, run shared diffusion and its two consumers, record decisions.'),
        ('Persistence', 'Write pad histories, axes, charge ledger, settings, provenance and ensemble reductions.'),
        ('Presentation', 'Python renders saved numerical outputs into PNGs, tables and generated LaTeX; compile to PDF.'),
        ('Documentation/tests', 'Extensive docstrings, worked examples, adaptation checks, conservation and convergence evidence.'),
    ]), 'reuse-integration')
    reuse += section('Repository organisation and root Julia project', code('''src/event_analysis/          # existing module; local Petit adaptation and event pipeline
tests/event_analysis/        # contracts, copied-algorithm regression and physics checks
scripts/                     # Julia execution; Python presentation and doc builder
config/detector.toml         # existing authoritative detector design inputs
config/event_analysis.toml   # proposed event-analysis card
results/event_analysis/      # proposed numerical HDF5 outputs
latex_generated/event_analysis/  # proposed generated report and figures
docs/event_analysis_review/index.html  # main design and review interface
docs/event_analysis_plan.md  # generated short roadmap
work/                       # scratch and development evidence''') + '''
<p>Keep one root <code>Project.toml</code>. Current dependencies include DataFrames,
HDF5, Dates, SHA, TOML and Test. The copied graph code requires Graphs and
NearestNeighbors; declare Random, LinearAlgebra, SparseArrays and Statistics
where directly used. Resolve and review compatibility bounds during extraction,
then record Project/Manifest and Julia versions with results. The present Julia
compatibility declaration is 1.11.</p>
<p>The existing <code>src/ion_signal</code> calculations remain the signal-model
boundary. Migrate shared event I/O deliberately with compatibility checks for its
current commands.</p>''', 'reuse-layout')

    contracts = section('How to review these contracts', '''
<p><strong>Implementation status:</strong> step-one through step-three interfaces are labelled implemented
and include their current source docstrings. Later stages retain draft labels. Petit’s exact original declarations
remain available in <a href="#interfaces">the source-interface appendix</a>.</p>
<p>Review one card at a time. Each card gives the data meaning, types, ownership,
relevant errors or edge cases, and tests that establish its contract. Functions
retain Petit names where their role is unchanged.</p>''', 'contract-reading')
    contracts += section('CoordinateFrame · implemented', code('''@enum CoordinateFrame source_frame vertex_frame''') + table(['Value', 'Meaning'], [
        ('source_frame', 'Source NEXUS axes and coordinates; used by the adopted event pipeline.'),
        ('vertex_frame', 'Coordinates after subtracting the common true generation vertex; source axis directions preserved.'),
    ]) + '<p>The frame labels the coordinates of EventHits. The implemented center_hits helper returns vertex-frame hits with independently owned arrays.</p>', 'contract-frame')
    contracts += section('EventInfo · implemented', code('''struct EventInfo
    event_id::Int64
    vertex_mm::NTuple{3,Float32}
    primary_particle_ids::Vector{Int32}
end''') + '''<p>One event’s source identifier, common primary generation vertex
and actual primary IDs. One electron has one primary ID; double beta has two.
Sample identity is supplied by the enclosing dataset. Construction stores a
sorted copy of nonempty distinct IDs and requires a finite vertex.</p>''', 'contract-info')
    contracts += section('EventHits · implemented', code('''struct EventHits
    event_id::Int64
    frame::CoordinateFrame
    x_mm::Vector{Float32}
    y_mm::Vector{Float32}
    z_mm::Vector{Float32}
    energy_MeV::Vector{Float32}
    time_ns::Vector{Float32}
    particle_id::Vector{Int32}
    hit_id::Vector{Int32}
end''') + '''<p>Entry j of each vector describes the same hit. Vectors have equal
positive length; physical values are finite and energy is nonnegative. The pair
(particle_id, hit_id) is unique within an event. Construction copies each array.
Current arrays are validated at computational boundaries because their contents
remain mutable. MC time supports the truth trajectory.</p>
<p>The v3 reduced reader supplies positions, energy, IDs, event vertices, time_ns and actual primary IDs. Both supplied samples have been regenerated and verified.</p>''', 'contract-hits')
    contracts += section('Petit table adapter · implemented; later tables · draft', code('''to_petit_dataframe(hits::EventHits)::DataFrame''') + table(['Table', 'Columns and concrete types'], [
        ('Input adapter', 'event_id::Int64; x/y/z/energy/time::Float32; particle_id/hit_id::Int32. Units mm/MeV/ns. One row per hit.'),
        ('Diffusion table D', 'event_id::Int64; x/y/z/energy::Float32; ions::Int64. Coordinates are histogram-bin centres; energy = ions × W_eV / 10⁶.'),
        ('Reconstruction bins', 'event_id::Int64; x/y/z/energy::Float32. One row per occupied graph cell; values are geometric centres and summed energies.'),
        ('Path', 'x/y/z/s::Float32, all mm. Separate Vector{Int} indexes the selected track rows in path order.'),
        ('Pad histories', 'pad_ix/pad_iy/slice_index::Int64; ions::Int64. Unique occupied-cell keys, ordered by pad and slice. Event identity comes from the enclosing event result.'),
    ]) + '''<p>The adapter constructs independent arrays and preserves row order,
values and IDs. Petit calls the integer
column <code>electrons</code>; the adapted positive-ion table calls it
<code>ions</code>. The implemented diffusion table uses ions::Int64.</p>''', 'contract-tables')
    contracts += section('Tracks, TrackCoords and DiffusionParams · implemented', code('''struct DiffusionParams
    drift_length_cm::Float32
    sigma_t_mm::Float32
    sigma_l_mm::Float32
    bin_width_mm::Float32
    max_distance_mm::Float32
    energy_threshold_keV::Float32
    histogram_nbins::Int
    histogram_nsigma::Float32
end

struct Tracks
    voxels::DataFrame
    graph::SimpleGraph{Int}
    components::Vector{Vector{Int}}
    diffusion::DiffusionParams
end

struct TrackCoords
    x::Vector{Float32}
    y::Vector{Float32}
    z::Vector{Float32}
end''') + '''<p>Retain Petit’s structures and index semantics. Local graph vertex i
refers to row i of track.voxels. The stored component vector refers to rows in the
energy-filtered parent table. Names in DiffusionParams make supplied units
explicit, including the runner’s drift distance in cm. The adapter constructs
new tables and graphs; callers preserve the returned row order when using indices.</p>''', 'contract-tracks')

    functions = [
        ('contract-production', 'Ion production · implemented',
         'ion_counts(hits::EventHits, ionisation_energy_eV::Float32)',
         'Returns a NamedTuple: expected::Vector{Float32}, sampled::Vector{Int64}, event_energy_MeV::Float32, total_expected::Float32, mean_per_hit::Float32 and rounding_difference::Float32. W must be finite and positive. Each count uses the approved per-hit rounding rule; input arrays retain their values.',
         'Known E/W examples, zero-energy hits, summed expectations, exact integer counts and event rounding difference.'),
        ('contract-widths', 'Thermal widths and reconstruction scales · implemented',
         'sigma_t_ion_mm(T_K::Float32, L_cm::Float32, E_V_cm::Float32)::Float32\nget_voxel_size_and_distance(sigma_t_mm::Float32;\n    bin_sigma_scale::Float32=3f0, edge_sigma_scale::Float32=6f0,\n    truth_bin_width_mm::Float32=0.5f0)',
         'Return the thermal width in mm. Use positive finite T/E and nonnegative L. The caller sets sigma_l_mm equal to sigma_t_mm for the baseline. The graph-scale helper returns (bin_width_mm, truth_bin_width_mm, max_distance_mm), each Float32, preserving Petit’s tuple order. Use b = 3σt and d = 6σt throughout the drift range, with scale factors from TOML. The interface applies these scales at every drift distance.',
         'Independent thermal expression, L/T/E scaling, zero L and unit conversion. Verify b/σt = 3 and d/σt = 6 at positive drift lengths, including both sides of 50 cm.'),
        ('contract-sampling', 'Event random choices · implemented',
         'sample_drift_cm(rng::AbstractRNG, min_cm::Float32, max_cm::Float32)::Float32\nsample_grid_offset_mm(rng::AbstractRNG, pitch_mm::Float32)::NTuple{2,Float32}',
         'Sample one L and two independent grid offsets. Equal drift bounds return the fixed value. Offsets are in [0, pitch). Positive finite pitch and ordered finite drift bounds are required. These functions advance the supplied RNG; stream assignment belongs to orchestration.',
         'Bounds, fixed-drift case, seeded reproducibility and ensemble distribution. Review per-sample/event/stage stream derivation to support order-independent runs.'),
        ('contract-diffuse', 'Monte Carlo diffusion · implemented',
         'diffuse_xyz_image_mc(rng::AbstractRNG, hits::EventHits, ions::Vector{Int64};\n    sigma_t_mm::Float32, sigma_l_mm::Float32, nbins::Int,\n    nsigma::Float32, ionisation_energy_eV::Float32)',
         'Returns a NamedTuple: bins::DataFrame, sampled_ions::Int64, retained_ions::Int64, outside_ions::Int64 and x/y/z_edges_mm::Vector{Float32}. Validate nonnegative widths/counts, matching lengths and positive histogram settings. Preserve hit inputs. Gaussian draws and the dense intermediate histogram follow Petit.',
         'Means/variances/cross-axis covariance, seeded reproducibility, count accounting, upper edges, coincident positions, zero longitudinal width and convergence with bin count/padding. A coincident axis with zero diffusion uses one bin centred at the input coordinate, with adjacent representable Float32 values as edges.'),
        ('contract-voxelize', 'Voxelisation · implemented',
         'voxelize_event(bins::DataFrame, bin_width_mm::Float32)::DataFrame',
         'Use floor(coordinate/bin_width) and geometric cell centres; accumulate energy per occupied cell. Width must be finite and positive. Return newly owned columns. Dictionary row order may be retained; compare results using integer cell keys.',
         'Positive/negative positions, exact boundaries, nearby Float32 values, known energy sums and unchanged input. Energy conservation is measured against Float64.'),
        ('contract-graph', 'Graph construction and selection · implemented',
         'make_tracks(bins::DataFrame; method::String="kNN", k::Int=10,\n    max_distance_mm::Float32, energy_threshold_kev::Float32=0f0,\n    diffusion::DiffusionParams=DiffusionParams())::Vector{Tracks}',
         'Retain Petit’s public method names and keV argument convention. Apply energy ≥ threshold, build the chosen undirected graph and return components in decreasing energy order. Finite maximum distance and k ≥ 1 are the baseline. Event acceptance is length(tracks) == 1.',
         'Require matching edges for fixtures with distinct neighbour distances. With ties, verify valid kNN choices and compare components, charge, acceptance and backbone observables; record consequential differences. Also test disconnected tracks, threshold equality, finite-distance boundaries, empty/single-bin input, local/parent indices and ownership.'),
        ('contract-backbone', 'Backbone extraction · implemented',
         'walk_track_from_extremes(track::Tracks;\n    use_energy_weighting::Bool=true, use_edge_energy_weighting::Bool=true,\n    use_mst_fallback::Bool=false, dense_track_threshold::Float32=6f0)\nget_raw_path(track::Tracks, path_indices::Vector{Int})::DataFrame',
         'Retain the walk NamedTuple: extremes, path_indices, path_voxels, total_length and confidence. Numeric values become Float32. Each endpoint is an owned NamedTuple copy. path_voxels is an independent copied table; path_indices uses local graph rows. get_raw_path returns x/y/z/s with cumulative Euclidean distance.',
         'Straight/curved/two-arm/branched examples, equivalent reversal, index correspondence, total length and single/empty paths. Check energy-weighted branch choices near thresholds.'),
        ('contract-truth', 'MC truth trajectory · implemented',
         'compute_mc_path(hits::EventHits, info::EventInfo, bin_width_mm::Float32)::DataFrame',
         'Use EventInfo’s actual primary IDs and the approved time-order/consecutive-bin algorithm. Return Float32 x/y/z/energy/s and actual particle_id::Int32. The one/two-primary paths preserve Petit’s joining procedure. Equal-time hits retain source row order; any missing primary gives a typed empty path, following Petit.',
         'Single electron, two arms with nonstandard particle IDs, revisited cells, zero-energy groups, endpoint comparisons and input preservation.'),
        ('contract-pads', 'Pad/slice counts · implemented',
         'compute_occupancy_analysis(diffused::DataFrame;\n    x_edges_mm::Vector{Float32}, y_edges_mm::Vector{Float32},\n    z_edges_mm::Vector{Float32}, drift_speed_mm_ms::Float32,\n    drift_direction::Int)',
         'Accept the full retained diffusion table of an accepted event. Return a NamedTuple with counts::DataFrame and independently owned x/y/z_edges_mm and time_edges_ms. Edges define collection cells; use the same shifted x/y edges for every slice. Positive speed and strictly increasing finite edges are required. Internal cells are half open and the final edge is included. Axes cover the full retained table; missing positive-count coverage raises an error. slice_index increases with relative time, using the explicit direction ±1 and t=0 at the first arrival edge.',
         'Known pad/slice assignment, cells below the graph threshold retained, sum(counts) = retained diffusion ions, integer-pitch relabelling, edges/corners, time conversion and reconstruction-input preservation.'),
    ]
    def live_doc(name, source):
        text = (root / source).read_text()
        pattern = r'(?:@doc )?"""((?:(?!""")[\s\S])*)"""\s*(?:function\s+|struct\s+)?' + re.escape(name) + r'\b'
        matches = re.findall(pattern, text)
        if not matches:
            raise ValueError(f"No source docstring for {name} in {source}")
        return '<details><summary>' + name + ' · implemented docstring</summary>' + code(matches[-1]) + '<p><a href="../../' + source + '">Open source</a></p></details>'

    implemented = {
        'contract-voxelize': [('voxelize_event','src/event_analysis/reconstruction.jl')],
        'contract-graph': [('select_diffused_event','src/event_analysis/collection.jl'), ('make_tracks','src/event_analysis/paths.jl'), ('is_single_track','src/event_analysis/reconstruction.jl')],
        'contract-backbone': [('walk_track_from_extremes','src/event_analysis/paths.jl'), ('get_raw_path','src/event_analysis/paths.jl')],
        'contract-truth': [('compute_mc_path','src/event_analysis/paths.jl')],
        'contract-widths': [('sigma_t_ion_mm','src/event_analysis/production.jl'), ('get_voxel_size_and_distance','src/event_analysis/reconstruction.jl')],
        'contract-production': [('ion_counts','src/event_analysis/production.jl')],
        'contract-sampling': [('sample_drift_cm','src/event_analysis/production.jl'), ('sample_grid_offset_mm','src/event_analysis/production.jl')],
        'contract-diffuse': [('diffuse_xyz_image_mc','src/event_analysis/diffusion.jl')],
        'contract-pads': [('compute_occupancy_analysis','src/event_analysis/collection.jl'), ('pad_axes','src/event_analysis/collection.jl'), ('pad_projections','src/event_analysis/collection.jl')],
    }
    contracts += section('Implemented reader, adapter and track records', ''.join(
        live_doc(name, source) for name, source in [
            ('read_event','src/event_analysis/io.jl'),
            ('to_petit_dataframe','src/event_analysis/reconstruction.jl'),
            ('DiffusionParams','src/event_analysis/reconstruction.jl'),
            ('Tracks','src/event_analysis/reconstruction.jl'),
            ('TrackCoords','src/event_analysis/petit/backbone.jl')]), 'contract-live-records')
    for anchor, title, sig, meaning, tests in functions:
        live = ''.join(live_doc(name, source) for name, source in implemented.get(anchor, []))
        contracts += section(title, code(sig) + f'<p>{meaning}</p><p><strong>Contract tests:</strong> {tests}</p>' + live, anchor)
    contracts += section('Ensemble alignment, reduction and persistence · implemented',
        '<p>These reductions extend Petit’s event-level collection to sample means. '
        'The source docstrings define every returned field, type, axis, denominator and ownership rule.</p>' +
        ''.join(live_doc(name, 'src/event_analysis/ensemble.jl') for name in
                ('align_pad_history', 'ensemble_pad_maps', 'write_ensemble', 'read_ensemble')),
        'contract-ensemble')
    contracts += section('Extensive docstrings as part of implementation', '''
<p>Each public function receives: purpose and physical role; complete signature;
every argument’s meaning, units, type and shape; all returned fields and indexing;
algorithm and Petit origin; assumptions and boundaries; ownership and RNG effects;
validation errors and edge cases; and a small worked example with expected output.</p>
<p>Document helper functions in proportion to their role, especially coordinate
conventions, weighted graph costs and histogram boundaries. Check examples with
the implementation. Generate or quote interface excerpts from code once implemented
so documentation and signatures evolve together.</p>''', 'contract-docstrings')

    contracts += section('read_event · implemented; writer and orchestration · draft', code('''read_event(path::AbstractString, event_id::Int64)::Tuple{EventHits,EventInfo}
write_event_results(path::AbstractString, results; provenance)
''') + '''<p><code>read_event</code> adapts the existing reduced-HDF5 reader into
the implemented event types. It validates schema version, units, event index,
primary identity and equal hit-column lengths, then closes owned file handles.
Missing required time/primary fields identify a file that needs regeneration.</p>
<p>The writer consumes the event and pad tables specified in the HDF5 schema.
Its concrete aggregate result type is reviewed after the computational returns;
the signature above identifies the I/O boundary. Round-trip tests compare
values, types, axes and sample/event identity.</p>
<p>The orchestration function resolves configuration, constructs stage RNGs,
runs the documented branches and passes results to the writer. Its final
signature follows the reviewed parameter/result contracts. Python presentation
consumes that persisted schema. This keeps filesystem and execution choices
outside the computational functions.</p>''', 'contract-io')

    outputs = section('Configuration ownership', '''
<p>The detector card is authoritative for field, temperature, speed and geometric
pitch. The event-analysis card supplies analysis settings and sample paths.
Resolve cards once at the application boundary, validate values, convert units
and pass explicit typed arguments into computation.</p>
<p>The example below gathers approved settings. Its key organisation and output
paths are draft contracts. The existing nominal W-value is mapped from
<code>config/ion_tracks.toml: ionization.w_eV</code> to the proposed explicit key;
implementation should establish one authoritative W input and save the resolved value.</p>'''
        + code(CONFIG), 'output-config')
    outputs += section('Reproducible event execution', '''
<p>Save the run seed, RNG algorithm/version, sampled L, grid offsets and evaluated
widths/thresholds. Proposed stream organisation: derive independent streams from
the run seed, sample identity, event ID and stage (drift, diffusion, grid phase).
Specify the stable derivation algorithm before implementation. This supports
repeatable event-order/threading changes and comparisons using the same diffusion.</p>
<p>Reuse the saved diffusion when scanning graph parameters. Compare histogram
settings with equivalent sampled ion positions where practical, so the comparison
measures binning rather than a different random event.</p>''', 'output-reproducibility')
    outputs += section('Executable event-analysis card · steps 2–3', code((root / 'config/event_analysis.toml').read_text()) +
        '<p>The validation runner reads this card and converts the detector SI temperature/field inputs at the thermal-width interface. '
        'It uses ionisation_energy_eV for both count production and equivalent-energy reconstruction. '
        'The integrated runner uses the graph settings, including the exact inclusive 10-ion threshold. Pad pitch and speed come from the detector card; slice width, minimum span and time direction come from this card. '
        'Seeds and resolved inputs are saved in the validation HDF5. The earlier full-run example above includes later-stage draft fields.</p>', 'step2-card')
    outputs += section('HDF5 layout · proposed for interface review', code('''/provenance/
    schema_version, source_files_and_hashes, source_revision
    resolved_cards, julia_version, project_manifest_hashes
    petit_source_hashes, rng_algorithm, rng_seed
/samples/<sample_id>/events/
    event_id, n_hits, n_components, accepted, event_energy_MeV
    total_expected_ions, mean_ions_per_hit, sampled_ions, rounding_difference
    diffusion_retained_ions, diffusion_outside_ions
    graph_retained_ions, graph_below_threshold_ions
    drift_length_cm, sigma_t_mm, sigma_l_mm
    reconstruction_bin_width_mm, max_distance_mm, threshold_ions, knn_k
    grid_offset_x_mm, grid_offset_y_mm
/samples/<sample_id>/pad_counts/
    event_id, pad_ix, pad_iy, slice_index, ions
/samples/<sample_id>/axes/<event_id>/
    x_edges_mm, y_edges_mm, z_edges_mm, time_edges_ms
/samples/<sample_id>/paths/<event_id>/
    truth/, reconstructed/, mst_comparison/
/samples/<sample_id>/reductions/
    acceptance, common_xy_axes, common_z_axis, mean_xy_ions, mean_z_ions
    denominators, settings_scan, validation_metrics''') + '''
<p>Use columnar datasets for event and pad-count tables. Proposed storage types:
physical arrays Float32, ion counts and table indices Int64, primary/hit IDs Int32,
accepted Bool. Paths carry their units and ordered-bin indices. Keep sample identity
with event identity throughout. Store axes so omitted sparse rows imply zero counts.</p>
<p>Record all event decisions; accepted events have pad histories. Optional debug
products can retain intermediate tables for selected events. The exact schema,
chunking/compression and event index are reviewed with the I/O contract and measured
on the target files.</p>''', 'output-hdf5')
    outputs += section('Axis and boundary conventions to settle', table(['Contract', 'Proposed convention / decision needed'], [
        ('Axis representation · approved', 'F12: retain explicit Float32 edge vectors in function interfaces and the initial HDF5 schema. Compact storage is a possible later writer/reader optimisation.'),
        ('Pad coverage · implemented', 'pad_axes preserves the selected-track energy-weighted x/y reference plus sampled phase and the selected minimum-z reference. Whole-cell extensions cover the full retained table. The initial transverse span is 160 mm; all axes remain explicit Float32 vectors.'),
        ('Boundary assignment', 'Use half-open internal bins and include the final upper edge in the final bin. Test exact edges. Diffusion uses one adjacent-Float32-edge bin centred on a coincident zero-diffusion axis. Pad axes use the same final-edge inclusion and cover all retained diffusion centres.'),
        ('Relative time · implemented', 'The caller supplies drift_direction=+1 or -1. Saved time edges increase from zero at the first arrival edge; slice indices follow time order. The checkpoint card declares increasing source z as later. The opposite sign reverses the slice sequence. Absolute event arrival offsets remain a separate contract.'),
        ('Ensemble axes · working convention implemented', 'Central-pad alignment: the pad containing the selected-track reference has centre (0,0). First-arrival alignment defines the longitudinal origin. Shared axes cover all cells, at least 81×81 transverse pads. Exact integer numerators and both accepted/processed denominators are saved. Review the convention with the step-four report.'),
        ('Ion threshold · implemented', 'voxelize_event sums an optional Int64 ions column alongside Float32 energy. select_diffused_event compares these counts inclusively, then calls Petit graph construction with the filtered rows. This implements the configured threshold exactly; backbone energy weights retain Petit arithmetic.'),
    ]), 'output-open-contracts')
    outputs += section('Ensemble HDF5 · implemented schema v1', '''
<p><code>results/event_analysis/step4_scans.h5</code> stores the configuration,
source hashes, stage seeds, per-event decisions and charge ledger. Each sample
and variant has <code>events</code>, <code>aligned_histories</code> and
<code>ensemble</code> groups; sample summaries store paired selection changes
and pad-response distances.</p>
<p>The ensemble schema <code>nausicaa.event.ensemble.v1</code> saves Float32
x/y, arrival-distance and time edges; Int64 x/y and arrival count sums;
accepted/processed event counts and total ions; and Float32 means for both
denominators. Alignment and matrix-axis metadata accompany the arrays. Julia
uses x as the first matrix dimension; h5py exposes the corresponding y,x array.
The reader checks shapes, axes, count conservation and each saved mean against
its numerator and denominator. Aligned histories retain event ID and integer
pad_dx, pad_dy, arrival_bin and ions columns.</p>''', 'output-ensemble-schema')
    outputs += section('Interactive event viewer · interface review', """
<p>The agreed MC energy / MC tracks workflow, Petit defaults and proposed
interfaces are in <a href="#viewer">chapter 7 · Interactive event viewer</a>.
Review <a href="#viewer-energy-types">energy-map types</a>,
<a href="#viewer-graph-types">graph records</a>,
<a href="#viewer-functions">function contracts</a>,
<a href="#viewer-hdf5">HDF5 schema</a> and
<a href="#viewer-state">session behaviour</a> before implementation.</p>""", 'output-event-viewer')
    outputs += section('Generated scientific report', '''
<p>Julia writes numerical results and reductions to HDF5. Python reads them to
produce final x/y ion maps, z/time profiles, pad histories, acceptance plots and
tables. The generated LaTeX includes those figures and generated numerical macros,
along with settings, assumptions, charge accounting, convergence and performance.</p>
<p>Provide commands for analysis, presentation-only regeneration and the complete
pipeline. Define exact CLI names during the interface review. Compile the report
with latexmk and check its log. Numerical physics and reductions remain in Julia;
presentation can be refreshed from saved HDF5.</p>''', 'output-report')

    validation = section('Validation follows the adaptation', '''
<p>Copy the required Petit functions and test the changes around them. Preserve
source identity and use a pinned reference for regression. Repair obsolete fixture
calls in the local tests where needed, documenting the repair. The source appendix
records the limitations of Petit’s current active test file.</p>
<p>This validation plan and the draft contracts are ready for joint review.
Measured agreement, performance and final numerical tolerances will be supplied
by the implementation runs.</p>''', 'validation-principles')
    checkpoint = root / 'work/event_analysis/step1_validation.toml'
    if checkpoint.exists():
        record = tomllib.loads(checkpoint.read_text())
        counts = record.get('counts', {})
        rows = [(g['name'], str(g['passed']), str(g['failed']), str(g['errors'])) for g in record.get('groups', [])]
        validation += section('Step 1 · execution checkpoint',
            '<p>Core status: <strong>' + record['status'] + '</strong>; ' + str(counts.get('passed', 0)) +
            ' passing assertions. The v3 conversion retains MC time and actual primary IDs. Both samples are verified against their original HDF5 values and hashes.</p>' +
            table(['Test group', 'Passed', 'Failed', 'Errors'], rows) +
            '<p><a href="../../latex_generated/event_analysis_step1/report.pdf">Generated validation PDF</a> · '
            '<a href="../../results/event_analysis/step1_validation.h5">Numerical validation HDF5</a> · '
            '<a href="../../src/event_analysis/PROVENANCE.toml">Petit copy provenance</a></p>' +
            '<p>Pinned regression fixtures use the original Petit graph/backbone source with the root dependency manifest. '
            'Tie-free graphs and paths match; tied fixtures check valid nearest choices, components, charge and path observables. '
            'Ten real events exercise the deterministic stages on raw hits. Ensemble selection efficiencies belong to the diffusion stage.</p>' +
            '<p>Implemented return detail: backbone endpoints are independent NamedTuple copies. The empty raw-path case returns a typed empty table. '
            'The new typed reader requires v3; the existing dataframe reader also accepts v1/v2.</p>' +
            '<p>This is the completed step-one record. The step-two checkpoint below adds ion production and diffusion.</p>', 'step1-checkpoint')
    step2_record = root / 'work/event_analysis/step2_validation.toml'
    step2_summary = root / 'latex_generated/event_analysis_step2/summary.json'
    if step2_record.exists():
        record = tomllib.loads(step2_record.read_text())
        count = record.get('counts', {}).get('passed', 0)
        measured = ''
        if step2_summary.exists():
            summary = json.loads(step2_summary.read_text())
            measured = table(['Sample', 'Production events', 'Mean rounding [ions]', 'Max relative rounding [%]', 'Mean diffusion [ms]'], [
                (item['sample'], str(int(item['events'])), f"{item['mean_rounding_difference']:.2f}",
                 f"{100*item['max_absolute_relative_rounding']:.3g}", f"{1000*item['mean_seconds']:.2f}")
                for item in summary['samples']])
        validation += section('Step 2 · ion-production and diffusion checkpoint',
            '<p>Cumulative core status: <strong>' + record['status'] + '</strong>; ' + str(count) + ' passing assertions. '
            'The kernel retains Petit’s per-hit rounding, three Gaussian draws per ion and dense intermediate histogram. '
            'Pure functions receive explicit RNG, widths, counts and W; physical arrays are Float32 and counts are Int64.</p>' + measured +
            '<p><a href="../../latex_generated/event_analysis_step2/report.pdf">Generated step-two PDF</a> · '
            '<a href="../../results/event_analysis/step2_validation.h5">Validation HDF5</a> · '
            '<a href="../../config/event_analysis.toml">Executable TOML card</a></p>' +
            '<p>The pinned Petit replay agrees away from boundaries. Separate tests cover exact edges, ties-to-even rounding, '
            'independent SI thermal scaling, Gaussian moments/covariances, uniform sampling, reproducibility, zero longitudinal '
            'diffusion and the exact sampled = retained + outside ledger. Histogram scans reuse identical Gaussian draws.</p>' +
            '<p>Numerical choices: assign cells using saved Float32 edges; include the final upper edge. '
            'A coincident axis with zero diffusion has one bin centred at its coordinate with adjacent representable edges. '
            'An unresolved positive-width axis raises a precision error before sampling. The thermal formula retains Petit’s rounded kB/e constant.</p>' +
            '<p>Ion production is checked on both complete samples; diffusion is measured for five events per sample at three fixed '
            'distances plus one sampled distance. Step 3 adds a selection/pad validation subset; full-ensemble selection and pad distributions remain subsequent work.</p>' +
            '<p>Reproduce with <code>make validate-event-step2 PYTHON=/path/to/python3.11-or-newer</code>. '
            'Regenerate figures and LaTeX from HDF5 with <code>make report-event-step2</code>.</p>', 'step2-checkpoint')
    step3_record = root / 'work/event_analysis/step3_validation.toml'
    step3_summary = root / 'latex_generated/event_analysis_step3/summary.json'
    if step3_record.exists():
        record = tomllib.loads(step3_record.read_text())
        measured = ''
        if step3_summary.exists():
            summary = json.loads(step3_summary.read_text())
            measured = table(['Sample', 'Accepted / tested', 'Accepted retained ions', 'Pad ions', 'Below graph cut, saved', 'Mean pipeline [ms]'], [
                (item['sample'], f"{item['accepted_events']} / {item['tested_events']}", str(item['accepted_retained_ions']),
                 str(item['pad_ions']), str(item['accepted_below_graph_threshold_ions']), f"{1000*item['mean_seconds']:.1f}")
                for item in summary['samples']])
        validation += section('Step 3 · selection and pad-history checkpoint',
            '<p>Cumulative core status: <strong>' + record['status'] + '</strong>; ' + str(record.get('counts', {}).get('passed', 0)) +
            ' passing assertions. The full chain reads hits, rounds ion counts, samples L, diffuses, reconstructs and selects exactly one component, '
            'then assigns the accepted event’s full retained diffusion charge to pads.</p>' + measured +
            '<p><a href="../../latex_generated/event_analysis_step3/report.pdf">Generated step-three PDF</a> · '
            '<a href="../../results/event_analysis/step3_validation.h5">Numerical event histories and validation</a> · '
            '<a href="../../src/event_analysis/collection.jl">Selection and collection code</a></p>' +
            '<p>The sample is a deterministic validation subset spread across each input file. These acceptance fractions describe that subset; '
            'cuts were left unchanged. Full-sample efficiency claims follow the planned binning, threshold and graph scans.</p>' +
            '<p>Every accepted history satisfies sum(pad counts) = retained diffusion ions exactly, including charge below the graph threshold. '
            'Tests also cover Petit occupancy parity, disconnected-event rejection, boundary inclusion, complete coverage, integer-pitch shifts, '
            'zero-count intervals and both time directions. Saved HDF5 histories are read back to check the ledger.</p>' +
            '<p>Grid reference: selected-track energy-weighted x/y and minimum z; random x/y edge phase in one pitch. '
            'Minimum transverse span 160 mm, extended by whole cells as needed. The checkpoint uses increasing-z relative time, with t=0 at the first saved arrival edge; '
            'the sign and origin are explicit in the output. Both signs are implemented.</p>' +
            '<p>RNG streams are keyed by run seed, sample, event ID and stage using length-prefixed SHA-256 inputs. '
            'Each realised UInt64 seed is saved; processing order preserves the event draws. '
            'Reproduce with <code>make validate-event-step3 PYTHON=/path/to/python3.11-or-newer</code>; '
            'presentation-only regeneration uses <code>make report-event-step3</code>.</p>', 'step3-checkpoint')
    step4_summary = root / 'latex_generated/event_analysis_step4/summary.json'
    if step4_summary.exists():
        summary = json.loads(step4_summary.read_text())
        scan_rows, resolution_rows = [], []
        for sample in summary['samples']:
            for row in sample['rows']:
                pct = lambda value: '—' if value is None else f'{100*value:.2f}'
                scan_rows.append((sample['sample'], row['variant'],
                    f"{row['accepted']}/{row['processed']}", str(row['lost']), str(row['gained']),
                    pct(row['mean_shared_pad_distance']), pct(row['mean_fixed_grid_pad_distance'])))
                if row['variant'] in ('bins_200', 'bins_400'):
                    resolution_rows.append((sample['sample'], row['resolution_reference'], row['variant'],
                        str(row['resolution_shared']), pct(row['mean_resolution_fixed_grid_distance'])))
        validation += section('Step 4 · ensembles and paired scans',
            '<p>Cumulative core suite: <strong>' + str(summary['passed']) + ' passing assertions</strong>. '
            'Twelve settings use the same 50 events per sample, drift draws, diffusion seeds and grid phases. '
            'Graph-only variations reuse the diffusion table. Baseline settings remain unchanged.</p>' +
            '<p><a href="../../latex_generated/event_analysis_step4/report.pdf">Generated step-four PDF</a> · '
            '<a href="../../results/event_analysis/step4_scans.h5">Numerical scans and ensembles</a> · '
            '<a href="../../src/event_analysis/ensemble.jl">Ensemble implementation</a></p>' +
            '<p>The working convention uses central-pad and first-arrival alignment, with both accepted and processed denominators. '
            'Tests cover translation, full coverage, exact charge conservation, zero acceptance, Float32 means, '
            'HDF5 round trips and corrupted counts/means. Every saved ensemble is read back and checked against its histories.</p>' +
            '<p>Pad distance is D = Σ|A−B|/(ΣA+ΣB), averaged over shared accepted events. '
            'Full uses each variant’s grid reference and arrival origin. Fixed uses the baseline physical lattice '
            'and signed geometric cell labels, separating histogram changes from grid-reference shifts. '
            'Lost/gained counts compare event identities with baseline.</p>' +
            table(['Sample', 'Setting', 'Accepted', 'Lost', 'Gained', 'Full D [%]', 'Fixed D [%]'], scan_rows) +
            '<h3>Successive histogram resolution</h3>' +
            table(['Sample', 'Reference', 'Setting', 'Shared', 'Fixed D [%]'], resolution_rows) +
            '<p>The 400-bin check follows the measurable 100→200 response. Successive refinement roughly halves '
            'the fixed-grid distance; identical selected events across these three resolutions establish selection stability '
            'on this subset. The required pad-history accuracy remains a scientific review decision before production.</p>' +
            '<p>Reproduce with <code>make validate-event-step4 PYTHON=/path/to/python3.11-or-newer</code>. '
            'Regenerate presentation with <code>make report-event-step4</code>.</p>', 'step4-checkpoint')
    validation += section('Evidence required by stage', table(['Stage', 'Checks', 'Evidence saved'], [
        ('Input', 'NEXUS/reduced-file round trip; units, time, IDs, one/two primaries, event indexing and ownership.', 'Schema version, fixtures and conversion checks.'),
        ('Copied bin/graph algorithms', 'Same inputs/settings; bin sums, exact edges on tie-free fixtures; valid neighbour choices and components/charge/acceptance on tied fixtures; threshold equality and local/parent indices.', 'Regression results, tie-dependent differences and Float32 effects.'),
        ('Backbone/truth', 'Known geometry, reversal, branch conditions, length, primary identity, equal-time and empty paths.', 'Endpoint/path differences and examples.'),
        ('Diffusion', 'Thermal scaling; Gaussian means/variances/covariances; seeded repeatability and count ledger.', 'Statistical checks, seeds and histogram convergence.'),
        ('Selection/pads', 'Exactly-one selection; full accepted charge retained; cell boundaries, grid phase, z/time axes.', 'Acceptance and conservation summaries.'),
        ('Ensemble', 'Both samples; scans of k, threshold, constant bin/edge scale factors and histogram padding. Baseline factors are 3 and 6 relative to σt for all L.', 'Acceptance versus L, pad distributions, runtime and peak memory.'),
        ('Presentation', 'Values read from HDF5; figures/macros regenerated; LaTeX compilation.', 'Saved settings/provenance, PNGs, report and build log.'),
    ]), 'validation-matrix')
    validation += section('Development increments', table(['Increment', 'Deliverable'], [(str(i) + '. ' + title, body) for i, (title, body) in enumerate(MILESTONES, 1)]), 'validation-roadmap')
    validation += section('Review queue for the next increments', '''
<ol><li>Review the working central-pad and first-arrival alignment, and the two explicit denominators.</li>
<li>Review histogram sensitivity and set the required pad-history accuracy before choosing production resolution.</li>
<li>Review the ensemble interfaces and saved schema, then scope the full-sample run and performance measurements.</li></ol>
<p>Use the HTML as the shared review surface. Record each approval beside its
contract; implementation docstrings and measured evidence then replace the draft labels.</p>''', 'validation-review-queue')
    validation += section('Maintaining these documents', '''
<p><code>scripts/event_design_content.py</code> is the source for this design and
the short generated roadmap. <code>scripts/petit_review_content.py</code> holds
the Petit evidence and decision record. The builder embeds source excerpts,
configuration-derived examples and provenance in an offline HTML file.</p>
<p>Run <code>make event-review PYTHON=/path/to/python3.11-or-newer</code>, then
refresh this tab. Navigation, deep links, search and printing work locally.
The <a href="archive/event_analysis_plan_before_consolidation.md">previous plan</a>
and <a href="archive/review_before_consolidation.html">previous review</a> preserve
the historical discussion. The current design is authoritative for new work.</p>''', 'validation-docs')

    validation += section('External review and current discussion', '''
<p>The <a href="review_report_2026-10-27.md">review dated 2026-10-27</a> checks
the source claims and proposes refinements. Findings are discussed before
changing approved scientific or algorithmic choices.</p>'''
        + table(['Finding', 'Disposition / discussion'], [
            ('F1 · common drift distance', 'Resolved: retain the common per-event drift approximation; per-hit refinement declined for this study.'),
            ('F2 · 50 cm reconstruction breakpoint', 'Resolved: use Petit’s low-drift prescription b = 3σt and d = 6σt = 2b throughout the drift range. Both factors are configurable in TOML. This supersedes the fixed 2 mm reconstruction-bin proposal and the earlier piecewise rule. Physical pad pitch remains 2 mm. The old discontinuous rule remains documented as Petit source behaviour.'),
            ('F7 · neighbour/path ties', 'Resolved: retain Petit’s kNN algorithm. Require matching edges for distinct neighbour distances. With ties, verify valid nearest-neighbour choices, compare components, component charge, acceptance and backbone observables, and record how often these change. Investigate acceptance differences. Treat equal-cost backbone alternatives explicitly.'),
            ('F12 · scalar grid description', 'Resolved: retain explicit Float32 edge vectors for the computational interface and initial HDF5 axes. This preserves inspectable grids and the existing Petit approach. Compact uniform-axis storage may be considered later at the writer/reader boundary. Grid reference, time origin/direction and ensemble alignment remain separate contracts.'),
        ]) + '<p>F1, F2, F7 and F12 are resolved. Other review findings have not been adopted automatically.</p>', 'validation-external-review')

    return [
        dict(id='design', title='1. Purpose and physics', summary='The scientific question, physical model and approved design.', status='Approved model · interface review next', body=guide),
        dict(id='flow', title='2. Complete event flow', summary='From simulated hits to selection, pad histories and ensemble results.', status='Approved data flow · axis contracts identified', body=flow),
        dict(id='reuse', title='3. Reuse and development', summary='A focused Petit adaptation, with explicit inputs, outputs and provenance.', status='Reuse strategy', body=reuse),
        dict(id='contracts', title='4. Types and functions', summary='Review types and functions one at a time, with units, ownership and tests.', status='Steps 1–4 implemented · contracts and docstrings', body=contracts),
        dict(id='outputs', title='5. Configuration and outputs', summary='Resolved cards, reproducibility, HDF5 and generated scientific reports.', status='Proposed schema · approved settings', body=outputs),
        dict(id='validation', title='6. Validation and roadmap', summary='Evidence, development increments and the remaining review queue.', status='Steps 1–4 tested · production accuracy for review', body=validation),
        viewer_page(section, table, code),
    ]
