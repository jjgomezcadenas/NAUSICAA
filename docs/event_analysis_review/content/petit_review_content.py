"""Source-based Petit assessment and adaptation choices for joint review.

Keep observed behaviour separate from proposals. Source excerpts are extracted
by the builder from the local Petit checkout.
"""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]



def assessment_pages(src, section, table, code):
    module = 'src/Petit.jl'
    runner = 'scripts/itaca_track_reco_mt.jl'
    simple = 'itacaScripts/itaca_simple_reco.jl'
    voxels = 'src/voxels.jl'
    graphs = 'src/track_building.jl'
    extremes = 'src/track_extreme_finding.jl'
    paths = 'src/track_path_reconstruction.jl'
    mc = 'src/mctrack.jl'
    diffusion = 'src/itaca_functions.jl'
    tests = 'test/test_voxels_tracks.jl'
    ex = src.excerpt

    active = section('What “active” means here', '''
<p>This assessment follows the definitions loaded by <code>src/Petit.jl</code>
and the calls made by the ITACA entry points in this checkout. The detailed
reference is <code>scripts/itaca_track_reco_mt.jl</code>, which runs both MC and
diffused-track reconstruction. The simpler reconstruction script is also traced
for its kNN default and ion-per-pixel occupancy calculation. These are source call paths;
execution and numerical equivalence checks come after this review.</p>'''
        + table(['Source', 'Role in the inspected implementation'], [
            ('<code>src/Petit.jl</code>', 'Loads voxels.jl, track_building.jl, track_extreme_finding.jl, track_path_reconstruction.jl, mctrack.jl and itaca_functions.jl.'),
            ('<code>scripts/itaca_track_reco_mt.jl</code>', 'Batch event loop: original MC, graph-derived MC path, diffusion, reconstructed path, MST comparisons and saved analysis.'),
            ('<code>scripts/itaca_single_track_analysis.jl</code>', 'Interactive/plotting companion with the same main MC and reconstruction stages.'),
            ('<code>itacaScripts/itaca_simple_reco.jl</code>', 'Simpler event reconstruction: optional diffusion, voxelisation and selectable graph construction and single-track pixel occupancy; kNN is its default.'),
            ('<code>tracks.jl</code>, <code>tracks_updated.jl</code>, <code>srcbkup/</code>', 'Alternative/older implementations outside the module include chain used by these entry points. The split files above define the reference behaviour.'),
            ('<code>src/itaca_aux.jl</code>', 'Companion helper file. The inspected module loads the helper definitions in itaca_functions.jl.'),
        ]) + ex(module, 'include("histos.jl")', 26)
        + ex(runner, '# Load Petit module', 5)
        + ex('scripts/itaca_single_track_analysis.jl', '        mcvoxels = Petit.voxelize_event', 8)
        + ex(simple, '    method = "kNN"', 5)
        + ex(simple, '        if method == "kNN"', 19), 'active-entry')

    active += section('The event pipeline', '''
<div class="flow"><div><small>INPUT</small>NEXUS hits</div><span>→</span>
<div><small>BINNING</small>Occupied bins</div><span>→</span>
<div><small>CONNECTIVITY</small>Tracks</div><span>→</span>
<div class="result"><small>TRAVERSAL</small>Ordered path</div></div>
<p>Petit runs the binning–graph–path sequence on undiffused MC hits and repeats
it on the diffused distribution. It also constructs a time-ordered primary-particle
path as a truth reference.</p>'''
        + table(['Stage', 'What the batch implementation actually does', 'Data produced'], [
            ('Read event', 'get_dataset_dfs reads datasets under /MC; get_event selects one event into a DataFrame.', 'Event rows with positions, energy, time and particle/hit IDs.'),
            ('MC reference', 'compute_mc_path selects primary IDs 1 or 1 and 2, sorts each particle by time, groups consecutive visits to the same bin and joins the two arms when requested.', 'Truth-reference path with x, y, z, energy, cumulative length s and primary_electron.'),
            ('MC graph', 'voxelize_event bins all event hits. make_tracks applies a voxel-energy threshold, creates a radius graph and returns components ordered by energy. The batch script chooses mctracks[1].', 'Highest-energy connected MC track.'),
            ('MC backbone', 'walk_track_from_extremes invokes the combined endpoint/path algorithm. get_raw_path extracts the ordered bin centres and cumulative length.', 'Graph-derived MC path; MST path computed separately for comparison.'),
            ('Diffusion', 'transform_hits_df converts deposited energy to integer carrier counts. diffuse_xyz_image_mc samples Gaussian displacements and accumulates them in a finite 3D histogram.', 'Histogram-bin centres, carrier counts and equivalent energy.'),
            ('Reconstructed graph/path', 'voxelize_event rebins the diffused histogram; make_tracks uses the requested graph method. The batch pipeline retains events with exactly one reconstructed component, then extracts its path.', 'Reconstructed path plus a separate MST comparison.'),
            ('Analysis and output', 'Compare path endpoints to MC truth; compute blob observables. ItacaResult collects event results; the runner writes CSV summaries and HDF5 tracks/paths.', 'Saved analysis with per-event settings.'),
        ]) + ex(runner, '    # Step 1: Get event hits', 49)
        + ex(runner, '    # Transform hits', 46), 'active-pipeline')

    active += section('Binning convention and settings', '''
<p><code>voxelize_event</code> uses a zero-anchored cubic grid in the coordinates
it receives. On each axis it computes <code>floor(Int, x / voxel_size)</code>.
The representative position is the geometric bin centre
<code>(index + 0.5) * voxel_size</code>. It adds all deposited energies in each
occupied bin and returns those bins as rows of a DataFrame.</p>
<p>The graph-MC branch passes the original event coordinates directly. Grid
alignment is implicit in this convention. The MC-reference helper uses the same
integer-bin rule but represents each consecutive group by an energy-weighted
position; it serves a different purpose from the graph input.</p>'''
        + table(['Setting', 'Petit batch value / rule', 'Interpretation'], [
            ('MC bin width', '0.5 mm from get_voxel_size_and_distance', 'Reference value observed in this entry point.'),
            ('MC connection radius', '2.5 × MC bin width = 1.25 mm', 'Radius used by the MC graph call.'),
            ('MC voxel threshold', '1.0 keV', 'Applied after energies have been summed per bin.'),
            ('Reconstructed bin width', '3 σt for drift ≤ 50 cm; 1.5 σt above 50 cm', 'Diffusion-dependent binning in the batch helper.'),
            ('Reconstructed connection distance', '2 × bin width for drift ≤ 50 cm; 1.5 × bin width above 50 cm', 'Passed to the selected graph builder.'),
            ('Graph method', 'MC: KDT; batch RECO: graph_method argument, default KDT', 'The simpler ITACA script instead defaults to kNN (k = 10, distance cutoff = Inf).'),
        ]) + ex(voxels, '    # Compute voxel indices', 28)
        + ex(diffusion, 'function get_voxel_size_and_distance', 17), 'active-binning')

    active += section('Checked graph defaults by reconstruction entry point', table(
        ['Entry point', 'Reconstructed-data default', 'Selection'], [
            ('itacaScripts/itaca_simple_reco.jl', 'kNN, k = 10; maximum distance defaults to Inf.', '--method selects kNN or KDT; the chosen value reaches make_tracks.'),
            ('scripts/itaca_track_reco_mt.jl', 'KDT.', '--graph-method selects KDT, kNN or kNN_mutual for RECO. Its MC branch calls make_tracks with the KDT default.'),
            ('scripts/itaca_single_track_analysis.jl', 'KDT.', '--graph-method also supports the two kNN variants for RECO.'),
        ]) + '<p>Petit therefore uses kNN for reconstructed data in the simple reconstruction entry point, while the other inspected entry points default to KDT. The graph method is a caller setting.</p>'
        + ex(simple, 'function process_events', 23)
        + ex(simple, '        if method == "kNN"', 14)
        + ex(runner, '    graph_method = "KDT"', 3)
        + ex('scripts/itaca_single_track_analysis.jl', '        "--graph-method", "-g"', 8), 'active-graph-defaults')

    active += section('Graph and backbone rules', '''
<p>The radius builder keeps bins with <code>energy ≥ threshold</code>, connects
nearby bin centres using a KDTree, splits the graph into connected components,
and constructs one <code>Tracks</code> object per component. <code>make_tracks</code>
sorts these by deposited energy. The kNN variants use the same component/result
structure with a different edge rule.</p>
<p>The combined backbone finder uses the following decision order. Its
confidence values are algorithmic scores used by these branches.</p>'''
        + table(['Step', 'Decision in find_extremes_combined'], [
            ('Topology first', 'Use degree-one vertices and BFS paths; if topology confidence ≥ 0.85, return that path.'),
            ('Density split', 'Dense means average degree 2E/V > 6.0 by default.'),
            ('Dense graph', 'Try energy-weighted edge costs, then energy-based endpoints. An early result is accepted for confidence ≥ 0.85 or length > 1.1 × topology length. Then compare spatial and curvature candidates; an optional MST fallback participates before final length comparison.'),
            ('Sparse graph', 'Accept topology at confidence ≥ 0.8; otherwise use the curvature path if it is more than 1.1 × the topology path length.'),
            ('Weighted edges', 'Use distance^α / (mean endpoint energy + ε)^β. Weighted Dijkstra double sweep selects endpoints and a weighted-shortest path.'),
            ('MST comparison', 'Build a minimum spanning tree using Euclidean edge lengths, then double BFS on that tree; report the physical length of the returned path.'),
        ]) + ex(extremes, 'function find_extremes_combined', 100), 'active-backbone')

    active += section('Diffusion and the later detector stage', '''
<p>In the ITACA batch runner, drift distance is a per-event parameter: either
fixed or sampled uniformly between lmin and lmax. The runner computes one pair
of diffusion widths from that distance and passes them to the event processor.
The sampled distance sets the widths; this code leaves the hit coordinates in
their original frame. <code>transform_hits_df</code> changes table columns and
adds carrier counts; it performs no spatial translation.</p>
<p>The active diffusion routine samples one displacement per integer carrier:
independent Gaussian x/y displacements with σt and a z displacement with σl.
It uses a histogram covering the hit extent plus <code>nsigma × sigma</code>
padding. Only samples inside that histogram are counted. The resulting histogram
is then rebinned by the ordinary voxelisation routine.</p>
<p>The current ion branch of <code>get_sigma</code> sets σt from the thermal
formula and σl to zero. Its drift argument is in cm; returned widths are in mm.
The other loaded routine, <code>select_events_itaca</code>, instead sets
σl = σt for ions. The pixel/slice occupancy calculation in the simple
reconstruction script is described below.</p>'''
        + ex(runner, '        # Determine drift length and diffusion parameters for this event', 22)
        + ex(diffusion, 'function get_sigma', 15)
        + ex(diffusion, '    # For each hit, generate diffused electron positions', 30), 'active-diffusion')

    active += section('Ion diffusion widths: two existing caller choices', '''
<p><code>sigma_t_ion_mm(Tk, Lcm, Evcm)</code> evaluates
σ = 10 √[2 (kB/e) T Lcm / Evcm] in mm. This is the one-coordinate thermal
width for a singly charged ion in a uniform bulk field. The batch runner
assigns the same widths to every hit of an event using its sampled L.</p>'''
        + table(['Caller', 'Transverse width', 'Longitudinal width'], [
            ('get_sigma("ion", L), used by itaca_track_reco_mt.jl', 'Thermal formula', '0'),
            ('select_events_itaca', 'Same thermal formula', 'σl = σt; marked as an approximation in the source'),
        ]) + '''<p>Our detector card specifies T = 293.15 K and a uniform bulk
field of 200 V/cm. Evaluating the thermal formula with these design inputs gives:</p>'''
        + table(['Drift L', 'One-coordinate σ', 'σ / v at v = 0.1 mm/ms'], [
            ('10 cm', '0.503 mm', '5.03 ms'),
            ('100 cm', '1.589 mm', '15.89 ms'),
            ('150 cm', '1.947 mm', '19.47 ms'),
        ]) + '''<p>These are uniform-field bulk estimates for charge +e and
thermal diffusion. The time widths additionally assume σl = σt and use the
configured constant drift speed. The NAUSICAA detector card already identifies
thermal isotropic diffusion as its transport approximation.</p>'''
        + ex(diffusion, 'function sigma_t_ion_mm', 5)
        + ex(diffusion, '    # Ion tracks', 7)
        + ex(diffusion, 'function get_sigma', 15), 'active-ion-widths')

    active += section('How Petit chooses the drift distance L', table(
        ['Code path', 'How L is obtained', 'Geometry interpretation'], [
            ('scripts/itaca_track_reco_mt.jl', 'lmin/lmax are caller-supplied drift bounds in cm. Equal bounds give fixed drift; otherwise L = lmin + rand() × (lmax − lmin). Defaults are both 100 cm.', 'The caller chooses a range appropriate to the detector. This runner does not derive L from event coordinates or a detector geometry object.'),
            ('src/itaca_functions.jl: select_events_itaca', 'Sample levt uniformly between lmin + lbuff and lmax − lbuff; set lion = lmax − levt and lele = levt. Defaults: lmin = 0, lmax = 200, lbuff = 10; lengths enter diffusion formulas in cm.', 'With lmin = 0, lmax acts as the longitudinal detector extent; ions and electrons drift to opposite ends. These are configurable keyword arguments.'),
        ]) + '''<p>The second routine is called by the module’s ITACA loop in
<code>src/hdt5_mt.jl</code>. For our requested ion drift range, the first runner’s
existing sampler can be used with <code>lmin = 10 cm</code> and
<code>lmax = 150 cm</code>, consistent with the 1.5 m bulk reference in our
detector configuration. Choosing this L requires a change of settings, not a
new placement algorithm.</p>'''
        + ex(runner, '    # Validate lmin/lmax', 14)
        + ex(diffusion, 'function select_events_itaca', 35)
        + ex('src/hdt5_mt.jl', '            mc_tracks, ion_tracks, ele_tracks = select_events_itaca', 13), 'active-drift-range')

    active += section('Ions per pixel and longitudinal slice', '''
<p><code>itacaScripts/itaca_simple_reco.jl</code> already calculates pixel
occupancy for events with exactly one reconstructed track. Its sequence is
optional diffusion → voxelisation → graph construction → single-track occupancy.</p>'''
        + table(['Operation', 'Actual implementation'], [
            ('Grid', 'Centre a square pixel grid on the energy-weighted track barycentre. The caller defaults to a 160 mm grid; pixel width is max_distance, or voxel_size when max_distance is infinite.'),
            ('Longitudinal slices', 'Start at the track’s minimum z and use 2 mm slices. The time conversion appears only in a comment assuming 1 mm/ms; the calculation itself bins z.'),
            ('Ion assignment', 'For each retained track voxel, add energy_MeV × 10⁶ / 22 to the pixel containing its x,y centre in the current z slice. This is square-cell collection with fractional expected counts.'),
            ('Results', 'Return total, average, maximum and nonzero minimum pixel-count matrices over slices, occupancy per slice and signal summaries above 10 ions. Each slice matrix is temporary; the full pad-by-slice distribution is not returned.'),
            ('Collection model', 'The event calculation assigns charge geometrically to pixels. A separate archived scripts/scriptBkup/mipa.jl prototype traces fields from a hole to a disk; it is independent of this event occupancy calculation.'),
        ]) + '''<p>The occupancy calculation uses the selected track’s energy after
the graph’s bin-energy threshold. Its average matrix is an average over slices
of one event. Sampling an ensemble requires an additional event accumulation.</p>
<p>Concrete adaptation checks: use one configured W-value consistently (the
diffusion routine defaults to 40,000 carriers/MeV, equivalent to 25 eV, whereas
occupancy uses 22 eV); test grid coverage and upper-bin boundaries. With the
current half-open z slices, an exact integer-multiple z extent puts its maximum
on the excluded final edge. Preserve the per-slice counts for our output.</p>'''
        + ex(simple, 'const TIME_SLOT_MM', 6)
        + ex(simple, 'function compute_occupancy_analysis', 132)
        + ex(simple, '            # Use max_distance as pixel size', 6)
        + ex('scripts/scriptBkup/mipa.jl', '# Solve Laplace', 2), 'active-occupancy')

    active += section('Existing event selection and distance diagnostics', table(
        ['Existing Petit operation', 'Actual behaviour', 'Where it is used'], [
            ('is_single_track', 'Bins raw MC hits, builds kNN with k = 10, a configurable maximum edge distance and a 1 keV bin threshold; accepts exactly one connected component.', 'itacaScripts/filter_single_tracks.jl applies it after an event-energy window and writes the selected events. CLI defaults: 2 mm bins and 2 mm maximum edge distance.'),
            ('Batch RECO single-track selection', 'Rejects events when length(tracks) != 1 after diffusion, binning and graph filtering.', 'scripts/itaca_track_reco_mt.jl; its MC branch instead chooses the highest-energy component when at least one exists.'),
            ('diagnose_skeleton_coverage', 'For the chosen track’s bins, finds the nearest discrete backbone vertex using a KDTree; returns the fraction of track energy within R_cover, plus energy inside and total energy.', 'Loaded from src/track_analysis.jl. Inspected diagnostic callers print the result; the batch selection does not impose a coverage cut.'),
            ('compute_track_distances', 'Measures the closest pair of bin centres between the highest-energy track and each satellite component.', 'itacaScripts/itaca_simple_reco.jl prints and saves these distances as event statistics.'),
        ]) + '''<p>These operations already address track multiplicity, satellite
separation and energy near the backbone. In the inspected code, the specific
maximum-distance cut over every original hit against continuous backbone
segments is not applied. The agreed NAUSICAA baseline adopts the kNN single-component criterion and
applies it after diffusion, as a reconstructed-event selection.
Skeleton coverage and satellite separation are retained as diagnostics.</p>'''
        + ex('itacaScripts/filter_single_tracks.jl', 'function is_single_track', 24)
        + ex('src/track_analysis.jl', 'function diagnose_skeleton_coverage', 49)
        + ex(simple, 'function closest_voxel_distance', 22)
        + ex(simple, 'function compute_track_distances', 29), 'active-selection')

    interfaces = section('How to read these interfaces', '''
<p>These are Petit’s existing types and function contracts. Concrete declarations
and source excerpts come directly from the inspected files. Table schemas describe
the columns used by the active routines. Float64 is Petit’s representation;
our Float32 adaptation is described in the design and decision record.</p>
<p>Review in order: <a href="#type-events">event tables</a> →
<a href="#type-voxels">binned tables</a> → <a href="#type-tracks">Tracks</a> →
<a href="#type-coords">TrackCoords</a> → <a href="#type-paths">path results</a> →
<a href="#type-truth">MC reference</a> → <a href="#type-diffusion">diffusion</a> →
<a href="#type-tests">existing tests</a>.</p>''')

    interfaces += section('Event tables and event selection', code('''get_dataset_dfs(filename::String)
get_event(hitsdf::DataFrame, event_id::Int)''')
        + table(['Data / interface', 'Existing contract', 'Ownership / consumer'], [
            ('Event hits', 'DataFrame columns event_id, x, y, z, energy; time and particle_id support the MC truth path. Original MC tables also retain hit_id and label.', 'Positions in mm, energy in MeV, MC time in the source convention (ns for our NEXUS samples).'),
            ('get_dataset_dfs', 'Returns Dict{String,DataFrame}; reads groups or compound datasets below /MC.', 'Owns file I/O. Loading closes the HDF5 file after constructing tables.'),
            ('get_event', 'Groups by event_id and returns the selected event as a DataFrame; raises an error for an absent ID.', 'Constructs an event table. The voxeliser receives it unchanged.'),
        ]) + ex('src/util.jl', 'function get_dataset_dfs', 30)
        + ex('src/event_and_hits.jl', 'function get_event', 22), 'type-events')

    interfaces += section('Occupied-bin DataFrame and voxelisation', code('''voxelize_event(evtdf::DataFrame, voxel_size::Float64; ef::Float64=1e5/2.5)
voxelize_event(hitsdf::DataFrame, event_number::Int, voxel_size::Float64;
               ef::Float64=1e5/2.5)
voxelize_hits(hitsdf::DataFrame, voxel_size::Float64; ef::Float64=1e5/2.5)''')
        + table(['Output column', 'Type / units', 'Meaning'], [
            ('event_id', 'Int64', 'Event identifier repeated in each occupied-bin row.'),
            ('x, y, z', 'Float64 · mm', 'Geometric centres of occupied cubic bins.'),
            ('energy', 'Float64 · MeV', 'Sum of all deposited energies assigned to that bin.'),
            ('electrons', 'Int64 · carriers', 'round(total_e × ef); default ef = 40,000 carriers/MeV.'),
        ]) + '''<p>The single-event function converts positions and energies to Float64,
uses a dictionary keyed by integer (ix, iy, iz), and creates a new result table.
Output row order follows dictionary iteration. The all-events wrapper concatenates
single-event results. Bin width is passed as an argument and carried separately
in track metadata; the table itself stores centres and energies.</p>'''
        + table(['Worked example: width 1 mm', 'Result'], [
            ('Hits: x = [0.25, −0.25, 0.75, 1.0], all y = z = 0.25 mm', 'Three occupied bins with x centres −0.5, 0.5, 1.5 mm; y = z = 0.5 mm.'),
            ('Energies = [0.001, 0.002, 0.003, 0.004] MeV', 'Energy by ascending x centre: [0.002, 0.004, 0.004] MeV.'),
        ]) + '<p class="caption">The example is derived from the source rule and displayed in coordinate order for readability.</p>'
        + ex(voxels, 'function voxelize_event(evtdf', 39), 'type-voxels')

    interfaces += section('Tracks and graph construction', src.declaration(graphs, 'Tracks')
        + table(['Field', 'Existing meaning'], [
            ('voxels::DataFrame', 'Rows belonging to one connected component after energy filtering.'),
            ('graph::SimpleGraph{Int}', 'Undirected graph; local vertex i corresponds to voxels row i.'),
            ('components::Vector{Vector{Int}}', 'The builder stores [comp], with comp containing indices in the filtered parent table. The induced graph is locally renumbered.'),
            ('diffusion::DiffusionParams', 'Diffusion and binning settings carried with the track.'),
        ]) + code('''make_tracks(event_data::DataFrame;
    max_distance_mm::Float64=1.0, energy_threshold_kev::Float64=0.0,
    diffusion::DiffusionParams=DiffusionParams(), method::String="KDT", k::Int=10)
# returns Vector{Tracks}, ordered by decreasing summed energy

build_tracks_kdtree(event_data::DataFrame;
    max_distance::Float64=1.5, energy_threshold::Float64=0.0,
    diffusion::DiffusionParams=DiffusionParams())

build_tracks_knn(event_data::DataFrame;
    k::Int=10, max_distance::Float64=Inf, energy_threshold::Float64=0.0,
    diffusion::DiffusionParams=DiffusionParams(), mutual::Bool=false)''')
        + table(['Behaviour', 'Existing rule'], [
            ('Threshold units', 'make_tracks accepts keV and multiplies by 10⁻³. Lower-level builders take MeV. Equality passes.'),
            ('KDT', 'Add an edge between distinct bin centres returned by the radius query.'),
            ('kNN', 'Query k nearest other bins, subject to the maximum edge distance. Ordinary kNN accepts either directed neighbour relation; mutual kNN requires both.'),
            ('Empty / single bin', 'Return Tracks[] after empty filtering, or a one-vertex Tracks object for a single surviving bin.'),
            ('Ownership', 'Build filtered/component DataFrames and new graphs; preserve the input table. Fields within Tracks are mutable objects.'),
            ('Consumers', 'The batch MC branch takes the highest-energy track; the RECO branch requires exactly one track.'),
        ]) + ex(graphs, 'function build_tracks_kdtree', 57)
        + ex(graphs, 'function make_tracks(event_data', 57), 'type-tracks')

    interfaces += section('TrackCoords and endpoint/path interfaces', src.declaration(extremes, 'TrackCoords')
        + '''<p><code>extract_coords(track)</code> collects each position column into
separate vectors for repeated access during graph traversal. Vertex indices
refer to these vectors and to the track’s bin rows.</p>'''
        + code('''find_track_extremes(track::Tracks;
    use_energy_weighting::Bool=true,
    use_edge_energy_weighting::Bool=true,
    use_mst_fallback::Bool=false,
    dense_track_threshold::Float64=6.0)
# returns (extreme1, extreme2, path::Vector{Int}, confidence)

find_extremes_mst_diameter(track::Tracks, coords::TrackCoords)
# returns (extreme1, extreme2, path, confidence, path_length)''')
        + '''<p>The public combined finder returns four values. Candidate helpers
also report physical path length as a fifth value. An empty graph yields missing
endpoints, an empty path and zero confidence; a single vertex returns that vertex
as both endpoints. The algorithm constructs temporary coordinate arrays and
paths while preserving the supplied graph and table.</p>'''
        + ex(extremes, 'function find_track_extremes(track::Tracks', 29)
        + ex(extremes, 'function extract_coords', 7)
        + ex(extremes, 'function find_extremes_mst_diameter', 25), 'type-coords')

    interfaces += section('Walk result and ordered path', code('''walk_track_from_extremes(track::Tracks;
    use_energy_weighting::Bool=true, use_edge_energy_weighting::Bool=true,
    use_mst_fallback::Bool=false, dense_track_threshold::Float64=6.0)

get_raw_path(track, path_indices::Vector{Int})
compute_extreme_distances(reco_path::DataFrame, mc_path::DataFrame)''')
        + table(['Return object', 'Fields and interpretation'], [
            ('Walk named tuple', 'extremes: two DataFrameRow views into track.voxels; path_indices: ordered local vertex IDs; path_voxels: selected rows in path order; total_length: mm; confidence: algorithm score.'),
            ('get_raw_path DataFrame', 'x, y, z in mm and s, the cumulative Euclidean path length in mm. The caller supplies a nonempty path; the batch runner checks this first.'),
            ('Endpoint comparison', 'Named tuple d1, d2, total in mm and pairing (:direct, :crossed or :none). Chooses the pairing with smaller summed endpoint distance; empty inputs yield NaN distances.'),
            ('project_voxels_to_path', 'Returns the arc length of the nearest discrete path point for each bin; used as an analysis projection.'),
        ]) + '''<p>The agreed selection uses exactly one connected component.
The earlier original-hit distance cut is deferred.</p>'''
        + ex(extremes, 'function walk_track_from_extremes(track::Tracks', 50)
        + ex(paths, 'function get_raw_path', 23)
        + ex(paths, 'function compute_extreme_distances', 50), 'type-paths')

    interfaces += section('Time-ordered MC reference', code('''compute_mc_path(event_df::DataFrame, mcvox_size::Float64;
                is_double_beta::Bool=false)
voxelize_particle_hits(hits_df::DataFrame, mcvox_size::Float64)''')
        + table(['Step', 'Existing behaviour'], [
            ('Primary selection', 'Single electron: particle_id == 1. Double beta: IDs 1 and 2. Missing requested primaries produce an empty path table.'),
            ('Ordering', 'Filter each primary into a new table and sort it in place by time. The original event table is preserved.'),
            ('Consecutive grouping', 'Group successive hits whose fixed-grid bin key is equal. A later revisit to a bin is a separate group.'),
            ('Position and energy', 'Within each group, sum energy and compute an energy-weighted position; use the mean position if energy is zero.'),
            ('Two-arm path', 'Reverse the first primary sequence, then append the second. Add cumulative Euclidean length s and primary_electron labels.'),
        ]) + '''<p>This is the truth-history comparator used by the batch analysis.
The graph-derived MC backbone is constructed from the whole event through
<code>voxelize_event → make_tracks → walk_track_from_extremes</code>.</p>'''
        + ex(mc, 'function voxelize_particle_hits', 58)
        + ex(mc, 'function compute_mc_path', 87), 'type-truth')

    interfaces += section('DiffusionParams and diffusion interfaces', src.declaration(graphs, 'DiffusionParams')
        + table(['Fields', 'Units / use in this pipeline'], [
            ('ldrift', 'The type docstring says mm, but the ITACA runner passes ldrft in cm. Its result field is explicitly ldrft_cm. Agreed adaptation: name and document the supplied unit explicitly, for example drift_length_cm.'),
            ('sigma_t, sigma_l, voxel_size, max_distance', 'mm; transport widths and graph/binning parameters.'),
            ('energy_threshold', 'keV in track metadata.'),
            ('nbins_df, nsigma_df', 'Histogram bins per axis and padding in units of sigma.'),
        ]) + code('''transform_hits_df(df::DataFrame; energy_to_electrons::Float64=1e5/2.5)
diffuse_xyz_image_mc(df::DataFrame;
    sigma_t_mm=0.6, sigma_l_mm=0.6, nbins=100, nsigma=3.0,
    ef::Float64=1e5/2.5)
get_sigma(particle_type, ldrft; dt=3.5, dl=0.9,
          tK=297.0, edrift=500.0, Pbar=15.0)''')
        + '''<p><code>transform_hits_df</code> returns a new table with integer carrier
counts derived from energy. <code>diffuse_xyz_image_mc</code> returns a new
DataFrame with event_id, x, y, z, energy and electrons; positions/energy are
Float64 and counts are integers. It consumes the default RNG through
<code>randn()</code>. The finite histogram records accepted counts and converts
them back to MeV using <code>ef</code>.</p>
<p>The convolution-based <code>diffuse_xyz_image_kernel</code> is inside a block
comment in this file. The inspected runner calls the Monte Carlo routine.</p>'''
        + ex(diffusion, 'function transform_hits_df', 5)
        + ex(diffusion, 'function diffuse_xyz_image_mc', 108)
        + ex(runner, '    dfpars = Petit.DiffusionParams', 1), 'type-diffusion')

    interfaces += section('Pixel occupancy interface', code('''compute_occupancy_analysis(track::Petit.Tracks;
    pixel_size_mm::Float64, nbins::Int=50)''') + '''
<p>The function is defined in <code>itacaScripts/itaca_simple_reco.jl</code>.
It reads <code>track.voxels</code> columns x, y, z and energy, and returns a
NamedTuple. Count matrices are Float64 arrays of size nbins × nbins;
occupancy_per_slot is a Float64 vector. The tuple also carries the barycentre,
z extent, number of slices, pixel width and scalar occupancy/signal summaries.</p>
<p>For our pad histories, retain the intermediate per-slice counts and their
axes as an explicit output. The existing <code>avg_counts</code> is the sum
divided by this event’s number of slices.</p>''', 'type-occupancy')

    interfaces += section('ItacaResult and output boundary', '''
<p><code>ItacaResult</code> is defined in the batch script. It groups the event ID,
thread ID, drift/diffusion/binning settings, MC and reconstructed Tracks objects,
five path tables (truth, MC, MC-MST, RECO, RECO-MST), lengths, blob measurements
and endpoint errors. It is the orchestration result consumed by the CSV/HDF5
writers.</p>
<p>Its path and track members show which computational outputs the existing
analysis needs. Our initial event-selection result should be scoped from those
outputs and the approved kNN single-track condition.</p>'''
        + ex(runner, 'struct ItacaResult', 74), 'type-output')

    interfaces += section('Existing tests: available fixtures and current limits', '''
<p><code>test/runtests.jl</code> currently includes only
<code>test_voxels_tracks.jl</code>. The other listed suites are commented out.
This assessment inspected the tests; it did not run the Petit suite.</p>'''
        + table(['Fixture family', 'What its assertions cover'], [
            ('Voxelisation', 'Energy aggregation, bin centres, negative coordinates, event separation and comparison of single/all-event wrappers.'),
            ('Graph construction', 'Components, threshold filtering and track-energy ordering.'),
            ('Endpoints / walks', 'Straight and branched tracks, turns, path membership, BFS and lengths.'),
            ('Blob measurements', 'Energy in endpoint spheres and overlapping/large-radius cases.'),
        ]) + '''<p>Several assertions still target older interfaces: the Make Tracks
block passes <code>max_dist</code>/<code>energy_thr</code>, while the active wrapper
expects <code>max_distance_mm</code>/<code>energy_threshold_kev</code>; endpoint
fixtures construct Tracks with three fields, while the loaded struct has four.
One voxel test expects 60,000 carriers for 0.6 MeV although the current default
factor gives 24,000. These specific mismatches must be reconciled when reusing
the fixtures as compatibility evidence.</p>'''
        + ex('test/runtests.jl', 'using Petit', 18)
        + ex(tests, '    voxels_combined =', 11)
        + ex(tests, '    tracks = make_tracks(test_data', 3)
        + ex(tests, '    linear_track = Tracks', 5), 'type-tests')

    proposal = section('Adaptation principle and review order', '''
<p>Start with the loaded Petit functions and their current data flow. Preserve
the binning convention, graph builders, component ordering and combined path
algorithm as the reference. Review each change below before writing the adapted
implementation.</p>
<p>The approved decisions below supply the current design. Continue with the
<a href="#contracts">adapted interface drafts</a> and
<a href="#validation">validation plan</a> before implementation.</p>''', 'proposal-order')

    proposal += section('Decision record and implementation status', '''
<p>We are reviewing the proposals one by one. Rows marked <strong>approved</strong>
record the agreed design. Steps 1–4 are implemented and tested; see the
<a href="#step4-checkpoint">execution checkpoint</a> for current evidence and
the pad-response accuracy decision before production.</p>
<ul><li>Explicit units in names and documentation: agreed.</li>
<li>Input boundary: approved; retain EventHits and adapt to Petit DataFrames.</li>
<li>Voxelisation and binning: approved; adapt Petit’s fixed
binning convention, return a DataFrame, and validate Float32 against the Float64 reference.
For event selection, bin the diffused distribution.</li>
<li>Carrier conversion: approved; calculate expected ion counts from the configured
ionisation energy in the ion-production stage, retaining fractional expectations.</li>
<li>Track representation: approved; retain Petit’s Tracks structure and energy
ordering, with Float32 numerical columns and documented index conventions.</li>
<li>Graph algorithms: approved; retain KDT, kNN and mutual kNN with TOML settings.
Use kNN with a finite edge-distance cutoff for the single-track selection baseline,
following filter_single_tracks.jl. Retain the other methods as comparison options.</li>
<li>Backbone extraction: approved; retain Petit’s combined endpoint/path algorithm,
its existing branch conditions and energy weighting, with MST as a separate comparison.
Validate Float32 decisions against the reference on shared fixtures.</li>
<li>MC truth comparator: approved; retain time ordering, consecutive-bin grouping
and two-arm joining, using primary IDs from EventInfo. Use it to assess the
graph-derived backbone of the diffused event, with both paths expressed in a common frame.</li>
<li>Whole-event selection: approved with corrected order; set drift parameters and diffuse the event,
then bin it, build kNN and accept exactly one connected component. Retain skeleton coverage and satellite distances
as diagnostics. The earlier maximum hit-to-backbone distance cut is deferred.</li>
<li>Pad histories: approved; retain every pad’s ion counts by longitudinal slice,
with pad and slice indices and axis information. Occupied entries can be stored
in HDF5; totals and other summaries are derived from these histories.</li>
<li>Pad pitch: approved; configure the physical pitch (2 mm for this study)
independently of the graph’s maximum connection distance in TOML.</li>
<li>Time bins: approved; convert longitudinal bin widths using
Δt = Δz / v_ion with the configured ion drift velocity. Save both z-bin edges
and corresponding time-bin edges in the output.</li>
<li>Common W-value: approved; use the configured ionisation_energy_eV for
both diffusion carrier conversion and pad counting.</li>
<li>Transverse grid offset: approved; keep event coordinates unchanged and shift
the pad grid by independent uniform x and y offsets in [0, pitch_mm). Use one
offset pair for every slice of an event and save it with the event.</li>
<li>Diffusion widths: approved; use isotropic thermal diffusion with configured
temperature and bulk field and one sampled drift distance per event. Retain
zero longitudinal diffusion as a configurable comparison.</li>
<li>Diffusion sampling: approved; retain Petit’s per-ion Gaussian Monte Carlo
as the initial reference, with an explicit RNG and the common configured W-value.
Benchmark runtime and memory on representative events.</li>
<li>Integer carrier conversion: approved; round each hit’s expected ion count
for Monte Carlo sampling, preserve fractional expectations and record the
event-level rounding difference.</li>
<li>Intermediate diffusion histogram: approved; retain Petit’s procedure with
nbins and nsigma in TOML. Start comparisons at 100 bins per axis and 3σ padding;
measure convergence, runtime, memory and generated/retained/out-of-range counts.</li>
<li>Reconstruction binning: F2 resolved with b = 3σt and maximum edge distance
d = 6σt = 2b for the full drift range, with factors in TOML. This supersedes the
fixed 2 mm reconstruction-bin proposal and Petit’s 50 cm switch. Physical pad
pitch remains independently fixed at 2 mm.</li>
<li>Graph-input threshold: approved; start at 10 ions per reconstruction bin,
configurable in TOML, using the common W-value for the equivalent energy threshold.
Retain bins at or above threshold and measure removed charge and acceptance.</li>
<li>kNN neighbour count: approved; start at k = 10, configurable in TOML,
preserving Petit’s small-input handling and undirected graph construction.
Check single-track acceptance versus k.</li>
<li>Backbone representation: approved; retain ordered local bin indices,
the x/y/z/s path DataFrame in mm with Float32 columns, total length and the
existing algorithm confidence score.</li>
<li>Pad-count input: approved; select events using the thresholded graph, then
fill pad histories from the accepted event’s full intermediate diffusion histogram.
The graph threshold is controlled by TOML and preserves subthreshold charge in
the pad-count input.</li>
<li>F7 validation: approved; retain Petit’s kNN. Require matching edges on
tie-free fixtures; with ties, check valid neighbour choices and compare components,
charge, acceptance and backbone observables. Record and investigate consequential
differences.</li>
<li>F12 axis representation: approved; retain explicit Float32 edge vectors in
computational interfaces and the initial HDF5 schema. Compact uniform-axis
storage remains a possible later I/O optimisation.</li></ul>''', 'proposal-decisions')

    proposal += section('Corrected event-selection order', '''
<p>The source reference for the ordering is the reconstructed-event branch in
<code>scripts/itaca_track_reco_mt.jl</code>: set event drift/diffusion parameters,
convert energy to carriers, diffuse, voxelise, build tracks, require exactly
one component, then extract the path. The graph method is selectable there;
our agreed choice is kNN. The raw-MC <code>filter_single_tracks.jl</code> is a
separate prefilter, not the event-preparation sequence being adopted.</p>
<p><strong>Event hits → set drift/diffusion parameters → diffusion → voxelisation →
kNN graph → exactly-one-component selection → backbone extraction.</strong></p>
<p>Diffusion changes bin occupancy and graph connectivity, so the single-track
condition is evaluated on the diffused event. The undiffused MC path provides
truth-reference information. The approved transverse alignment keeps event
coordinates unchanged and shifts the pad grid within one pitch.</p>''', 'proposal-pipeline')

    proposal += section('What to take; what to change', table(
        ['Component', 'Take from Petit', 'Proposed change and reason'], [
            ('Input boundary · approved', 'Event DataFrame with coordinates, energy, time and particle identity; original NEXUS as source.', '<strong>Agreed:</strong> keep the existing Julia/HDF5 reader and reduced hit files; add MC time and primary-particle IDs to the reduced schema. Keep EventHits as the event input and adapt it to Petit’s DataFrame columns (x_mm → x, energy_MeV → energy, time_ns → time), preserving values and units.'),
            ('Voxelisation · approved', 'floor(coordinate/bin_width), geometric bin centres, dictionary energy sums, occupied rows.', '<strong>Agreed, with corrected pipeline order:</strong> adapt voxelize_event’s fixed binning rule to the coordinates supplied. For event selection, its input is the diffused distribution. Return a DataFrame with event ID, bin-centre coordinates and summed energy, matching Petit’s graph input. Use Float32 and validate numerical differences against Petit’s Float64 results. Document units explicitly, preserving mm, MeV and ns.'),
            ('Binning settings · approved', 'One width supplied to voxelize_event; implicit zero-anchored boundaries.', 'Accept the bin width in mm and use Petit’s fixed binning convention. Pass bin_width_mm as a parameter.'),
            ('Carrier conversion · approved', 'Petit derives carriers from energy using ef.', '<strong>Agreed:</strong> keep deposited energy as the input to binning and graph construction. In the ion-production stage calculate mean_ions = energy_MeV × 10⁶ / ionisation_energy_eV using the configured ionisation energy. Retain fractional expected counts; any later integer sampling is a separate explicit step.'),
            ('Track type / graph · approved', 'Tracks groups a bin table, SimpleGraph, component information and diffusion/binning settings; make_tracks returns components sorted by energy.', '<strong>Agreed:</strong> retain this structure and decreasing-energy ordering, with Float32 numerical columns and documented units. Local graph vertex i refers to row i of the track’s voxels; stored component indices refer to the filtered parent table. Preserve and document both conventions. Read graph settings from TOML cards and pass them explicitly to functions.'),
            ('Graph algorithms · approved', 'KDTree radius graph and the existing kNN/mutual-kNN variants.', '<strong>Agreed:</strong> retain all three implementations with method and parameters selected through TOML. Use kNN as the default for single-track selection, following filter_single_tracks.jl, with its finite maximum edge distance. Retain KDT and mutual kNN as comparison options. This updates the earlier KDT baseline choice after tracing the dedicated selection script.'),
            ('Backbone · approved', 'Combined topology/energy/spatial/curvature decision tree; walk result and raw path; separate MST comparison.', '<strong>Agreed:</strong> retain Petit’s combined endpoint/path algorithm, including its existing branch conditions and energy weighting, with MST as a separate comparison. Validate whether the Float32 adaptation preserves endpoint and path decisions on shared fixtures.'),
            ('Truth comparator · approved', 'Time ordering, consecutive-bin grouping and two-arm joining in compute_mc_path.', '<strong>Agreed:</strong> retain the procedure, including energy-weighted positions for consecutive groups, and read primary IDs from EventInfo instead of assuming IDs 1/2. Use this trajectory as the truth reference for assessing the graph-derived backbone, whose input is the diffused event. Express both paths in a common frame for comparison.'),
            ('Whole-event selection · approved', 'The ITACA batch RECO branch diffuses, voxelises, builds tracks using the selected method, and requires exactly one component before path extraction. The separate raw-MC filter demonstrates the kNN single-track criterion.', '<strong>Agreed:</strong> follow the existing RECO ordering with kNN. Carry bin width, k, finite maximum edge distance and energy threshold in TOML. The raw-MC filter’s 2 mm bins, 2 mm edge cutoff and 1 keV threshold describe that prefilter; the approved diffused-event settings are recorded below. Retain skeleton coverage and main-to-satellite distances as diagnostics. The earlier maximum hit-to-backbone distance cut is deferred.'),
            ('Drift sampling and grid offset', 'Petit samples or fixes one drift distance per event and uses it to set common diffusion widths. Its drift bounds are configurable; another ITACA routine derives ion/electron drift distances from a sampled longitudinal position and detector extent.', 'Reuse the existing uniform drift sampler with lmin = 10 cm and lmax = 150 cm for our requested range. Take the range from configuration, use an explicit RNG and save sampled values. Approved transverse alignment: keep event coordinates unchanged and shift the pixel-grid edges by independent uniform x and y offsets in [0, pitch_mm). Use the same pair for all slices and save it with the event. See the pad-count decision below.'),
            ('Diffusion · approved', 'Existing thermal width formula and per-carrier Monte Carlo diffusion.', 'Use isotropic thermal widths from configured detector inputs, an explicit RNG and the common W-value. Retain per-hit rounding and the intermediate histogram with charge accounting. Benchmark runtime/memory and histogram convergence before considering another algorithm.'),
            ('Presentation', 'HDF5 persistence already separates saved numerical results from plots.', 'Use the agreed Julia → HDF5 → Python → generated LaTeX/PDF workflow.'),
        ]), 'proposal-components')

    proposal += section('Drift sampling and approved pad-count decisions', '''
<p>Reuse Petit’s drift-distance sampler and pixel occupancy calculation.
For common per-event diffusion widths, sampling L supplies the transport input;
the hit coordinates can remain in their source frame. The pixel grid can be
defined around the selected track as Petit already does.</p>
<p><strong>Approved: retain every pad’s z-slice ion counts.</strong> Save the
counts with their event, pad and slice indices and axis information in HDF5.
Occupied entries suffice to represent the full sequence, including zero-count
slices through the saved axes. Derive totals, averages and extrema from these
histories. The sequence supplies the charge distribution needed for pad signals.</p>
<p><strong>Approved: configure pad pitch independently of the graph connection
distance.</strong> The physical pitch is 2 mm for this study. TOML carries the
pad pitch and maximum graph edge distance as separate parameters.</p>
<p><strong>Approved: derive time bins from the configured ion drift velocity.</strong>
Use Δt = Δz / v_ion and save both longitudinal and corresponding time-bin edges.
For example, 2 mm corresponds to 20 ms at 0.1 mm/ms.</p>
<p><strong>Approved: use a common configured W-value throughout diffusion and
pad counting.</strong> Both stages use <code>ionisation_energy_eV</code> for
energy-to-ion conversion.</p>
<p><strong>Approved: shift the pad grid within one pitch.</strong> Petit centres
each pixel grid on the track barycentre. Shift those grid edges by independent
uniform offsets δx and δy in [0, pitch_mm), with pitch_mm = 2 mm for this study.
Keep event coordinates unchanged. Use the same offset pair for every slice of
the event and save the offsets with its results. An integer-pitch shift
reproduces the pattern with different pad labels.</p>''', 'proposal-pad-counts')

    proposal += section('Diffusion widths: approved', '''
<p>Retain Petit’s thermal formula with temperature and bulk field from
<code>config/detector.toml</code>, and one sampled drift distance per event.
Use σx = σy = σz = σ, following <code>select_events_itaca</code> and the
existing detector card’s isotropic thermal approximation. Preserve the batch
runner’s σz = 0 choice as a configurable comparison.</p>
<p><strong>Review F1 resolved:</strong> retain one common drift distance and width
per event. The proposed per-hit drift refinement is declined for this study.</p>
<p>Width tests: verify σ² = 2(kB/e)TL/E independently; check square-root scaling
with L and T, inverse-square-root scaling with E, zero width at zero drift,
unit conversion and both longitudinal settings.</p>''', 'proposal-diffusion-widths')

    proposal += section('Diffusion sampling: approved', '''
<p>Petit’s active <code>diffuse_xyz_image_mc</code> takes an integer carrier
count per hit, samples independent Gaussian displacements in x, y and z for
each carrier, and accumulates the resulting positions in a 3D histogram.
When deriving counts from energy, Petit rounds each hit’s expected count.</p>
<p><strong>Approved:</strong> retain this per-ion Monte Carlo procedure as the initial reference,
with an explicit RNG and the approved common W-value. Benchmark representative
events for runtime and memory before deciding whether an alternative is needed.
Fractional expected counts remain available as the charge reference. The approved
integer conversion and histogram procedure are detailed below.</p>
<p>Sampling validation should compare displacement means, variances and cross-axis
covariances with the configured Gaussian distribution, check seeded reproducibility,
and account for generated counts and counts outside the histogram.</p>''', 'proposal-diffusion-sampling')

    proposal += section('Integer carrier conversion: approved', '''
<p>Retain Petit’s per-hit rounding rule: compute mean_ions = energy_MeV × 10⁶ / W_eV,
then n_ions = round(Int, mean_ions) for Monte Carlo sampling. For example,
120.4 expected ions gives 120 sampled ions and 120.7 gives 121.</p>
<p>Preserve the fractional expectation alongside the integer count. Save the
event-level difference between summed integer counts and summed expectations,
and evaluate its size on the input samples. This makes the effect of rounding
explicit while retaining Petit’s algorithm.</p>''', 'proposal-ion-rounding')

    proposal += section('Intermediate diffusion histogram: approved', '''
<p>Petit extends the event’s hit bounds on each axis by nsigma times the diffusion
width, divides each padded axis into nbins bins, and counts sampled ions in the
resulting dense 3D array. It returns occupied bin centres, counts and equivalent
energy; the graph voxeliser subsequently bins that table.</p>
<p>The batch runner defaults to nbins = 100 and nsigma = 3. The simple
reconstruction script uses 300 and 3. These settings define the intermediate
diffusion resolution; physical pad pitch is configured separately.</p>
<p><strong>Approved:</strong> retain this histogram procedure as the initial reference, with
nbins and nsigma in TOML. Start comparisons from the batch defaults, measure
runtime and memory, and check convergence of track selection and pad histories
as bin count and padding increase. Record generated, retained and out-of-range
ion counts so their accounting is explicit. Review final settings using those
results.</p>''', 'proposal-diffusion-histogram')

    proposal += section('Reconstruction binning: revised after F2', '''
<p>After diffusion, the batch runner uses <code>get_voxel_size_and_distance</code>
to choose the graph-input bin width b and maximum graph edge distance d:</p>'''
        + table(['Sampled drift distance', 'Bin width b', 'Maximum edge distance d'], [
            ('L ≤ 50 cm', '3σt', '2b'),
            ('L > 50 cm', '1.5σt', '1.5b'),
        ]) + '''<p>The table records Petit’s source behaviour. <strong>The revised
baseline applies the low-drift prescription throughout: b = 3σt and d = 6σt = 2b.</strong>
Both constant factors are configured in TOML. This supersedes the fixed 2 mm
reconstruction-bin proposal and the earlier approval of the piecewise rule.</p>
<p>Save evaluated reconstruction settings per event. The physical pad pitch
remains independently configured at 2 mm.</p>''', 'proposal-reco-binning')

    proposal += section('Graph-input ion threshold: approved', '''
<p>The batch runner defaults to a threshold of 10 ions per reconstruction bin.
<code>get_energy_threshold</code> converts that to keV and the graph builder
retains bins whose energy is greater than or equal to the threshold. Petit’s
hard-coded conversion of 40 ions/keV makes this 0.25 keV.</p>
<p><strong>Approved:</strong> retain the 10-ion starting threshold, configured as
<code>[graph] min_ions_per_bin = 10</code> in TOML. Use the approved common W-value for conversion:
E_threshold_keV = threshold_ions × ionisation_energy_eV / 1000.
For W = 22 eV this is 0.22 keV. Apply the threshold after reconstruction
voxelisation and before kNN construction. Measure the removed charge and the
single-track acceptance as the threshold is varied.</p>'''
        + ex(diffusion, 'function get_energy_threshold', 12), 'proposal-graph-threshold')

    proposal += section('kNN neighbour count: approved', '''
<p>Petit’s batch runner and simple reconstruction default to k = 10.
<strong>Approved:</strong> retain this starting value, configurable in TOML. Each retained bin considers
up to ten nearest other bins and adds connections within the approved maximum
edge distance. With fewer than eleven bins, Petit uses k_actual = min(k, N − 1).</p>
<p>Ordinary kNN adds an undirected edge when either endpoint selects the other;
the final degree can therefore exceed k. Use this existing rule for the approved
single-component selection, and assess acceptance as k varies.</p>''', 'proposal-knn-count')

    proposal += section('Backbone path representation: approved', '''
<p>Retain Petit’s ordered local vertex indices, which identify rows in the
selected track’s voxel table, and its path DataFrame with x, y, z and cumulative
distance s. Distances are in mm; adapt numerical columns to the agreed Float32
representation. The sequence follows the backbone from one endpoint to the other.</p>
<p>Retain the walk result’s total length and confidence score. The score is the
existing algorithm’s diagnostic. Validate the row/index correspondence,
cumulative Euclidean distances and total length on known paths, including
single-bin and empty-path handling.</p>''', 'proposal-path-representation')

    proposal += section('Input to pad counting: approved', '''
<p>Petit passes the selected <code>Tracks</code> object to
<code>compute_occupancy_analysis</code>, which reads <code>track.voxels</code>.
Consequently its pixel counts use reconstruction-bin centres and include only
the energy that survives the graph-input threshold.</p>
<p><strong>Approved change for physical pad histories:</strong> use the graph and its threshold
to accept or reject the event, then fill pad/slice counts from the accepted
event’s full intermediate diffusion histogram, before reconstruction rebinning
and thresholding. Retain the selected track to define the grid reference and
calculate topology diagnostics. This preserves charge below the graph threshold
and the intermediate histogram’s spatial resolution in the pad calculation.</p>
<p>The graph threshold is supplied through TOML. The existing charge accounting
records losses outside the finite diffusion histogram separately.</p>''', 'proposal-pad-input')

    proposal += section('Data representation: approved decisions and remaining review', '''
<p>Petit’s bin table already stores separate x, y, z and energy columns. Its
<code>TrackCoords</code> already uses separate coordinate vectors. These match
the representation you requested without requiring a new graph abstraction.</p>
<p>The approved input approach uses a DataFrame adapter from EventHits to the
existing functions. Voxelisation returns a DataFrame of occupied bins and their
energies; the bin width is passed as a parameter. This supersedes the earlier
VoxelGrid/VoxelCloud type proposals. Retaining Petit’s Tracks representation is
approved, together with its ordered path indices and x/y/z/s path table.</p>
<p>Float32 changes must reach allocations, keyword annotations, coordinate
helpers, reductions and edge weights. Simply supplying Float32 columns to Petit
would still trigger its explicit Float64 conversions. The comparison baseline
therefore uses original Petit behaviour, followed by the proposed precision
changes as a separately measured step.</p>''', 'proposal-data')

    proposal += section('Compatibility assessment of implemented group one', table(
        ['Implemented item', 'How it fits Petit', 'Proposed treatment'], [
            ('CoordinateFrame', 'Petit uses the input coordinates without a frame tag.', 'Retain the explicit tag as a useful NAUSICAA addition. The first topology baseline uses source_frame.'),
            ('EventInfo', 'Petit’s truth comparator selects hard-coded primary IDs; event information can supply actual IDs and the generation vertex.', 'Retain event-level information for truth comparison. Confirm the reader supplies the agreed fields.'),
            ('EventHits', 'Its x_mm/y_mm/z_mm/energy_MeV/time_ns vectors map directly to Petit’s x/y/z/energy/time columns; identifiers are retained.', 'Approved: keep it as the event input and adapt to Petit’s DataFrame layout.'),
            ('center_hits', 'Petit’s MC graph path bins original coordinates. Arbitrary translation before fixed-grid binning can change bin membership.', 'The helper remains available. The proposed Petit-based pipeline uses source coordinates, sampled drift and track-relative pixel-grid edges; it requires no call to center_hits.'),
            ('Copying and validation', 'Petit generally builds new tables/graphs, but the walk endpoints are DataFrameRow views.', 'Review ownership at each boundary. Preserve the agreed pure-function behaviour; explicitly address views when defining the adapted return type.'),
        ]) + '''<p>The four items remain implemented in <code>src/event_analysis/</code>.
Their existing contract tests concern that implementation; they do not establish
compatibility with Petit. Their adapted integration is documented in the
<a href="#contracts">interface drafts</a>.</p>''', 'proposal-group-one')

    proposal += section('Reference tests against Petit: proposal for review', '''
<p>Compare each adapted stage with the loaded Petit implementation using the
same input and settings. Use fixed small fixtures with known geometry: a straight
track, curved track, two arms and separated components. Compare bin energies,
graph components, backbone observables and lengths. Require matching graph edges
for distinct neighbour distances. With ties, validate the neighbour choices and
record effects on components, charge, acceptance and paths, following approved F7.
Account for row ordering and equivalent path reversal when comparing outputs.</p>
<p>Supply identical diffused tables to graph tests so graph comparisons isolate
the reconstruction. Test the stochastic diffusion stage through seeded
reproducibility and its Gaussian/count statistics. Measure Float32 differences
against the Float64 reference, including changes near graph thresholds.</p>
<p>Approved changes such as isotropic widths, a common W-value and full-diffusion
pad input have their own physical expectations and charge-accounting tests.
After these small tests, validate acceptance, pad histories and performance on
representative events from both HDF5 samples. This validation approach awaits review.</p>''', 'proposal-reference-tests')

    proposal += section('Validation strategy for the adaptation', table(
        ['Step', 'Comparison / evidence'], [
            ('Reference fixtures', 'Repair the specific obsolete calls in copied/adapted Petit fixtures. Exercise the loaded functions in isolation from plotting and the batch command-line runner.'),
            ('Binning parity', 'Compare cell keys, geometric centres and per-bin energy sums after ordering rows by cell key. Include negative positions, exact boundaries and neighbouring Float32 values. Supply identical coordinates to the reference and adapted voxelisers.'),
            ('Graph parity', 'Use identical bin tables/settings and account for vertex labels. Require matching kNN edges for distinct neighbour distances. With ties, validate neighbour choices and compare components, component charge and event acceptance; quantify and investigate consequential differences. Radius-graph edge checks also cover distance-cutoff precision.'),
            ('Backbone parity', 'Compare endpoints, ordered vertex paths and physical lengths on straight, curved, two-arm, branched and dense fixtures. Treat path reversal as equivalent where orientation is arbitrary; record genuine tie-dependent differences.'),
            ('Float32 effects', 'Compare energy errors against Float64 references and the 0.0125 MeV resolution scale. Separately check graph edges and path branch decisions near thresholds.'),
            ('Single-track selection', 'Compare the kNN single-component criterion using identical diffused input tables and settings. Include one connected track, two separated components, empty/fully filtered events, threshold equality and finite-edge-cutoff cases. Compare coverage and satellite distances as diagnostics.'),
            ('Representative samples', 'Run both original NEXUS samples through the reference and adapted stages; report graph/path differences and single-track acceptance versus bin width, neighbour count, edge cutoff and energy threshold. Keep settings and source provenance with HDF5 results.'),
        ]) + '''<p>Numerical compatibility and acceptance results will be generated
by execution after interface approval. This page records the proposed evidence,
not measured agreement.</p>''', 'proposal-validation')

    proposal += section('Software scope and source provenance', '''
<p>Keep the new work in <code>src/event_analysis/</code>, using the root Julia
project. Reuse only the required loaded Petit functions and their MIT notice.
The root project now includes Graphs, NearestNeighbors, LinearAlgebra,
SparseArrays and Statistics alongside the existing table/I/O dependencies.
Manifest.toml pins the tested package versions.</p>
<p>Document each adapted function with its originating Petit file, retained
algorithm, changes, units and tests. Keep computational functions separate from
HDF5 access, TOML parsing and presentation. Random sampling receives an explicit
RNG. Adapted types and signatures are drafted in <a href="#contracts">chapter 4</a>.</p>
<p><a href="previous_interfaces.html">Earlier interface pages</a> are retained
as <strong>superseded pending review</strong>. The earlier plan’s API sketches
are historical until reconciled with this assessment.</p>''', 'proposal-scope')

    return [
        dict(id='active', title='1. Active Petit implementation',
             summary='The loaded code, its call paths and the behaviour we can reuse.',
             status='Source assessment · ready for joint review', body=active),
        dict(id='interfaces', title='2. Petit types and interfaces',
             summary='Existing data structures, arguments, results and test fixtures.',
             status='Existing Petit contracts · source excerpts', body=interfaces),
        dict(id='proposal', title='3. Adaptation proposal',
             summary='What we take, what we change and how we establish compatibility.',
             status='Proposal · awaiting prerequisite review', body=proposal),
    ]
