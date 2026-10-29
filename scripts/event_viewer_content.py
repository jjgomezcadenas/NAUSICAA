"""Viewer implementation status and reviewed interfaces for the HTML guide."""


VIEWER_CARD = '''# config/event_viewer.toml — applied MC defaults
[samples]
bb0nu = "events/0nubb_hitdf.h5"
bb0nu_100um = "events/0nubb_100um_hitdf.h5"
electrons = "events/electron_0nubb_energy_hitdf.h5"

[display]
nbins = 100                     # Petit histogram edge-selection target
extent_scale = 2.0              # Frame spans twice the track extent
colormap = "viridis"            # Petit palette
colour_scale = "linear"
voxel_opacity = 1.0             # Viewer setting, adjustable
auto_rotate = false

[mc_graph]
voxel_width_mm = 0.5
method = "KDT"
distance_mode = "voxel_multiple" # or "fixed_mm"
distance_voxel_multiple = 2.5
max_distance_mm = 1.25          # Used in fixed_mm mode
energy_threshold_keV = 1.0
k = 10                         # Used by kNN / kNN_mutual

[session]
file_directory = "events"
initial_sample = "bb0nu"
restore = true
work_dir = "work/event_viewer"

[server]
host = "127.0.0.1"
port = 0                       # Select an available local port
open_browser = true
'''


def viewer_page(section, table, code):
    body = section('Step 3 · MC graph inspection', '''
<p><strong>The MC energy viewer is implemented and accepted after browser review.</strong>
Select double beta or single electrons; use Return, Next, Previous or an event
ID. The four views show XY/XZ/YZ summed-energy heatmaps and occupied 3D energy
voxels. Step 3 adds the MC tracks tab, component table, graph layers, backbones
and endpoints. Settings are read from TOML and displayed with their applied
values; interactive parameter editing is step 4.</p>
<p>The user confirmed event 2 after switching to independent axis scales in
all four views. Projection zoom and pan are independent.</p>
<p>The user approved the step-one interfaces as implementation guidance and
requested application-level validation on real events. Existing reconstruction
library tests remain the numerical reference.</p>
<p>Launch from the repository root: <code>make view-events</code>. It selects
the installed local environment in <code>work/event-viewer-venv</code>, opens
the browser and prints its local URL. Stop with Ctrl-C. For an explicit port:
<code>work/event-viewer-venv/bin/python scripts/view_events.py --port 8766</code>.
Use <code>--no-browser</code> to open the printed URL yourself, or
<code>--work-dir work/another-viewer</code> for an independent session.</p>
<p>The first event warms Julia’s reader and map calculations; subsequent events
reuse the process. Session state and the current numerical HDF5 snapshot are in
<code>work/event_viewer/</code>. Reloading the page restores the current view;
restarting the process restores the sample, event, tab, bin target, camera,
opacity, selected component and layer choices.</p>
<p><a href="../../src/event_analysis/display.jl">Energy reduction and docstrings</a> ·
<a href="../../src/event_analysis/mc_inspection.jl">MC graph adapter and docstrings</a> ·
<a href="../../src/event_viewer/server.py">Viewer server</a> ·
<a href="../../config/event_viewer.toml">Executable viewer card</a></p>''', 'viewer-review')

    body += section('Location and responsibility', code('''src/event_analysis/
  display.jl                   # Implemented pure energy-map reduction
  mc_inspection.jl             # Pure component/backbone extraction adapter
src/event_viewer/
  worker.jl                    # Persistent Julia reader, cache and HDF5 writer
  server.py                    # Local server and navigation/session state
  presentation.py              # HDF5-to-Plotly figures
  web/index.html
  web/viewer.js
  web/viewer.css
scripts/view_events.py         # Launch command
config/event_viewer.toml       # Resolved viewer defaults
tests/event_viewer/            # Worker/presentation/navigation integration
tests/event_analysis/          # Numerical helper tests
work/event_viewer/             # Request files, snapshots and session.toml''') + '''
<p>Entry point: <code>python scripts/view_events.py</code>, with optional
<code>--config</code>. The server starts one Julia process and opens one browser
page. Julia owns reading, binning, summaries, graphs and paths. Python handles
navigation and plotting. Browser code controls layers, camera and pending edits.
All project paths resolve against the repository root.</p>
<p>Dependencies: StatsBase 0.34.13 is resolved in Manifest.toml; Plotly 6.9.0 is
pinned in requirements.txt and installed in the local viewer environment. The
Plotly JavaScript bundle is served locally. A fresh Python environment can be
prepared with <code>python -m pip install -r requirements.txt</code>; Julia uses
the root Project/Manifest. Existing graph and HDF5 dependencies are reused.</p>''', 'viewer-location')

    body += section('Defaults and their origin', table(['Setting', 'Default', 'Origin'], [
        ('Display bin target', '100 per axis', 'Petit src/plots.jl: plot_event → plot_hits → StatsBase fit(Histogram, …; nbins). Actual edges/counts are saved.'),
        ('Palette / axis extent', 'Viridis / 1.3 × coordinate span', 'Petit src/plots.jl.'),
        ('MC voxel width', '0.5 mm', 'Petit src/itaca_functions.jl: get_voxel_size_and_distance.'),
        ('MC graph / distance', 'KDT / 2.5 × voxel width = 1.25 mm', 'Petit scripts/itaca_single_track_analysis.jl: make_tracks call and its default method.'),
        ('MC voxel-energy threshold', '1 keV, inclusive', 'Same active MC graph call.'),
        ('kNN neighbour count', '10 when a kNN method is selected', 'Petit make_tracks default.'),
        ('Initial component focus', 'Highest deposited energy', 'Petit MC branch selects the first energy-sorted component.'),
    ]) + '''
<p>Approved presentation changes: sum deposited energy in the projection pixels;
represent occupied 3D display voxels by their energy. Petit’s original 2D
histograms count entries, and its 3D panel is an energy-coloured hit scatter.
The new display uses mm, MeV, independently scaled axes and common data bounds
for each coordinate across views. The existing diffused-event kNN card retains
its own settings.</p>''' + code(VIEWER_CARD) + '''
<p>Parse physical settings to Float32 at the Julia boundary. Voxel width and
resolved distance must be finite and positive; threshold is finite and
nonnegative; nbins and k are positive integers. Methods are KDT, kNN and
kNN_mutual. In multiplier mode, resolve distance = width × multiplier; in fixed
mode, use max_distance_mm. Save the mode, entered values and resolved distance.
Reset to Petit defaults restores the numeric MC defaults shown above.</p>''', 'viewer-config')

    body += section('Data types · event and energy maps', '''
<p>The reader returns the existing <code>EventHits</code> and
<code>EventInfo</code>: source-frame Float32 x/y/z, MeV energy and ns time vectors;
Int32 particle/hit IDs; Int64 event ID, generation vertex and primary IDs.
The enclosing sample key identifies the file. Event ID and sample position are
separate: position is one-based in the file’s stored event index.</p>
<p>Implemented pure reduction output:</p>''' + code('''struct MCEnergyMaps
    x_edges_mm::Vector{Float32}
    y_edges_mm::Vector{Float32}
    z_edges_mm::Vector{Float32}
    xy_energy_MeV::Matrix{Float32}  # nx × ny; sum over z
    xz_energy_MeV::Matrix{Float32}  # nx × nz; sum over y
    yz_energy_MeV::Matrix{Float32}  # ny × nz; sum over x
    xy_hit_count::Matrix{Int64}    # nx × ny; contributing hits
    xz_hit_count::Matrix{Int64}    # nx × nz
    yz_hit_count::Matrix{Int64}    # ny × nz
    occupied::DataFrame
    total_energy_MeV::Float32
    hit_count::Int64
end

# occupied columns, one row per cell containing at least one hit:
# ix, iy, iz::Int64              one-based indices into the saved axes
# energy_MeV::Float32           summed energy in that cell
# hit_count::Int64              contributing hits''') + '''
<p>Each axis is selected once using Petit’s StatsBase bin-target convention.
The resulting axes define both projections and sparse occupied 3D cells.
Final edges cover every hit; the Float32 conversion must preserve increasing
edges and coverage. A coincident-coordinate axis uses one finite cell of
width 0.5 mm, recorded explicitly. Display padding changes plot limits.</p>
<p>Sum occupied voxel energies to obtain the three projections. Internal bins
are half open; the last bin includes its upper edge. Empty cells display as
background, while occupied zero-energy cells remain in the saved table. Hover
reports bin bounds, summed energy and contributing hit count. Original hit
identifiers remain in the snapshot’s hit columns.</p>
<p>Output arrays are newly owned and inputs retain their values. Float32
reductions are checked against independent Float64 references. The sum of each
projection and of occupied energies agrees with the event energy to the tested
rounding tolerance. 3D data remain sparse to keep event-sized storage.</p>''', 'viewer-energy-types')

    body += section('Data types · implemented MC graph settings and components', code('''struct MCGraphSettings
    voxel_width_mm::Float32
    method::String
    distance_mode::String
    distance_voxel_multiple::Float32
    max_distance_mm::Float32       # fixed-distance entry
    energy_threshold_keV::Float32
    k::Int
end

struct MCComponent
    component_id::Int64            # rank by decreasing component energy
    nodes::DataFrame               # x_mm/y_mm/z_mm/energy_MeV::Float32
    edge_src::Vector{Int64}        # local node indices, one-based
    edge_dst::Vector{Int64}        # one undirected edge per pair; src < dst
    path_indices::Vector{Int64}    # ordered local node indices
    endpoint_indices::Union{Nothing,NTuple{2,Int64}}
    path_s_mm::Vector{Float32}     # same length as path_indices
    energy_MeV::Float32
    length_mm::Float32
    confidence::Float32           # existing backbone algorithm score
end''') + '''
<p>MC graph input is the existing voxelize_event DataFrame, with Float32
x/y/z/energy and Int64 event_id. The settings describe its voxel width. Graph
construction applies the inclusive bin-energy threshold and returns every
connected component in Petit’s order. The viewer initially focuses component 1.</p>
<p>Adapt existing Tracks and backbone results into owned MCComponent arrays.
Local node i is row i of nodes; edges, path and endpoint indices all use this
same numbering. Path distance comes from get_raw_path. An empty backbone has
empty path arrays and endpoint_indices = nothing; a single-node path has both
endpoints at that node and length zero. An event with no surviving components
returns an empty vector and its energy maps remain available.</p>
<p>MCComponent is a presentation/persistence record derived from the existing
graph objects. The graph algorithm and endpoint/path extraction are the local
Petit adaptation. Properties include energy, node/edge counts, endpoints,
backbone length and the named algorithm score.</p>''', 'viewer-graph-types')

    body += section('Function contracts and tests', table(['Interface', 'Behaviour and output', 'Required evidence'], [
        ('<code>bin_mc_energy(hits::EventHits; nbins::Int=100)::MCEnergyMaps</code>',
         'Select shared display axes; sum energy and hit counts into occupied cells and projections; calculate the event summary. Preserve the input frame and arrays.',
         'Petit edge parity on representative ranges; negative/coincident coordinates and final edges; independently calculated cell assignments, energy sums and Float32 errors; input ownership.'),
        ('<code>resolved_mc_distance(settings::MCGraphSettings)::Float32</code>',
         'Validate settings and resolve multiplier or fixed-mm distance. Default result: 1.25f0.',
         'Changing width updates the multiplier distance; fixed distance remains constant; reject invalid modes and values.'),
        ('<code>voxelize_event(to_petit_dataframe(hits), settings.voxel_width_mm)::DataFrame</code>',
         'Reuse the existing pure voxelisation. The worker caches this table with the event identity and voxel width.',
         'Existing voxelisation tests plus cache invalidation when event or width changes.'),
        ('<code>build_mc_components(voxels::DataFrame; method::String, max_distance_mm::Float32, energy_threshold_keV::Float32, k::Int)::Vector{MCComponent}</code>',
         'Call make_tracks with explicit settings, then walk_track_from_extremes and get_raw_path per component. Rename columns at the presentation boundary and preserve local node numbering. Preserve the voxel table.',
         'Match existing graph edges, components, paths and statistics; empty/multiple/single-node components; threshold equality and valid indices.'),
        ('<code>write_viewer_snapshot(path::AbstractString, hits::EventHits, info::EventInfo, maps::MCEnergyMaps; context, graph=nothing)::Nothing</code>',
         'Worker-side I/O writes a completed snapshot atomically. Context and optional graph records have the explicit fields listed below. Close owned handles before announcing completion.',
         'Julia-write/Python-read field, type, shape and index parity; successful request identity; failed writes preserve the previous snapshot.'),
    ]) + '''
<p>Docstrings specify complete inputs and returned fields, units, shapes,
ownership, errors, Petit provenance and worked examples. Pure helpers live in
event_analysis; filesystem and session orchestration live in event_viewer.</p>
<p>Snapshot <code>context</code> fields: request_id::String, sample::String,
source_path::String, source_schema::String, source_size_bytes::Int64,
source_mtime_ns::Int64, event_position::Int64, event_count::Int64,
display_nbins::Int, resolved_card_toml::String and code_sha256::String.
The source fingerprint and any available upstream hash are recorded once on
opening the sample. The optional <code>graph</code> record contains
settings::MCGraphSettings, resolved_distance_mm::Float32, voxels::DataFrame,
components::Vector{MCComponent}, retained_energy_MeV::Float32 and
below_threshold_energy_MeV::Float32, calculated in Julia.</p>''', 'viewer-functions')

    body += section('Worker requests and completion', table(['Request', 'Inputs', 'Completed result'], [
        ('Open sample', 'request_id; configured sample key.', 'Sample metadata and ordered Int64 event-ID vector in an HDF5 catalogue. Python uses this index for navigation.'),
        ('Load event', 'request_id; sample; event_id; display nbins.', 'Fresh EventHits/EventInfo, energy maps and summary in the snapshot; graph status pending. Commit the selected event when complete.'),
        ('Reconstruct', 'request_id; sample; event_id; complete MCGraphSettings.', 'Graph bins, all component records, properties, applied settings and resolved distance in a completed snapshot.'),
        ('Update display binning', 'request_id; sample; event_id; display nbins.', 'New energy maps from the loaded hits; existing graph result remains valid.'),
        ('Close', 'request_id.', 'Close input handles and end the Julia worker.'),
    ]) + '''
<p>Use TOML request/response files in the session work directory, with a
line-based stdin/stdout notification carrying the request identifier. Julia’s
existing TOML library parses the files; diagnostic logs go to stderr. Each
response contains request_id, status, snapshot/catalogue path or an error code
and message. Python sends one request at a time and only publishes a completed
result matching the requested sample/event. The browser shows loading/error
state and the most recent successful result.</p>
<p>The browser talks to the local Python server. JSON carries browser commands
and plotted arrays; numerical persistence is HDF5. Encode Int64 event IDs as
decimal strings in browser messages to preserve their full range. The browser
uses the saved local graph indices when resolving line endpoints.</p>''', 'viewer-protocol')

    body += section('HDF5 snapshot · energy and graph inspection', code('''attributes: schema = "nausicaa.event.viewer.v1", status = "complete"
/context/                  request, source and sample metadata
/configuration/            resolved card; applied display/graph settings
/event/                    event_id, position, count, frame
/event/info/               vertex_mm, primary_particle_ids
/event/hits/               x_mm, y_mm, z_mm, energy_MeV, time_ns,
                           particle_id, hit_id
/event/summary/            hit_count, total_energy_MeV
/display/axes/             x_edges_mm, y_edges_mm, z_edges_mm
/display/projections/      xy_energy_MeV, xz_energy_MeV, yz_energy_MeV,
                           xy_hit_count, xz_hit_count, yz_hit_count
/display/occupied/         ix, iy, iz, energy_MeV, hit_count
/reconstruction/           status = pending | complete
                           build_request_id (when complete)
/reconstruction/voxels/    event_id, x_mm, y_mm, z_mm, energy_MeV
/reconstruction/summary/   component_count, retained_energy_MeV,
                           below_threshold_energy_MeV
/reconstruction/components/1/, /2/, ...
    nodes/                x_mm, y_mm, z_mm, energy_MeV
    edges/                src, dst
    backbone/             path_indices, s_mm, endpoint_indices
    summary/              energy_MeV, energy_fraction, node_count, edge_count,
                           length_mm, confidence''') + '''
<p>The snapshot contains context, configuration, event, display and source-code
provenance hashes. Reconstruction is pending until the MC tracks tab is first
opened, then includes all graph groups above. The worker caches that result
across tab changes and display rebins; build_request_id identifies its build.
Its writer is <code>write_snapshot</code> in worker.jl, with <code>save_graph</code>
for reconstruction. Applied settings and resolved distance are in
/configuration/mc_graph. Energy fractions use the full original event energy.
Sample position/count are in /event; the card and applied bin target are in
/configuration. Core library and worker source hashes are
recorded along with the source file size, modification time and schema.</p>
<p>Floats are Float32; counts, event IDs and persisted indices are Int64;
particle/hit IDs are Int32. Frame is the CoordinateFrame name. Every index
dataset declares one-based indexing. Endpoint storage is a length-two Int64
vector, or an empty vector for an absent path. Completed reconstruction with
zero components has component_count = 0 and empty component groups.</p>
<p>In Julia the projection shapes are (nx,ny), (nx,nz), (ny,nz). HDF5 as read by
h5py exposes (ny,nx), (nz,nx), (nz,ny); save axis-order attributes. Python verifies
these shapes and supplies them directly to the corresponding heatmaps. Cell
bounds derive from the saved edge vectors. Every panel labels energy units.</p>
<p>The writer finishes a temporary HDF5 file, closes it, then replaces the
current snapshot. Pending edits live in session state; the snapshot records
the settings that produced its arrays. Original hits, maps and graph records
share the snapshot’s event identity. Read-back tests verify that relationship.</p>''', 'viewer-hdf5')

    body += section('Session state and parameter updates', '''
<p>Implemented in step 3: active tab, component highlight, layer choices,
independent projection ranges, camera and opacity are saved with the selected
event. Entering MC tracks builds the graph once for the current event. Event
changes invalidate the graph; display bin changes reuse it. The parameter-edit
behaviour below describes step 4.</p>
<p>Save <code>work/event_viewer/session.toml</code> after successful navigation
or settings application. Fields: sample key, event ID, active tab, selected
component, display settings, last applied graph settings, pending graph edits,
and the current camera/zoom. Restore the event and settings on restart and
regenerate its numerical snapshot; identify the restored source by its saved
fingerprint. Runtime request IDs and loading flags describe the live process.</p>
<p>The graph panel displays both applied settings and a pending-edits badge.
Apply validates all edits together and reconstructs once. On success the new
result and applied settings become current; on failure the previous result
and pending edits remain visible. Reset stages Petit’s defaults; Apply executes
them. The distance control shows the resolved value alongside its mode.</p>''' + table(['Change', 'Work triggered', 'View behaviour'], [
        ('Voxel width + Apply', 'Rebin original hits, rebuild components and backbones.', 'Preserve event, camera and zoom; focus highest-energy component in the new result.'),
        ('Graph method/distance/threshold/k + Apply', 'Reuse current graph voxels; rebuild components and backbones.', 'Preserve event, camera and zoom; refresh the component table.'),
        ('Display nbins + Apply', 'Recompute display maps using original hits.', 'Update pixels/voxels; keep graph geometry and camera.'),
        ('Palette, opacity, visible layers, component focus, camera', 'Redraw browser figures.', 'Numerical arrays and graph settings retain their values.'),
        ('Next/Previous/event ID/sample', 'Load energy maps; reconstruct on entering the MC tracks tab using last applied graph settings.', 'Fit the new event extent; retain chosen layers and graph settings. Pending edits remain marked until Apply.'),
        ('Tab switch within an event', 'Build the graph on first entry if needed; reuse completed data thereafter.', 'Retain the event and available view state.'),
    ]) + '''
<p>Return advances outside text/number inputs; in the event-ID input it submits
that ID. At the sample boundaries disable the corresponding navigation button.
Controls show the sample position and actual event ID. The 3D camera supports
rotation, zoom, pan, standard axis views, reset and an optional Start/Stop orbit.
The component table highlights the same component across all four panels.</p>''', 'viewer-state')

    body += section('Using and checking the energy viewer', '''
<p>Enter a filename and click Open file, or use Browse to choose an HDF5 file.
The default directory is <code>events/</code>, configurable through
<code>session.file_directory</code>. Bare filenames resolve there; absolute
paths are also accepted. Browse supports directory navigation and starts in
the current file’s directory. The Samples / opened files menu retains the
configured samples and files opened during the session.</p>
<p>Reduced <code>_hitdf.h5</code> files load directly. Opening a NEXUS file
uses the existing Julia converter to create and verify its reduced table on
first use. Raw files live in <code>~/NAUSICAA/</code>; their reduced tables
are written to the project’s <code>events/</code> directory. Other inputs write
beside the source. Existing companions are checked against the source SHA-256
before reuse. A failed open keeps the current display; an incompatible or stale
companion reports regeneration guidance. Session restore reopens the selected
reduced file and event.</p>
<p>All hits are included. The bin target defaults to Petit’s 100; the header
shows realised cell counts and widths. Edit the target and click Apply to
change display binning. This is independent of the MC graph width.
XY/XZ/YZ use a shared linear energy scale; 3D uses MeV per occupied voxel.</p>
<p>Each projection has independent zoom and pan controls. Its ranges are saved
separately. Reset views fits the whole event. Drag the 3D scene to rotate,
use the wheel to zoom, select a standard axis view, or use Rotate/Stop for a
camera orbit. The opacity slider exposes overlapping voxels. Navigation fits
the new event; display rebinning preserves the current view.</p>
<p>Full-event framing includes distant clusters. Each 2D axis uses its own
padded coordinate bounds and fills its panel. The 3D scene maps the three
coordinate ranges to a cube. Axis labels and hover bounds remain in mm;
the footer identifies the independent axis scales. Each projection can be
zoomed separately to inspect the main track or a distant cluster.</p>
<p>Validation: the Julia core suite includes independent cell-assignment and
Float64 energy-sum checks. <code>make check-event-viewer</code> reads real events
from both files, verifies the HDF5 hit columns against their source rows,
compares projections with independent NumPy histograms, validates Plotly
figures and checks navigation, missing-ID handling and session restart.
Application acceptance is the hands-on browser review.</p>''', 'viewer-usage')
    body += section('Using and checking MC tracks', '''
<p>Open MC tracks to reconstruct the current event. The applied-settings line
shows voxel width, graph method, connection distance and energy threshold.
The table lists every connected component in decreasing energy order, with
event-energy fraction, node/edge counts, backbone length, endpoint coordinates
and Petit’s backbone score. Select a component number to highlight it across
all four panels; Highlight all restores equal emphasis.</p>
<p>Energy, Nodes, Edges, Backbone and Endpoints are independent layer switches.
Energy maps retain their energy colour scale. Graph colours identify components;
the selected component is emphasised while the others remain visible. Zoom and
pan remain independent between projections. Switching tabs retains the view.</p>
<p>Validation compares the adapter directly with make_tracks,
walk_track_from_extremes and get_raw_path on synthetic and real events. The
viewer integration check verifies one-based HDF5 indices, projected backbone
coordinates, energy accounting, graph reuse across display rebins, event changes,
session restore and a completed reconstruction with zero surviving components.</p>''', 'viewer-tracks-usage')
    body += section('Transported hits: Gaussian ion transport', '''
<p>The existing viewer’s Transported hits tab applies the implemented Petit
ion-production and Gaussian diffusion functions to each source event. One
uniform event drift in 10–150 cm sets the thermal width, using the detector
card’s temperature and bulk field. The baseline uses the same transverse and
longitudinal width. Fixed mode accepts a distance within those limits. The
prominent label shows distance in cm and both widths in mm.</p>
<p>XY/XZ/YZ and rotatable 3D show ion counts on the intermediate diffusion
histogram, before reconstruction cuts. Its 100 cells per axis and 3σ padding
are configured under <code>transport</code> in the viewer card. The visible
ledger lists produced, retained and outside-histogram ions. The MC display-bin
control applies to the two MC tabs. Histogram clipping and per-hit ion rounding
follow the existing transport calculation.</p>
<p>Apply transport evaluates the selected random/fixed mode and fixed distance.
New realisation advances the event’s realisation counter; in random mode it
draws a new distance too. A SHA-256-derived seed identifies the canonical input
path, event ID, configured base seed and realisation. Separate random streams
set the distance and ion displacements. Returning to an event or restarting
the viewer reproduces its realisation. Tab switches and view controls retain
the transport. New transport calculations use the original source hits.</p>
<p><code>src/event_viewer/transport.jl</code> adapts the tested Julia analysis
functions to the application. The snapshot’s <code>/transport</code> group saves
evaluated distance, widths, detector card, seeds, ion ledger, occupied bins and
exact Int64 projections. Python reads those arrays for colour maps and mesh
geometry. <code>tests/event_viewer/check_transport.py</code> checks both species,
projection counts, charge accounting, deterministic navigation, fixed-distance
width scaling, new realisations, invalid distances and session restart.</p>
<p>Review this transported-hit display before adding the reconstructed-track
tab with its separate voxelisation, graph and single-component selection.</p>''', 'viewer-transport')
    body += section('Pad charge and the three reconstruction stages', '''
<p>The Pad charge tab collects the same transported event on the detector’s
2 mm regular hole/pad lattice. Its reference is the full transported cloud’s
ion-weighted x/y centre and minimum z. One independent uniform x/y phase within
a pitch is fixed per event realisation. Axes extend to cover every retained
diffusion-bin centre. This uses the existing ideal geometric assignment to the
containing hole/pad cell.</p>
<p>XY displays each pad’s total collected charge in fC. X–time and Y–time sum
the opposite transverse coordinate. The 3D heatmap displays individual
pad/time cells. The saved history preserves integer ions per pad and slice;
charge uses +e per ion. Source-z slices are 2 mm, corresponding to 20 ms at
the detector card’s 0.1 m/s speed. Time is relative to the first arrival edge;
the drift sign and source-z edges are saved. The visible drift label refers
to the same Gaussian transport used in Transported hits.</p>
<p>All retained ions are collected, including low-count cells. Collection and
the transport realisation are shared across tabs, navigation and restart.
The existing pad/time history and covering-axis routines perform the calculation;
the new explicit-reference pad_axes method supports collection before track
selection. Python reads saved charge arrays to construct the plots. Framing
uses twice the occupied track extent, keeping empty grid margins outside the
initial view.</p>
<p>Keep three independent reconstruction observations: A on original MC,
B on the diffused cloud before pads, and C on pad/time discretised charge.
Selections at A/B do not gate transport or collection. The MC-track display
implements A; B/C graph displays remain the next increment. Existing selected
batch analyses retain their saved reference and selection conventions.</p>
<p><code>tests/event_viewer/check_pads.py</code> independently assigns saved
diffusion centres to pad/time cells, checks exact total-ion conservation and
charge projection values, 2 mm pitch, 20 ms slices, shared transport, random
realisations and restart for both species. The core collection tests check
covering axes, phase periodicity and charge conservation.</p>''', 'viewer-pads')
    body += section('Batch pixel charge per electronics frame', '''
<p><code>make pixel-charge</code> extends the batch runner through
<code>scripts/analyze_mc_selection.jl --pixel-charge</code>, using
<code>config/pixel_charge.toml</code>. It transports and collects the first
requested events of each 100 μm sample, limited by the available catalogue.
The request is set by <code>run.n_events</code>. All events reach collection. The shared
<code>collect_ion_cloud</code> function supplies the same full-cloud reference,
phase convention and geometric assignment as the interactive pad display.</p>
<p>A pulse is a positive pixel–frame cell. The inclusive histogram contains
all pulses. Event mean and maximum use all pulses with no threshold. Compute
an event’s true minimum over all pulses, then include the event in the minimum
histogram only if that minimum &gt; n_min when a minimum cut is requested.
The current cards use n_min = 0: all four distributions include every positive
cell/event with charge, with no charge cut. An event with a smaller
pulse is excluded from that histogram, rather than discarding its smaller pulses.
An empty minimum distribution is labelled explicitly.</p>
<p>The 2×2 figure shows inclusive, mean, maximum and minimum ions per pulse,
with logarithmic count and ion axes. Statistics use the same event populations
as their histograms. Exact event pulse sums and counts are saved. All-pixel
retention compares collected ions with produced ions; the minimum display cut
introduces no threshold loss. The reducer refreshes these quantities from saved
sparse histories without repeating transport or collection.</p>
<p>Save complete sparse event/pixel/frame counts, per-event axes, ion ledgers,
drift distances, seeds, grid phases, histograms, event pulse means, provenance and
timings in <code>results/pixel_charge/bb0nu.h5</code> and
<code>electrons.h5</code>. <code>scripts/present_pixel_charge.py</code> generates
one 2×2 PNG per species, <code>results/pixel_charge/pixel_charge.html</code>,
and <code>latex_generated/pixel_charge/pixel_charge.tex</code>/<code>.pdf</code>.
<code>make present-pixel-charge</code> refreshes presentation from saved results.</p>
<p><code>make pixel-charge-1mm</code> applies the same study to the original
1 mm-step inputs, using <code>config/pixel_charge_1mm.toml</code>. Its HDF5,
figures and HTML are separate under <code>results/pixel_charge_1mm/</code>,
with LaTeX/PDF under <code>latex_generated/pixel_charge_1mm/</code>.
Both cards record the MC generation step, and presentation accepts explicit
<code>--results-dir</code> and <code>--latex-dir</code> paths.</p>
<p>Diffused ion clouds are cached separately in <code>results/diffused_tracks/</code>.
Each event stores losslessly compressed Float32 diffusion-bin centres/energies,
Int64 ions, explicit diffusion edges, distance, width, seed and ion ledger.
Cache identity includes the source content/path, physical transport parameters,
histogram resolution/padding, random seed, Julia version and transport code.
Pitch, longitudinal bucket size and requested event count are independent.
The optional <code>--repixelize-only</code> flag refuses to generate diffusion
when an event is missing from the cache, guaranteeing a collection-only run.
Completed events are reused and additional events can extend the cache.</p>
<p>Set <code>pads.pitch_mm</code> and <code>pads.z_bin_width_mm</code> in the card
to re-pixelise. <code>config/pixel_charge_1mm_3mm.toml</code> and
<code>make pixel-charge-1mm-3mm</code> use 3 mm pitch and 3 mm longitudinal buckets
(30 ms). Outputs go to <code>results/pixel_charge_1mm_3mm/</code> and
<code>latex_generated/pixel_charge_1mm_3mm/</code>. Per-sample timings record cache
hits/misses and read/write time separately from newly computed diffusion.
<code>config/pixel_charge_1mm_4mm.toml</code> and
<code>make pixel-charge-1mm-4mm</code> use the same cache for 4 mm pitch and
40 ms frames, with separate <code>pixel_charge_1mm_4mm/</code> output directories.
The cache preserves the chosen diffusion-histogram resolution; grid changes
use those saved occupied cells directly.</p>
<p>Tests independently reconstruct event extrema, all histogram counts and
event pulse means from the sparse histories. Core tests cover strict
threshold equality, empty extrema, ownership and charge conservation. Next
interactive work: fix the random-distance control and add separate B/C track
reconstruction displays.</p>''', 'pixel-charge-batch')
    body += section('Electron readout: Petit diffusion prescription', '''
<p>Scope added 2026-10-29: each single-electron or double-electron ionisation
track supplies both an electron trail and a positive-ion trail. The existing
studies propagate positive ions. Extend the simulation to electrons drifting
to the anode, preserving the connection to the same MC energy deposits.</p>
<p>Petit’s loaded <code>src/itaca_functions.jl:17–25</code> defines
<code>sigma_t_mm(Lcm,Pbar;dtmm) = dtmm × sqrt(Lcm/Pbar)</code> and
<code>sigma_l_mm(Lcm,Pbar;dlmm) = dlmm × sqrt(Lcm/Pbar)</code>.
Widths are mm, L is cm, P is bar. The pressure-normalised coefficients have
units mm √(bar/cm) under this numeric convention.</p>
<table><thead><tr><th>Gas prescription</th><th>Transverse coefficient</th>
<th>Longitudinal coefficient</th><th>Petit source</th></tr></thead><tbody>
<tr><td>Pure Xe</td><td>3.5</td><td>0.9</td><td>Loaded function defaults and
active scripts/itaca_single_track_analysis.jl CLI defaults</td></tr>
<tr><td>Xe–He</td><td>1.6</td><td>0.75</td><td>Explicit Xe–He variables in
pluto/archive/itaca_st.jl:372–379; loaded select_events_itaca defaults</td></tr>
</tbody></table>
<p>The active caller passes its dt/dl and pressure to
<code>get_sigma("electron",...)</code>, then passes the resulting widths to
<code>diffuse_xyz_image_mc</code>: independent x/y Gaussian draws with sigma_t
and z draws with sigma_l. Gas selection enters through the coefficient pair.
The electron branch uses L, P, dt and dl; the electric-field argument is used
in the ion branch.</p>
<p>These source prescriptions identify the starting model. The inspected Xe–He
calls specify the coefficient pair. The user identifies the operating pressure
as 15 bar and the Xe–He mixture as 15% helium; these design inputs and both
coefficient pairs are in <code>config/electron_transport.toml</code>. The
operating-pressure coefficients are the tabulated values divided by √15.
<code>scripts/calculate_electron_diffusion.jl</code> saves the evaluated widths
and Float64 precision checks; <code>scripts/present_electron_diffusion.py</code>
generates <code>results/electron_diffusion/electron_diffusion.html</code> and
<code>latex_generated/electron_diffusion/electron_diffusion.pdf</code>. Carrier and
gas identity must enter the transport settings, cache identity and outputs.</p>''', 'electron-readout-petit')
    body += section('Paired drift and ELCC electron readout', '''
<p>User design inputs, 2026-10-29: sample one event location along the
1.5 m TPC. Ion drift is L towards the cathode; electron drift is 1.5 m − L
towards the anode. Save both diffused clouds before mapping them through
their grids, retaining their common event and ionisation-deposit identity.</p>
<p>The electron grid has cylindrical holes of 5 mm diameter and 5 mm depth,
on 10 mm pitch. Geometric collection assigns electrons to holes. The EL
field E/P = 3 kV/cm/bar at 15 bar gives E = 45 kV/cm and a 22.5 kV
potential rise along electron passage. Entry-grid and exit-anode potentials
belong in TOML; the absolute voltage reference is an outstanding input.</p>
<p>Electrons produce VUV photons along the hole. Calculate the direct optical
acceptance analytically and sample that fraction. Save photoelectron
histories per SiPM in configurable time buckets. User inputs:
1 μs per mm of drift, initially 2 μs buckets, and SiPM PDE = 0.30 in TOML.
For isotropic emission uniformly along the hole axis, an exit-plane SiPM
covering the full hole, and absorbing walls, the depth-averaged acceptance
is (h + R − √(h² + R²))/(2h), or 0.190983 for R=2.5 mm and h=5 mm.
This benchmark specifies an axial electron path. Sample incident photons
and photoelectrons with successive binomial draws.</p>
<p>Source: J. Hikida et al., NIM A 1080 (2025) 170706,
<code>papers/ELCC/Axel-ELCC.pdf</code>, Fig. 2 and sections 2 and 4.
The paper describes focusing, EL and SiPM readout; its upgraded unit has
4.5 mm holes covered by 4.7 mm round SiPMs. The 5/10/5 mm geometry above
is the user's design input. The paper supplies no quantitative EL yield
equation or numerical SiPM PDE.</p>
<p>New yield sources: <code>papers/ELCC/NEXT-EL-Xe-He.pdf</code>,
Table 1, uses Y = A max(E/P − X0,0) P h. At 15 bar, E/P=3 kV/cm/bar
and h=0.5 cm: pure Xe (A=137, X0=0.73) gives 2332.425 photons/electron;
Xe–15% He (A=139, X0=0.80) gives 2293.5. These matched fits are the
proposed starting prescriptions. Their measurements were near 1.2 bar;
application at 15 bar assumes reduced-field scaling.</p>
<p><code>papers/ELCC/EL-Next-DEMO.pdf</code>, Eq. (2), quotes the conventional
pure-Xe A=140, X0=0.83 prescription (2278.5 photons/electron here). Its
2.0–9.4 bar measurements show about 5% relative slope increase above 5 bar,
with about 30% absolute calibration uncertainty. It supplies no measured
15 bar or Xe–He correction. See <code>docs/decisions.md</code> for the source
audit, uncertainties and the internal percentage discrepancy in the Xe–He
paper. Yield prescription and EL emission fluctuations await review.</p>''', 'paired-drift-elcc')
    body += section('Cloud viewer II: paired ion and electron events', '''
<p><code>scripts/serve_cloud_viewer.py</code> opens a second local viewer at
<code>http://127.0.0.1:8767/</code>. It reads the six completed HDF5 streams in
<code>results/three_clouds_8000/</code>. The event number selects the same
ordinal event from the double-beta and single-electron samples (0–7999).</p>
<p>The first row shows ions, pure-Xe electrons and Xe–15% He electrons from
the double-beta event; the second row shows the corresponding three clouds
from the single-electron event. Every panel displays its event ID, ion drift
length, electron drift length and retained carrier count, followed by XY,
XZ, YZ and a rotatable 3D view. Colour encodes carrier counts. The 2D
projections retain the saved 100-bin resolution. The 3D view groups 4 × 4 × 4
saved cells for display and preserves the total carrier count. Event navigation
supports Previous, Next, an event-number field and left/right arrow keys.</p>
<p>The browser receives one event at a time from
<code>src/cloud_viewer/server.py</code>. The original interactive viewer
remains available for MC hits, reconstruction, transport and pad studies.</p>''', 'cloud-viewer-ii')
    body += section('Batch MC single-component selection', '''
<p><code>make mc-selection</code> runs <code>scripts/analyze_mc_selection.jl</code>
with <code>config/mc_selection.toml</code>, then generates presentation through
<code>scripts/present_mc_selection.py</code>. The run processes the first
requested events from each 100 μm sample, up to the available catalogue size.
The count is set by <code>run.n_events</code>. Event track length is the longest component
backbone. Selection accepts exactly one connected component after the configured
voxel threshold.</p>
<p>Per-event summaries, common-bin before/after histograms, efficiency estimates,
energy checks, provenance and serial/two-process benchmarks are saved to
<code>results/mc_selection/bb0nu.h5</code> and <code>electrons.h5</code>.
Plots and tables are presented in <code>results/mc_selection/mc_selection.html</code>
and <code>latex_generated/mc_selection/mc_selection.tex</code>/<code>.pdf</code>.
<code>make present-mc-selection</code> refreshes presentation from saved HDF5.</p>''', 'mc-selection-batch')
    body += section('Implementation checkpoints and review', table(['Step', 'Deliverable', 'Review/checks'], [
        ('1 · complete', 'Types, signatures, schema, defaults, state and repository layout in this chapter.', 'User accepted the interfaces as implementation guidance.'),
        ('2 · accepted', 'MC energy viewer: navigation, energy maps, 3D controls and session restore.', 'User accepted the viewer including event 2 and independent axes/zoom.'),
        ('3 · implemented', 'MC graph tab, overlays, component properties and backbone.', 'Parity with existing Julia outputs; inspect track construction together.'),
        ('4', 'Parameter editing, Apply, distance modes and reset.', 'Cache invalidation, stable event/view state, reproducibility on returning to settings.'),
        ('5', 'Session integration checks, responsiveness, usage documentation and selected-event report export.', 'Saved HDF5/settings, Python figures, generated LaTeX/PDF and build-log checks.'),
    ]) + '''
<p>Review now: inspect graph components and backbones on both samples, including
event 2 with its separated deposits. Next: interactive parameter editing and
Apply controls.</p>''', 'viewer-next')

    return dict(id='viewer', title='7. Interactive event viewer',
                summary='MC energy and graph inspection, component properties and usage.',
                status='MC tracks implemented · hands-on review', body=body)
