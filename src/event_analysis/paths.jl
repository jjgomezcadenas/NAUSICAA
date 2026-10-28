"""
    compute_mc_path(hits::EventHits, info::EventInfo, bin_width_mm::Float32) -> DataFrame

Construct Petit's time-ordered primary trajectory in the supplied hit frame.
Use actual sorted primary IDs from EventInfo (one or two), group consecutive
visits to the same floor-defined cell and use energy-weighted coordinates in mm
(arithmetic means for zero-energy groups). For two primaries reverse the first
arm and append the second. Stable equal-time ordering follows source row order.

Return newly owned Float32 x/y/z/energy/s columns (mm/MeV/mm) and particle_id::Int32
identifying the actual primary for each path row. A primary with no deposited
hits gives a typed empty path, as in Petit. Revisited cells separated in time
remain distinct path entries. Inputs retain their values. Validate matching IDs,
finite positive width, current hit values and one/two distinct primary IDs;
invalid arguments raise `ArgumentError` or the EventHits validation error.
The MC trajectory is a truth comparator; graph construction uses all event hits.

# Example
```julia
h = EventHits(42, source_frame, Float32[2,0], zeros(Float32,2), zeros(Float32,2),
              Float32[1,1], Float32[1,0], Int32[7,7], Int32[0,1])
p = compute_mc_path(h, EventInfo(42,(0f0,0f0,0f0),Int32[7]), 1f0)
@assert p.x == Float32[0,2] && p.s == Float32[0,2]
```
"""
function compute_mc_path(hits::EventHits, info::EventInfo, bin_width_mm::Float32)::DataFrame
    hits.event_id == info.event_id || throw(ArgumentError("Event IDs must match"))
    _validate_event_info(info.vertex_mm, info.primary_particle_ids)
    length(info.primary_particle_ids) in (1,2) || throw(ArgumentError("Expected one or two primaries"))
    isfinite(bin_width_mm) && bin_width_mm > 0f0 || throw(ArgumentError("bin_width_mm must be finite and positive"))
    return _compute_mc_path(to_petit_dataframe(hits), bin_width_mm, sort(info.primary_particle_ids))
end

@doc """
    get_raw_path(track::Tracks, path_indices::Vector{Int}) -> DataFrame

Extract selected track rows in the specified order and calculate cumulative
Euclidean distance s starting at zero. Output x/y/z/s columns are Float32 in mm;
arrays are newly allocated. Empty indices give a typed empty table; a single
point has s = 0. Invalid row indices raise BoundsError. Coordinates and graph
are preserved. Adapted directly from Petit src/track_path_reconstruction.jl.

# Example
```julia
using DataFrames
d = DataFrame(x=Float32[0,3],y=Float32[0,4],z=Float32[0,0],energy=Float32[1,1])
t = only(make_tracks(d; max_distance_mm=5f0))
@assert get_raw_path(t,[1,2]).s == Float32[0,5]
```
""" get_raw_path

@doc """
    make_tracks(event_data::DataFrame; max_distance_mm::Float32,
        energy_threshold_kev::Float32=0f0, diffusion=DiffusionParams(),
        method="kNN", k=10) -> Vector{Tracks}

Build Petit's graph from one event's Float32 x/y/z/energy table (mm/MeV). Keep
bins with energy ≥ energy_threshold_kev × 1f-3; threshold is in keV. Require an
explicit finite positive maximum edge distance in mm and k ≥ 1. Supported
methods are KDT (radius), kNN (undirected union of nearest-neighbour choices)
and kNN_mutual (both endpoints must select one another). Distance equality
passes. kNN queries min(k, N-1) other points and applies the distance cutoff.
Equal-distance neighbour choices follow NearestNeighbors and the supplied row
order. The maximum degree can exceed k through incoming choices.

Return components sorted by decreasing total Float32 energy. Each owns a new
voxel table, graph and component-index vector; local vertex i denotes row i,
while components stores parent indices after energy filtering. Empty filtered
input returns Tracks[]; a single surviving bin returns a one-vertex track.
Invalid columns, physical values, method, cutoff, k or threshold raise
ArgumentError. Input values and row order remain unchanged.

Source: Petit src/track_building.jl. Default method is adapted to the approved
kNN baseline; the required cutoff avoids the original wrapper's implicit 1 mm
default. Explicit energy thresholds keep this deterministic stage usable with
raw MC and future diffused-energy tables.

# Example
```julia
using DataFrames
d = DataFrame(x=Float32[0,1,20],y=zeros(Float32,3),z=zeros(Float32,3),energy=Float32[1,1,3])
t = make_tracks(d;max_distance_mm=2f0,k=1)
@assert length(t) == 2 && t[1].voxels.x == Float32[20]
@assert !is_single_track(t)
```
""" make_tracks

@doc """
    walk_track_from_extremes(track::Tracks; use_energy_weighting=true,
        use_edge_energy_weighting=true, use_mst_fallback=false,
        dense_track_threshold::Float32=6f0)

Run Petit's combined topology/energy/spatial/curvature path algorithm on an
aligned connected Tracks object. Preserve its branch conditions and diagnostic
confidence calculation; MST remains an optional fallback and separately callable
comparison. Numerical weights, lengths and confidence use Float32.

Return a NamedTuple with extremes (two owned NamedTuple endpoint rows),
path_indices::Vector{Int} (local track row indices), path_voxels::DataFrame
(copied rows in path order), total_length::Float32 (mm) and confidence::Float32
(algorithm score). Empty tracks return (nothing,nothing) endpoints and an empty
path; a single vertex is both endpoints with zero length. Results preserve the
input table and graph. Endpoint copies remain independent of both input and
returned path table. The algorithm assumes finite Float32 coordinates/energies
and matching graph vertices, as supplied by make_tracks.

Source: Petit src/track_extreme_finding.jl. Endpoint ownership is the deliberate
adaptation from Petit's DataFrameRow views. Equal-cost path alternatives are
handled in validation by comparing path validity and physical observables.

# Example
```julia
using DataFrames
d = DataFrame(x=Float32[0,1,2],y=zeros(Float32,3),z=zeros(Float32,3),energy=ones(Float32,3))
t = only(make_tracks(d;max_distance_mm=1f0))
p = walk_track_from_extremes(t)
@assert p.total_length === 2f0 && length(p.path_indices) == 3
```
""" walk_track_from_extremes
