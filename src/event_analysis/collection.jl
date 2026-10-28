# Petit occupancy adaptation: itacaScripts/itaca_simple_reco.jl (MIT).
# Origins and targeted changes are recorded in PROVENANCE.toml.

"""Validate the full single-event diffusion table and sum its Int64 ion counts."""
function _diffused_total(diffused::DataFrame)::Int64
    validate_bin_table(diffused)
    hasproperty(diffused,:ions) && eltype(diffused.ions) == Int64 && all(>=(0),diffused.ions) ||
        throw(ArgumentError("Diffusion input requires nonnegative ions::Int64"))
    return foldl(Base.Checked.checked_add,diffused.ions;init=Int64(0))
end

"""
    select_diffused_event(diffused::DataFrame; bin_width_mm::Float32,
        max_distance_mm::Float32, min_ions_per_bin::Int64=10,
        method::String="kNN", k::Int=10, diffusion::DiffusionParams=DiffusionParams())

Apply Petit's reconstructed-event sequence to the full retained diffusion table:
voxelise, threshold, build components and accept exactly one component. Input
columns are event_id::Int64 (optional), x/y/z/energy::Float32 in mm/MeV, and
ions::Int64. Both energy and ions are summed in each reconstruction cell. Compare
the integer count with min_ions_per_bin inclusively, then call Petit's graph
builder with the already filtered table and zero additional energy threshold.
This implements the configured ion cut exactly while preserving energy weights.

Return a NamedTuple: accepted::Bool, tracks::Vector{Tracks}, bins::DataFrame
(all reconstruction cells), retained_ions::Int64, graph_retained_ions::Int64 and
graph_below_threshold_ions::Int64. Tracks are sorted by energy as in Petit;
their component indices refer to the thresholded parent table. Empty or multiple
components reject the whole event. retained = graph_retained + below_threshold.
Pad counting uses the original diffused table for accepted events.

Inputs retain their values. Returned tables/graphs own their arrays. Require
nonnegative threshold and counts; inherit the finite bin-width, cutoff, method
and k validation from voxelize_event and make_tracks. Count overflow raises
OverflowError. Diffusion metadata is supplied by the caller and preserved.

# Example
```julia
using DataFrames
d = DataFrame(x=Float32[0,1,20],y=zeros(Float32,3),z=zeros(Float32,3),
    energy=Float32[14,14,6].*21.9f-6,ions=Int64[14,14,6])
s = select_diffused_event(d;bin_width_mm=1f0,max_distance_mm=2f0)
@assert s.accepted && s.graph_retained_ions == 28 && s.retained_ions == 34
```
"""
function select_diffused_event(diffused::DataFrame;bin_width_mm::Float32,
        max_distance_mm::Float32,min_ions_per_bin::Int64=10,method::String="kNN",
        k::Int=10,diffusion::DiffusionParams=DiffusionParams())
    retained_ions = _diffused_total(diffused)
    min_ions_per_bin >= 0 || throw(ArgumentError("min_ions_per_bin must be nonnegative"))
    bins = voxelize_event(diffused,bin_width_mm)
    keep = (bins.ions .>= min_ions_per_bin) .& (bins.ions .> 0)
    graph_retained_ions = sum(bins.ions[keep])
    tracks = make_tracks(bins[keep,:];max_distance_mm,method,k,diffusion,energy_threshold_kev=0f0)
    return (;accepted=is_single_track(tracks),tracks,bins,retained_ions,graph_retained_ions,
            graph_below_threshold_ions=retained_ions-graph_retained_ions)
end

"""Require at least two finite, strictly increasing Float32 edge values."""
function _check_edges(edges::Vector{Float32})
    length(edges) >= 2 && all(isfinite,edges) && all(diff(edges) .> 0f0) ||
        throw(ArgumentError("Edges must be finite and strictly increasing"))
    return nothing
end

"""
Extend a regular grid by whole cells around base, retaining its phase and covering
all values. Return edges and the integer index of the first edge relative to base.
minimum_half_cells supplies a symmetric initial transverse extent; zero gives a
minimal covering axis. Float32 edge precision is checked after construction.
"""
function _covering_axis(values::Vector{Float32},base::Float32,width::Float32,minimum_half_cells::Int)
    lo,hi = extrema(values)
    first_index = min(-minimum_half_cells,floor(Int,(lo-base)/width))
    last_index = max(minimum_half_cells,ceil(Int,(hi-base)/width),first_index+1)
    # Include any coordinate just outside an endpoint after Float32 arithmetic.
    while base+Float32(first_index)*width > lo; first_index -= 1; end
    while base+Float32(last_index)*width < hi; last_index += 1; end
    edges = Float32[base+Float32(i)*width for i in first_index:last_index]
    _check_edges(edges)
    return edges,Int64(first_index)
end

"""
    pad_axes(diffused::DataFrame, track::Tracks; pitch_mm::Float32,
        z_bin_width_mm::Float32, grid_offset_mm::NTuple{2,Float32},
        minimum_xy_span_mm::Float32=160f0)

Build source-frame collection edges around an accepted reconstructed track,
following Petit's energy-weighted x/y centre and minimum-z slice reference.
Compute the centre from track.voxels, then shift the transverse edge lattice by
the supplied x/y offsets. The same grid serves all z slices. The normal caller
samples each offset in [0,pitch); any finite offset is accepted to support
integer-pitch relabelling tests. Pad pitch is independent of graph distance.

Start with an even number of transverse cells covering minimum_xy_span_mm
(Petit runner default 160 mm). Extend each axis by whole cells to cover every
retained diffusion-bin centre, including charge below the graph threshold.
Z starts from the selected track's minimum z and extends both ways in whole
z_bin_width_mm cells as needed. A zero-length track gets at least one slice.

Return owned Float32 x_edges_mm, y_edges_mm, z_edges_mm; reference_mm is the
NTuple{3,Float32} (energy-weighted x,y, selected minimum z); grid_offset_mm is
the supplied pair; first_edge_indices::NTuple{3,Int64} locates each lower edge
relative to its shifted reference lattice. Local pad indices are one-based in
these saved arrays. A centre reference is an edge, as in Petit's even-cell grid.

Require nonempty tables in the same frame, positive finite widths and track
energy, finite offsets and nonnegative minimum span. Validate Float32 edge
precision. Invalid arguments raise ArgumentError. Preserve all input data.
Collection assigns all ions at each diffusion-bin centre to its containing cell.
This is the approved geometric hole/pad assignment approximation.

# Example
```julia
using DataFrames
d = DataFrame(x=Float32[0,1,100],y=zeros(Float32,3),z=Float32[0,2,4],
    energy=Float32[20,20,1].*21.9f-6,ions=Int64[20,20,1])
s = select_diffused_event(d;bin_width_mm=1f0,max_distance_mm=3f0)
a = pad_axes(d,only(s.tracks);pitch_mm=2f0,z_bin_width_mm=2f0,
             grid_offset_mm=(0.25f0,0.5f0))
@assert first(a.x_edges_mm) <= 0f0 && last(a.x_edges_mm) >= 100f0
```
"""
function pad_axes(diffused::DataFrame,track::Tracks;pitch_mm::Float32,z_bin_width_mm::Float32,
        grid_offset_mm::NTuple{2,Float32},minimum_xy_span_mm::Float32=160f0)
    _diffused_total(diffused)
    validate_bin_table(track.voxels)
    nrow(diffused) > 0 && nrow(track.voxels) > 0 || throw(ArgumentError("Pad axes require nonempty diffusion and track tables"))
    all(v -> isfinite(v) && v > 0f0,(pitch_mm,z_bin_width_mm)) || throw(ArgumentError("Cell widths must be finite and positive"))
    all(isfinite,grid_offset_mm) && isfinite(minimum_xy_span_mm) && minimum_xy_span_mm >= 0f0 ||
        throw(ArgumentError("Offsets must be finite; minimum span must be finite and nonnegative"))
    e = track.voxels.energy
    total = sum(e)
    isfinite(total) && total > 0f0 || throw(ArgumentError("Selected track must have positive finite energy"))
    reference_mm = (sum(track.voxels.x.*e)/total,sum(track.voxels.y.*e)/total,minimum(track.voxels.z))
    return pad_axes(diffused;reference_mm,pitch_mm,z_bin_width_mm,grid_offset_mm,minimum_xy_span_mm)
end

"""
    pad_axes(diffused; reference_mm::NTuple{3,Float32}, pitch_mm::Float32,
        z_bin_width_mm::Float32, grid_offset_mm::NTuple{2,Float32},
        minimum_xy_span_mm::Float32=160f0)

Construct the same covering pad lattice from an explicit source-frame reference
in mm. This supports collection before reconstruction. The caller supplies the
x/y reference and minimum-z reference; all retained diffusion charge is covered.
Return owned axes, reference, phase and first-edge indices as in the track method.
Validate positive widths, finite references/offsets, nonnegative minimum span and
nonempty diffusion data. Preserve input arrays and exact integer charge counts.
"""
function pad_axes(diffused::DataFrame;reference_mm::NTuple{3,Float32},pitch_mm::Float32,
        z_bin_width_mm::Float32,grid_offset_mm::NTuple{2,Float32},minimum_xy_span_mm::Float32=160f0)
    _diffused_total(diffused)
    nrow(diffused)>0 || throw(ArgumentError("Pad axes require nonempty diffusion"))
    all(v->isfinite(v)&&v>0f0,(pitch_mm,z_bin_width_mm)) || throw(ArgumentError("Cell widths must be finite and positive"))
    all(isfinite,grid_offset_mm) && isfinite(minimum_xy_span_mm) && minimum_xy_span_mm>=0f0 || throw(ArgumentError("Invalid grid phase or span"))
    all(isfinite,reference_mm) || throw(ArgumentError("Track reference overflows Float32"))
    half = ceil(Int,(minimum_xy_span_mm/2f0)/pitch_mm)
    x,ix = _covering_axis(diffused.x,reference_mm[1]+grid_offset_mm[1],pitch_mm,half)
    y,iy = _covering_axis(diffused.y,reference_mm[2]+grid_offset_mm[2],pitch_mm,half)
    z,iz = _covering_axis(diffused.z,reference_mm[3],z_bin_width_mm,0)
    return (;x_edges_mm=x,y_edges_mm=y,z_edges_mm=z,reference_mm,grid_offset_mm,
             first_edge_indices=(ix,iy,iz))
end

"""
    compute_occupancy_analysis(diffused::DataFrame;
        x_edges_mm::Vector{Float32}, y_edges_mm::Vector{Float32},
        z_edges_mm::Vector{Float32}, drift_speed_mm_ms::Float32, drift_direction::Int)

Retain the full pad-by-slice sequence from Petit's pixel occupancy calculation.
Assign each diffusion-bin centre to its x/y cell and z slice using the saved
edges. Sum ions::Int64 exactly, including the full accepted event's charge below
the reconstruction threshold. Internal cells are half open; final upper edges
are included. Every positive-count centre must be covered, otherwise throw
ArgumentError rather than returning a truncated charge history.

Return a NamedTuple with counts::DataFrame (pad_ix, pad_iy, slice_index, ions,
all Int64), copied Float32 x/y/z_edges_mm and time_edges_ms, drift_direction::Int,
time_origin_z_mm::Float32 and total_ions::Int64. Store occupied entries only,
sorted by pad_ix, pad_iy, slice_index. Missing entries in the saved axis ranges
represent zero counts; local pad indices are one-based.

slice_index is increasing relative arrival time. For drift_direction=+1,
increasing source z is later and t=(z-z_first)/v. For -1, decreasing source z
is later and t=(z_last-z)/v; slice i maps to ascending-z bin n_z-i+1. time_edges
always increase from zero; time_origin_z_mm records the corresponding end of
the saved z range. Thus 2 mm gives 20 ms at 0.1 mm/ms for either direction.
The caller supplies the sign explicitly. These are relative event times.

Require finite increasing edges, positive finite speed and sign ±1. Validate
the single-event diffusion columns and checked total counts. Outputs own their
arrays and preserve inputs. Empty diffusion gives a typed empty count table.
The function consumes no RNG. Source: Petit compute_occupancy_analysis in
itacaScripts/itaca_simple_reco.jl; full-diffusion input, integer counts, complete
coverage, explicit axes/time and retained sequences are the approved adaptations.

# Example
```julia
using DataFrames
d = DataFrame(x=Float32[0,0],y=Float32[0,0],z=Float32[0,2],
    energy=Float32[14,6].*21.9f-6,ions=Int64[14,6])
p = compute_occupancy_analysis(d;x_edges_mm=Float32[-1,1],y_edges_mm=Float32[-1,1],
    z_edges_mm=Float32[0,2,4],drift_speed_mm_ms=0.1f0,drift_direction=1)
@assert p.counts.ions == Int64[14,6] && p.time_edges_ms == Float32[0,20,40]
```
"""
function compute_occupancy_analysis(diffused::DataFrame;x_edges_mm::Vector{Float32},
        y_edges_mm::Vector{Float32},z_edges_mm::Vector{Float32},
        drift_speed_mm_ms::Float32,drift_direction::Int)
    total_ions = _diffused_total(diffused)
    foreach(_check_edges,(x_edges_mm,y_edges_mm,z_edges_mm))
    isfinite(drift_speed_mm_ms) && drift_speed_mm_ms > 0f0 || throw(ArgumentError("Drift speed must be finite and positive"))
    drift_direction in (-1,1) || throw(ArgumentError("drift_direction must be +1 or -1"))
    nz = length(z_edges_mm)-1
    time_origin_z_mm = drift_direction == 1 ? first(z_edges_mm) : last(z_edges_mm)
    time_edges_ms = drift_direction == 1 ? (z_edges_mm.-time_origin_z_mm)./drift_speed_mm_ms :
                                          (time_origin_z_mm.-reverse(z_edges_mm))./drift_speed_mm_ms
    _check_edges(time_edges_ms)
    sums = Dict{NTuple{3,Int64},Int64}()
    for row in eachrow(diffused)
        row.ions == 0 && continue
        i = _histogram_index(row.x,x_edges_mm)
        j = _histogram_index(row.y,y_edges_mm)
        z = _histogram_index(row.z,z_edges_mm)
        min(i,j,z) > 0 || throw(ArgumentError("Pad axes do not cover all retained diffusion charge"))
        k = drift_direction == 1 ? z : nz-z+1
        key = (Int64(i),Int64(j),Int64(k))
        sums[key] = Base.Checked.checked_add(get(sums,key,Int64(0)),row.ions)
    end
    counts = DataFrame(pad_ix=Int64[],pad_iy=Int64[],slice_index=Int64[],ions=Int64[])
    for key in sort(collect(keys(sums)));push!(counts,(key...,sums[key]));end
    return (;counts,x_edges_mm=copy(x_edges_mm),y_edges_mm=copy(y_edges_mm),
             z_edges_mm=copy(z_edges_mm),time_edges_ms,drift_direction,time_origin_z_mm,total_ions)
end

"""
    pad_projections(history) -> NamedTuple

Reduce a compute_occupancy_analysis result into xy_ions::Matrix{Int64}
(x index first, y second), slice_ions::Vector{Int64} in time order, and
total_ions::Int64. All slices and transverse cells are present, including zeros.
Each projection sums the entire event charge, preserving the sparse input.
Allocate independent arrays. The caller controls whether dense projections are
needed; the primary history remains sparse. These are per-event totals.

# Example
```julia
using DataFrames
d = DataFrame(x=Float32[0],y=Float32[0],z=Float32[0],energy=Float32[0.001],ions=Int64[10])
h = compute_occupancy_analysis(d;x_edges_mm=Float32[-1,1],y_edges_mm=Float32[-1,1],
    z_edges_mm=Float32[0,2,4],drift_speed_mm_ms=0.1f0,drift_direction=1)
p = pad_projections(h)
@assert p.xy_ions == reshape(Int64[10],1,1) && p.slice_ions == Int64[10,0]
```
"""
function pad_projections(history)
    xy_ions = zeros(Int64,length(history.x_edges_mm)-1,length(history.y_edges_mm)-1)
    slice_ions = zeros(Int64,length(history.time_edges_ms)-1)
    for row in eachrow(history.counts)
        xy_ions[row.pad_ix,row.pad_iy] += row.ions
        slice_ions[row.slice_index] += row.ions
    end
    return (;xy_ions,slice_ions,total_ions=history.total_ions)
end
