"""
    MCGraphSettings(; voxel_width_mm=0.5f0, method="KDT",
        distance_mode="voxel_multiple", distance_voxel_multiple=2.5f0,
        max_distance_mm=1.25f0, energy_threshold_keV=1f0, k=10)

Settings for graph inspection of original MC hits. Physical scalars are
Float32, distances are mm, and the inclusive voxel-energy threshold is keV.
`method` selects KDT, kNN or kNN_mutual. `distance_mode` selects either
voxel_width_mm × distance_voxel_multiple or the fixed max_distance_mm value.
`k::Int` is the nearest-neighbour count for the kNN methods.

Defaults follow Petit's active MC branch in
scripts/itaca_single_track_analysis.jl and get_voxel_size_and_distance.
Construction converts scalar inputs, owns the settings, and validates them
through `resolved_mc_distance`. Invalid modes, nonfinite/nonpositive widths,
distances or multipliers, a negative/nonfinite threshold, or k < 1 raise
ArgumentError. Both distance entries are retained for provenance.
"""
struct MCGraphSettings
    voxel_width_mm::Float32
    method::String
    distance_mode::String
    distance_voxel_multiple::Float32
    max_distance_mm::Float32
    energy_threshold_keV::Float32
    k::Int
end
function MCGraphSettings(;voxel_width_mm=0.5f0,method="KDT",distance_mode="voxel_multiple",
        distance_voxel_multiple=2.5f0,max_distance_mm=1.25f0,energy_threshold_keV=1f0,k=10)
    s=MCGraphSettings(Float32(voxel_width_mm),String(method),String(distance_mode),
        Float32(distance_voxel_multiple),Float32(max_distance_mm),Float32(energy_threshold_keV),Int(k))
    resolved_mc_distance(s)
    return s
end

"""
    resolved_mc_distance(settings::MCGraphSettings) -> Float32

Validate all MC graph settings and return the maximum connection distance in
mm. Multiplier mode follows voxel width; fixed_mm mode uses max_distance_mm.
Distance equality is accepted by the existing graph builder. Overflow of the
resolved Float32 distance raises ArgumentError, as do invalid settings listed
in `MCGraphSettings`. Input settings retain their values.

# Example
```julia
@assert resolved_mc_distance(MCGraphSettings()) === 1.25f0
@assert resolved_mc_distance(MCGraphSettings(voxel_width_mm=1)) === 2.5f0
```
"""
function resolved_mc_distance(s::MCGraphSettings)::Float32
    for x in (s.voxel_width_mm,s.distance_voxel_multiple,s.max_distance_mm)
        isfinite(x) && x>0f0 || throw(ArgumentError("Graph widths, distances and multipliers must be finite and positive"))
    end
    s.method in ("KDT","kNN","kNN_mutual") || throw(ArgumentError("Unknown MC graph method"))
    s.distance_mode in ("voxel_multiple","fixed_mm") || throw(ArgumentError("Unknown MC distance mode"))
    isfinite(s.energy_threshold_keV) && s.energy_threshold_keV>=0f0 || throw(ArgumentError("MC threshold must be finite and nonnegative"))
    s.k>=1 || throw(ArgumentError("k must be positive"))
    distance=s.distance_mode=="voxel_multiple" ? s.voxel_width_mm*s.distance_voxel_multiple : s.max_distance_mm
    isfinite(distance) && distance>0f0 || throw(ArgumentError("Resolved MC distance must be finite and positive"))
    return distance
end

"""
    MCComponent

Owned presentation record for one connected MC graph component. component_id
is its one-based rank by decreasing deposited energy, following make_tracks.
`nodes` owns Float32 x_mm/y_mm/z_mm/energy_MeV columns. Local node index i
identifies row i; edge_src/edge_dst contain each undirected edge once, src < dst.
path_indices lists ordered backbone nodes; path_s_mm gives cumulative distance
in mm. All indices are Int64 and one-based. endpoint_indices contains the first
and last path nodes, or nothing for an empty path. A one-node path repeats that
node as both endpoints and has zero length. energy_MeV, length_mm and confidence
are Float32; confidence is Petit's backbone algorithm score.

`build_mc_components` creates independent arrays and tables. Their contents are
mutable and retain the coordinate frame of the input voxel table.
"""
struct MCComponent
    component_id::Int64
    nodes::DataFrame
    edge_src::Vector{Int64}
    edge_dst::Vector{Int64}
    path_indices::Vector{Int64}
    endpoint_indices::Union{Nothing,NTuple{2,Int64}}
    path_s_mm::Vector{Float32}
    energy_MeV::Float32
    length_mm::Float32
    confidence::Float32
end

"""
    build_mc_components(voxels::DataFrame; method::String,
        max_distance_mm::Float32, energy_threshold_keV::Float32, k::Int)
        -> Vector{MCComponent}

Inspect a voxel table using the existing Petit adaptation: make_tracks,
walk_track_from_extremes and get_raw_path. Input columns are Float32 x/y/z in
mm and energy in MeV, with optional event_id::Int64. The inclusive threshold
applies to voxel energy. Return every component in make_tracks energy order,
with owned nodes, edges, backbone indices, endpoints and cumulative distance.
The walker uses its existing energy-weighting defaults. An empty filtered
table returns MCComponent[]. Inputs and their row order are preserved.

Validation errors follow make_tracks: invalid physical values, columns,
method, threshold, distance or k raise ArgumentError. This adapter converts
column names and records graph-local numbering for HDF5 and plotting.

# Example
```julia
using DataFrames
v=DataFrame(x=Float32[0,1,20],y=zeros(Float32,3),z=zeros(Float32,3),energy=Float32[1,1,0.5])
c=build_mc_components(v;method="KDT",max_distance_mm=1.25f0,energy_threshold_keV=1f0,k=10)
@assert length(c)==2 && c[1].energy_MeV==2f0
@assert c[2].endpoint_indices==(1,1) && c[2].length_mm==0f0
```
"""
function build_mc_components(voxels::DataFrame;method::String,max_distance_mm::Float32,
        energy_threshold_keV::Float32,k::Int)::Vector{MCComponent}
    tracks=make_tracks(voxels;method,max_distance_mm,energy_threshold_kev=energy_threshold_keV,k)
    result=MCComponent[]
    for (i,track) in enumerate(tracks)
        walk=walk_track_from_extremes(track)
        path=Int64.(walk.path_indices)
        raw=get_raw_path(track,walk.path_indices)
        nodes=DataFrame(x_mm=track.voxels.x,y_mm=track.voxels.y,z_mm=track.voxels.z,energy_MeV=track.voxels.energy)
        es=collect(edges(track.graph))
        endpoints=isempty(path) ? nothing : (first(path),last(path))
        push!(result,MCComponent(Int64(i),nodes,Int64[src(e) for e in es],Int64[dst(e) for e in es],
            path,endpoints,copy(raw.s),sum(track.voxels.energy),walk.total_length,walk.confidence))
    end
    return result
end
