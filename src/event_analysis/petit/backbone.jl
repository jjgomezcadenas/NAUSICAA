# Adapted from Petit src/track_extreme_finding.jl (MIT; see LICENSE).
# Source SHA-256: b35c21e21cb9504f789e33b9985a860b06a57600b21ef9b35111e8cb89820a12
# Float32 adaptation; algorithm branches retained. See ../PROVENANCE.toml.

# Track Extreme Finding Functions
# Functions for finding track endpoints (extremes)

using LinearAlgebra
using Statistics
using DataFrames
using Graphs
using SparseArrays

"""
    calculate_path_length_from_coords(coords::AbstractMatrix{<:Real},
                                      path::AbstractVector{Int}) -> Float32

Compute the total Euclidean length of a path defined by vertex indices.

# Arguments
- `coords`: 3×N matrix of voxel coordinates
- `path`: ordered vector of vertex indices along the path

# Returns
- Total path length (Float32)
"""
function calculate_path_length_from_coords(coords::AbstractMatrix{<:Real},
                                           path::AbstractVector{Int})

    n = length(path)
    n < 2 && return 0.0f0

    L = 0.0f0
    @inbounds for i in 1:(n-1)
        u = path[i]
        v = path[i+1]

        dx = coords[1,u] - coords[1,v]
        dy = coords[2,u] - coords[2,v]
        dz = coords[3,u] - coords[3,v]

        L += sqrt(dx*dx + dy*dy + dz*dz)
    end

    return L
end


"""
    find_track_extremes(track::Tracks; kwargs...)

Find track endpoints using combined topology/curvature/energy analysis.

# Arguments
- `track::Tracks`: Track object
- `use_energy_weighting::Bool=true`: Enable energy-weighted extreme finding
- `use_edge_energy_weighting::Bool=true`: Enable edge-energy-weighted method
- `use_mst_fallback::Bool=false`: Enable MST-based fallback method
- `dense_track_threshold::Float32=6.0`: Average degree threshold for dense tracks

# Returns
- `(extreme1_idx, extreme2_idx, path, confidence)`: Endpoint indices, path between them, confidence score (0-1)
"""
function find_track_extremes(track::Tracks;
                             use_energy_weighting::Bool=true,
                             use_edge_energy_weighting::Bool=true,
                             use_mst_fallback::Bool=false,
                             dense_track_threshold::Float32=6.0f0)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0f0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0f0)
    end

    # Pre-extract coordinates once
    coords = extract_coords(track)

    return find_extremes_combined(track, coords;
                                  use_energy_weighting=use_energy_weighting,
                                  use_edge_energy_weighting=use_edge_energy_weighting,
                                  use_mst_fallback=use_mst_fallback,
                                  dense_track_threshold=dense_track_threshold)
end


"""
    walk_track_from_extremes(track::Tracks; kwargs...)

Walk through a track from one endpoint to the other.

# Arguments
- `track::Tracks`: Track object
- `use_energy_weighting::Bool=true`: Enable energy-weighted extreme finding
- `use_edge_energy_weighting::Bool=true`: Enable edge-energy-weighted method
- `use_mst_fallback::Bool=false`: Enable MST-based fallback method
- `dense_track_threshold::Float32=6.0`: Average degree threshold for dense tracks

# Returns
NamedTuple with:
- `extremes`: (start_voxel, end_voxel) owned NamedTuple rows
- `path_indices`: Vertex indices along the path
- `path_voxels`: DataFrame of voxels in path order
- `total_length`: Path length in mm
- `confidence`: Confidence score (0-1)
"""
function walk_track_from_extremes(track::Tracks;
                                  use_energy_weighting::Bool=true,
                                  use_edge_energy_weighting::Bool=true,
                                  use_mst_fallback::Bool=false,
                                  dense_track_threshold::Float32=6.0f0)
    extreme1, extreme2, path, confidence = find_track_extremes(track;
                                            use_energy_weighting=use_energy_weighting,
                                            use_edge_energy_weighting=use_edge_energy_weighting,
                                            use_mst_fallback=use_mst_fallback,
                                            dense_track_threshold=dense_track_threshold)

    if isnothing(extreme1)
        return (extremes = (nothing, nothing),
                path_indices = Int[],
                path_voxels = DataFrame(),
                total_length = 0.0f0,
                confidence = 0.0f0)
    end

    # Get voxel data for extremes
    start_voxel = copy(track.voxels[extreme1, :])
    end_voxel = copy(track.voxels[extreme2, :])

    # Get voxels along the path
    path_voxels = track.voxels[path, :]

    # Calculate total path length
    total_length = 0.0f0
    for i in 1:length(path)-1
        v1, v2 = path[i], path[i+1]
        total_length += euclidean_distance(
            track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
            track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
        )
    end

    return (extremes = (start_voxel, end_voxel),
            path_indices = path,
            path_voxels = path_voxels,
            total_length = total_length,
            confidence = confidence)
end


"""
    find_track_extremes(trk; i=1, kwargs...)

Convenience function: find extremes for i-th track in a vector.

# Arguments
- `trk`: Vector of Tracks
- `i::Int=1`: Index of track to analyze
- `use_energy_weighting::Bool=true`: Enable energy-weighted extreme finding
- `use_edge_energy_weighting::Bool=true`: Enable edge-energy-weighted method
- `use_mst_fallback::Bool=false`: Enable MST-based fallback method
- `dense_track_threshold::Float32=6.0`: Average degree threshold for dense tracks
"""
function find_track_extremes(trk::Vector{Tracks}; i=1,
                             use_energy_weighting::Bool=true,
                             use_edge_energy_weighting::Bool=true,
                             use_mst_fallback::Bool=false,
                             dense_track_threshold::Float32=6.0f0)
    return walk_track_from_extremes(trk[i];
                                    use_energy_weighting=use_energy_weighting,
                                    use_edge_energy_weighting=use_edge_energy_weighting,
                                    use_mst_fallback=use_mst_fallback,
                                    dense_track_threshold=dense_track_threshold)
end


# =============================================================================
# TRACK EXTREMES ALGORITHMS
# =============================================================================

"""
    TrackCoords

Float32 x/y/z vectors in mm, indexed by the local graph vertex. The direct
constructor stores supplied arrays. `extract_coords(track)` makes independent
copies for backbone calculations; it preserves the track table and row order.
"""
struct TrackCoords
    x::Vector{Float32}
    y::Vector{Float32}
    z::Vector{Float32}
end

"""Extract coordinates from track voxels once."""
function extract_coords(track::Tracks)
    return TrackCoords(
        collect(track.voxels.x),
        collect(track.voxels.y),
        collect(track.voxels.z)
    )
end

"""
    find_path_bfs(g, start_vertex, end_vertex)

Optimized BFS using parent pointers instead of copying paths.
Returns the ordered vertex indices as Vector{Int}; an unreachable target gives Int[].
"""
function find_path_bfs(g::SimpleGraph, start_vertex::Int, end_vertex::Int)
    if start_vertex == end_vertex
        return [start_vertex]
    end

    n = nv(g)
    parent = zeros(Int, n)
    visited = falses(n)

    # Use a simple queue (circular buffer would be faster but this is clearer)
    queue = Int[start_vertex]
    visited[start_vertex] = true

    while !isempty(queue)
        current = popfirst!(queue)

        for neighbor in neighbors(g, current)
            if !visited[neighbor]
                visited[neighbor] = true
                parent[neighbor] = current

                if neighbor == end_vertex
                    # Reconstruct path from parent pointers
                    path = Int[end_vertex]
                    node = end_vertex
                    while parent[node] != 0
                        node = parent[node]
                        pushfirst!(path, node)
                    end
                    return path
                end

                push!(queue, neighbor)
            end
        end
    end

    return Int[]  # No path found
end

"""
    calculate_path_length(coords, path)

Calculate path length using pre-extracted coordinates.
"""
function calculate_path_length(coords::TrackCoords, path::Vector{Int})
    length(path) < 2 && return 0.0f0

    total = 0.0f0
    @inbounds for i in 1:length(path)-1
        v1, v2 = path[i], path[i+1]
        dx = coords.x[v2] - coords.x[v1]
        dy = coords.y[v2] - coords.y[v1]
        dz = coords.z[v2] - coords.z[v1]
        total += sqrt(dx*dx + dy*dy + dz*dz)
    end
    return total
end

"""
    euclidean_distance_coords(coords, v1, v2)

Euclidean distance using pre-extracted coordinates.
"""
@inline function euclidean_distance_coords(coords::TrackCoords, v1::Int, v2::Int)
    dx = coords.x[v2] - coords.x[v1]
    dy = coords.y[v2] - coords.y[v1]
    dz = coords.z[v2] - coords.z[v1]
    return sqrt(dx*dx + dy*dy + dz*dz)
end

"""
    find_track_extremes_internal(track)

Optimized version of find_track_extremes using pre-extracted coordinates
and early-exit strategies.
"""
function find_track_extremes_internal(track::Tracks)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0f0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0f0)
    end

    # Pre-extract coordinates once
    coords = extract_coords(track)

    return find_extremes_combined(track, coords)
end

"""
    find_extremes_topology(track, coords)

Optimized topology-based extreme finding.
Now includes path coverage check to avoid false positives from adjacent degree-1 vertices.
"""
function find_extremes_topology(track::Tracks, coords::TrackCoords)
    g = track.graph
    n_vertices = nv(g)

    # Compute track extent for coverage check
    x_range = maximum(coords.x) - minimum(coords.x)
    y_range = maximum(coords.y) - minimum(coords.y)
    z_range = maximum(coords.z) - minimum(coords.z)
    track_extent = sqrt(x_range^2 + y_range^2 + z_range^2)

    # Find degree-1 vertices efficiently
    endpoints = Int[]
    sizehint!(endpoints, 10)
    for v in 1:n_vertices
        if degree(g, v) == 1
            push!(endpoints, v)
        end
    end

    # Exactly 2 endpoints - check if they give good coverage
    if length(endpoints) == 2
        path = find_path_bfs(g, endpoints[1], endpoints[2])
        path_length = calculate_path_length(coords, path)
        coverage = path_length / max(track_extent, 1.0f0)

        # High confidence only if coverage is good
        if coverage >= 0.5f0
            confidence = 0.85f0 + 0.1f0 * min(1.0f0, coverage)  # 0.85-0.95
        else
            confidence = 0.5f0 + 0.3f0 * coverage  # 0.5-0.65 for poor coverage
        end

        return (endpoints[1], endpoints[2], path, confidence, path_length)
    end

    # More than 2 endpoints - find pair with maximum path length
    if length(endpoints) > 2
        max_path_length = 0.0f0
        best_pair = (endpoints[1], endpoints[2])
        best_path = Int[]

        for i in 1:length(endpoints)-1
            for j in i+1:length(endpoints)
                path = find_path_bfs(g, endpoints[i], endpoints[j])
                path_length = calculate_path_length(coords, path)

                if path_length > max_path_length
                    max_path_length = path_length
                    best_pair = (endpoints[i], endpoints[j])
                    best_path = path
                end
            end
        end

        coverage = max_path_length / max(track_extent, 1.0f0)
        if coverage >= 0.5f0
            confidence = 0.75f0 + 0.15f0 * min(1.0f0, coverage)  # 0.75-0.90
        else
            confidence = 0.4f0 + 0.3f0 * coverage  # 0.4-0.55 for poor coverage
        end

        return (best_pair[1], best_pair[2], best_path, confidence, max_path_length)
    end

    # Fallback: no clear endpoints
    return find_extremes_distance_fallback_opt(track, coords, 0.6f0)
end

"""
    calculate_vertex_curvatures(track, coords)

Optimized curvature calculation avoiding temporary array allocations.
"""
function calculate_vertex_curvatures(track::Tracks, coords::TrackCoords)
    g = track.graph
    n_vertices = nv(g)
    curvatures = zeros(Float32, n_vertices)

    @inbounds for v in 1:n_vertices
        nbrs = neighbors(g, v)
        n_nbrs = length(nbrs)

        if n_nbrs < 2
            curvatures[v] = 0.0f0  # Endpoint
        elseif n_nbrs == 2
            v1, v2 = nbrs[1], nbrs[2]

            # Compute vectors inline (no array allocation)
            dx1 = coords.x[v1] - coords.x[v]
            dy1 = coords.y[v1] - coords.y[v]
            dz1 = coords.z[v1] - coords.z[v]

            dx2 = coords.x[v2] - coords.x[v]
            dy2 = coords.y[v2] - coords.y[v]
            dz2 = coords.z[v2] - coords.z[v]

            norm1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
            norm2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2)

            if norm1 > 1f-6 && norm2 > 1f-6
                # Normalize and compute dot product
                inv_norm1 = 1.0f0 / norm1
                inv_norm2 = 1.0f0 / norm2
                cos_angle = (dx1*dx2 + dy1*dy2 + dz1*dz2) * inv_norm1 * inv_norm2
                cos_angle = clamp(cos_angle, -1.0f0, 1.0f0)
                curvatures[v] = 1.0f0 - abs(cos_angle)
            end
        else
            # Multiple neighbors - compute average curvature
            total_curv = 0.0f0
            count = 0

            for i in 1:n_nbrs-1
                for j in i+1:n_nbrs
                    vi, vj = nbrs[i], nbrs[j]

                    dx1 = coords.x[vi] - coords.x[v]
                    dy1 = coords.y[vi] - coords.y[v]
                    dz1 = coords.z[vi] - coords.z[v]

                    dx2 = coords.x[vj] - coords.x[v]
                    dy2 = coords.y[vj] - coords.y[v]
                    dz2 = coords.z[vj] - coords.z[v]

                    norm1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
                    norm2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2)

                    if norm1 > 1f-6 && norm2 > 1f-6
                        inv_norm1 = 1.0f0 / norm1
                        inv_norm2 = 1.0f0 / norm2
                        cos_angle = (dx1*dx2 + dy1*dy2 + dz1*dz2) * inv_norm1 * inv_norm2
                        cos_angle = clamp(cos_angle, -1.0f0, 1.0f0)
                        total_curv += 1.0f0 - abs(cos_angle)
                        count += 1
                    end
                end
            end

            curvatures[v] = count > 0 ? total_curv / count : 1.0f0
        end
    end

    return curvatures
end

"""
    find_extremes_curvature(track, coords)

Optimized curvature-based extreme finding.
"""
function find_extremes_curvature(track::Tracks, coords::TrackCoords)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices < 3
        return find_extremes_distance_fallback_opt(track, coords, 0.7f0)
    end

    curvatures = calculate_vertex_curvatures(track, coords)
    sorted_indices = sortperm(curvatures)

    max_score = 0.0f0
    best_pair = (sorted_indices[1], sorted_indices[2])
    best_path = Int[]
    best_path_length = 0.0f0
    confidence = 0.0f0

    n_candidates = min(6, n_vertices)

    for i in 1:n_candidates-1
        for j in i+1:n_candidates
            v1, v2 = sorted_indices[i], sorted_indices[j]

            # Skip if too close
            if euclidean_distance_coords(coords, v1, v2) < 2.0f0
                continue
            end

            path = find_path_bfs(g, v1, v2)
            if !isempty(path)
                path_length = calculate_path_length(coords, path)
                curvature_bonus = 2.0f0 / (curvatures[v1] + curvatures[v2] + 1f-6)
                score = path_length + curvature_bonus

                if score > max_score
                    max_score = score
                    best_pair = (v1, v2)
                    best_path = path
                    best_path_length = path_length
                    confidence = 0.7f0 + 0.2f0 * curvature_bonus / 10.0f0
                end
            end
        end
    end

    return (best_pair[1], best_pair[2], best_path, min(confidence, 0.9f0), best_path_length)
end

"""
    find_extremes_spatial(track, coords)

Optimized spatial-based extreme finding for dense tracks.
"""
function find_extremes_spatial(track::Tracks, coords::TrackCoords)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0f0, 0.0f0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0f0, 0.0f0)
    end

    # Spatial extreme candidates
    spatial_candidates = [
        argmin(coords.x), argmax(coords.x),
        argmin(coords.y), argmax(coords.y),
        argmin(coords.z), argmax(coords.z)
    ]

    # Add minimum degree vertices
    min_deg = minimum(degree(g, v) for v in 1:n_vertices)
    for v in 1:n_vertices
        if degree(g, v) == min_deg
            push!(spatial_candidates, v)
        end
    end

    all_candidates = unique(spatial_candidates)

    max_path_length = 0.0f0
    best_pair = (all_candidates[1], all_candidates[min(2, length(all_candidates))])
    best_path = Int[]

    for i in 1:length(all_candidates)-1
        for j in i+1:length(all_candidates)
            v1, v2 = all_candidates[i], all_candidates[j]
            path = find_path_bfs(g, v1, v2)
            if !isempty(path)
                path_length = calculate_path_length(coords, path)
                if path_length > max_path_length
                    max_path_length = path_length
                    best_pair = (v1, v2)
                    best_path = path
                end
            end
        end
    end

    # Calculate confidence
    confidence = 0.75f0
    if !isempty(best_path) && max_path_length > 0
        straight_dist = euclidean_distance_coords(coords, best_pair[1], best_pair[2])
        efficiency = straight_dist / max_path_length

        if efficiency > 0.8f0
            confidence += 0.1f0
        end
        confidence += min(0.1f0, max_path_length / 500.0f0)

        deg1 = degree(g, best_pair[1])
        deg2 = degree(g, best_pair[2])
        avg_endpoint_degree = (deg1 + deg2) / 2.0f0
        avg_track_degree = 2f0 * ne(g) / n_vertices

        if avg_endpoint_degree < avg_track_degree * 0.7f0
            confidence += 0.05f0
        end
        confidence = min(confidence, 0.95f0)
    end

    return (best_pair[1], best_pair[2], best_path, confidence, max_path_length)
end

"""
    find_extremes_combined(track, coords; use_energy_weighting=true, dense_track_threshold=6.0)

Optimized combined approach with early-exit strategies.

# Key optimizations
1. Pre-extracted coordinates passed through
2. Path lengths returned with results (no recomputation)
3. Early exit when topology gives high confidence
4. Energy-weighted method for dense tracks (Bragg peak detection)

# Arguments
- `track::Tracks`: Track to analyze
- `coords::TrackCoords`: Pre-extracted coordinates
- `use_energy_weighting::Bool`: Enable energy-weighted extreme finding (default: true)
- `use_edge_energy_weighting::Bool`: Enable edge-energy-weighted method (default: true)
- `use_mst_fallback::Bool`: Enable MST-based fallback method (default: false)
- `dense_track_threshold::Float32`: Average degree above which track is considered "dense"
  and energy-weighted methods are preferred. Default 6.0 corresponds to 3D lattice connectivity.
"""
function find_extremes_combined(track::Tracks, coords::TrackCoords;
                                    use_energy_weighting::Bool=true,
                                    use_edge_energy_weighting::Bool=true,
                                    use_mst_fallback::Bool=false,
                                    dense_track_threshold::Float32=6.0f0)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0f0)
    end

    # Check track density
    avg_degree = 2f0 * ne(g) / n_vertices
    is_dense = avg_degree > dense_track_threshold

    # Always try topology first - it's fast and often sufficient
    topo_result = find_extremes_topology(track, coords)
    topo_confidence = topo_result[4]
    topo_length = topo_result[5]

    # Early exit for high-confidence topology result (2 degree-1 vertices found)
    if topo_confidence >= 0.85f0
        return (topo_result[1], topo_result[2], topo_result[3], topo_confidence)
    end

    if is_dense
        # Dense track: use energy-weighted method (best for Bragg peak detection)

        # Dense track: energy-aware traversal can suppress geometric shortcuts
        if use_edge_energy_weighting
            edge_result = find_extremes_edge_energy_weighted_opt(track, coords)
            edge_confidence = edge_result[4]
            edge_length = edge_result[5]
            if edge_confidence >= 0.85f0 || edge_length > topo_length * 1.1f0
                return (edge_result[1], edge_result[2], edge_result[3], edge_confidence)
            end
        end

        # Compute energy result once if enabled
        energy_result = nothing
        energy_length = 0.0f0
        if use_energy_weighting
            energy_result = find_extremes_energy_weighted(track, coords)
            energy_length = energy_result[5]

            # Energy method is preferred for dense tracks - early exit if high confidence
            if energy_result[4] >= 0.85f0 || energy_length > topo_length * 1.1f0
                return (energy_result[1], energy_result[2], energy_result[3], energy_result[4])
            end
        end

        # Fallback to spatial method
        spatial_result = find_extremes_spatial(track, coords)
        spatial_length = spatial_result[5]

        if spatial_length > topo_length * 1.2f0
            return (spatial_result[1], spatial_result[2], spatial_result[3], spatial_result[4])
        end

        # Try curvature as final fallback
        curv_result = find_extremes_curvature(track, coords)
        curv_length = curv_result[5]

        # MST fallback if enabled
        if use_mst_fallback
            mst_result = find_extremes_mst_diameter(track, coords)
            mst_length = mst_result[5]
            max_other = use_energy_weighting ? max(spatial_length, curv_length, topo_length, energy_length) :
                                               max(spatial_length, curv_length, topo_length)
            if mst_length >= max_other
                return (mst_result[1], mst_result[2], mst_result[3], 1.0f0)
            end
        end

        # Pick the best by path length
        if use_energy_weighting && energy_length >= spatial_length && energy_length >= curv_length && energy_length >= topo_length
            return (energy_result[1], energy_result[2], energy_result[3], energy_result[4])
        elseif spatial_length >= topo_length && spatial_length >= curv_length
            return (spatial_result[1], spatial_result[2], spatial_result[3], spatial_result[4])
        elseif curv_length > topo_length
            return (curv_result[1], curv_result[2], curv_result[3], curv_result[4])
        else
            return (topo_result[1], topo_result[2], topo_result[3], topo_confidence)
        end
    else
        # Sparse track: topology + curvature comparison
        if topo_confidence >= 0.8f0
            return (topo_result[1], topo_result[2], topo_result[3], topo_confidence)
        end

        curv_result = find_extremes_curvature(track, coords)
        curv_length = curv_result[5]

        if curv_length > topo_length * 1.1f0
            return (curv_result[1], curv_result[2], curv_result[3], curv_result[4])
        else
            return (topo_result[1], topo_result[2], topo_result[3], topo_confidence)
        end
    end
end

"""
    find_extremes_energy_weighted(track, coords; min_coverage=0.6)

Energy-weighted extreme finding. Bragg peaks at track endpoints have high energy.
Combines spatial extent with energy to find true endpoints.

Confidence calibration based on:
1. Energy at endpoints (Bragg peak detection)
2. Path coverage (path_length / track_extent)
3. Endpoint separation quality
"""
function find_extremes_energy_weighted(track::Tracks, coords::TrackCoords; min_coverage::Float32=0.6f0)
    g = track.graph
    n_vertices = nv(g)
    vox = track.voxels

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0f0, 0.0f0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0f0, 0.0f0)
    end

    # Get energy values
    energies = collect(vox.energy)
    mean_energy = mean(energies)
    max_energy = maximum(energies)

    # Track spatial extent (diagonal of bounding box)
    x_range = maximum(coords.x) - minimum(coords.x)
    y_range = maximum(coords.y) - minimum(coords.y)
    z_range = maximum(coords.z) - minimum(coords.z)
    track_extent = sqrt(x_range^2 + y_range^2 + z_range^2)

    # Find high-energy voxels (potential Bragg peaks)
    # Threshold: voxels with energy > mean + 0.5*(max-mean)
    energy_threshold = mean_energy + 0.5f0 * (max_energy - mean_energy)
    high_energy_indices = findall(e -> e >= energy_threshold, energies)

    # Also include spatial extremes
    spatial_candidates = [
        argmin(coords.x), argmax(coords.x),
        argmin(coords.y), argmax(coords.y),
        argmin(coords.z), argmax(coords.z)
    ]

    # Combine candidates: high energy + spatial extremes
    all_candidates = unique(vcat(high_energy_indices, spatial_candidates))

    # Limit candidates to avoid O(n²) explosion
    if length(all_candidates) > 20
        sorted_by_energy = sortperm(energies[high_energy_indices], rev=true)
        top_energy = high_energy_indices[sorted_by_energy[1:min(10, length(sorted_by_energy))]]
        all_candidates = unique(vcat(top_energy, spatial_candidates))
    end

    # Score each pair by: path_length * (1 + energy_factor) * coverage_bonus
    best_score = -Inf32
    best_pair = (all_candidates[1], all_candidates[min(2, length(all_candidates))])
    best_path = Int[]
    best_path_length = 0.0f0
    best_coverage = 0.0f0

    for i in 1:length(all_candidates)-1
        for j in i+1:length(all_candidates)
            v1, v2 = all_candidates[i], all_candidates[j]

            # Skip if too close
            dist = euclidean_distance_coords(coords, v1, v2)
            if dist < 5.0f0  # mm, minimum separation
                continue
            end

            path = find_path_bfs(g, v1, v2)
            if isempty(path)
                continue
            end

            path_length = calculate_path_length(coords, path)

            # Coverage: how much of the track extent does this path cover?
            coverage = path_length / max(track_extent, 1.0f0)

            # Skip pairs with poor coverage
            if coverage < min_coverage
                continue
            end

            # Energy bonus: normalized energy at both endpoints
            e1_norm = energies[v1] / max_energy
            e2_norm = energies[v2] / max_energy
            energy_factor = 0.5f0 * (e1_norm + e2_norm)

            # Coverage bonus: reward paths that cover more of the track
            coverage_bonus = min(coverage, 1.5f0)  # cap at 1.5

            # Combined score
            score = path_length * (1.0f0 + energy_factor) * coverage_bonus

            if score > best_score
                best_score = score
                best_pair = (v1, v2)
                best_path = path
                best_path_length = path_length
                best_coverage = coverage
            end
        end
    end

    # If no valid pair found with min_coverage, relax and try again
    if best_score == -Inf32
        return find_extremes_energy_weighted_relaxed(track, coords, energies, energy_threshold,
                                                      max_energy, track_extent, all_candidates)
    end

    # Calibrated confidence based on multiple factors
    e1 = energies[best_pair[1]]
    e2 = energies[best_pair[2]]
    both_high_energy = (e1 >= energy_threshold) && (e2 >= energy_threshold)
    one_high_energy = (e1 >= energy_threshold) || (e2 >= energy_threshold)

    # Base confidence from energy
    if both_high_energy
        conf_energy = 0.4f0
    elseif one_high_energy
        conf_energy = 0.25f0
    else
        conf_energy = 0.1f0
    end

    # Confidence from coverage (0.0 to 0.4)
    conf_coverage = min(0.4f0, 0.4f0 * best_coverage / 1.0f0)

    # Confidence from path quality (0.0 to 0.2)
    # Straight-line efficiency: how direct is the path?
    straight_dist = euclidean_distance_coords(coords, best_pair[1], best_pair[2])
    efficiency = straight_dist / max(best_path_length, 1.0f0)
    conf_efficiency = 0.2f0 * min(1.0f0, efficiency / 0.8f0)  # max at 80% efficiency

    # Total confidence (0.0 to 1.0)
    confidence = min(0.95f0, conf_energy + conf_coverage + conf_efficiency)

    return (best_pair[1], best_pair[2], best_path, confidence, best_path_length)
end


"""
Fallback when min_coverage cannot be satisfied.
Returns results with LOW confidence to flag unreliable extremes.
"""
function find_extremes_energy_weighted_relaxed(track::Tracks, coords::TrackCoords,
                                                energies::Vector{Float32}, energy_threshold::Float32,
                                                max_energy::Float32, track_extent::Float32,
                                                all_candidates::Vector{Int})
    g = track.graph

    best_score = -Inf32
    best_pair = (all_candidates[1], all_candidates[min(2, length(all_candidates))])
    best_path = Int[]
    best_path_length = 0.0f0

    for i in 1:length(all_candidates)-1
        for j in i+1:length(all_candidates)
            v1, v2 = all_candidates[i], all_candidates[j]

            dist = euclidean_distance_coords(coords, v1, v2)
            if dist < 3.0f0
                continue
            end

            path = find_path_bfs(g, v1, v2)
            if isempty(path)
                continue
            end

            path_length = calculate_path_length(coords, path)
            e1_norm = energies[v1] / max_energy
            e2_norm = energies[v2] / max_energy
            energy_factor = 0.5f0 * (e1_norm + e2_norm)
            score = path_length * (1.0f0 + energy_factor)

            if score > best_score
                best_score = score
                best_pair = (v1, v2)
                best_path = path
                best_path_length = path_length
            end
        end
    end

    # LOW confidence for relaxed results - coverage-dependent
    # This flags events where we couldn't find good extremes
    coverage = best_path_length / max(track_extent, 1.0f0)
    confidence = min(0.5f0, 0.2f0 + 0.3f0 * coverage)  # max 0.5 for relaxed

    return (best_pair[1], best_pair[2], best_path, confidence, best_path_length)
end


"""
    find_extremes_distance_fallback_opt(track, coords, confidence)

Optimized fallback using maximum Euclidean distance.
"""
function find_extremes_distance_fallback_opt(track::Tracks, coords::TrackCoords, confidence::Float32)
    g = track.graph
    n_vertices = nv(g)

    max_dist = 0.0f0
    best_pair = (1, min(2, n_vertices))

    # For small tracks, check all pairs
    # For larger tracks, use spatial extremes as candidates
    if n_vertices <= 50
        for i in 1:n_vertices-1
            for j in i+1:n_vertices
                dist = euclidean_distance_coords(coords, i, j)
                if dist > max_dist
                    max_dist = dist
                    best_pair = (i, j)
                end
            end
        end
    else
        # Use spatial extremes for large tracks
        candidates = unique([
            argmin(coords.x), argmax(coords.x),
            argmin(coords.y), argmax(coords.y),
            argmin(coords.z), argmax(coords.z)
        ])

        for i in 1:length(candidates)-1
            for j in i+1:length(candidates)
                v1, v2 = candidates[i], candidates[j]
                dist = euclidean_distance_coords(coords, v1, v2)
                if dist > max_dist
                    max_dist = dist
                    best_pair = (v1, v2)
                end
            end
        end
    end

    path = find_path_bfs(g, best_pair[1], best_pair[2])
    path_length = calculate_path_length(coords, path)
    return (best_pair[1], best_pair[2], path, confidence, path_length)
end

############################
# Energy-weighted traversal #
############################

"""
    energy_weight_matrix(g, coords, energy; epsE=1e-6, α=1.0, β=1.0) -> SparseMatrixCSC{Float32,Int}

Build a sparse weight matrix for `dijkstra_shortest_paths` on an unweighted graph `g`.

Edge cost model (default):
    w_ij = (d_ij^α) / ((0.5*(E_i + E_j) + epsE)^β)

- `coords` must be 3×N (Float32) with columns aligned to graph vertices.
- `energy` is a length-N vector (Float32), typically `track.voxels.energy` or `track.voxels.electrons`.
- `epsE` prevents blow-up when energies are small.
- `α` controls geometric emphasis (α=1 is linear distance).
- `β` controls energy emphasis (β=1 is inverse-energy).
"""
function energy_weight_matrix(g::SimpleGraph,
                              coords::AbstractMatrix{<:Real},
                              energy::AbstractVector{<:Real};
                              epsE::Float32 = 1f-6,
                              α::Float32 = 1.0f0,
                              β::Float32 = 1.0f0)

    n = nv(g)
    I = Int[]
    J = Int[]
    V = Float32[]

    @inbounds for e in edges(g)
        i = src(e); j = dst(e)

        dx = float(coords[1,i] - coords[1,j])
        dy = float(coords[2,i] - coords[2,j])
        dz = float(coords[3,i] - coords[3,j])
        d  = sqrt(dx*dx + dy*dy + dz*dz)

        Ej = 0.5f0*(float(energy[i]) + float(energy[j])) + epsE
        w  = (d^α) / (Ej^β)

        # symmetric
        push!(I,i); push!(J,j); push!(V,w)
        push!(I,j); push!(J,i); push!(V,w)
    end

    return sparse(I, J, V, n, n)
end


"""
    dijkstra_path(g, src, dst, W) -> Vector{Int}

Reconstruct a shortest path from `src` to `dst` using Dijkstra parents.
Returns an empty vector if `dst` is unreachable.
"""
function dijkstra_path(g::SimpleGraph, srcv::Int, dstv::Int, W)
    sp = dijkstra_shortest_paths(g, srcv, W)
    parents = sp.parents

    dstv > length(parents) && return Int[]
    isinf(sp.dists[dstv]) && return Int[]

    path = Int[dstv]
    v = dstv
    while v != srcv
        p = parents[v]
        (p == 0 || p == v) && return Int[]  # unreachable / broken parent chain
        push!(path, p)
        v = p
    end
    reverse!(path)
    return path
end


"""
    find_extremes_edge_energy_weighted_opt(track, coords; epsE=1e-6, α=1.0, β=1.0)
        -> (extreme1, extreme2, path, confidence, path_length)

Energy-weighted *edge-cost* extreme finding:
1) Build sparse weight matrix from geometry + voxel energies
2) Double-sweep using weighted Dijkstra distances
3) Return endpoints and the weighted-shortest path between them

This is specifically designed to suppress "geometric shortcuts" that run through
low-support (low-energy) regions while preserving the connected manifold.
"""
function find_extremes_edge_energy_weighted_opt(track::Tracks, coords::TrackCoords;
                                               epsE::Float32 = 1f-6,
                                               α::Float32 = 1.0f0,
                                               β::Float32 = 1.0f0)

    g = track.graph
    n = nv(g)
    n == 0 && return (nothing, nothing, Int[], 0.0f0, 0.0f0)
    n == 1 && return (1, 1, [1], 1.0f0, 0.0f0)

    # 3×N coordinate matrix aligned to graph vertices
    coord_matrix = hcat(coords.x, coords.y, coords.z)'

    # energy vector aligned to vertices (use energy by default)
    E = Float32.(track.voxels.energy)

    # robust epsilon scaled to typical energy
    epsE_eff = max(epsE, 1f-6 * max(maximum(E), 1f-12))

    W = energy_weight_matrix(g, coord_matrix, E; epsE=epsE_eff, α=α, β=β)

    # Double sweep in weighted metric
    d1 = dijkstra_shortest_paths(g, 1, W).dists
    u = argmax(d1)
    d2 = dijkstra_shortest_paths(g, u, W).dists
    v = argmax(d2)

    path = dijkstra_path(g, u, v, W)
    isempty(path) && return (u, v, Int[], 0.0f0, 0.0f0)

    # Physical (Euclidean) path length along the returned vertex sequence
    path_length = calculate_path_length_from_coords(coord_matrix, path)

    # Confidence: combine coverage proxy + endpoint separation proxy
    # (keep it simple and diagnostic-driven; you can refine later)
    D_end = euclidean_distance(coords.x[u], coords.y[u], coords.z[u],
                               coords.x[v], coords.y[v], coords.z[v])
    η = (D_end > 0) ? (path_length / D_end) : Inf32
    # η>1 indicates non-trivial curvature; extremely small η indicates shortcutting
    confidence = clamp(0.5f0 + 0.5f0 * tanh((η - 1.05f0) / 0.25f0), 0.0f0, 1.0f0)

    return (u, v, path, confidence, path_length)
end


#################
# MST + diameter #
#################

"""
    compute_mst_graph(g, coords) -> SimpleGraph

Compute an MST (minimum spanning tree) as an unweighted `SimpleGraph` using
Euclidean distances for edges already present in `g`.

Note: MST can be used as a fallback/topology-cleaner for pathological cases.
"""
function compute_mst_graph(g::SimpleGraph, coords::AbstractMatrix{<:Real})
    n = nv(g)
    (n <= 1 || ne(g) == 0) && return copy(g)

    # Build a dense distance matrix only for existing edges; others are ignored by `kruskal_mst`.
    distmx = fill(Inf32, n, n)
    @inbounds for e in edges(g)
        i = src(e); j = dst(e)
        dx = float(coords[1,i] - coords[1,j])
        dy = float(coords[2,i] - coords[2,j])
        dz = float(coords[3,i] - coords[3,j])
        d  = sqrt(dx*dx + dy*dy + dz*dz)
        distmx[i,j] = d
        distmx[j,i] = d
    end

    mst_edges = kruskal_mst(g, distmx)  # returns an edge iterator/collection
    mst = SimpleGraph(n)
    for e in mst_edges
        add_edge!(mst, src(e), dst(e))
    end
    return mst
end


"""
    find_extremes_mst_diameter(track, coords) -> (extreme1, extreme2, path, confidence, path_length)

Fallback method:
1) Build MST from `track.graph`
2) Find diameter endpoints by double BFS (unweighted)
3) Return BFS path and Euclidean path length
"""
function find_extremes_mst_diameter(track::Tracks, coords::TrackCoords)
    g = track.graph
    n = nv(g)
    n == 0 && return (nothing, nothing, Int[], 0.0f0, 0.0f0)
    n == 1 && return (1, 1, [1], 1.0f0, 0.0f0)

    coord_matrix = hcat(coords.x, coords.y, coords.z)'
    mst = compute_mst_graph(g, coord_matrix)

    d1 = gdistances(mst, 1)
    u = argmax(d1)
    d2 = gdistances(mst, u)
    v = argmax(d2)

    path = find_path_bfs(mst, u, v)
    path_length = calculate_path_length_from_coords(coord_matrix, path)

    return (u, v, path, 1.0f0, path_length)
end
