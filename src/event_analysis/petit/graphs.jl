# Adapted from Petit src/track_building.jl (MIT; see LICENSE).
# Source SHA-256: 1510e17ddce2f75b85820ff1c80828355d7981447e4e596f04ff5873f818a577
# Float32 adaptation; algorithm branches retained. See ../PROVENANCE.toml.

# helper: coord matrix 3×N
function _coords_from_df(df::DataFrame)
    n = nrow(df)
    coords = Matrix{Float32}(undef, 3, n)
    @inbounds for i in 1:n
        coords[1, i] = df.x[i]
        coords[2, i] = df.y[i]
        coords[3, i] = df.z[i]
    end
    return coords
end

function build_tracks_kdtree(event_data::DataFrame;
                             max_distance::Float32 = 1.5f0,
                             energy_threshold::Float32 = 0.0f0,
                             diffusion::DiffusionParams = DiffusionParams())

    n_voxels = nrow(event_data)
    n_voxels == 0 && return Tracks[]

    # filter by energy
    valid_indices = findall(e -> e >= energy_threshold, event_data.energy)
    isempty(valid_indices) && return Tracks[]

    filtered_data = event_data[valid_indices, :]
    n_filtered = nrow(filtered_data)

    # trivial case
    n_filtered == 1 && return [Tracks(filtered_data, SimpleGraph(1), [collect(1:1)], diffusion)]

    # coords and KDTree
    coords = _coords_from_df(filtered_data)
    tree   = KDTree(coords)

    g = SimpleGraph(n_filtered)
    r  = max_distance

    # build graph using neighbors in range
    @inbounds for i in 1:n_filtered
        idxs = inrange(tree, coords[:, i], r, true)
        for j in idxs
            j == i && continue
            j <  i && continue   # avoid double edges; we'll see (i,j) once with i<j
            add_edge!(g, i, j)
        end
    end

    # components → tracks
    components = connected_components(g)
    tracks = Tracks[]

    for comp in components
        isempty(comp) && continue

        comp_data = filtered_data[comp, :]

        # induced subgraph for this component (returns tuple: subgraph, vertex_map)
        subg, _ = induced_subgraph(g, comp)

        # Tracks constructor: data, graph, list of components, diffusion params
        push!(tracks, Tracks(comp_data, subg, [comp], diffusion))
    end

    return tracks
end


"""
    build_tracks_knn(event_data::DataFrame; k=10, max_distance=Inf,
                     energy_threshold=0.0, diffusion=DiffusionParams(),
                     mutual=false)

Build tracks using k-Nearest Neighbor graph instead of radius graph.

# Arguments
- `event_data::DataFrame`: Voxelized event data with x, y, z, energy columns
- `k::Int`: Number of nearest neighbors per voxel (default: 10)
- `max_distance::Float32`: Maximum edge length allowed (default: Inf, no limit)
- `energy_threshold::Float32`: Minimum voxel energy in MeV (default: 0.0)
- `diffusion::DiffusionParams`: Diffusion parameters for the tracks
- `mutual::Bool`: If true, keep only mutual kNN edges (default: false)

# Returns
- `Vector{Tracks}`: Vector of tracks (connected components)

# Notes
kNN limits outgoing neighbour choices to the nearest points; the finite distance
cutoff additionally bounds every accepted edge. Incoming choices can increase
a vertex degree above k.
"""
function build_tracks_knn(event_data::DataFrame;
                          k::Int = 10,
                          max_distance::Float32 = Inf32,
                          energy_threshold::Float32 = 0.0f0,
                          diffusion::DiffusionParams = DiffusionParams(),
                          mutual::Bool = false)

    n_voxels = nrow(event_data)
    n_voxels == 0 && return Tracks[]

    # Filter by energy
    valid_indices = findall(e -> e >= energy_threshold, event_data.energy)
    isempty(valid_indices) && return Tracks[]

    filtered_data = event_data[valid_indices, :]
    n_filtered = nrow(filtered_data)

    # Trivial case: single voxel
    n_filtered == 1 && return [Tracks(filtered_data, SimpleGraph(1), [collect(1:1)], diffusion)]

    # Adjust k if we have fewer points
    k_actual = min(k, n_filtered - 1)
    k_actual < 1 && return [Tracks(filtered_data, SimpleGraph(n_filtered), [collect(1:n_filtered)], diffusion)]

    # Build KDTree
    coords = _coords_from_df(filtered_data)
    tree = KDTree(coords)

    # Build kNN graph
    g = SimpleGraph(n_filtered)

    if mutual
        # Mutual kNN: keep edge only if both endpoints have each other in their kNN
        # First, collect all kNN sets
        knn_sets = Vector{Set{Int}}(undef, n_filtered)
        @inbounds for i in 1:n_filtered
            idxs, dists = knn(tree, coords[:, i], k_actual + 1, true)
            # Filter by max_distance and exclude self
            neighbors = Set{Int}()
            for (idx, d) in zip(idxs, dists)
                idx == i && continue
                if d <= max_distance
                    push!(neighbors, idx)
                end
            end
            knn_sets[i] = neighbors
        end

        # Add edge only if mutual
        @inbounds for i in 1:n_filtered
            for j in knn_sets[i]
                j <= i && continue  # Avoid duplicate edges
                if i in knn_sets[j]
                    add_edge!(g, i, j)
                end
            end
        end
    else
        # Standard kNN: add edge if j is in kNN(i)
        # Since we process all i, edges become symmetric
        @inbounds for i in 1:n_filtered
            idxs, dists = knn(tree, coords[:, i], k_actual + 1, true)
            for (idx, d) in zip(idxs, dists)
                idx == i && continue
                if d <= max_distance
                    add_edge!(g, i, idx)  # SimpleGraph handles duplicates
                end
            end
        end
    end

    # Connected components → tracks
    components = connected_components(g)
    tracks = Tracks[]

    for comp in components
        isempty(comp) && continue
        comp_data = filtered_data[comp, :]
        subg, _ = induced_subgraph(g, comp)
        push!(tracks, Tracks(comp_data, subg, [comp], diffusion))
    end

    return tracks
end


"""
    make_tracks(event_data; max_distance_mm, energy_threshold_kev=0f0,
                diffusion=DiffusionParams(), method="kNN", k=10)

Build and sort tracks by total energy (descending).

# Arguments
- `event_data::DataFrame`: Voxelized event data
- `max_distance_mm::Float32`: Max voxel connection distance in mm
- `energy_threshold_kev::Float32`: Min voxel energy in keV
- `diffusion::DiffusionParams`: Diffusion parameters for the tracks
- `method::String`: Graph construction method:
  - "KDT": Radius graph (all neighbors within max_distance)
  - "kNN": k-Nearest Neighbor graph (k closest neighbors)
  - "kNN_mutual": Mutual kNN (edge only if both are in each other's kNN)
- `k::Int`: Number of neighbors for kNN methods (default: 10)

# Returns
- `Vector{Tracks}`: Tracks sorted by total energy (highest first)
"""
function make_tracks(event_data::DataFrame;
                    max_distance_mm::Float32,
                    energy_threshold_kev::Float32=0.0f0,
                    diffusion::DiffusionParams=DiffusionParams(),
                    method::String="kNN",
                    k::Int=10)

    validate_bin_table(event_data)
    isfinite(max_distance_mm) && max_distance_mm > 0f0 || throw(ArgumentError("max_distance_mm must be finite and positive"))
    isfinite(energy_threshold_kev) && energy_threshold_kev >= 0f0 || throw(ArgumentError("energy_threshold_kev must be finite and nonnegative"))
    k >= 1 || throw(ArgumentError("k must be positive"))
    energy_threshold = energy_threshold_kev * 1f-3  # Convert keV to MeV

    valid_methods = ("KDT", "kNN", "kNN_mutual")
    if !(method in valid_methods)
        throw(ArgumentError("Invalid method '$method'. Valid options: $(join(valid_methods, ", "))"))
    end

    if method == "KDT"
        tracks = build_tracks_kdtree(event_data;
                                     max_distance=max_distance_mm,
                                     energy_threshold=energy_threshold,
                                     diffusion=diffusion)
    elseif method == "kNN"
        tracks = build_tracks_knn(event_data;
                                  k=k,
                                  max_distance=max_distance_mm,
                                  energy_threshold=energy_threshold,
                                  diffusion=diffusion,
                                  mutual=false)
    else  # kNN_mutual
        tracks = build_tracks_knn(event_data;
                                  k=k,
                                  max_distance=max_distance_mm,
                                  energy_threshold=energy_threshold,
                                  diffusion=diffusion,
                                  mutual=true)
    end

    # Sort by total energy (descending)
    if length(tracks) > 0
        track_energies = [sum(track.voxels.energy) for track in tracks]
        sorted_indices = sortperm(track_energies, rev=true)
        tracks = tracks[sorted_indices]
    end

    return tracks
end
