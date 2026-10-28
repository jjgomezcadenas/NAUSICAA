# Adapted from Petit src/mctrack.jl (MIT; see LICENSE).
# Source SHA-256: 9e2e481adfe159011af571839d24014287cabe0f86a0716e87d7938cfa4dec81
# Float32 adaptation; algorithm branches retained. See ../PROVENANCE.toml.

"""
    voxelize_particle_hits(hits_df, mcvox_size)

Helper function to voxelize hits from a single particle.
Returns vectors (xs, ys, zs, Es) of voxelized positions and energies.
"""
function voxelize_particle_hits(hits_df::DataFrame, mcvox_size::Float32)
    if nrow(hits_df) == 0
        return Float32[], Float32[], Float32[], Float32[]
    end

    hits = copy(hits_df)
    hits.vox_ix = floor.(Int, hits.x ./ mcvox_size)
    hits.vox_iy = floor.(Int, hits.y ./ mcvox_size)
    hits.vox_iz = floor.(Int, hits.z ./ mcvox_size)

    # Group consecutive hits in same voxel
    n_hits = nrow(hits)
    voxel_groups = Vector{UnitRange{Int}}()

    i = 1
    while i <= n_hits
        start_i = i
        current_vox = (hits.vox_ix[i], hits.vox_iy[i], hits.vox_iz[i])

        while i <= n_hits &&
              (hits.vox_ix[i], hits.vox_iy[i], hits.vox_iz[i]) == current_vox
            i += 1
        end

        push!(voxel_groups, start_i:(i-1))
    end

    # Compute voxel positions and energies
    n_voxels = length(voxel_groups)
    xs = Vector{Float32}(undef, n_voxels)
    ys = Vector{Float32}(undef, n_voxels)
    zs = Vector{Float32}(undef, n_voxels)
    Es = Vector{Float32}(undef, n_voxels)

    for (k, group) in enumerate(voxel_groups)
        hits_in_group = @view hits[group, :]
        total_E = sum(hits_in_group.energy)

        if total_E > 0
            xs[k] = sum(hits_in_group.x .* hits_in_group.energy) / total_E
            ys[k] = sum(hits_in_group.y .* hits_in_group.energy) / total_E
            zs[k] = sum(hits_in_group.z .* hits_in_group.energy) / total_E
        else
            xs[k] = mean(hits_in_group.x)
            ys[k] = mean(hits_in_group.y)
            zs[k] = mean(hits_in_group.z)
        end
        Es[k] = total_E
    end

    return xs, ys, zs, Es
end

"""
    _compute_mc_path(event_df, mcvox_size, primary_ids)

Petit truth-path kernel using one or two supplied Int32 primary IDs. Input uses
x/y/z in mm, energy in MeV and time in ns. Sort each primary by time; group only
consecutive visits to a cell. Reverse the first of two arms before concatenating.
Return owned Float32 x/y/z/energy/s columns and actual particle_id::Int32.
The public EventHits wrapper validates arguments; absent primary hits return a
typed empty table. All input tables remain unchanged.
"""
function _compute_mc_path(event_df::DataFrame, mcvox_size::Float32, primary_ids::Vector{Int32})
    is_double_beta = length(primary_ids) == 2
    if is_double_beta
        # Double-beta decay: combine both electrons
        hits1 = filter(row -> row.particle_id == primary_ids[1], event_df)
        hits2 = filter(row -> row.particle_id == primary_ids[2], event_df)

        if nrow(hits1) == 0 || nrow(hits2) == 0
            return DataFrame(x=Float32[], y=Float32[], z=Float32[],
                            energy=Float32[], s=Float32[], particle_id=Int32[])
        end

        # Sort each by time
        sort!(hits1, :time)
        sort!(hits2, :time)

        # Voxelize each particle separately
        x1, y1, z1, E1 = voxelize_particle_hits(hits1, mcvox_size)
        x2, y2, z2, E2 = voxelize_particle_hits(hits2, mcvox_size)

        if isempty(x1) || isempty(x2)
            return DataFrame(x=Float32[], y=Float32[], z=Float32[],
                            energy=Float32[], s=Float32[], particle_id=Int32[])
        end

        # Create continuous path: electron1 reversed (Bragg→vertex) + electron2 (vertex→Bragg)
        # Reverse electron 1
        x1_rev = reverse(x1)
        y1_rev = reverse(y1)
        z1_rev = reverse(z1)
        E1_rev = reverse(E1)

        # Concatenate: electron1_reversed + electron2
        xs = vcat(x1_rev, x2)
        ys = vcat(y1_rev, y2)
        zs = vcat(z1_rev, z2)
        Es = vcat(E1_rev, E2)

        # Track which electron each voxel belongs to
        pe = vcat(fill(primary_ids[1], length(x1)), fill(primary_ids[2], length(x2)))
    else
        # Single electron: use its supplied primary ID
        primary_hits = filter(row -> row.particle_id == primary_ids[1], event_df)

        if nrow(primary_hits) == 0
            return DataFrame(x=Float32[], y=Float32[], z=Float32[],
                            energy=Float32[], s=Float32[], particle_id=Int32[])
        end

        sort!(primary_hits, :time)
        xs, ys, zs, Es = voxelize_particle_hits(primary_hits, mcvox_size)

        if isempty(xs)
            return DataFrame(x=Float32[], y=Float32[], z=Float32[],
                            energy=Float32[], s=Float32[], particle_id=Int32[])
        end

        # All path rows belong to the supplied primary
        pe = fill(primary_ids[1], length(xs))
    end

    # Compute cumulative arc-length
    n_voxels = length(xs)
    ss = Vector{Float32}(undef, n_voxels)
    ss[1] = 0.0f0
    for k in 2:n_voxels
        dx = xs[k] - xs[k-1]
        dy = ys[k] - ys[k-1]
        dz = zs[k] - zs[k-1]
        ss[k] = ss[k-1] + sqrt(dx^2 + dy^2 + dz^2)
    end

    return DataFrame(x=xs, y=ys, z=zs, energy=Es, s=ss, particle_id=pe)
end
