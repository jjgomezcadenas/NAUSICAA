"""
    energy_in_spheres_around_extremes(track, walk_result, radius)

Sum all component hits within each endpoint sphere, including its boundary.
Return centres, counts and energies in the input energy units. Blob 1 has
the larger energy. Hits in overlapping spheres contribute to both sums.
Copied from Petit/src/track_blob_analysis.jl; calculation unchanged.
"""
function energy_in_spheres_around_extremes(track::Tracks, walk_result, radius::Float64)

    # Check if we have valid extremes
    if isnothing(walk_result.extremes[1])
        return (start_sphere_energy = 0.0,
                end_sphere_energy = 0.0,
                start_voxel_count = 0,
                end_voxel_count = 0,
                start_center = (NaN, NaN, NaN),
                end_center = (NaN, NaN, NaN))
    end

    start_voxel, end_voxel = walk_result.extremes

    # Get coordinates of extreme voxels
    start_center = (start_voxel.x, start_voxel.y, start_voxel.z)
    end_center = (end_voxel.x, end_voxel.y, end_voxel.z)

    # Initialize counters
    start_sphere_energy = 0.0
    end_sphere_energy = 0.0
    start_voxel_count = 0
    end_voxel_count = 0

    # Check each voxel in the track
    for i in 1:nrow(track.voxels)
        voxel_x = track.voxels.x[i]
        voxel_y = track.voxels.y[i]
        voxel_z = track.voxels.z[i]
        voxel_energy = track.voxels.energy[i]

        # Calculate distance to start voxel
        dist_to_start = euclidean_distance(voxel_x, voxel_y, voxel_z,
                                          start_center[1], start_center[2], start_center[3])

        # Calculate distance to end voxel
        dist_to_end = euclidean_distance(voxel_x, voxel_y, voxel_z,
                                        end_center[1], end_center[2], end_center[3])

        # Check if voxel is within start sphere
        if dist_to_start <= radius
            start_sphere_energy += voxel_energy
            start_voxel_count += 1
        end

        # Check if voxel is within end sphere
        if dist_to_end <= radius
            end_sphere_energy += voxel_energy
            end_voxel_count += 1
        end
    end

    # Determine which sphere has more energy
    # blob1 is the sphere with larger energy
    if start_sphere_energy >= end_sphere_energy
        blob1_energy = start_sphere_energy
        blob1_voxel_count = start_voxel_count
        blob1_center = start_center
        blob2_energy = end_sphere_energy
        blob2_voxel_count = end_voxel_count
        blob2_center = end_center
    else
        blob1_energy = end_sphere_energy
        blob1_voxel_count = end_voxel_count
        blob1_center = end_center
        blob2_energy = start_sphere_energy
        blob2_voxel_count = start_voxel_count
        blob2_center = start_center
    end

    # Return both old format (for compatibility) and new format
    return (start_sphere_energy = start_sphere_energy,
            end_sphere_energy = end_sphere_energy,
            start_voxel_count = start_voxel_count,
            end_voxel_count = end_voxel_count,
            start_center = start_center,
            end_center = end_center,
            # New format: blob1 has higher energy than blob2
            blob1_energy = blob1_energy,
            blob1_voxel_count = blob1_voxel_count,
            blob1_center = blob1_center,
            blob2_energy = blob2_energy,
            blob2_voxel_count = blob2_voxel_count,
            blob2_center = blob2_center)
end
