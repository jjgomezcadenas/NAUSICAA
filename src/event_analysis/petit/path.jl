# Adapted from Petit src/track_path_reconstruction.jl (MIT; see LICENSE).
# Source SHA-256: 9be8b1a1dccf555c4711665724f631543f5b41e08770fe1dfd27ee536983caa3
# Float32 adaptation; algorithm branches retained. See ../PROVENANCE.toml.

function get_raw_path(track::Tracks, path_indices::Vector{Int})
    voxels = track.voxels
    n = length(path_indices)

    xs = voxels.x[path_indices]
    ys = voxels.y[path_indices]
    zs = voxels.z[path_indices]

    # Compute cumulative arc-length
    ss = Vector{Float32}(undef, n)
    n > 0 && (ss[1] = 0f0)
    for i in 2:n
        dx = xs[i] - xs[i-1]
        dy = ys[i] - ys[i-1]
        dz = zs[i] - zs[i-1]
        ss[i] = ss[i-1] + sqrt(dx^2 + dy^2 + dz^2)
    end

    return DataFrame(x=xs, y=ys, z=zs, s=ss)
end
