#=
Regenerate deterministic reference data directly from a local Petit checkout.
Usage: julia --project=. scripts/validation/generate_petit_reference.jl [PETIT_DIRECTORY]

Loads the original graph/backbone functions in an isolated module with the root
project's pinned dependencies. No detector simulation or plotting packages are
loaded. Normal tests read the saved HDF5 and need no Petit checkout. Source hashes
and dependency manifest hash identify the exact reference calculation.
=#
using DataFrames, Graphs, HDF5, SHA
const ROOT = normpath(joinpath(@__DIR__,"../.."))
const PETIT = isempty(ARGS) ? normpath(joinpath(ROOT,"../Petit")) : abspath(ARGS[1])
module PetitReference
using DataFrames, Graphs, NearestNeighbors, LinearAlgebra, SparseArrays, Statistics
euclidean_distance(x1,y1,z1,x2,y2,z2) = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
end
for name in ("track_building.jl","track_extreme_finding.jl")
    Base.include(PetitReference,joinpath(PETIT,"src",name))
end
include(joinpath(ROOT,"tests/event_analysis/reconstruction_fixtures.jl"))
path = joinpath(ROOT,"tests/event_analysis/fixtures/petit_reference.h5")
mkpath(dirname(path))
h5open(path,"w") do f
    attributes(f)["schema"] = "nausicaa.petit.reference.v1"
    attributes(f)["manifest_sha256"] = bytes2hex(sha256(read(joinpath(ROOT,"Manifest.toml"))))
    for name in ("track_building.jl","track_extreme_finding.jl")
        attributes(f)[name*"_sha256"] = bytes2hex(sha256(read(joinpath(PETIT,"src",name))))
    end
    for (name,case) in reconstruction_fixtures(Float64), method in ("KDT","kNN","kNN_mutual")
        tracks = PetitReference.make_tracks(case.data; max_distance_mm=case.cutoff,method,k=case.k)
        group = create_group(f,name*"/"*method)
        group["edges"] = parent_edges(tracks)
        group["component_energy_MeV"] = [sum(t.voxels.energy) for t in tracks]
        group["component_sizes"] = [nv(t.graph) for t in tracks]
        for (i,t) in enumerate(tracks)
            r = PetitReference.walk_track_from_extremes(t)
            g = create_group(group,string(i))
            g["path"] = r.path_indices
            g["length_mm"] = r.total_length
            g["confidence"] = r.confidence
        end
    end
end
println(path)
