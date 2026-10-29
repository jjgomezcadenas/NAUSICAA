# Generate a pinned reference from the active Petit diffusion body.
# Only the RNG injection is changed; all arithmetic/binning remains original.
using DataFrames, Random, HDF5, SHA
const ROOT = normpath(joinpath(@__DIR__,"../.."))
const PETIT = isempty(ARGS) ? normpath(joinpath(ROOT,"../Petit")) : abspath(ARGS[1])
module PetitDiffusionReference
using DataFrames, Random
mutable struct ReplayRNG <: AbstractRNG
    draws::Vector{Float64}
    cursor::Int
end
function Random.randn(rng::ReplayRNG)
    rng.cursor += 1
    return rng.draws[rng.cursor]
end
end
source = read(joinpath(PETIT,"src/itaca_functions.jl"),String)
first = findfirst("function diffuse_xyz_image_mc(",source).start
last = findnext("\nend",source,first).stop
body = source[first:last]
body = replace(body,"diffuse_xyz_image_mc(df::DataFrame;"=>"diffuse_xyz_image_mc(rng,df::DataFrame;",
                    "randn()"=>"randn(rng)")
Base.include_string(PetitDiffusionReference,body,"Petit/src/itaca_functions.jl:diffuse_xyz_image_mc")
draws = Float64[-0.3,0.2,0.1, 0.8,-0.7,0.4, 1.3,0.9,-1.1, -5,0,0,
                -0.9,0.8,1.2, 0.2,-0.3,0.6, -1.4,0.1,-0.8, 0,5,0]
input = DataFrame(event_id=Int64[42,42],x=[-1.,1.],y=[-2.,2.],z=[-3.,3.],electrons=Int64[4,4])
result = PetitDiffusionReference.diffuse_xyz_image_mc(PetitDiffusionReference.ReplayRNG(draws,0),input;
    sigma_t_mm=0.5,sigma_l_mm=0.5,nbins=8,nsigma=2.0,ef=50000.0)
path = joinpath(ROOT,"tests/event_analysis/fixtures/petit_diffusion_reference.h5")
h5open(path,"w") do f
    attributes(f)["source_sha256"] = bytes2hex(sha256(source))
    attributes(f)["adaptation"] = "Explicit replay RNG injection only; original Float64 diffusion body"
    f["draws"] = draws
    for name in names(input); f["input/"*name] = input[!,name]; end
    for name in names(result); f["output/"*name] = result[!,name]; end
end
println(path)
