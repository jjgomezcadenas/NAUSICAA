# Pinned comparison against the original active simple-RECO occupancy functions.
using DataFrames, HDF5, SHA
const ROOT = dirname(@__DIR__)
const PETIT = isempty(ARGS) ? normpath(joinpath(ROOT,"../Petit")) : abspath(ARGS[1])
module OccupancyReference
using DataFrames, Statistics, Graphs
module Petit
using DataFrames, Graphs, NearestNeighbors
end
const TIME_SLOT_MM = 2.0
const W_VALUE_EV = 22.0
const ION_THRESHOLD = 10.0
end
Base.include(OccupancyReference.Petit,joinpath(PETIT,"src/track_building.jl"))
source = read(joinpath(PETIT,"itacaScripts/itaca_simple_reco.jl"),String)
for name in ("compute_track_z_extension","compute_barycenter","compute_occupancy_analysis")
    a = findfirst("function "*name*"(",source).start
    b = findnext("\nend",source,a).stop
    Base.include_string(OccupancyReference,source[a:b],"Petit/itacaScripts/itaca_simple_reco.jl:"*name)
end
counts = Int64[20,10,30,40]
d = DataFrame(event_id=fill(Int64(42),4),x=[-3.25,-0.75,1.25,3.75],
              y=[-0.5,0.25,0.5,-0.25],z=[0.25,1.25,2.75,3.5],energy=counts.*22e-6)
t = OccupancyReference.Petit.Tracks(d,OccupancyReference.SimpleGraph(4),[collect(1:4)],OccupancyReference.Petit.DiffusionParams())
r = OccupancyReference.compute_occupancy_analysis(t;pixel_size_mm=2.0,nbins=8)
path = joinpath(ROOT,"tests/event_analysis/fixtures/petit_collection_reference.h5")
h5open(path,"w") do f
    attributes(f)["source_sha256"] = bytes2hex(sha256(source))
    attributes(f)["scope"] = "Original occupancy body, helpers, W=22 eV, 2 mm pitch/slices, 8x8 grid; fixture avoids final-edge exclusion"
    for name in names(d); f["input/"*name] = d[!,name]; end
    f["input/ions"] = counts
    f["total_counts"] = r.total_counts
    f["avg_counts"] = r.avg_counts
    f["n_time_slots"] = r.n_time_slots
    f["occupancy_per_slot"] = r.occupancy_per_slot
end
println(path)
