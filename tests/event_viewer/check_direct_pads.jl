using HDF5, TOML, Random, DataFrames, SHA
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CONFIG=TOML.parsefile(joinpath(ROOT,"config/event_viewer.toml"))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis
include(joinpath(ROOT,"src/event_viewer/transport.jl"))
include(joinpath(ROOT,"src/event_viewer/pads.jl"))

source=realpath(joinpath(ROOT,CONFIG["samples"]["bb0nu"]))
hits,_=h5open(f->read_event(f,Int64(0)),source,"r")
transport=transport_event(hits,source,"fixed",50f0,0)
pads=collect_pads(transport)
@assert sum(pads.history.counts.ions)==transport.d.retained_ions
@assert pads.projections.total_ions==transport.d.retained_ions
@assert pads.history.sampled_ions==transport.d.sampled_ions
@assert pads.history.outside_ions==transport.d.outside_ions
@assert length(pads.history.time_edges_ms)>1
println("Original viewer direct pad collection: ",pads.history.total_ions,
        " ions in ",nrow(pads.history.counts)," occupied pad/time cells")
