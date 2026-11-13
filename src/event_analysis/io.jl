"""
    read_event(path::AbstractString, event_id::Int64) -> Tuple{EventHits,EventInfo}
    read_event(file::HDF5.File, event_id::Int64) -> Tuple{EventHits,EventInfo}

Read an indexed event from `nausicaa.hitdf.v3` into the implemented event types.
Hits are in `source_frame`, with Float32 mm/MeV/ns arrays; IDs retain their Int64/
Int32 types. EventInfo carries the common primary vertex and actual primary IDs,
read through the flat primary index. The result owns its arrays.

The path method opens and closes the HDF5 file. The open-file method preserves
the caller's handle. Missing IDs raise `KeyError`; invalid units/types/schema or
index entries raise an error. Older reduced schemas require regeneration; the
error supplies the existing converter command. EventHits/EventInfo validation
checks current values before returning. Sample identity belongs to the enclosing
dataset and is kept separate from the event ID.
"""
function read_event(f::HDF5.File, event_id::Int64)::Tuple{EventHits,EventInfo}
    schema = read(attributes(f)["schema"])
    schema == EventIO.SCHEMA || error("read_event requires $(EventIO.SCHEMA), got $schema. Regenerate: julia --project=. scripts/events/extract_hits.jl --overwrite SOURCE.next.h5")
    df = EventIO.read_hitdf(f; event_id)
    metadata = EventIO.read_event_metadata(f)
    i = findfirst(==(event_id), metadata.event_id)
    i === nothing && throw(KeyError(event_id))
    primaries = EventIO.read_primary_ids(f; event_id)
    info = EventInfo(event_id, (metadata.vertex_x_mm[i], metadata.vertex_y_mm[i],
                               metadata.vertex_z_mm[i]), primaries)
    hits = EventHits(event_id, source_frame, df.x_mm, df.y_mm, df.z_mm,
                    df.energy_MeV, df.time_ns, df.particle_id, df.hit_id)
    return (hits, info)
end
read_event(path::AbstractString, event_id::Int64) = h5open(f -> read_event(f,event_id), path, "r")

"""
    read_mc_event(file, event_id; event_row=nothing) -> (MCEventHits, MCEventInfo)

Read one indexed reduced HDF5 event. Positions remain Float32 mm. Convert
original Float32 MeV to Float32 keV once; read actual primary IDs and vertex.
`event_row` optionally supplies the known 1-based event-index row for batch
reading; it is checked against event_id. G4 time remains in the input file.
The file method leaves the handle open; the path method closes its own handle.
"""
function read_mc_event(f::HDF5.File,id::Int64;event_row::Union{Nothing,Int}=nothing)
    read(attributes(f)["schema"])==EventIO.SCHEMA || throw(ArgumentError("Expected reduced hitdf v3"))
    i=event_row===nothing ? findfirst(==(id),read(f["events/event_id"])) : event_row
    i===nothing && throw(KeyError(id))
    f["events/event_id"][i]==id || throw(ArgumentError("Event index does not match event number"))
    start=f["events/first_row"][i];n=f["events/n_hits"][i];ix=start:start+n-1
    ps=f["events/first_primary_row"][i];np=f["events/n_primaries"][i]
    ids=Int32.(f["events/primary_particle_ids"][ps:ps+np-1])
    v=Tuple(Float32(f["events/vertex_$(a)_mm"][i]) for a in ("x","y","z"))
    info=MCEventInfo(id,v,ids)
    h=MCEventHits(id,mc_source_frame,f["hits/x_mm"][ix],f["hits/y_mm"][ix],f["hits/z_mm"][ix],f["hits/energy_MeV"][ix].*1000f0,f["hits/particle_id"][ix],f["hits/hit_id"][ix])
    (h,info)
end
read_mc_event(path::AbstractString,id::Int64;kwargs...)=h5open(f->read_mc_event(f,id;kwargs...),path,"r")

"""
    read_mc_particles(file, event_ids; block_rows=8192)

Read MCParticle vectors for the requested event numbers from raw NEXUS
MC/particles. Scan the table in bounded blocks and group by event number;
particle rows may be unsorted. Missing requested events raise KeyError.
"""
function read_mc_particles(f::HDF5.File,ids::AbstractVector{Int64};block_rows::Int=8192)
    block_rows>0 || throw(ArgumentError("block_rows must be positive"))
    result=Dict(id=>MCParticle[] for id in ids);ds=f["MC/particles"]
    for first in 1:block_rows:length(ds)
        for p in ds[first:min(first+block_rows-1,length(ds))]
            haskey(result,p.event_id) || continue
            push!(result[p.event_id],MCParticle(p.particle_id,p.mother_id,String(p.particle_name)))
        end
    end
    for id in ids;isempty(result[id]) && throw(KeyError(id));end
    result
end
read_mc_particles(path::AbstractString,ids::AbstractVector{Int64};kwargs...)=h5open(f->read_mc_particles(f,ids;kwargs...),path,"r")
