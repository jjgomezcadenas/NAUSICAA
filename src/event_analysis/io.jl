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
    schema == EventIO.SCHEMA || error("read_event requires $(EventIO.SCHEMA), got $schema. Regenerate: julia --project=. scripts/extract_hits.jl --overwrite SOURCE.next.h5")
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
