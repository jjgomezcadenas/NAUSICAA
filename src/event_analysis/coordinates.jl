"""
    center_hits(hits::EventHits, info::EventInfo) -> EventHits

Express each hit relative to the event's true generation vertex, preserving the
source axis directions and hit ordering.

# Arguments
- `hits::EventHits`: one event in `source_frame`. Its `x_mm`, `y_mm`, `z_mm`,
  `energy_MeV` and `time_ns` fields contain Float32 vectors in mm, MeV and ns;
  `particle_id` and `hit_id` contain Int32 vectors with the same length.
- `info::EventInfo`: matching event ID, finite Float32 generation vertex in
  source coordinates, and distinct primary particle IDs.

# Returns and arithmetic
A new `EventHits` in `vertex_frame`. Each axis is translated in Float32:
`x_mm = hits.x_mm .- info.vertex_mm[1]`, and likewise for y and z. Coordinates
retain Float32 precision, so comparisons under arbitrary translations use a
tolerance consistent with the source coordinate spacing.

Energy, MC time and identifiers retain their exact values and ordering. Hit
identity remains `(event_id, particle_id, hit_id)` within its source sample.

# Ownership and determinism
The result depends on the supplied arguments. Inputs retain their values and
every returned array is independently owned, including the preserved fields.

# Errors
Throw `ArgumentError` for mismatched event IDs, an input frame other than
`source_frame`, an empty event, invalid physical values, duplicate hit identities
or empty/duplicate primary IDs. Translated coordinates must remain finite.
Throw `DimensionMismatch` for unequal hit-vector lengths. Current array contents
are validated on entry because array fields are mutable.

# Example
```julia
info = EventInfo(Int64(42), (10f0, 20f0, 100f0), Int32[1, 2])
hits = EventHits(Int64(42), source_frame,
    Float32[10, 13], Float32[20, 24], Float32[100, 102],
    Float32[0.001, 0.002], Float32[0, 0.5], Int32[1, 2], Int32[0, 0])
local_hits = center_hits(hits, info)
@assert local_hits.x_mm == Float32[0, 3]
@assert local_hits.y_mm == Float32[0, 4]
@assert local_hits.z_mm == Float32[0, 2]
@assert local_hits.frame == vertex_frame
```
"""
function center_hits(hits::EventHits, info::EventInfo)::EventHits
    hits.event_id == info.event_id || throw(ArgumentError("hits and info must have the same event_id"))
    hits.frame == source_frame || throw(ArgumentError("center_hits requires source_frame"))
    _validate_event_info(info.vertex_mm, info.primary_particle_ids)
    _validate_hit_vectors(hits.x_mm, hits.y_mm, hits.z_mm, hits.energy_MeV,
                          hits.time_ns, hits.particle_id, hits.hit_id)
    return EventHits(hits.event_id, vertex_frame,
        hits.x_mm .- info.vertex_mm[1],
        hits.y_mm .- info.vertex_mm[2],
        hits.z_mm .- info.vertex_mm[3],
        hits.energy_MeV, hits.time_ns, hits.particle_id, hits.hit_id)
end
