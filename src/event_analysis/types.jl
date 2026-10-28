"""
    CoordinateFrame

Coordinate system used by the `x_mm`, `y_mm` and `z_mm` vectors in [`EventHits`](@ref).
The `@enum` declaration defines this type and its two values:

- `source_frame`: positions in the original simulation coordinate system.
- `vertex_frame`: positions relative to the event's true generation vertex,
  obtained by subtracting that vertex from each source position. The vertex is
  at `(0, 0, 0)` in this frame, with the source axis directions preserved.

The field `hits.frame` records the coordinate system. [`center_hits`](@ref)
maps a source-frame event into the vertex frame. Distances are expressed in mm
in both frames.

# Example
```julia
frame = source_frame
@assert frame isa CoordinateFrame
@assert frame != vertex_frame
```
"""
@enum CoordinateFrame source_frame vertex_frame

"""Validate a finite vertex and a nonempty set of distinct primary IDs."""
function _validate_event_info(vertex_mm::NTuple{3,Float32}, primary_particle_ids::Vector{Int32})
    all(isfinite, vertex_mm) || throw(ArgumentError("vertex_mm must contain finite coordinates"))
    isempty(primary_particle_ids) && throw(ArgumentError("primary_particle_ids must be nonempty"))
    allunique(primary_particle_ids) || throw(ArgumentError("primary_particle_ids must be unique"))
    return nothing
end

"""
    EventInfo(event_id::Int64, vertex_mm::NTuple{3,Float32},
              primary_particle_ids::Vector{Int32})

Event-level information: source event ID, common primary generation vertex and
primary particle IDs. Electron-background events have one primary electron;
double-beta events have two.

# Fields and arguments
- `event_id::Int64`: event identifier within its source sample. Zero and gaps
  are valid; sample identity is supplied by the enclosing dataset.
- `vertex_mm::NTuple{3,Float32}`: finite source-frame coordinates `(x, y, z)` of
  the common generation vertex, in mm. The reader obtains the common vertex
  from the MC particle table.
- `primary_particle_ids::Vector{Int32}`: nonempty, distinct particle IDs read
  from the MC particle table. Construction stores a sorted copy. Its length
  gives the number of primaries.

# Ownership and validation
The constructor preserves its arguments and owns its primary-ID vector.
The tuple is immutable; array contents remain mutable. Computational functions
validate those contents when called.

Throw `ArgumentError` for a nonfinite vertex or an empty/duplicate primary-ID
vector. Arguments use the declared field types; conversion from file or table
representations belongs to the reader.

# Examples
```julia
electron = EventInfo(Int64(43), (10f0, 20f0, 100f0), Int32[7])
double_beta = EventInfo(Int64(42), (10f0, 20f0, 100f0), Int32[8, 3])
@assert electron.primary_particle_ids == Int32[7]
@assert double_beta.primary_particle_ids == Int32[3, 8]
```
"""
struct EventInfo
    event_id::Int64
    vertex_mm::NTuple{3,Float32}
    primary_particle_ids::Vector{Int32}

    function EventInfo(event_id::Int64, vertex_mm::NTuple{3,Float32},
                       primary_particle_ids::Vector{Int32})
        _validate_event_info(vertex_mm, primary_particle_ids)
        return new(event_id, vertex_mm, sort(primary_particle_ids))
    end
end

"""
Validate equal positive vector lengths, finite physical values, nonnegative
energies and unique `(particle_id, hit_id)` pairs within one event.
"""
function _validate_hit_vectors(x_mm::Vector{Float32}, y_mm::Vector{Float32},
        z_mm::Vector{Float32}, energy_MeV::Vector{Float32}, time_ns::Vector{Float32},
        particle_id::Vector{Int32}, hit_id::Vector{Int32})
    n = length(x_mm)
    all(v -> length(v) == n, (y_mm, z_mm, energy_MeV, time_ns, particle_id, hit_id)) ||
        throw(DimensionMismatch("all EventHits vectors must have equal lengths"))
    n > 0 || throw(ArgumentError("an event must contain at least one hit"))
    for (name, values) in ((:x_mm, x_mm), (:y_mm, y_mm), (:z_mm, z_mm),
                            (:energy_MeV, energy_MeV), (:time_ns, time_ns))
        all(isfinite, values) || throw(ArgumentError("$name must contain finite values"))
    end
    all(>=(0f0), energy_MeV) || throw(ArgumentError("energy_MeV must be nonnegative"))
    seen = Set{Tuple{Int32,Int32}}()
    for key in zip(particle_id, hit_id)
        key in seen && throw(ArgumentError("duplicate (particle_id, hit_id) pair: $key"))
        push!(seen, key)
    end
    return nothing
end

"""
    EventHits(event_id::Int64, frame::CoordinateFrame,
              x_mm::Vector{Float32}, y_mm::Vector{Float32}, z_mm::Vector{Float32},
              energy_MeV::Vector{Float32}, time_ns::Vector{Float32},
              particle_id::Vector{Int32}, hit_id::Vector{Int32})

One event's hits stored as parallel, independently owned vectors. Entry `j` in
each vector describes the same hit. The same concrete type represents source
and vertex-centred coordinates.

# Fields and arguments
- `event_id::Int64`: event identifier within the enclosing source sample.
- `frame::CoordinateFrame`: coordinate system of the three position vectors;
  see [`CoordinateFrame`](@ref).
- `x_mm`, `y_mm`, `z_mm`: Float32 coordinates in mm, each with length `N`.
- `energy_MeV`: length-`N` Float32 vector of deposited energies in MeV.
- `time_ns`: length-`N` Float32 vector of MC deposition times in ns.
- `particle_id`, `hit_id`: length-`N` Int32 identifier vectors. The pair is
  unique within the event; `(event_id, particle_id, hit_id)` identifies a hit
  within its sample. Hit IDs may repeat across particles.

# Ownership and validation
Construction preserves hit ordering and copies every vector separately, even
when several arguments refer to the same array. All vectors have the same
positive length. Physical values are finite; energy deposits are nonnegative.
The constructor accepts the declared field types.

Throw `DimensionMismatch` for unequal vector lengths. Throw `ArgumentError`
for an empty event, nonfinite physical values, negative energy or duplicate hit
identities. Array contents remain mutable and are validated at computational
boundaries.

# Example
```julia
hits = EventHits(Int64(42), source_frame,
    Float32[10, 13], Float32[20, 24], Float32[100, 102],
    Float32[0.001, 0.002], Float32[0, 0.5], Int32[1, 2], Int32[0, 0])
@assert hits.x_mm == Float32[10, 13]
@assert hits.hit_id == Int32[0, 0]
@assert length(hits.energy_MeV) == 2
```
"""
struct EventHits
    event_id::Int64
    frame::CoordinateFrame
    x_mm::Vector{Float32}
    y_mm::Vector{Float32}
    z_mm::Vector{Float32}
    energy_MeV::Vector{Float32}
    time_ns::Vector{Float32}
    particle_id::Vector{Int32}
    hit_id::Vector{Int32}

    function EventHits(event_id::Int64, frame::CoordinateFrame,
            x_mm::Vector{Float32}, y_mm::Vector{Float32}, z_mm::Vector{Float32},
            energy_MeV::Vector{Float32}, time_ns::Vector{Float32},
            particle_id::Vector{Int32}, hit_id::Vector{Int32})
        _validate_hit_vectors(x_mm, y_mm, z_mm, energy_MeV, time_ns, particle_id, hit_id)
        return new(event_id, frame, copy(x_mm), copy(y_mm), copy(z_mm),
                   copy(energy_MeV), copy(time_ns), copy(particle_id), copy(hit_id))
    end
end
