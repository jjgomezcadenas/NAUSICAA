module EventIO

using DataFrames, Dates, HDF5, SHA, TOML

export read_nexus_hits, read_hitdf, read_event_metadata, read_primary_ids, convert_hits, verify_hitdf, default_output

const ROOT = normpath(joinpath(@__DIR__, "../.."))
const SCHEMA = "nausicaa.hitdf.v3"
const COLUMNS = (:event_id, :particle_id, :hit_id, :x_mm, :y_mm, :z_mm, :energy_MeV, :time_ns)
const SOURCE_COLUMNS = (:event_id, :particle_id, :hit_id, :x, :y, :z, :energy, :time)
const TYPES = (Int64, Int32, Int32, Float32, Float32, Float32, Float32, Float32)
const UNITS = ("1", "1", "1", "mm", "mm", "mm", "MeV", "ns")
const DEFAULT_BLOCK_ROWS = 65_536
const DEFAULT_CHUNK_ROWS = 4_096

file_sha256(path) = open(io -> bytes2hex(sha256(io)), path)
empty_hits() = DataFrame([name => T[] for (name, T) in zip(COLUMNS, TYPES)])

function hit_dataframe(rows)
    # Reject unexpected input types rather than silently rounding source data.
    columns = Pair{Symbol,AbstractVector}[]
    for (out, src, T) in zip(COLUMNS, SOURCE_COLUMNS, TYPES)
        hasfield(eltype(rows), src) || error("Missing NEXUS hit field: $src")
        fieldtype(eltype(rows), src) == T || error("Expected $src::$T in NEXUS input")
        push!(columns, out => getproperty.(rows, src))
    end
    df = DataFrame(columns; copycols=false)
    for name in (:x_mm, :y_mm, :z_mm, :energy_MeV, :time_ns)
        all(isfinite, df[!, name]) || error("Non-finite hit values in $name")
    end
    all(>=(0), df.energy_MeV) || error("Negative deposited energy")
    return df
end

"""
    read_nexus_hits(path; event_id=nothing, block_rows=65536) -> DataFrame

Read original /MC/hits in blocks, preserving row order, all particles, mm, MeV, ns
and Float32. With event_id, return only that event (KeyError if absent). The
original file has no event index, so this selection scans its hit table.
"""
function read_nexus_hits(path::AbstractString; event_id=nothing,
                         block_rows::Int=DEFAULT_BLOCK_ROWS)
    block_rows > 0 || throw(ArgumentError("block_rows must be positive"))
    h5open(path, "r") do f
        d = f["MC/hits"]
        df = empty_hits()
        for start in 1:block_rows:length(d)
            rows = d[start:min(start + block_rows - 1, length(d))]
            event_id === nothing || (rows = filter(h -> h.event_id == event_id, rows))
            isempty(rows) || append!(df, hit_dataframe(rows))
        end
        event_id !== nothing && isempty(df) && throw(KeyError(event_id))
        return df
    end
end

"""
    default_output(path) -> String

Map name.next.h5 (or name.h5) to name_hitdf.h5. Sources in the raw archive
`~/NAUSICAA` produce reduced tables in the project's events/ directory.
NAUSICAA_RAW_DIRECTORY overrides the archive location. Other inputs produce
their reduced table beside the source.
"""
function default_output(path::AbstractString)
    endswith(path, ".h5") || throw(ArgumentError("Input must end in .h5"))
    stem = replace(basename(path), r"(?:\.next)?\.h5$" => "")
    archive=abspath(expanduser(get(ENV,"NAUSICAA_RAW_DIRECTORY",joinpath(homedir(),"NAUSICAA"))))
    directory=dirname(abspath(path))==archive ? joinpath(ROOT,"events") : dirname(path)
    return joinpath(directory, stem * "_hitdf.h5")
end

function hit_columns(f)
    schema = read(attributes(f)["schema"])
    schema in (SCHEMA, "nausicaa.hitdf.v2", "nausicaa.hitdf.v1") || error("Unsupported hitdf schema: $schema")
    return schema == SCHEMA ? COLUMNS : COLUMNS[1:end-1]
end

function check_schema(f)
    columns = hit_columns(f)
    read(attributes(f["events"])["index_base"]) == 1 || error("Expected one-based row index")
    lengths = [length(f["hits/$name"]) for name in columns]
    allequal(lengths) || error("Hit columns have different lengths")
    for (name, T, unit) in zip(columns, TYPES, UNITS)
        d = f["hits/$name"]
        eltype(d) == T || error("Unexpected type of $name")
        read(attributes(d)["units"]) == unit || error("Unexpected units of $name")
    end
    return first(lengths)
end

"""
    read_hitdf(path_or_open_file; event_id=nothing) -> DataFrame

Read all reduced hits or only an indexed event. An open HDF5 file can be reused
for repeated calls. Event IDs are identifiers, not row numbers. Missing IDs
raise KeyError. v1/v2 files return the original seven columns; v3 adds time_ns.
Source identity is stored in /metadata/source_file.
"""
function read_hitdf(f::HDF5.File; event_id=nothing)
    n = check_schema(f)
    rows = 1:n
    if event_id !== nothing
        ids = read(f["events/event_id"])
        i = findfirst(==(event_id), ids)
        i === nothing && throw(KeyError(event_id))
        start = only(f["events/first_row"][i:i])
        count = only(f["events/n_hits"][i:i])
        1 <= start <= n && 0 < count <= n-start+1 || error("Invalid event index")
        rows = start:(start+count-1)
    end
    df = DataFrame([name => f["hits/$name"][rows] for name in hit_columns(f)]; copycols=false)
    event_id === nothing || all(==(event_id), df.event_id) || error("Incorrect event index")
    return df
end

read_hitdf(path::AbstractString; kwargs...) = h5open(f -> read_hitdf(f; kwargs...), path, "r")

"""Read one row per event, including the shared primary generation vertex (mm)."""
function read_event_metadata(f::HDF5.File)
    check_schema(f)
    names = (:event_id,:first_row,:n_hits,:energy_MeV,:vertex_x_mm,:vertex_y_mm,:vertex_z_mm,:n_primaries)
    all(name -> haskey(f,"events/$name"),names) || error("Missing generation vertices; regenerate with scripts/events/extract_hits.jl --overwrite")
    df = DataFrame([name => read(f["events/$name"]) for name in names]; copycols=false)
    allunique(df.event_id) || error("Repeated event IDs in metadata")
    all(>(0),df.n_primaries) || error("Missing primary particle")
    for name in (:vertex_x_mm,:vertex_y_mm,:vertex_z_mm)
        eltype(df[!,name]) == Float32 && all(isfinite,df[!,name]) || error("Invalid generation vertex")
        read(attributes(f["events/$name"])["units"]) == "mm" || error("Invalid vertex units")
    end
    return df
end
read_event_metadata(path::AbstractString) = h5open(read_event_metadata,path,"r")

function primary_vertices(source, ids)
    vertices = Dict{Int64,NTuple{3,Float32}}()
    counts = Dict{Int64,Int32}()
    primaries = Dict(id => Int32[] for id in ids)
    wanted = Set(ids)
    isempty(ids) && return (x=Float32[],y=Float32[],z=Float32[],counts=Int32[],primary_ids=Int32[],primary_starts=Int64[])
    haskey(source,"MC/particles") || error("Missing /MC/particles: cannot identify the generation vertex")
    for p in read(source["MC/particles"])
        p.event_id in wanted && p.primary == 1 || continue
        v = (p.initial_x,p.initial_y,p.initial_z)
        v isa NTuple{3,Float32} && all(isfinite,v) || error("Invalid primary vertex in event $(p.event_id)")
        haskey(vertices,p.event_id) && vertices[p.event_id] != v &&
            error("Primaries in event $(p.event_id) have different generation vertices")
        get!(vertices,p.event_id,v)
        p.particle_id isa Int32 || error("Expected Int32 primary particle ID")
        p.particle_id in primaries[p.event_id] && error("Duplicate primary ID in event $(p.event_id)")
        push!(primaries[p.event_id], p.particle_id)
        counts[p.event_id] = get(counts,p.event_id,Int32(0))+Int32(1)
    end
    all(id -> haskey(vertices,id),ids) || error("An event has no primary generation vertex")
    flat = Int32[]
    starts = Int64[]
    for id in ids
        push!(starts, length(flat)+1)
        append!(flat, sort(primaries[id]))
    end
    return (x=Float32[vertices[id][1] for id in ids],y=Float32[vertices[id][2] for id in ids],
            z=Float32[vertices[id][3] for id in ids],counts=Int32[counts[id] for id in ids],
            primary_ids=flat, primary_starts=starts)
end

"""
    read_primary_ids(file::HDF5.File; event_id) -> Vector{Int32}

Read the sorted, actual primary particle IDs for one event from a v3 reduced
file. The flat `/events/primary_particle_ids` column is indexed by one-based
`first_primary_row` and `n_primaries`. Return an independently owned vector.
One-electron and two-primary events use the same representation. Missing event
IDs raise `KeyError`; missing v3 fields or invalid ranges raise an error with
regeneration guidance. The supplied file remains open.
"""
function read_primary_ids(f::HDF5.File; event_id)
    read(attributes(f)["schema"]) == SCHEMA || error("Primary IDs require $SCHEMA; regenerate with julia --project=. scripts/events/extract_hits.jl --overwrite SOURCE.next.h5")
    ids = read(f["events/event_id"])
    i = findfirst(==(event_id), ids)
    i === nothing && throw(KeyError(event_id))
    start = only(f["events/first_primary_row"][i:i])
    count = only(f["events/n_primaries"][i:i])
    d = f["events/primary_particle_ids"]
    eltype(d) == Int32 || error("Primary IDs must be Int32")
    1 <= start <= length(d) && 0 < count <= length(d)-start+1 || error("Invalid primary index")
    result = d[start:start+count-1]
    allunique(result) && issorted(result) || error("Primary IDs must be sorted and distinct")
    return result
end

mutable struct EventIndex
    ids::Vector{Int64}
    starts::Vector{Int64}
    counts::Vector{Int64}
    energies::Vector{Float32}
    seen::Set{Int64}
    hit_keys::Set{Tuple{Int32,Int32}}
    labels::Set{String}
end
EventIndex() = EventIndex(Int64[], Int64[], Int64[], Float32[], Set{Int64}(),
                          Set{Tuple{Int32,Int32}}(), Set{String}())

function index_rows!(idx, rows, first_row)
    for (i, h) in enumerate(rows)
        push!(idx.labels, h.label)
        length(idx.labels) <= 1 || error("Multiple hit labels: cannot omit the label column losslessly")
        if isempty(idx.ids) || h.event_id != last(idx.ids)
            h.event_id in idx.seen && error("Event $(h.event_id) is not contiguous; input order is preserved")
            push!(idx.seen, h.event_id)
            push!(idx.ids, h.event_id); push!(idx.starts, first_row+i-1)
            push!(idx.counts, 0); push!(idx.energies, 0f0)
            empty!(idx.hit_keys)
        end
        key = (h.particle_id, h.hit_id)
        key in idx.hit_keys && error("Duplicate particle_id/hit_id in event $(h.event_id)")
        push!(idx.hit_keys, key)
        idx.counts[end] += 1
        idx.energies[end] += h.energy
    end
end

function write_metadata(out, source, input, idx; chunk_rows, compression)
    events = create_group(out, "events")
    for (name, vals) in zip(("event_id", "first_row", "n_hits", "energy_MeV"),
                            (idx.ids, idx.starts, idx.counts, idx.energies))
        events[name] = vals
    end
    attributes(events)["index_base"] = 1
    vertices = primary_vertices(source,idx.ids)
    for (name,values) in zip(("vertex_x_mm","vertex_y_mm","vertex_z_mm"),(vertices.x,vertices.y,vertices.z))
        events[name] = values
        attributes(events[name])["units"] = "mm"
        attributes(events[name])["source"] = "/MC/particles: shared initial position of all primary=1 rows in event"
    end
    events["n_primaries"] = vertices.counts
    events["first_primary_row"] = vertices.primary_starts
    events["primary_particle_ids"] = vertices.primary_ids
    attributes(events["primary_particle_ids"])["source"] = "/MC/particles: particle_id of primary=1 rows, sorted within event"
    attributes(events["energy_MeV"])["units"] = "MeV"
    attributes(events["energy_MeV"])["summation"] = "Float32 accumulation in source row order"
    meta = create_group(out, "metadata")
    meta["source_file"] = basename(input)
    meta["source_sha256"] = file_sha256(input)
    meta["source_dataset"] = "/MC/hits"
    meta["label"] = isempty(idx.labels) ? "" : only(idx.labels)
    meta["n_hits"] = sum(idx.counts)
    meta["n_events"] = length(idx.ids)
    meta["created_utc"] = string(now(UTC)) * "Z"
    meta["julia_version"] = string(VERSION)
    meta["hdf5_jl_version"] = string(pkgversion(HDF5))
    meta["dataframes_jl_version"] = string(pkgversion(DataFrames))
    meta["reader_sha256"] = file_sha256(@__FILE__)
    manifest = joinpath(ROOT, "Manifest.toml")
    isfile(manifest) && (meta["manifest_sha256"] = file_sha256(manifest))
    meta["chunk_rows"] = chunk_rows
    meta["deflate_level"] = compression
    meta["source_columns"] = string.(collect(SOURCE_COLUMNS))
    meta["output_columns"] = string.(collect(COLUMNS))
    meta["units"] = collect(UNITS)
    meta["unit_conversion_factors"] = ones(Float64, length(COLUMNS))
    meta["selection"] = "All hits, all particles; original row order; no energy or topology cuts"
    meta["omitted_fields"] = ["label (constant, retained in metadata)"]
    meta["time_assumption"] = "Original MC deposition time retained in ns for truth ordering"
    if haskey(source, "MC/configuration")
        cfg = read(source["MC/configuration"])
        g = create_group(meta, "configuration")
        g["param_key"] = String[h.param_key for h in cfg]
        g["param_value"] = String[h.param_value for h in cfg]
        attributes(g)["ordering"] = "Original row order; duplicate keys retained"
    end
end

"""
    verify_hitdf(input, output; block_rows=65536)

Independently compare each reduced value and type against the original HDF5,
then check event ranges, Float32 source-order energy sums, primary IDs, labels,
configuration and SHA-256. Only v3 outputs have the full verification contract.
No tolerance or renormalization is used for hit values or event energy sums.
"""
function verify_hitdf(input, output; block_rows::Int=DEFAULT_BLOCK_ROWS)
    block_rows > 0 || throw(ArgumentError("block_rows must be positive"))
    h5open(input, "r") do src
        h5open(output, "r") do dst
            read(attributes(dst)["schema"]) == SCHEMA || error("Verification requires $SCHEMA output")
            hits = src["MC/hits"]
            n = check_schema(dst)
            length(hits) == n || error("Hit count differs from source")
            totals = Dict{Int64,Float32}()
            counts = Dict{Int64,Int64}()
            first_rows = Dict{Int64,Int64}()
            last_rows = Dict{Int64,Int64}()
            label = read(dst["metadata/label"])
            for start in 1:block_rows:n
                range = start:min(start+block_rows-1, n)
                original = hits[range]
                for (out, field, T) in zip(COLUMNS, SOURCE_COLUMNS, TYPES)
                    actual = dst["hits/$out"][range]
                    expected = getproperty.(original, field)
                    eltype(actual) == eltype(expected) == T || error("Type mismatch: $out")
                    isequal(actual, expected) || error("Hit values differ: $out, rows $range")
                end
                for (offset, h) in enumerate(original)
                    h.label == label || error("Source label mismatch")
                    row = start+offset-1
                    get!(first_rows, h.event_id, row)
                    last_rows[h.event_id] = row
                    counts[h.event_id] = get(counts, h.event_id, 0) + 1
                    totals[h.event_id] = get(totals, h.event_id, 0f0) + h.energy
                end
            end
            ids = read(dst["events/event_id"])
            starts = read(dst["events/first_row"])
            sizes = read(dst["events/n_hits"])
            energies = read(dst["events/energy_MeV"])
            length(ids) == length(starts) == length(sizes) == length(energies) == length(totals) || error("Invalid event index length")
            length(unique(ids)) == length(ids) && Set(ids) == Set(keys(totals)) || error("Event IDs differ")
            for i in eachindex(ids)
                id = ids[i]
                starts[i] == first_rows[id] || error("Incorrect event start: $id")
                sizes[i] == counts[id] == last_rows[id]-first_rows[id]+1 || error("Incorrect event size: $id")
                isequal(energies[i], totals[id]) || error("Event energy differs: $id")
            end
            vertices = primary_vertices(src,ids)
            for (name,expected) in zip(("vertex_x_mm","vertex_y_mm","vertex_z_mm","n_primaries"),
                                       (vertices.x,vertices.y,vertices.z,vertices.counts))
                actual = read(dst["events/$name"])
                eltype(actual) == eltype(expected) && isequal(actual,expected) || error("Primary vertex metadata differs: $name")
            end
            isequal(read(dst["events/first_primary_row"]), vertices.primary_starts) || error("Primary index differs")
            isequal(read(dst["events/primary_particle_ids"]), vertices.primary_ids) || error("Primary IDs differ")
            read_event_metadata(dst)
            read(dst["metadata/n_hits"]) == n || error("Metadata hit count differs")
            read(dst["metadata/n_events"]) == length(ids) || error("Metadata event count differs")
            read(dst["metadata/source_sha256"]) == file_sha256(input) || error("Source SHA-256 differs")
            if haskey(src, "MC/configuration")
                cfg = read(src["MC/configuration"])
                read(dst["metadata/configuration/param_key"]) == [r.param_key for r in cfg] || error("Configuration keys differ")
                read(dst["metadata/configuration/param_value"]) == [r.param_value for r in cfg] || error("Configuration values differ")
            end
            return (n_hits=n, n_events=length(ids))
        end
    end
end

"""
    convert_hits(input; output=default_output(input), overwrite=false, ...)

Write a losslessly compressed columnar HDF5 file in blocks. Input events must
occupy contiguous rows and label must be constant. Keep the original untouched;
verify a temporary file under work/ before renaming it to the destination.
"""
function convert_hits(input::AbstractString; output::AbstractString=default_output(input),
                      overwrite::Bool=false, block_rows::Int=DEFAULT_BLOCK_ROWS,
                      chunk_rows::Int=DEFAULT_CHUNK_ROWS, compression::Int=3)
    block_rows > 0 && chunk_rows > 0 || throw(ArgumentError("Block and chunk sizes must be positive"))
    0 <= compression <= 9 || throw(ArgumentError("Compression must be 0:9"))
    isfile(input) || throw(ArgumentError("Input does not exist: $input"))
    # Reject aliases/symlinks of the source as well as the same spelling.
    target = ispath(output) ? realpath(output) : joinpath(realpath(dirname(abspath(output))), basename(output))
    realpath(input) != target || throw(ArgumentError("Output cannot replace the source"))
    ispath(output) && !overwrite && throw(ArgumentError("Output exists; use --overwrite: $output"))
    scratch = joinpath(ROOT, "work")
    mkpath(scratch)
    mktempdir(scratch; prefix="hitdf-") do tempdir
        temporary = joinpath(tempdir, "hits.h5")
        h5open(input, "r") do source
            hits = source["MC/hits"]
            h5open(temporary, "w") do out
                attributes(out)["schema"] = SCHEMA
                group = create_group(out, "hits")
                n = length(hits)
                for (name, T, unit) in zip(COLUMNS, TYPES, UNITS)
                    d = n == 0 ? create_dataset(group, string(name), T, (0,)) :
                        create_dataset(group, string(name), T, (n,);
                                       chunk=(min(n, chunk_rows),), shuffle=true, deflate=compression)
                    attributes(d)["units"] = unit
                    close(d)
                end
                idx = EventIndex()
                for start in 1:block_rows:n
                    range = start:min(start+block_rows-1, n)
                    rows = hits[range]
                    df = hit_dataframe(rows)
                    index_rows!(idx, rows, start)
                    for name in COLUMNS
                        group[string(name)][range] = df[!, name]
                    end
                end
                write_metadata(out, source, input, idx; chunk_rows, compression)
            end
        end
        checked = verify_hitdf(input, temporary; block_rows)
        h5open(temporary, "r+") do f
            f["metadata/validation"] = "All hit values/types, event sums, MC times, actual primary IDs and common primary generation vertices verified exactly against source"
        end
        ispath(output) && !overwrite && error("Output appeared during conversion: $output")
        Base.Filesystem.rename(temporary, output)
        return (; input=String(input), output=String(output), checked...,
                input_bytes=filesize(input), output_bytes=filesize(output))
    end
end

end # module
