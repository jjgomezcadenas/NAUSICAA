using Distributions
export primary_code_hash, legacy_realisation_identities, check_primary_code

"""Canonical representation of parsed TOML values, independent of key order/comments."""
function canonical(value)
    if value isa AbstractDict
        return "{"*join((repr(k)*":"*canonical(value[k]) for k in sort!(collect(keys(value)))),",")*"}"
    elseif value isa AbstractVector
        return "["*join(canonical.(value),",")*"]"
    elseif value isa AbstractString || value isa Number || value isa Bool
        return repr(value)
    end
    throw(ArgumentError("Unsupported TOML value $(typeof(value))"))
end

"""SHA-256 hex digest of a string or file, used for source and configuration IDs."""
content_hash(value::AbstractString;file::Bool=false)=
    bytes2hex(sha256(file ? read(value) : codeunits(value)))

"""Derive one UInt64 RNG seed from an identity, event ID and stage label."""
function event_seed(identity::AbstractString,event_id::Int64,label::AbstractString)::UInt64
    digest=content_hash(join((identity,string(event_id),label),"|"))
    return parse(UInt64,digest[1:16];base=16)
end

"""Hash ordered Int64 pair counts for branch-by-branch identity checks."""
pair_checksum(counts::Vector{Int64})=
    bytes2hex(sha256(Vector{UInt8}(reinterpret(UInt8,counts))))

"""Resolve primary RNG and common-realization identities independently of readout.

A `primary_reference` card anchors the primary-pair stream to an existing
realization. Its seed, W and Fano factor must match the new card; placement,
chamber and event-limit changes then alter only the common-realization identity.
Cards without a reference preserve the original stream hashing exactly.
"""
function legacy_realisation_identities(root,sample,source_sha,r)
    reference=haskey(r,"primary_reference") ? TOML.parsefile(joinpath(root,r["primary_reference"])) : r
    if reference!==r
        haskey(reference,"primary_reference") && error("Primary references must point to an original card")
        all(reference[k]==r[k] for k in ("seed","w_eV","fano_factor")) || error("Primary reference physics differs")
        reference["samples"][sample]==r["samples"][sample] || error("Primary reference sample differs")
    end
    stream_values=Dict(k=>v for (k,v) in reference if !(k in ("samples","max_events_per_sample")))
    stream_identity=content_hash(join((sample,source_sha,canonical(stream_values)),"|"))
    values=Dict(k=>v for (k,v) in r if k!="samples")
    realisation_id=content_hash(join((stream_identity,canonical(values)),"|"))
    (;stream_identity,realisation_id)
end

"""Hash the code and runtime dependencies that determine shared primary ledgers.

Pair sampling, event placement, seed derivation and ordered hit reading determine
this identity. Ion/EL collection code and detector cards have separate identities.
The runtime versions pin the RNG and distribution implementations.
"""
function primary_code_hash(root)
    files=("src/event_analysis/primary_sampling.jl","src/event_analysis/types.jl",
        "src/event_analysis/io.jl","src/ion_signal/event_io.jl",
        "scripts/production/primary_identity.jl")
    runtime=string(VERSION,"|",Base.pkgversion(Distributions))
    content_hash(join((content_hash(joinpath(root,p);file=true) for p in files),"|")*runtime)
end

"""Identify a v2 shared ledger while retaining the original random-stream seeds.

The source contents, card and sampling implementation identify the ledger. Code
changes create a new directory; existing ledger files remain immutable.
"""
function realisation_identities(root,sample,source_sha,r)
    legacy=legacy_realisation_identities(root,sample,source_sha,r)
    code=primary_code_hash(root)
    (;stream_identity=legacy.stream_identity,
      realisation_id=content_hash(join((legacy.realisation_id,code,"primary.v2"),"|")),
      primary_code_hash=code,legacy_realisation_id=legacy.realisation_id)
end

"""Require the same sampler identity in either collection arm before replay."""
function check_primary_code(root,metadata)
    haskey(metadata,"primary_code_hash") ||
        error("Legacy ledger: run run_common_realisation.jl to create a verified v2 ledger")
    String(read(metadata["primary_code_hash"]))==primary_code_hash(root) ||
        error("Shared primary sampling implementation differs from the saved ledger")
    nothing
end
