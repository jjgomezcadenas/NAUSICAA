"""Original MC coordinates and coordinates relative to the generation vertex."""
const mc_source_frame = source_frame
const mc_vertex_frame = vertex_frame

"""Event number, source generation vertex in mm and actual primary particle IDs.

Owns a sorted copy of the primary IDs. Coordinates are Float32. The event number
is unique within its sample; single electrons have one primary, bb0nu has two.
"""
struct MCEventInfo
    event_id::Int64
    generation_vertex_mm::NTuple{3,Float32}
    primary_particle_ids::Vector{Int32}
    function MCEventInfo(id::Int64,v::NTuple{3,Float32},ids::Vector{Int32})
        _validate_event_info(v,ids)
        new(id,v,sort(ids))
    end
end

"""An MC particle identified by its ID, mother ID and original particle name."""
struct MCParticle
    particle_id::Int32
    mother_id::Int32
    particle_name::String
end

"""One event's MC hits: Float32 positions in mm and deposited energies in keV.

Each vector entry describes one hit. IDs retain input order and precision.
Construction copies the vectors, checks equal lengths, finite coordinates,
nonnegative energies and unique (particle_id, hit_id) pairs. G4 time remains
in the source file. The reader converts original MeV energies to keV once.
"""
struct MCEventHits
    event_id::Int64
    frame::CoordinateFrame
    x_mm::Vector{Float32}
    y_mm::Vector{Float32}
    z_mm::Vector{Float32}
    energy_keV::Vector{Float32}
    particle_id::Vector{Int32}
    hit_id::Vector{Int32}
    function MCEventHits(id::Int64,f::CoordinateFrame,x::Vector{Float32},y::Vector{Float32},z::Vector{Float32},e::Vector{Float32},p::Vector{Int32},h::Vector{Int32})
        n=length(x)
        n>0 || throw(ArgumentError("Empty MC event"))
        all(v->length(v)==n,(y,z,e,p,h)) || throw(DimensionMismatch("MC hit column lengths differ"))
        all(v->all(isfinite,v),(x,y,z,e)) && all(>=(0f0),e) || throw(ArgumentError("Invalid MC hit values"))
        allunique(zip(p,h)) || throw(ArgumentError("Duplicate MC hit IDs"))
        new(id,f,copy(x),copy(y),copy(z),copy(e),copy(p),copy(h))
    end
end

"""Inclusive XYZ bounds, in mm, in the stated coordinate frame."""
struct EventFrame
    frame::CoordinateFrame
    lower_mm::NTuple{3,Float32}
    upper_mm::NTuple{3,Float32}
    function EventFrame(f::CoordinateFrame,lo::NTuple{3,Float32},hi::NTuple{3,Float32})
        all(isfinite,lo) && all(isfinite,hi) && all(lo[i]<=hi[i] for i in 1:3) || throw(ArgumentError("Invalid frame bounds"))
        new(f,lo,hi)
    end
end

"""Total deposited energy in keV and hit bounds for a gamma and its descendants."""
struct GammaDeposit
    gamma_particle_id::Int32
    energy_keV::Float32
    bounds::EventFrame
end

"""Detector length, vertex interval and frame margins in mm; gamma veto in keV.

Margins are applied to both sides of X, Y and Z. The detector begins at Z=0.
Settings are supplied by the caller; TOML parsing belongs in scripts.
"""
struct FastTriggerSettings
    chamber_length_mm::Float32
    vertex_z_min_mm::Float32
    vertex_z_max_mm::Float32
    margin_mm::NTuple{3,Float32}
    gamma_energy_cut_keV::Float32
    function FastTriggerSettings(L::Float32,lo::Float32,hi::Float32,m::NTuple{3,Float32},e::Float32)
        all(isfinite,(L,lo,hi,m...,e)) && L>0 && 0<=lo<=hi<=L && all(>=(0f0),m) && e>=0 || throw(ArgumentError("Invalid fast-trigger settings"))
        new(L,lo,hi,m,e)
    end
end

"""First failed cut, or acceptance. Numeric values are stable in saved files."""
@enum FastTriggerDecision::UInt8 trigger_accepted=0 outside_detector=1 outside_vertex_fiducial=2 outside_frame_gamma=3

"""The selection decision and the event's expanded acquisition frame."""
struct FastTriggerOutput
    decision::FastTriggerDecision
    frame::EventFrame
end

"""Disjoint counts for one sample. Normalization is input minus outside detector."""
struct FastTriggerStatistics
    input_events::Int64
    outside_detector_events::Int64
    outside_vertex_events::Int64
    outside_frame_gamma_events::Int64
    accepted_events::Int64
end
FastTriggerStatistics()=FastTriggerStatistics(0,0,0,0,0)
