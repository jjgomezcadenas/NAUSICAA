"""Subtract the generation vertex from every MC hit. Preserve energies and IDs."""
function center_hits(h::MCEventHits,info::MCEventInfo)::MCEventHits
    h.event_id==info.event_id && h.frame==mc_source_frame || throw(ArgumentError("Expected matching source-frame event"))
    v=info.generation_vertex_mm
    MCEventHits(h.event_id,mc_vertex_frame,h.x_mm.-v[1],h.y_mm.-v[2],h.z_mm.-v[3],h.energy_keV,h.particle_id,h.hit_id)
end

"""Sample a uniform vertex along the chamber; return complementary drift lengths in cm.

The ion end is Z=0. Detector Z equals vertex-relative Z plus ion drift in mm.
The caller owns the RNG, records its seed, and saves the returned placement.
"""
function place_event(rng::AbstractRNG,id::Int64,L::Float32)::EventPlacement
    isfinite(L) && L>0 || throw(ArgumentError("Invalid chamber length"))
    z=Float32(rand(rng)*L)
    EventPlacement(id,z/10f0,(L-z)/10f0)
end

"""Map particles to the first gamma in their chain of mothers; zero if none.

Primary IDs come from MCEventInfo. Reject duplicate IDs, missing mothers,
cycles or non-electron primaries. Input order does not affect the map.
"""
function gamma_map(particles::Vector{MCParticle},info::MCEventInfo)::Dict{Int32,Int32}
    byid=Dict(p.particle_id=>p for p in particles);prim=Set(info.primary_particle_ids)
    length(byid)==length(particles) || throw(ArgumentError("Duplicate particle IDs"))
    all(id->haskey(byid,id)&&byid[id].particle_name=="e-",prim) || throw(ArgumentError("Missing or non-electron primary"))
    roots=Dict{Int32,Int32}();visiting=Set{Int32}()
    function walk(id)
        haskey(roots,id) && return roots[id]
        id in visiting && throw(ArgumentError("Particle mother cycle"))
        haskey(byid,id) || throw(ArgumentError("Missing mother particle $id"))
        push!(visiting,id);p=byid[id]
        parent=id in prim ? Int32(0) : walk(p.mother_id)
        root=parent!=0 ? parent : (p.particle_name=="gamma" ? id : Int32(0))
        roots[id]=root;delete!(visiting,id);root
    end
    for id in keys(byid);walk(id);end
    roots
end

function hit_bounds(h::MCEventHits,ix)::EventFrame
    isempty(ix) && throw(ArgumentError("No hits for bounds"))
    EventFrame(h.frame,ntuple(i->minimum((h.x_mm,h.y_mm,h.z_mm)[i][ix]),3),ntuple(i->maximum((h.x_mm,h.y_mm,h.z_mm)[i][ix]),3))
end

"""Find the XYZ extrema of direct primary-electron hits using their actual IDs."""
function primary_box(h::MCEventHits,info::MCEventInfo)::EventFrame
    h.event_id==info.event_id || throw(ArgumentError("Mismatched event numbers"))
    ids=Set(info.primary_particle_ids)
    hit_bounds(h,findall(id->id in ids,h.particle_id))
end

"""For each gamma, sum deposited energy in Float32 keV and calculate hit bounds.

Descendant gammas remain with their first gamma. Return positive-energy
GammaDeposit records sorted by gamma ID; each hit is counted once.
"""
function gamma_deposits(h::MCEventHits,roots::Dict{Int32,Int32})::Vector{GammaDeposit}
    all(id->haskey(roots,id),h.particle_id) || throw(ArgumentError("Hit particle missing from MC table"))
    groups=Dict{Int32,Vector{Int}}()
    for (j,id) in enumerate(h.particle_id)
        root=roots[id];root==0 && continue
        push!(get!(groups,root,Int[]),j)
    end
    result=GammaDeposit[]
    for id in sort(collect(keys(groups)))
        ix=groups[id];e=sum(h.energy_keV[ix]);e>0 && push!(result,GammaDeposit(id,e,hit_bounds(h,ix)))
    end
    result
end

"""Construct the event frame by adding the margins to each side of the box."""
function make_frame(b::EventFrame,m::NTuple{3,Float32})::EventFrame
    all(x->isfinite(x)&&x>=0,m) || throw(ArgumentError("Invalid margins"))
    EventFrame(b.frame,ntuple(i->b.lower_mm[i]-m[i],3),ntuple(i->b.upper_mm[i]+m[i],3))
end

"""Return true when a point or box lies inside the frame, including its edges."""
Base.contains(b::EventFrame,p::NTuple{3,Float32})=all(b.lower_mm[i]<=p[i]<=b.upper_mm[i] for i in 1:3)
function Base.contains(a::EventFrame,b::EventFrame)
    a.frame==b.frame || throw(ArgumentError("Coordinate frames differ"))
    contains(a,b.lower_mm)&&contains(a,b.upper_mm)
end

"""Return true if a gamma outside the frame deposits more energy than the cut."""
gamma_veto(gs::Vector{GammaDeposit},frame::EventFrame,cut::Float32)::Bool=any(g->g.energy_keV>cut&&!contains(frame,g.bounds),gs)

"""Define the event frame and apply detector, vertex and gamma cuts.

Input hits must be vertex-centred and match the info and placement event numbers.
Detector coordinates add 10*ion_drift_cm to hit Z. All hits must lie within
0..chamber_length_mm. The vertex must then lie in the fiducial interval.
The final cut rejects any gamma above the threshold outside the expanded frame.
Return the first failed decision and frame. Inputs are preserved.
"""
function mc_fast_trigger(h::MCEventHits,info::MCEventInfo,particles::Vector{MCParticle},p::EventPlacement,s::FastTriggerSettings)::FastTriggerOutput
    h.event_id==info.event_id==p.event_id || throw(ArgumentError("Mismatched event numbers"))
    h.frame==mc_vertex_frame || throw(ArgumentError("Expected vertex-centred hits"))
    isapprox(10f0*(p.ion_drift_cm+p.electron_drift_cm),s.chamber_length_mm;rtol=4eps(Float32)) || throw(ArgumentError("Drift lengths do not sum to chamber length"))
    gs=gamma_deposits(h,gamma_map(particles,info));frame=make_frame(primary_box(h,info),s.margin_mm)
    z=10f0*p.ion_drift_cm
    decision=if !all(v->0<=v+z<=s.chamber_length_mm,h.z_mm)
        outside_detector
    elseif !(s.vertex_z_min_mm<=z<=s.vertex_z_max_mm)
        outside_vertex_fiducial
    elseif gamma_veto(gs,frame,s.gamma_energy_cut_keV)
        outside_frame_gamma
    else
        trigger_accepted
    end
    FastTriggerOutput(decision,frame)
end

"""Add one event to its selection count; return the updated statistics."""
function update_statistics(s::FastTriggerStatistics,o::FastTriggerOutput)::FastTriggerStatistics
    d=o.decision
    FastTriggerStatistics(s.input_events+1,s.outside_detector_events+(d==outside_detector),s.outside_vertex_events+(d==outside_vertex_fiducial),s.outside_frame_gamma_events+(d==outside_frame_gamma),s.accepted_events+(d==trigger_accepted))
end
"""Add the event counts from two batches of the same sample."""
merge_statistics(a::FastTriggerStatistics,b::FastTriggerStatistics)=FastTriggerStatistics((getfield(a,k)+getfield(b,k) for k in fieldnames(FastTriggerStatistics))...)
