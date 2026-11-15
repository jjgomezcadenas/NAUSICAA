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

"""Take the event where it was generated; return its drift lengths in cm.

The source coordinates are detector coordinates: Z = 0 at the ion readout
plane and Z = L at the EL entry. Ion drift is the vertex Z, electron drift
L minus it. The vertex must lie within 0..L.
"""
function placement_from_vertex(info::MCEventInfo,L::Float32)::EventPlacement
    isfinite(L) && L>0 || throw(ArgumentError("Invalid chamber length"))
    z=info.generation_vertex_mm[3]
    0f0<=z<=L || throw(ArgumentError("Vertex of event $(info.event_id) outside 0..$L mm: $z"))
    EventPlacement(info.event_id,z/10f0,(L-z)/10f0)
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

"""
    main_track_box(hits::MCEventHits, info::MCEventInfo,
                   particles::Vector{MCParticle})::EventFrame

Find XYZ bounds of the primary electrons and their electron descendants.
Follow particle mothers back to the primary IDs. An electron whose chain passes
through any gamma belongs to that gamma deposit and does not set the main-track
bounds. Validate the particle table with gamma_map and preserve hit order/values.
"""
function main_track_box(h::MCEventHits,info::MCEventInfo,particles::Vector{MCParticle})::EventFrame
    h.event_id==info.event_id || throw(ArgumentError("Mismatched event numbers"))
    roots=gamma_map(particles,info)
    all(id->haskey(roots,id),h.particle_id) || throw(ArgumentError("Hit particle missing from MC table"))
    ids=Set(p.particle_id for p in particles if p.particle_name=="e-" && roots[p.particle_id]==0)
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
    gs=gamma_deposits(h,gamma_map(particles,info));frame=make_frame(main_track_box(h,info,particles),s.margin_mm)
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

"""Cube of the given side centred on the centre of a box, in the box's frame."""
function fixed_frame(b::EventFrame,side::Float32)::EventFrame
    isfinite(side)&&side>0 || throw(ArgumentError("Invalid frame side"))
    c=ntuple(i->(b.lower_mm[i]+b.upper_mm[i])/2f0,3)
    EventFrame(b.frame,ntuple(i->c[i]-side/2f0,3),ntuple(i->c[i]+side/2f0,3))
end

"""Largest side of a box in mm."""
box_side(b::EventFrame)::Float32=maximum(b.upper_mm[i]-b.lower_mm[i] for i in 1:3)

"""True when a vertex-relative box lies inside the detector cylinder."""
function _inside_detector(b::EventFrame,v::NTuple{3,Float32},L::Float32,R::Float32)::Bool
    0f0<=b.lower_mm[3]+v[3] && b.upper_mm[3]+v[3]<=L || return false
    all(hypot(x+v[1],y+v[2])<=R for x in (b.lower_mm[1],b.upper_mm[1]), y in (b.lower_mm[2],b.upper_mm[2]))
end

"""Frame the event in a fixed cube and apply detector, frame and gamma cuts.

Input hits are vertex-centred; detector coordinates add the generation vertex.
The placement must be the one given by `placement_from_vertex`. Cuts, in order:

1. `outside_detector`: a main-track hit, or a hit of a gamma deposit above the
   energy cut, lies outside the detector cylinder.
2. `track_larger_than_frame`: a side of the main-track box exceeds the frame side.
3. `frame_outside_detector`: the cube centred on the main-track box is not
   fully inside the detector.
4. `outside_frame_gamma`: a gamma deposit above the energy cut is not fully
   inside the cube.

Return the first failed decision and the cube. Inputs are preserved.
"""
function mc_fast_trigger(h::MCEventHits,info::MCEventInfo,particles::Vector{MCParticle},p::EventPlacement,s::FixedFrameTriggerSettings)::FastTriggerOutput
    h.event_id==info.event_id==p.event_id || throw(ArgumentError("Mismatched event numbers"))
    h.frame==mc_vertex_frame || throw(ArgumentError("Expected vertex-centred hits"))
    v=info.generation_vertex_mm;L=s.chamber_length_mm;R=s.detector_radius_mm
    isapprox(10f0*p.ion_drift_cm,v[3];atol=1f-2) && isapprox(10f0*(p.ion_drift_cm+p.electron_drift_cm),L;rtol=4eps(Float32)) ||
        throw(ArgumentError("Placement does not match the generation vertex"))
    roots=gamma_map(particles,info);gs=gamma_deposits(h,roots)
    box=main_track_box(h,info,particles);frame=fixed_frame(box,s.frame_side_mm)
    energetic=Set(g.gamma_particle_id for g in gs if g.energy_keV>s.gamma_energy_cut_keV)
    outside(j)=!(0f0<=h.z_mm[j]+v[3]<=L && hypot(h.x_mm[j]+v[1],h.y_mm[j]+v[2])<=R)
    escaped=any(j->(roots[h.particle_id[j]]==0 || roots[h.particle_id[j]] in energetic) && outside(j),eachindex(h.x_mm))
    decision=if escaped
        outside_detector
    elseif box_side(box)>s.frame_side_mm
        track_larger_than_frame
    elseif !_inside_detector(frame,v,L,R)
        frame_outside_detector
    elseif gamma_veto(gs,frame,s.gamma_energy_cut_keV)
        outside_frame_gamma
    else
        trigger_accepted
    end
    FastTriggerOutput(decision,frame)
end

"""Box to simulate inside a fixed frame: main track and the gamma deposits
fully inside the frame, enlarged by the margin on each side and limited to the
frame. With no margin given, the whole frame."""
function simulation_box(h::MCEventHits,info::MCEventInfo,particles::Vector{MCParticle},frame::EventFrame,margin_mm::Union{Nothing,Float32})::EventFrame
    margin_mm===nothing && return frame
    isfinite(margin_mm)&&margin_mm>=0 || throw(ArgumentError("Invalid margin"))
    boxes=[main_track_box(h,info,particles)]
    for g in gamma_deposits(h,gamma_map(particles,info));contains(frame,g.bounds)&&push!(boxes,g.bounds);end
    lo=ntuple(i->max(frame.lower_mm[i],minimum(b.lower_mm[i] for b in boxes)-margin_mm),3)
    hi=ntuple(i->min(frame.upper_mm[i],maximum(b.upper_mm[i] for b in boxes)+margin_mm),3)
    EventFrame(frame.frame,lo,hi)
end

"""Add one event to its selection count; return the updated statistics."""
function update_statistics(s::FixedFrameStatistics,o::FastTriggerOutput)::FixedFrameStatistics
    d=o.decision
    d in (trigger_accepted,outside_detector,track_larger_than_frame,frame_outside_detector,outside_frame_gamma) || throw(ArgumentError("Decision $d does not belong to the fixed-frame selection"))
    FixedFrameStatistics(s.input_events+1,s.outside_detector_events+(d==outside_detector),s.track_larger_than_frame_events+(d==track_larger_than_frame),
        s.frame_outside_detector_events+(d==frame_outside_detector),s.outside_frame_gamma_events+(d==outside_frame_gamma),s.accepted_events+(d==trigger_accepted))
end
merge_statistics(a::FixedFrameStatistics,b::FixedFrameStatistics)=FixedFrameStatistics((getfield(a,k)+getfield(b,k) for k in fieldnames(FixedFrameStatistics))...)

"""Add one event to its selection count; return the updated statistics."""
function update_statistics(s::FastTriggerStatistics,o::FastTriggerOutput)::FastTriggerStatistics
    d=o.decision
    FastTriggerStatistics(s.input_events+1,s.outside_detector_events+(d==outside_detector),s.outside_vertex_events+(d==outside_vertex_fiducial),s.outside_frame_gamma_events+(d==outside_frame_gamma),s.accepted_events+(d==trigger_accepted))
end
"""Add the event counts from two batches of the same sample."""
merge_statistics(a::FastTriggerStatistics,b::FastTriggerStatistics)=FastTriggerStatistics((getfield(a,k)+getfield(b,k) for k in fieldnames(FastTriggerStatistics))...)
