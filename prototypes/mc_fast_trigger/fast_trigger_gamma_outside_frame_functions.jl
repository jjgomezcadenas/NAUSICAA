module MCFastTrigger
using Statistics, LinearAlgebra
include(joinpath(@__DIR__,"../../src/event_analysis/EventAnalysis.jl"))
const EA=EventAnalysis
export passes_gamma_cut, mc_fast_trigger, gamma_roots, summary, detector_selection

"""
    gamma_roots(particles) -> Dict{Int32,Int32}

Map each MC particle number to its first gamma ancestor, walking from a primary
outwards. Zero means the particle has no gamma ancestor. A gamma's electron
and later photon descendants belong to that first gamma, so deposited energy
is counted once. Direct electron-ionisation secondaries without a gamma
ancestor remain part of the gamma-free track. Particle names are the original
NEXUS names (`gamma`, `e-`). Missing parents and ancestry cycles raise errors.
"""
function gamma_roots(particles)
    prim=[p for p in particles if p.primary==1]
    info=EA.MCEventInfo(Int64(0),(0f0,0f0,0f0),Int32[p.particle_id for p in prim])
    EA.gamma_map(EA.MCParticle[EA.MCParticle(p.particle_id,p.mother_id,String(p.particle_name)) for p in particles],info)
end

"""Summary values for plot annotations. RMS is sqrt(mean(x^2)); std is population SD."""
function summary(x)
    isempty(x) && return (;n=Int64(0),mean=NaN,std=NaN,rms=NaN,median=NaN,p95=NaN,max=NaN)
    y=Float64.(x)
    (;n=Int64(length(y)),mean=mean(y),std=std(y;corrected=false),rms=sqrt(mean(abs2,y)),median=median(y),p95=quantile(y,.95),max=maximum(y))
end
"""
    detector_selection(z_relative_mm, vertex_z_mm; length_mm=1500, fiducial=(50,1450))

Translate all hit Z coordinates by the sampled primary-vertex Z position.
Return translated coordinates and sequential containment and vertex-fiducial
flags. Containment includes every deposited hit, including gamma descendants.
The fiducial requirement acts on the primary vertex. Boundaries are inclusive.
"""
function detector_selection(z_relative_mm,vertex_z_mm;length_mm=1500f0,fiducial=(50f0,1450f0))
    z=Float32.(z_relative_mm).+Float32(vertex_z_mm)
    contained=all(v->0<=v<=length_mm,z)
    (;z,contained,accepted=contained && fiducial[1]<=vertex_z_mm<=fiducial[2])
end
"""
    passes_gamma_cut(gammas, energy_cut_keV) -> Bool

Accept the event when every gamma deposit extending outside the primary frame
has total deposited energy at most `energy_cut_keV`. An outside deposit with
energy strictly greater than the threshold rejects the whole event. Deposits
entirely inside the frame impose no energy restriction.
"""
passes_gamma_cut(gammas,energy_cut_keV)=!any(g->g.inside_frame==0 && g.energy_keV>energy_cut_keV,gammas)
extent(x)=vec(maximum(x;dims=2)-minimum(x;dims=2))

"""
    mc_fast_trigger(hits, particles; margin_mm=20f0, expected_primaries=2)

Descriptive prototype for one event: `hits` supplies x_mm/y_mm/z_mm,
energy_MeV, particle_id; `particles` supplies NEXUS ancestry. Returns event
measurements, gamma-family measurements, and classified hits. No visibility
cut or acceptance decision is applied in this pilot.

The reference track consists of hits from the primary electron(s). A second
extent includes all hits without a gamma ancestor, retaining electron
ionisation secondaries. All deposited-energy sums use MeV->keV only.
Each gamma-family energy is the sum of its descendant hits, not initial
photon energy. For each hit compute distance outside the primary XYZ box on
each axis and Euclidean distance to the nearest primary hit. Family minima,
maxima and energy-weighted centroids are all saved so split deposits are not
hidden by a centroid alone. The minimum-distance pair also supplies absolute
Z separation and transverse R separation, satisfying D² = Z² + R². Ties
use the first pair in input hit order. `margin_mm` draws an illustrative frame and
measures excluded energy, without deciding acceptance.

Input coordinates/energies remain Float32. Distances and reductions use a
Float64 reference here; a Float32 sum comparison is saved. Coordinates in the
output are relative to the common primary vertex, preserving axis directions.
"""
function mc_fast_trigger(hits, particles; margin_mm::Float32=20f0,expected_primaries::Int=2)
    prim=[p for p in particles if p.primary==1]
    length(prim)==expected_primaries || error("Unexpected primary count")
    all(p->String(p.particle_name)=="e-",prim) || error("Non-electron primary")
    vertex=Float32[prim[1].initial_x,prim[1].initial_y,prim[1].initial_z]
    all(p->Float32[p.initial_x,p.initial_y,p.initial_z]==vertex,prim) || error("Primary vertices differ")
    roots=gamma_roots(particles);primary_ids=Set(Int32(p.particle_id) for p in prim)
    all(id->haskey(roots,id),hits.particle_id) || error("Hit particle missing from table")
    xyz=permutedims(hcat(hits.x_mm,hits.y_mm,hits.z_mm));xyz .-= vertex
    energy=Float64.(hits.energy_MeV.*1000f0)
    all(>=(0),energy) || error("Negative deposited energy")
    root=Int32[roots[id] for id in hits.particle_id]
    direct=map(id->id in primary_ids,hits.particle_id);free=root.==0
    any(direct) || error("No primary-electron hits")
    core=Float64.(xyz[:,direct]);lo=vec(minimum(core;dims=2));hi=vec(maximum(core;dims=2))
    outside=max.(lo.-xyz,xyz.-hi,0.0)
    framed=vec(any((xyz.<(lo.-margin_mm)).|(xyz.>(hi.+margin_mm));dims=1))
    nearest=zeros(Float64,length(energy));nearest_z=similar(nearest);nearest_r=similar(nearest)
    for j in findall(root.!=0)
        delta=core.-xyz[:,j];squared=vec(sum(abs2,delta;dims=1));k=argmin(squared)
        nearest[j]=sqrt(squared[k]);nearest_z[j]=abs(delta[3,k]);nearest_r[j]=hypot(delta[1,k],delta[2,k])
    end
    mh=EA.MCEventHits(Int64(0),EA.mc_vertex_frame,xyz[1,:],xyz[2,:],xyz[3,:],Float32.(energy),Int32.(hits.particle_id),Int32.(1:length(energy)))
    mi=EA.MCEventInfo(Int64(0),(0f0,0f0,0f0),collect(primary_ids))
    box=EA.primary_box(mh,mi);frame=EA.make_frame(box,(margin_mm,margin_mm,margin_mm))
    deposits=Dict(g.gamma_particle_id=>g for g in EA.gamma_deposits(mh,roots))
    families=NamedTuple[]
    byid=Dict(Int32(p.particle_id)=>p for p in particles)
    for id in sort(unique(root[root.!=0]))
        ix=findall(==(id),root);e=energy[ix];sum(e)>0 || continue
        c=vec(sum(Float64.(xyz[:,ix]).*permutedims(e);dims=2))./sum(e)
        d=vec(maximum(outside[:,ix];dims=2));dm=vec(minimum(outside[:,ix];dims=2))
        closest=ix[argmin(nearest[ix])]
        push!(families,(gamma_particle_id=id,energy_keV=deposits[id].energy_keV,initial_energy_keV=Float64(byid[id].kin_energy)*1000,
          n_hits=Int64(length(ix)),centroid_x_mm=c[1],centroid_y_mm=c[2],centroid_z_mm=c[3],
          outside_x_mm=d[1],outside_y_mm=d[2],outside_z_mm=d[3],
          outside_frame_r_mm=maximum(hypot.(max.(outside[1,ix].-margin_mm,0),max.(outside[2,ix].-margin_mm,0))),
          inside_frame=UInt8(contains(frame,deposits[id].bounds)),
          outside_frame_x_mm=max(d[1]-margin_mm,0),outside_frame_y_mm=max(d[2]-margin_mm,0),outside_frame_z_mm=max(d[3]-margin_mm,0),
          min_outside_x_mm=dm[1],min_outside_y_mm=dm[2],min_outside_z_mm=dm[3],
          nearest_hit_min_mm=nearest[closest],nearest_pair_z_mm=nearest_z[closest],nearest_pair_r_mm=nearest_r[closest],nearest_hit_max_mm=maximum(nearest[ix]),
          outside_primary_energy_keV=sum(energy[ix[vec(any(outside[:,ix].>0;dims=1))]]),
          outside_frame_energy_keV=sum(energy[ix[framed[ix]]]),
          creator_process=String(byid[id].creator_proc)))
    end
    full=extent(xyz);pr=extent(xyz[:,direct]);gf=extent(xyz[:,free])
    ev=(total_energy_keV=sum(energy),primary_energy_keV=sum(energy[direct]),
      gamma_free_energy_keV=sum(energy[free]),gamma_energy_keV=sum(energy[.!free]),
      primary_dx_mm=pr[1],primary_dy_mm=pr[2],primary_dz_mm=pr[3],
      gamma_free_dx_mm=gf[1],gamma_free_dy_mm=gf[2],gamma_free_dz_mm=gf[3],
      full_dx_mm=full[1],full_dy_mm=full[2],full_dz_mm=full[3],
      frame_dx_mm=pr[1]+2margin_mm,frame_dy_mm=pr[2]+2margin_mm,frame_dz_mm=pr[3]+2margin_mm,
      n_gamma_families=Int64(length(families)),n_primary_hits=Int64(count(direct)),
      n_nongamma_secondary_hits=Int64(count(free .& .!direct)),
      outside_frame_energy_keV=sum(energy[framed]),
      float32_sum_error_keV=Float64(sum(hits.energy_MeV))*1000-sum(energy))
    isapprox(ev.gamma_free_energy_keV+sum(f.energy_keV for f in families;init=0.0),ev.total_energy_keV;atol=0.01) || error("Energy partition failed")
    all(lo.<=core.<=hi) || error("Primary containment failed")
    (;event=ev,gammas=families,hits=(x_mm=xyz[1,:],y_mm=xyz[2,:],z_mm=xyz[3,:],energy_keV=Float32.(energy),
       particle_id=hits.particle_id,gamma_particle_id=root,primary=UInt8.(direct)),lo=Float32.(lo),hi=Float32.(hi))
end
end
