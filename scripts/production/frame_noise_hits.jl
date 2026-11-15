# Ion noise hits in the part of the fixed frame that is not simulated.
using Distributions: Normal, Poisson, ccdf, quantile, pdf

"""Draw the ion noise hits of the frame cells outside the simulated vectors.

A pad with no ions gives, after reconstruction, Gaussian noise of standard
deviation `sigma_ions` per pad and bucket. Only cells above the hit cut become
hits, so their number is Poisson with mean (free cells) x P(noise > cut) and
each value follows the Gaussian tail above the cut. Cells are the pads with
positive overlap with the frame and the buckets of its Z range, on the same
pad lattice and time lattice as the simulated vectors; the simulated pad and
bucket cells are excluded. Correlations between buckets of a pad are neglected.
Return hits in vertex-relative mm and keV, as `reconstruct_hits` does.
"""
function frame_noise_hits(rng,frame::EventFrame,ri,placement,detector,sigma_ions::Real,keV_per_count::Real,cut_keV::Real)
    sigma_ions>0&&keV_per_count>0&&cut_keV>0 || throw(ArgumentError("Noise, scale and cut must be positive"))
    pads=readout_channels(frame,detector.ion_grid);v=Float64(detector.ion_transport.ion_drift_velocity_mm_ms)
    time=ion_readout_time(frame,placement,detector.chamber_length_mm,detector.ion_transport.ion_drift_velocity_mm_ms,ri.time.bucket_ms)
    dt=Float64(time.bucket_ms);first=round(Int,Float64(time.start_ms)/dt);count_k=Int(time.bucket_count)
    sim_first=round(Int,Float64(ri.time.start_ms)/dt);sim=sim_first:(sim_first+Int(ri.time.bucket_count)-1)
    simulated=Set(zip(ri.channels.ix,ri.channels.iy))
    pad_simulated=[(pads.ix[j],pads.iy[j]) in simulated for j in eachindex(pads.ix)]
    overlap=length(intersect(first:(first+count_k-1),sim))
    free=length(pads.ix)*count_k-count(pad_simulated)*overlap
    threshold=Float64(cut_keV)/Float64(keV_per_count)/Float64(sigma_ions);tail=ccdf(Normal(),threshold)
    x=Float32[];y=Float32[];z=Float32[];q=Float32[];energy=Float32[]
    (free>0&&tail>0) || return ReconstructedHits(x,y,z,q,energy)
    wanted=rand(rng,Poisson(free*tail));vertex_z=10Float64(placement.ion_drift_cm)
    while length(x)<wanted
        j=rand(rng,eachindex(pads.ix));k=first+rand(rng,0:count_k-1)
        pad_simulated[j]&&(k in sim)&&continue
        charge=Float64(sigma_ions)*quantile(Normal(),1-tail*rand(rng))
        e=Float32(charge*keV_per_count)
        e>Float32(cut_keV)||(e=nextfloat(Float32(cut_keV)))
        push!(x,pads.x_mm[j]);push!(y,pads.y_mm[j]);push!(z,(k+.5)*dt*v-vertex_z);push!(q,charge);push!(energy,e)
    end
    ReconstructedHits(x,y,z,q,energy)
end

"""Append the second set of hits to the first."""
join_hits(a::ReconstructedHits,b::ReconstructedHits)=ReconstructedHits(vcat(a.x_mm,b.x_mm),vcat(a.y_mm,b.y_mm),vcat(a.z_mm,b.z_mm),vcat(a.charge,b.charge),vcat(a.energy_keV,b.energy_keV))
