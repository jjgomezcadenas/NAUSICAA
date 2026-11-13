"""Limit a vertex-relative frame to detector Z=0…L, using the event placement."""
function acquisition_frame(f::EventFrame,p::EventPlacement,L::Float32)::EventFrame
    f.frame==mc_vertex_frame || throw(ArgumentError("Use vertex-relative coordinates"))
    isfinite(L)&&L>0&&isapprox(10f0*(p.ion_drift_cm+p.electron_drift_cm),L;atol=1f-3) ||
        throw(ArgumentError("Inconsistent chamber length"))
    z=10f0*p.ion_drift_cm
    lo=(f.lower_mm[1],f.lower_mm[2],max(f.lower_mm[3],-z))
    hi=(f.upper_mm[1],f.upper_mm[2],min(f.upper_mm[3],L-z))
    EventFrame(mc_vertex_frame,lo,hi)
end
"""Select hits in the frame without reordering; boundary hits remain included."""
function hits_in_frame(h::MCEventHits,f::EventFrame)::MCEventHits
    h.frame==f.frame || throw(ArgumentError("Coordinate mismatch"))
    k=findall(i->contains(f,(h.x_mm[i],h.y_mm[i],h.z_mm[i])),eachindex(h.x_mm))
    MCEventHits(h.event_id,h.frame,h.x_mm[k],h.y_mm[k],h.z_mm[k],h.energy_keV[k],h.particle_id[k],h.hit_id[k])
end
"""Enumerate every channel cell with positive XY overlap with the frame.

An acquired channel is read in full, even when its cell extends beyond the frame.
Zero-signal channels are included. Grid indexing uses the production lattice rule.
"""
function readout_channels(f::EventFrame,g::ReadoutGrid)::ReadoutChannels
    idx(x,d)=EventAnalysis._direct_lattice_index(x,g.edge_mm[d],g.pitch_mm)
    low=ntuple(d->idx(f.lower_mm[d],d),2)
    high=ntuple(d->begin
        i=idx(f.upper_mm[d],d)
        g.edge_mm[d]+Float32(i)*g.pitch_mm==f.upper_mm[d] ? i-1 : i
    end,2)
    all(d->high[d]>=low[d],1:2) || throw(ArgumentError("Frame needs positive XY area"))
    n=Base.Checked.checked_mul(high[1]-low[1]+1,high[2]-low[2]+1)
    n<=1_000_000 || throw(ArgumentError("Channel list exceeds 1,000,000 channels"))
    ix=Int32[];iy=Int32[]
    for x in low[1]:high[1],y in low[2]:high[2];push!(ix,x);push!(iy,y);end
    ReadoutChannels(ix,iy,g.edge_mm[1].+(Float32.(ix).+.5f0).*g.pitch_mm,
        g.edge_mm[2].+(Float32.(iy).+.5f0).*g.pitch_mm)
end
"""Return bucket-aligned start and count using a common t0 lattice, in native time units."""
function _interval(f,p,L,v,dt,ion)
    all(isfinite,(L,v,dt))&&min(L,v,dt)>0 || throw(ArgumentError("Positive length, velocity and bucket required"))
    z0=Float64(p.ion_drift_cm)*10
    lo=max(0.,Float64(f.lower_mm[3])+z0);hi=min(Float64(L),Float64(f.upper_mm[3])+z0)
    lo<=hi || throw(ArgumentError("Frame outside chamber"))
    a,b=ion ? (lo/v,hi/v) : ((L-hi)/v,(L-lo)/v)
    first=floor(Int64,a/Float64(dt));last=ceil(Int64,b/Float64(dt))
    Float32(first*Float64(dt)),Int32(max(1,last-first))
end
"""Ion interval from frame Z bounds; start and bucket width are ms after event t0."""
function ion_readout_time(f::EventFrame,p::EventPlacement,L::Float32,
        ion_drift_velocity_mm_ms::Float32,bucket_ms::Float32)::IonReadoutTime
    a,n=_interval(f,p,L,ion_drift_velocity_mm_ms,bucket_ms,true)
    IonReadoutTime(a,bucket_ms,n)
end
"""Electron interval from frame Z bounds; start and bucket width are μs after event t0."""
function electron_readout_time(f::EventFrame,p::EventPlacement,L::Float32,
        electron_drift_velocity_mm_us::Float32,bucket_us::Float32)::ElectronReadoutTime
    a,n=_interval(f,p,L,electron_drift_velocity_mm_us,bucket_us,false)
    ElectronReadoutTime(a,bucket_us,n)
end
"""Sample Fano pair counts from keV hits, preserving the production integer sampler.

The mean for each hit is 1000*energy_keV/W_eV, with target variance F*mean.
At small counts, the existing adjacent-integer minimum variance applies. The
returned counts are supplied unchanged to both response functions.
"""
function sample_fano_pairs(rng::AbstractRNG,h::MCEventHits,W_eV::Float32,F::Float32)::Vector{Int64}
    isfinite(W_eV)&&W_eV>0&&isfinite(F)&&0<=F<1 || throw(ArgumentError("Invalid W or F"))
    [sample_integer_moments(rng,Float64(e)*1000/Float64(W_eV),Float64(F)) for e in h.energy_keV]
end
