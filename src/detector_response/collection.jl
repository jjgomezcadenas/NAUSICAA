_time(t::IonReadoutTime)=(t.start_ms,t.bucket_ms,t.bucket_count)
_time(t::ElectronReadoutTime)=(t.start_us,t.bucket_us,t.bucket_count)

"""Check event inputs before allocation or random sampling. Counts refer to supplied hits."""
function _check(h,pairs,p,f,s)
    h.event_id==p.event_id&&h.frame==f.frame==mc_vertex_frame || throw(ArgumentError("Event/coordinate mismatch"))
    length(pairs)==length(h.x_mm)&&all(>=(0),pairs) || throw(ArgumentError("Invalid pair counts"))
    all(i->contains(f,(h.x_mm[i],h.y_mm[i],h.z_mm[i])),eachindex(pairs)) || throw(ArgumentError("Select hits inside the frame first"))
    bounded=acquisition_frame(f,p,s.chamber_length_mm)
    bounded.lower_mm==f.lower_mm && bounded.upper_mm==f.upper_mm || throw(ArgumentError("Frame exceeds detector Z"))
    for tr in (s.ion_transport,s.electron_transport)
        all(x->isfinite(x)&&x>=0,(tr.sigma_t_mm,tr.sigma_l_mm)) || throw(ArgumentError("Invalid diffusion widths"))
    end
end

"""Diffuse counts directly into an unbounded physical channel/time dictionary.

Sample Float32 XYZ offsets per carrier, as in production_readout.jl. The common
event t0 determines bucket indices. Do not clip XY before physical channel
assignment. Aggregate excluded cells too, so losses and examples are measurable.
"""
function _collect(rng,h,counts,p,g,time,v,sigmat,sigmal,L,ion)
    start,dt,nb=_time(time);first=round(Int64,Float64(start)/dt)
    cells=Dict{NTuple{3,Int64},Int64}()
    for i in eachindex(counts),_ in 1:counts[i]
        x=h.x_mm[i]+randn(rng,Float32)*sigmat
        y=h.y_mm[i]+randn(rng,Float32)*sigmat
        z=h.z_mm[i]+randn(rng,Float32)*sigmal
        zd=Float64(z)+10*Float64(p.ion_drift_cm)
        t=(ion ? zd : Float64(L)-zd)/Float64(v)
        k=floor(Int64,t/Float64(dt))-first+1
        # The final acquisition endpoint belongs to the last bucket.
        t==(first+nb)*Float64(dt) && (k=Int64(nb))
        # Whole buckets must not admit diffusion beyond physical detector Z.
        if zd<0 || zd>Float64(L)
            k=t<0 ? Int64(0) : Int64(nb)+1
        end
        key=(EventAnalysis._direct_lattice_index(x,g.edge_mm[1],g.pitch_mm),
             EventAnalysis._direct_lattice_index(y,g.edge_mm[2],g.pitch_mm),k)
        cells[key]=get(cells,key,Int64(0))+1
    end
    cells
end

"""Build all-channel vectors, record XY/time losses, and retain excluded aggregate cells.

The dense allocation is checked before creating the matrix. max_vector_bytes is
an explicit memory bound, based on 8*channels*buckets, independent of input size.
"""
function _vectors(cells,channels,time;max_vector_bytes::Int=256*1024^2)
    n=length(channels.ix);nb=Int(time.bucket_count)
    bytes=Base.Checked.checked_mul(8,Base.Checked.checked_mul(n,nb))
    bytes<=max_vector_bytes || throw(ArgumentError("Detector vector needs $bytes bytes; limit is $max_vector_bytes"))
    data=zeros(Int64,nb,n)
    lookup=Dict((Int64(x),Int64(y))=>i for (i,(x,y)) in enumerate(zip(channels.ix,channels.iy)))
    allch=Set{Tuple{Int64,Int64}}();outch=Set{Tuple{Int64,Int64}}()
    excluded=Dict{NTuple{3,Int64},Int64}()
    total=Int64(0);xy=Int64(0);tm=Int64(0);out=Int64(0)
    for (key,c) in cells
        c>0 || continue
        x,y,k=key;total+=c;push!(allch,(x,y));j=get(lookup,(x,y),0)
        j==0 && (xy+=c;push!(outch,(x,y)))
        intime=1<=k<=nb
        !intime && (tm+=c)
        if j!=0&&intime;data[k,j]=c
        else;out+=c;excluded[key]=c;end
    end
    vectors=DetectorVectors(channels,time,data)
    losses=ReadoutLosses(total,xy,tm,out,length(allch),length(outch))
    @assert sum(data)+out==total
    vectors,losses,EventAnalysis._sparse_cells(excluded)
end

"""One-pass positive-cell statistics: event mean/max, total and occupied channels.

Zero cells are excluded only from mean/max. No threshold, sorting or event
reordering is used. Empty vectors return zero mean/max and positive-cell count.
"""
function vector_statistics(v::DetectorVectors)
    n=Int64(0);total=Int64(0);mx=Int64(0);channels=Int64(0)
    for j in axes(v.counts,2)
        active=false
        for i in axes(v.counts,1)
            c=v.counts[i,j]
            if c>0;n+=1;total+=c;mx=max(mx,c);active=true;end
        end
        channels+=active
    end
    (;positive_cells=n,total,mean=n==0 ? 0. : total/n,maximum=mx,channels_with_signal=channels)
end
