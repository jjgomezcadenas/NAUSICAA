"""
    rebin_electron_signal(signal, factor) -> ReconstructedPE

Sum consecutive corrected PE buckets in each SiPM. Preserve signed charge,
channel order and the common event-t0 start. Keep the last incomplete group
and the original acquisition end. A group is valid only if every contributing
bucket is valid; clipping flags use logical OR. Invalid sums are NaN.
"""
function rebin_electron_signal(s::ReconstructedPE,factor::Int64)
 factor>0||throw(ArgumentError("Positive rebin factor required"))
 n=cld(size(s.pe,1),factor);nc=size(s.pe,2)
 out=fill(NaN32,n,nc);valid=falses(n,nc);clip=falses(n,nc)
 for j in 1:nc,k in 1:n
  rows=(k-1)*factor+1:min(k*factor,size(s.pe,1))
  valid[k,j]=all(view(s.valid,rows,j));clip[k,j]=any(view(s.asic_clipped,rows,j))
  valid[k,j]&&(out[k,j]=sum(Float64,view(s.pe,rows,j)))
 end
 time=ElectronReadoutTime(s.time.start_us,s.time.bucket_us*factor,Int32(n))
 ReconstructedPE(s.event_id,s.channels,time,s.acquired_end_us,out,valid,clip)
end
