"""Numerical controls for reconstructed signal and proposed zero-suppression cuts.
Truth labels are used only here, never inside signal reconstruction."""
function control_histogram(x,edges)
 values=Float64.(filter(isfinite,x));counts=zeros(Int64,length(edges)-1)
 under=0;over=0
 for v in values
  if v<first(edges);under+=1
  elseif v>last(edges);over+=1
  else;counts[clamp(searchsortedlast(edges,v),1,length(counts))]+=1;end
 end
 Dict{String,Any}("edges"=>edges,"counts"=>counts,"entries"=>length(values),"underflow"=>under,"overflow"=>over,
 "mean"=>isempty(values) ? NaN : mean(values),"std"=>isempty(values) ? NaN : std(values;corrected=false),
 "median"=>isempty(values) ? NaN : median(values),"rms"=>isempty(values) ? NaN : sqrt(mean(abs2,values)))
end
function histogram_edges(values,bins)
 lo,hi=extrema(values);lo==hi&&(hi=lo+1.)
 collect(range(lo,hi;length=bins+1))
end
"""Group channel–bucket entries by input signal in that bucket.
Histogram signed reconstructed values and record invalid entries.
"""
function channel_histograms(values,valid,input,bins)
 size(values)==size(input)==size(valid)||throw(DimensionMismatch("Control vectors"))
 good=valid .& isfinite.(values)
 has_signal=input.>0
 hist=Dict{String,Any}()
 for (name,mask) in (("buckets_with_signal",has_signal),("buckets_without_signal",.!has_signal))
  selected=values[good .& mask]
  edges=isempty(selected) ? collect(range(-1.,1.;length=bins+1)) : histogram_edges(selected,bins)
  h=control_histogram(selected,edges)
  h["channels"]=count(vec(any(mask;dims=1)))
  h["invalid_buckets"]=count(mask .& .!good)
  hist[name]=h
 end
 Dict{String,Any}("histograms"=>hist)
end
"""Sum truth counts in the same groups used for electron rebinning."""
function grouped_input(x,factor)
 n=cld(size(x,1),factor);out=zeros(Float64,n,size(x,2))
 for j in axes(x,2),k in 1:n;out[k,j]=sum(Float64,view(x,(k-1)*factor+1:min(k*factor,size(x,1)),j));end
 out
end
function write_control_tree(group,data)
 for (name,value) in data
  if value isa AbstractDict;write_control_tree(create_group(group,name),value)
  else;group[name]=value;end
 end
end
