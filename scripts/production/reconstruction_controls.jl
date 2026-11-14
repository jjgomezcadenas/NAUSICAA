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
"""Histograms for all valid cells, cells with input, zero-input cells and wholly
zero-input channels. Threshold loss uses input-weighted charge; invalid signal
is reported separately. No suppression is applied to reconstructed vectors."""
function zs_controls(values,valid,input,thresholds,bins)
 size(values)==size(input)==size(valid)||throw(DimensionMismatch("Control vectors"))
 good=valid .& isfinite.(values);signal=input.>0;empty=.!signal
 quiet=repeat(reshape(vec(all(iszero,input;dims=1)),1,:),size(input,1),1)
 allvalues=values[good];isempty(allvalues)&&error("No valid reconstruction cells")
 edges=histogram_edges(allvalues,bins)
 hist=Dict{String,Any}()
 for (name,mask) in [("all",good),("signal",good.&signal),("zero_input",good.&empty),("quiet_channels",good.&quiet)]
  hist[name]=control_histogram(values[mask],edges)
 end
 zeros=values[good.&empty]
 noise=Dict{String,Any}()
 if !isempty(zeros)
  # A noise-only range avoids losing the pedestal in the full signal range.
  nedges=histogram_edges(zeros,bins)
  noise["zero_input"]=control_histogram(zeros,nedges)
  noise["quiet_channels"]=control_histogram(values[good.&quiet],nedges)
 end
 total=sum(Float64,input);invalid=sum(Float64,input[.!good];init=0.)
 noise_count=count(good.&empty);quiet_count=count(good.&quiet)
 ns=Int64[];qs=Int64[];loss=Float64[];cells=Int64[]
 for cut in thresholds
  survive=good .& (values.>cut)
  push!(ns,count(survive.&empty));push!(qs,count(survive.&quiet));push!(cells,count(survive))
  push!(loss,sum(Float64,input[good .& .!survive];init=0.))
 end
 Dict{String,Any}("histograms"=>hist,"noise_histograms"=>noise,"thresholds"=>thresholds,
 "noise_above"=>ns,"noise_cells"=>noise_count,"quiet_above"=>qs,"quiet_cells"=>quiet_count,
 "input_below"=>loss,"input_total"=>total,"input_invalid"=>invalid,"cells_above"=>cells,
 "valid_cells"=>count(good),"invalid_cells"=>count(.!good),"input_sum"=>total,
 "reconstructed_sum"=>sum(Float64,values[good];init=0.),"input_valid_sum"=>sum(Float64,input[good];init=0.),
 "signal_residual"=>control_histogram(Float64.(values[good.&signal]).-input[good.&signal],
 histogram_edges(vcat(Float64.(values[good.&signal]).-input[good.&signal],0.),bins)))
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
