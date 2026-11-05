module PadHistoryReadout
using HDF5,TOML,SHA,Statistics,Test
include("discharge_calibration.jl")
const D=DischargeCalibration
const B=D.B
sha(p)=open(io->bytes2hex(sha256(io)),p)
settings(root)=TOML.parsefile(joinpath(root,"config/pad_history_readout.toml"))
function kernel(root,c,frame)
 key="$(round(Int,1000frame))ms"
 h5open(joinpath(root,c["circuit_analysis"]),"r") do a
  at=attributes(a[c["circuit_name"]][key]);tau=read(at["tau_s"]);alpha=exp(-frame/tau)
  h5open(joinpath(root,c["circuit_results"]),"r") do h
   t,v=D.waveform(h[key],100,2);t.-=c["prehistory_s"];v./=-100
   if first(t)>-c["averaging_s"]/2;t=vcat(-c["averaging_s"]/2-.0001,t);v=vcat(0.,v);end
   raw=[B.interp(t,v,j*frame) for j in 0:2]
   avg=[D.average(t,v,j*frame,c["averaging_s"]) for j in 0:2]
   gain=avg[2]-alpha*avg[1]
   @assert abs(B.interp(t,v,3frame)/(raw[3]*alpha)-1)<1e-4
   baseline=only(read(h["op/V(out)"]))
   (;frame,tau,alpha,raw,avg,gain,baseline)
  end
 end
end
"""Propagate all buckets with a calibrated packet kernel and exponential late tail."""
function response(counts,k)
 n=length(counts);raw=zeros(n+3);avg=zeros(n+3);oldr=0.;olda=0.
 for boundary in 0:n+2
  if boundary>=2
   q=boundary-1<=n ? counts[boundary-1] : 0
   oldr=k.alpha*oldr+k.raw[3]*q;olda=k.alpha*olda+k.avg[3]*q
  end
  now=1<=boundary<=n ? counts[boundary] : 0
  future=boundary<n ? counts[boundary+1] : 0
  raw[boundary+1]=k.baseline-(oldr+k.raw[2]*now+k.raw[1]*future)
  avg[boundary+1]=k.baseline-(olda+k.avg[2]*now+k.avg[1]*future)
 end
 measured=-((avg[2:end].-k.baseline).-k.alpha.*(avg[1:end-1].-k.baseline))
 raw,avg,measured
end
function columns(h)
 (px=read(h["cells/pixel_x"]),py=read(h["cells/pixel_y"]),slice=read(h["cells/time_slice"]),ions=read(h["cells/n_ions"]))
end
function screen(root)
 c=settings(root);dir=joinpath(root,c["output"]);mkpath(dir)
 @testset "Complete pad-history accounting and kernel" begin
  for sample in c["samples"]
   source=joinpath(root,sample["path"]);k=kernel(root,c,sample["frame_s"])
   h5open(source,"r") do input
    @test read(input["settings/frame_ms"])/1000≈sample["frame_s"]
    @test read(input["settings/pitch_mm"])/1000≈sample["pitch_m"]
    ids=read(input["events/event_id"]);sizes=read(input["events/cells_with_charge"]);collected=read(input["events/collected_ions"])
    x=columns(input);@test sum(sizes)==length(x.ions)
    names=(:event_id,:pixel_x,:pixel_y,:total_ions,:max_bucket_ions,:positive_buckets,:n_buckets,:minimum_output_V,:minimum_average_V,:maximum_frame_signal_V,:first_guard_slice,:guard_buckets,:guard_positive_buckets)
    values=Dict{Symbol,Any}(key=>Float64[] for key in names)
    for key in names[1:7];values[key]=Int64[];end
    for key in names[11:13];values[key]=Int64[];end
    candidates=Dict{String,Any}();event_pads=zeros(Int64,length(ids));event_ions=zeros(Int64,length(ids));offset=0
    for (event_index,id) in enumerate(ids)
     lo=offset+1;hi=offset+sizes[event_index];offset=hi
     edges=Float64.(read(input["event_axes/$id/time_edges_ms"])).*.001
     nb=length(edges)-1
     @assert maximum(abs.(diff(edges).-k.frame))<1e-6
     order=sortperm(lo:hi;by=i->(x.px[i],x.py[i],x.slice[i]));order .+=lo-1
     # sortperm(range) returns positions, shifted to global indices.
     cursor=1
     while cursor<=length(order)
      row=order[cursor];px=x.px[row];py=x.py[row];n=zeros(Int64,nb)
      while cursor<=length(order) && x.px[order[cursor]]==px && x.py[order[cursor]]==py
       j=order[cursor];bin=x.slice[j]+1
       @assert 1<=bin<=nb && x.ions[j]>0
       n[bin]+=x.ions[j];cursor+=1
      end
      rv,av,dv=response(n,k);minimumv=minimum(rv)
      guard=rv[2:nb+1].<c["range_guard_low_V"]
      firstguard=findfirst(guard);total=sum(n);peak=maximum(n)
      rowvalues=(id,px,py,total,peak,count(>(0),n),nb,minimumv,minimum(av),maximum(dv),isnothing(firstguard) ? -1 : firstguard-1,count(guard),count(guard .& (n.>0)))
      for (key,value) in zip(names,rowvalues);push!(values[key],value);end
      event_pads[event_index]+=1;event_ions[event_index]+=total
      for (label,score) in (("worst_voltage",minimumv),("largest_bucket",-Float64(peak)),("near_guard",abs(minimumv-c["range_guard_low_V"])))
       if !haskey(candidates,label) || score<candidates[label].score
        candidates[label]=(;score,id,px,py,counts=copy(n),edges=copy(edges),raw=rv,avg=av,signal=dv)
       end
      end
     end
     @assert event_ions[event_index]==collected[event_index]
    end
    path=joinpath(dir,sample["name"]*".h5")
    h5open(path,"w") do h
     a=attributes(h)
     for (key,value) in ("source"=>sample["path"],"source_sha256"=>sha(source),"config_sha256"=>sha(joinpath(root,"config/pad_history_readout.toml")),"code_sha256"=>sha(@__FILE__),"frame_s"=>k.frame,"pitch_m"=>sample["pitch_m"],"baseline_V"=>k.baseline,"tau_s"=>k.tau,"alpha"=>k.alpha,"gain_V_per_ion"=>k.gain,"range_guard_low_V"=>c["range_guard_low_V"]);a[key]=value;end
     a["assumptions"]="Uniform arrivals inside stored buckets; calibrated 100-ion finite-gap SPICE packet response with exponential late tail. Full zero buckets included. Independent event starts at quiescent baseline. Noise-free linear screening; the voltage guard is not an established transistor saturation threshold."
     for key in names;h["pads/"*string(key),chunk=(min(length(values[key]),16384),),compress=6]=values[key];end
     h["events/event_id"]=ids;h["events/pad_count"]=event_pads;h["events/collected_ions"]=event_ions
     for (key,value) in ("events"=>length(ids),"pads"=>length(values[:event_id]),"collected_ions"=>sum(event_ions),"pads_below_guard"=>count(<(c["range_guard_low_V"]),values[:minimum_output_V]),"pads_below_ground"=>count(<(0.),values[:minimum_output_V]),"minimum_linear_output_V"=>minimum(values[:minimum_output_V]),"maximum_bucket_ions"=>maximum(values[:max_bucket_ions]),"events_with_guard_risk"=>length(unique(values[:event_id][values[:minimum_output_V].<c["range_guard_low_V"]])),"guard_positive_cells"=>sum(values[:guard_positive_buckets]))
      attributes(h)[string(key)]=value
     end
     h["kernel/raw_V_per_ion"]=k.raw;h["kernel/average_V_per_ion"]=k.avg
     for (name,case) in candidates
      for (key,value) in ("event_id"=>case.id,"pixel_x"=>case.px,"pixel_y"=>case.py);h["candidates/$name/$key"]=value;end
      h["candidates/$name/counts"]=case.counts;h["candidates/$name/edges_s"]=case.edges
      h["candidates/$name/linear_output_V"]=case.raw;h["candidates/$name/linear_average_V"]=case.avg;h["candidates/$name/linear_frame_signal_V"]=case.signal
     end
    end
    @test event_ions==collected
    @test sum(values[:total_ions])==sum(x.ions)
    @test all(values[:positive_buckets].<=values[:n_buckets])
    @test response(zeros(Int,10),k)[1]==fill(k.baseline,13)
    @test maximum(abs.(response(zeros(Int,10),k)[3]))==0.
    println(sample["name"]," events=",length(ids)," pads=",length(values[:event_id])," minV=",minimum(values[:minimum_output_V])," max bucket=",maximum(values[:max_bucket_ions])," below guard=",count(<(c["range_guard_low_V"]),values[:minimum_output_V]));flush(stdout)
   end
   GC.gc()
  end
 end
end
function spice_inputs(root)
 c=settings(root);dir=joinpath(root,c["output"]);file=joinpath(dir,"spice_inputs.h5")
 isfile(file) && error("Preserving existing SPICE inputs")
 h5open(file,"w") do out
  attributes(out)["generator_sha256"]=sha(@__FILE__)
  for sample in c["samples"]
   source=joinpath(dir,sample["name"]*".h5")
   h5open(source,"r") do h
    for label in ("worst_voltage","largest_bucket","near_guard")
     orig=h["candidates/$label"];counts=read(orig["counts"]);edges=read(orig["edges_s"])
     firstbin=findfirst(>(0),counts);lastbin=findlast(>(0),counts)
     key=sample["name"]*"_"*label;g=create_group(out,key);a=attributes(g)
     for id in ("event_id","pixel_x","pixel_y");a[id]=read(orig[id]);end
     a["sample"]=sample["name"];a["label"]=label;a["frame_s"]=sample["frame_s"]
     a["source_sha256"]=sha(source);a["first_nonzero_bin"]=firstbin-1
     a["event_time_offset_s"]=edges[firstbin];a["prehistory_s"]=.003
     g["counts"]=counts;g["edges_s"]=edges
     g["linear_output_V"]=read(orig["linear_output_V"])
     g["linear_average_V"]=read(orig["linear_average_V"])
     g["linear_frame_signal_V"]=read(orig["linear_frame_signal_V"])
     for (r,step) in enumerate(c["max_steps_s"])
      frame=sample["frame_s"];gap=.005;v=3.;half=.00025
      u=collect(range(-gap/v,0.,length=ceil(Int,gap/v/step)+1))
      ut=vcat(u,frame.+u);ui=B.current.(ut,frame,half,gap,v)
      @assert abs(B.integral(ut,ui)/(B.E*(1-B.phi(gap,half)))-1)<1e-10
      active=findall(>(0),counts)
      times=sort!(unique(round.(vcat(0.,[x+.003+(k-firstbin)*frame for k in active for x in ut],last(edges)-edges[firstbin]+.003+c["recovery_tail_s"]);digits=12)))
      current=[sum(counts[k]*((first(ut)<=x-.003-(k-firstbin)*frame<=last(ut)) ? B.interp(ut,ui,x-.003-(k-firstbin)*frame) : 0.) for k in active) for x in times]
      @assert abs(B.integral(times,current)/(sum(counts)*B.E*(1-B.phi(gap,half)))-1)<1e-8
      group=create_group(g,"r$r");group["time_s"]=times;group["current_A"]=current
      attributes(group)["charge_C"]=B.integral(times,current)
     end
    end
   end
  end
 end
end

end
