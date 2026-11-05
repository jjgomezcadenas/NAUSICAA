module PadHistorySpice
using HDF5,TOML,SHA,Test,Statistics
include("pad_history_readout.jl")
const P=PadHistoryReadout
const D=P.D
const B=P.B
function analyze(root)
 c=P.settings(root);dir=joinpath(root,c["output"])
 h5open(joinpath(dir,"spice_inputs.h5"),"r") do inp
  h5open(joinpath(dir,"spice_analysis.h5"),"w") do out
   attributes(out)["analysis_sha256"]=P.sha(@__FILE__);attributes(out)["inputs_sha256"]=P.sha(joinpath(dir,"spice_inputs.h5"))
   @testset "Complete-history SPICE validation" begin
    for name in keys(inp)
     source=joinpath(dir,"spice",name*".h5");isfile(source) || continue
     h5open(source,"r") do raw
      g=create_group(out,name);a=attributes(g)
      for key in keys(attributes(raw));a[key]=read(attributes(raw)[key]);end
      a["raw_sha256"]=P.sha(source)
      if haskey(a,"acquisition_failure");println(name," acquisition failed");return;end
      case=inp[name];counts=read(case["counts"]);edges=read(case["edges_s"])
      pre=read(a["prehistory_s"]);offset=read(a["event_time_offset_s"]);firstbin=read(a["first_nonzero_bin"])
      frame=read(a["frame_s"]);k=P.kernel(root,c,frame);width=c["averaging_s"]
      times=edges[firstbin+1:end].-offset.+pre
      stoptime=read(raw["history_r2/time"])[end]-width/2
      nextra=max(0,floor(Int,(stoptime-last(times))/frame))
      times=vcat(times,last(times).+frame.*collect(1:nextra))
      truth=vcat(counts[firstbin+1:end],zeros(Int64,nextra))
      # Baseline before the first retained nonzero bucket is exactly the quiescent pad.
      function reduce(r)
       d=raw["history_r$r"];b=raw["baseline_r$r"]
       t=read(d["time"]);v=read(d["V(out)"]);tb=read(b["time"]);vb=read(b["V(out)"])
       delta=v.-[B.interp(tb,vb,x) for x in t]
       point=[B.interp(t,v,x) for x in times]
       cumulative=vcat(0.,cumsum(diff(t).*(delta[1:end-1].+delta[2:end])./2))
       primitive(x)=begin
        j=clamp(searchsortedlast(t,x),1,length(t)-1);dx=x-t[j]
        cumulative[j]+delta[j]*dx+(delta[j+1]-delta[j])/(t[j+1]-t[j])*dx^2/2
       end
       avg=[(primitive(x+width/2)-primitive(x-width/2))/width for x in times]
       signal=-(avg[2:end].-k.alpha.*avg[1:end-1])
       (;t,v,delta,point,avg,signal,tb,vb)
      end
      fine=reduce(2);coarse=reduce(1)
      scale=max(maximum(abs.(fine.signal)),1e-8)
      error=maximum(abs.(fine.signal.-coarse.signal))/scale
      charge=B.integral(fine.t,read(raw["history_r2/I(Is1)"]))
      expected=sum(counts)*B.E*(1-B.phi(.005,.00025))
      chargeerror=abs(charge/expected-1)
      predicted_full,_,linear_full=P.response(truth,k)
      predicted=predicted_full[1:length(times)]
      delta_predicted=predicted.-k.baseline
      point_error=maximum(abs.((fine.point.-[B.interp(fine.tb,fine.vb,x) for x in times]).-delta_predicted))
      linear_signal=linear_full[1:length(truth)]
      signal_error=maximum(abs.(fine.signal.-linear_signal))
      below=fine.v.<c["range_guard_low_V"];nearrail=fine.v.<.05
      lastcollection=edges[findlast(>(0),counts)+1]-offset+pre
      after=findall(x->x>=lastcollection,fine.t)
      recovery=NaN
      if any(below)
       lastbelow=findlast(below)
       if lastbelow<length(fine.t);recovery=max(0.,fine.t[lastbelow+1]-lastcollection);end
      else;recovery=0.;end
      modelnoise=h5open(joinpath(root,c["circuit_analysis"]),"r") do ref
       sourcegroup=ref[c["circuit_name"]]["$(round(Int,frame*1000))ms"]
       j=read(attributes(sourcegroup)["best_index"]);read(sourcegroup["noise_rms_V"])[j]
      end
      post=findall(times[2:end].>=lastcollection+frame)
      falsepeak=maximum(abs.(fine.signal[post].-linear_signal[post]))/k.gain
      badpost=filter(j->abs(fine.signal[j]-linear_signal[j])>modelnoise,post)
      useful_recovery=isempty(badpost) ? 0. : max(0.,times[maximum(badpost)+1]-lastcollection)
      a["post_collection_false_peak_e"]=falsepeak
      residual=fine.signal[post].-linear_signal[post]
      a["post_collection_positive_peak_e"]=max(0.,maximum(residual)/k.gain)
      a["post_collection_negative_peak_e"]=min(0.,minimum(residual)/k.gain)
      a["recovery_to_one_ENC_s"]=useful_recovery
      a["quiet_frames_above_one_ENC"]=count(j->abs(fine.signal[j]-linear_signal[j])>modelnoise,post)
      a["baseline_noise_rms_V"]=modelnoise
      for (key,value) in ("minimum_output_V"=>minimum(fine.v),"maximum_output_V"=>maximum(fine.v),"baseline_minimum_V"=>minimum(fine.vb),"baseline_maximum_V"=>maximum(fine.vb),"crossed_guard"=>any(below),"below_50mV"=>any(nearrail),"recovery_after_collection_s"=>recovery,"transient_refinement_error"=>error,"charge_error"=>chargeerror,"maximum_boundary_error_V"=>point_error,"maximum_frame_signal_error_V"=>signal_error,"maximum_frame_signal_V"=>maximum(fine.signal),"input_charge_C"=>charge,"range_guard_low_V"=>c["range_guard_low_V"])
       a[key]=value
      end
      g["edge_time_s"]=times.-pre.+offset;g["counts"]=truth
      g["output_at_edges_V"]=fine.point;g["linear_output_at_edges_V"]=predicted
      g["frame_signal_V"]=fine.signal;g["linear_frame_signal_V"]=linear_signal
      # Exact solver extrema and dense samples around input activity remain available.
      keep=unique(vcat(collect(1:20:length(fine.t)),findmin(fine.v)[2],length(fine.t)))
      sort!(keep);g["time_s"]=fine.t[keep].-pre.+offset;g["output_V"]=fine.v[keep]
      @test chargeerror<1e-6
      @test error<.005
      @test minimum(fine.v)>-.1
      @test maximum(fine.v)<c["supply_V"]+.1
      println(name," min=",minimum(fine.v)," guard=",any(below)," rail=",any(nearrail)," recovery=",recovery," refinement=",error," charge=",chargeerror," frameerr=",signal_error," quiet false peak=",falsepeak," recovery1ENC=",useful_recovery);flush(stdout)
     end
    end
   end
  end
 end
end
function summarize(root)
 c=P.settings(root);dir=joinpath(root,c["output"])
 h5open(joinpath(dir,"summary.h5"),"w") do h
  attributes(h)["analysis_sha256"]=P.sha(@__FILE__)
  attributes(h)["config_sha256"]=P.sha(joinpath(root,"config/pad_history_readout.toml"))
  for sample in c["samples"]
   source=joinpath(dir,sample["name"]*".h5")
   h5open(source,"r") do f
    g=create_group(h,sample["name"]);a=attributes(g)
    for key in keys(attributes(f));a[key]=read(attributes(f)[key]);end
    a["screen_sha256"]=P.sha(source)
    a["pad_guard_fraction"]=read(a["pads_below_guard"])/read(a["pads"])
    a["event_guard_fraction"]=read(a["events_with_guard_risk"])/read(a["events"])
    thresholds=collect(-1.2:.005:.9);v=sort!(read(f["pads/minimum_output_V"]))
    g["voltage_threshold_V"]=thresholds
    g["pad_fraction_below_threshold"]=[searchsortedlast(v,x)/length(v) for x in thresholds]
   end
  end
  h5open(joinpath(dir,"spice_analysis.h5"),"r") do f
   rootgroup=create_group(h,"spice")
   for name in keys(f)
    g=create_group(rootgroup,name)
    for key in keys(attributes(f[name]));attributes(g)[key]=read(attributes(f[name])[key]);end
   end
   attributes(rootgroup)["source_sha256"]=P.sha(joinpath(dir,"spice_analysis.h5"))
  end
  adc=create_group(h,"adc")
  h5open(joinpath(root,c["circuit_analysis"]),"r") do f
   for frame in (.02,.04)
    key="$(round(Int,frame*1000))ms";a=attributes(f[c["circuit_name"]][key]);j=read(a["best_index"])
    gg=read(f[c["circuit_name"]][key]["gain_V_per_ion"])[j];sigma=read(f[c["circuit_name"]][key]["noise_rms_V"])[j];alpha=read(a["alpha"])
    rows=[(span,bits,span/2.0^bits) for span in (1.,2.) for bits in (10,12,13,16)]
    g=create_group(adc,key)
    g["span_V"]=[r[1] for r in rows];g["bits"]=[r[2] for r in rows];g["lsb_V"]=[r[3] for r in rows]
    qnoise=[r[3]*sqrt(1+alpha^2)/sqrt(12) for r in rows]
    g["frame_quantization_rms_V"]=qnoise;g["quantization_ENC_e"]=qnoise./gg
    g["total_noise_relative_increase"]=[sqrt(sigma^2+x^2)/sigma-1 for x in qnoise]
    attributes(g)["assumptions"]="Independent uniform quantization errors in the two averaged readings; no extra reduction credited for ADC oversampling. ADC input/buffer noise not included."
   end
  end
 end
end

end
