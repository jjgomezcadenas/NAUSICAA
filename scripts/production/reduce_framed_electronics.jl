include("run_framed_detector_response.jl")
"""Reduce saved electronics example vectors for noise plots and timing summaries."""
function reduce_electronics(path)
 target=joinpath(dirname(path),"plot_controls.h5");tmp=target*".tmp";budget=4*1024^2
 GuardedFiles.space_check(dirname(path),budget,1024^3)
 try
 h5open(path,"r") do f
  h5open(tmp,"w") do out
   if haskey(f["events"],"peak_current_exceedances")
    counts=read(f["events/peak_current_exceedances"])
    out["summary/events_with_current_exceedances"]=count(>(0),counts)
    out["summary/pads_with_current_exceedances"]=sum(counts)
    out["summary/maximum_current_ratio"]=maximum(read(f["events/maximum_current_ratio"]))
   end
   out["summary/event_count"]=length(read(f["events/event_number"]))
   out["summary/adc_out_of_range"]=sum(read(f["events/adc_below"]))+sum(read(f["events/adc_above"]))
   out["summary/asic_clipped_buckets"]=sum(read(f["events/asic_clipped_buckets"]))
   out["summary/nonlinear_channels"]=sum(read(f["events/nonlinear_channels"]))
   for name in filter(n->haskey(f["events"],n),["ion_s","sipm_s","ion_total_s","setup_s","current_check_s","ion_allocated_bytes","ion_prepared_bytes"])
    values=read(f["events/$name"]);out["timing/$name/mean_excluding_first"]=mean(values[2:end]);out["timing/$name/median_excluding_first"]=median(values[2:end])
   end
   sizes=create_group(out,"vector_sizes")
   for (name,field) in [("adc","ion_vector_bytes"),("asic","sipm_vector_bytes")]
    values=Float64.(read(f["events/$field"]))./1024^2
    group=create_group(sizes,name)
    histograms(group,[(;value=x) for x in values],20)
    group["value/median"]=median(values);group["value/maximum"]=maximum(values)
   end
   for id in keys(f["examples"]),name in keys(f["examples/$id"])
    g=f["examples/$id/$name"]
    for field in ("input", startswith(name,"ions") ? "adc_midpoint_V" : "response_pe")
     values=Float64.(read(g[field]));stats=create_group(out,"statistics/$id/$name/$field")
     stats["entries"]=length(values);stats["mean"]=mean(values);stats["std"]=std(values;corrected=false);stats["sum"]=sum(values)
    end
    if startswith(name,"ions")
      x=read(g["adc_midpoint_V"])-read(g["voltage_V"])
    else
      all(iszero,read(g["input"]))||continue
      x=Float64.(read(g["response_pe"]))
    end
    h=create_group(out,"examples/$id/$name");rows=[(;value=y) for y in x];histograms(h,rows,50)
    h["lag1_covariance"]=mean((x[1:end-1].-mean(x)).*(x[2:end].-mean(x)))
   end
   reference=joinpath(dirname(dirname(path)),"electronics_noise_reference","controls.h5")
   if isfile(reference) && h5open(h->read(h["events/event_number"])==read(f["events/event_number"]),reference,"r")
    h5open(reference,"r") do old
     read(old["events/event_number"])==read(f["events/event_number"])||error("Comparison events differ")
     for field in ["ion_s","sipm_s"]
      before=mean(read(old["events/$field"])[2:end]);after=mean(read(f["events/$field"])[2:end])
      out["comparison/$field/before_s"]=before;out["comparison/$field/after_s"]=after;out["comparison/$field/speedup"]=before/after
     end
     # A new noise sampler changes individual ion noise vectors, not the physical signal or SiPMs.
     for event in keys(f["examples"]),name in keys(f["examples/$event"])
      g=f["examples/$event/$name"];r=old["examples/$event/$name"]
      fields=startswith(name,"ions") ? ["input","voltage_V","current_A"] : ["input","response_pe","clipped"]
      for field in fields;read(g[field])==read(r[field])||error("Changed signal: event $event $name $field");end
     end
     out["comparison/signals_unchanged"]=true
     out["comparison/reference_sha256"]=bytes2hex(open(sha256,reference))
    end
   end
   previous=joinpath(dirname(dirname(path)),"electronics_100events","controls.h5")
   if haskey(f["events"],"ion_total_s") && isfile(previous)
    h5open(previous,"r") do old
     read(old["events/event_number"])==read(f["events/event_number"])||error("Comparison events differ")
     before=read(old["events/setup_s"])+read(old["events/current_check_s"])+read(old["events/ion_s"])
     after=read(f["events/ion_total_s"])
     out["efficiency/before_s"]=mean(before[2:end]);out["efficiency/after_s"]=mean(after[2:end])
     out["efficiency/speedup"]=mean(before[2:end])/mean(after[2:end])
     out["efficiency/predicted_8000_s"]=8000mean(after[2:end])
     out["efficiency/max_prepared_MiB"]=maximum(read(f["events/ion_prepared_bytes"]))/1024^2
     residual=0.;adc_changes=0;adc_max_difference=0
     for event in keys(f["examples"]),name in keys(f["examples/$event"])
      g=f["examples/$event/$name"];r=old["examples/$event/$name"]
      read(g["input"])==read(r["input"])||error("Input changed")
      if startswith(name,"ions")
       read(g["current_A"])==read(r["current_A"])||error("Current changed")
       residual=max(residual,maximum(abs.(read(g["voltage_V"])-read(r["voltage_V"]))))
       delta=abs.(Int.(read(g["adc_counts"]))-Int.(read(r["adc_counts"])))
       adc_changes+=count(!iszero,delta);adc_max_difference=max(adc_max_difference,maximum(delta))
      else
       read(g["response_pe"])==read(r["response_pe"])||error("SiPM response changed")
       read(g["clipped"])==read(r["clipped"])||error("SiPM clipping changed")
      end
     end
     residual<1e-5||error("Linear recurrence disagrees with prior examples")
     adc_max_difference<=1||error("ADC change exceeds one count")
     for field in ["linear_channels","nonlinear_channels","peak_current_exceedances","adc_below","adc_above"]
      read(old["events/$field"])==read(f["events/$field"])||error("Changed selection: $field")
     end
     out["efficiency/maximum_voltage_change_V"]=residual
     out["efficiency/changed_adc_samples"]=adc_changes
     out["efficiency/maximum_adc_change"]=adc_max_difference
     out["efficiency/reference_sha256"]=bytes2hex(open(sha256,previous))
    end
   end
   out["source_sha256"]=bytes2hex(open(sha256,path));out["reduction_code_sha256"]=bytes2hex(open(sha256,@__FILE__))
   flush(out);GuardedFiles.check_size(tmp,budget)
  end
 end
 mv(tmp,target;force=true)
 catch
 isfile(tmp)&&rm(tmp);rethrow()
 end
end
if abspath(PROGRAM_FILE)==@__FILE__
 p=isempty(ARGS) ? joinpath(ROOT,"results/readout/bb0nu_full_ion3mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42/electronics_100events_efficient/controls.h5") : abspath(only(ARGS))
 reduce_electronics(p)
end
