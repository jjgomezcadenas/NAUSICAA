using HDF5,Statistics
include("../../src/ion_signal/guarded_files.jl")
include("reconstruction_controls.jl")
path=ARGS[1];target=joinpath(dirname(path),"hit_cut_summary.h5");isfile(target)&&error("Preserve existing summary")
budget=2*1024^2;GuardedFiles.space_check(dirname(path),budget,1024^3);tmp=target*".tmp"
try
 h5open(path,"r") do f
  ids=read(f["events/event_number"])
  h5open(tmp,"w") do out
   out["schema"]="itaca.hit_cut_summary.v1"
   for arm in ("ions","pe")
    groups=[f["hit_histograms/$id/cut_statistics/$arm"] for id in ids]
    v(key)=[read(g[key]) for g in groups]
    result=Dict{String,Any}()
    for key in ("positive_hits","kept_hits","noise_hits_before","noise_hits_after","signal_hits_before","signal_hits_after","input_energy_valid_keV","input_energy_before_keV","input_energy_after_keV","input_energy_invalid_keV")
     result[key]=sum(v(key))
    end
    result["noise_retained_percent"]=100result["noise_hits_after"]/max(1,result["noise_hits_before"])
    result["signal_hits_retained_percent"]=100result["signal_hits_after"]/max(1,result["signal_hits_before"])
    result["additional_input_energy_lost_percent"]=100*(result["input_energy_before_keV"]-result["input_energy_after_keV"])/result["input_energy_valid_keV"]
    result["total_input_energy_lost_percent"]=100*(result["input_energy_valid_keV"]-result["input_energy_after_keV"])/result["input_energy_valid_keV"]
    for key in ("components_before","components_after","signal_components_before","signal_components_after","signal_largest_fraction_before","signal_largest_fraction_after","length_before_mm","length_after_mm","endpoint_mean_shift_mm","endpoint_max_shift_mm","max_distance_mm")
     values=v("topology/"*key);finite=filter(isfinite,values)
     result[key*"_mean"]=isempty(finite) ? NaN : mean(finite)
     result[key*"_median"]=isempty(finite) ? NaN : median(finite)
     result[key*"_p90"]=isempty(finite) ? NaN : quantile(finite,.9)
    end
    result["single_component_before"]=count(==(1),v("topology/components_before"))
    result["single_component_after"]=count(==(1),v("topology/components_after"))
    sb=v("topology/signal_components_before");sa=v("topology/signal_components_after")
    result["signal_single_before"]=count(==(1),sb);result["signal_single_after"]=count(==(1),sa)
    result["signal_single_to_multiple"]=count((sb.==1).&(sa.>1))
    result["signal_components_increased"]=count(sa.>sb)
    losses=100 .* (v("input_energy_before_keV")-v("input_energy_after_keV"))./v("input_energy_valid_keV")
    result["energy_loss_histogram"]=control_histogram(losses,histogram_edges(losses,60))
    shifts=filter(isfinite,v("topology/endpoint_mean_shift_mm"));result["endpoint_shift_histogram"]=control_histogram(shifts,histogram_edges(shifts,60))
    write_control_tree(create_group(out,arm),result)
    println(arm,": noise retained ",result["noise_retained_percent"],"%; input energy lost (additional) ",result["additional_input_energy_lost_percent"],"%; total ",result["total_input_energy_lost_percent"],"%; signal connected -> split ",result["signal_single_to_multiple"])
   end
   out["events"]=length(ids);flush(out);GuardedFiles.check_size(tmp,budget)
  end
 end
 mv(tmp,target)
catch
 isfile(tmp)&&rm(tmp);rethrow()
end
