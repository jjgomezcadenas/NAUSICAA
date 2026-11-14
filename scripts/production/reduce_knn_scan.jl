using HDF5,TOML
include("../../src/ion_signal/guarded_files.jl")
path=ARGS[1];target=length(ARGS)>1 ? ARGS[2] : joinpath(dirname(path),"knn_scan.h5");isfile(target)&&error("Preserve existing scan")
budget=2*1024^2;GuardedFiles.space_check(dirname(path),budget,1024^3);tmp=target*".tmp"
try
 h5open(path,"r") do f
  ids=read(f["events/event_number"]);card=TOML.parse(read(f["card"],String))
  h5open(tmp,"w") do out
   out["schema"]="itaca.knn_scan.v1";out["events"]=length(ids);out["card"]=read(f["card"],String)
   for arm in ("ions","pe"),k in card["track_analysis"]["graph_k_scan"]
    records=[f["hit_histograms/$id/cut_statistics/$arm/topology_scan/k_$k"] for id in ids]
    before=[read(g["signal_components_before"]) for g in records]
    after=[read(g["signal_components_after"]) for g in records]
    connected=before.==1;broken=connected.&(after.!=1)
    allafter=[read(g["components_after"]) for g in records]
    dest=create_group(out,"$arm/k_$k")
    dest["k"]=k;dest["baseline_connected"]=count(connected);dest["broken"]=count(broken)
    dest["broken_fraction"]=count(broken)/max(1,count(connected))
    dest["signal_disconnected_after"]=count(!=(1),after)
    dest["all_disconnected_after"]=count(!=(1),allafter)
    reference=[read(f["hit_histograms/$id/cut_statistics/$arm/topology_scan/k_10/signal_components_before"])==1 for id in ids]
    fixed_broken=count(reference.&(after.!=1))
    dest["reference_connected"]=count(reference);dest["reference_broken"]=fixed_broken
    dest["reference_broken_fraction"]=fixed_broken/max(1,count(reference))
    dest["broken_event_numbers"]=ids[broken]
    dest["signal_components_before"]=before;dest["signal_components_after"]=after
    println("$arm k=$k: $(count(broken))/$(count(connected)) broken ($(100count(broken)/max(1,count(connected)))%); all-hit disconnected $(count(!=(1),allafter))/$(length(ids))")
   end
   flush(out);GuardedFiles.check_size(tmp,budget)
  end
 end
 mv(tmp,target)
catch
 isfile(tmp)&&rm(tmp);rethrow()
end
