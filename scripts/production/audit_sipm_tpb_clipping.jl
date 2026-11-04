#!/usr/bin/env julia
# Audit event-level clipping tails in the completed TPB pilot.
using HDF5,TOML,Statistics
const ROOT=normpath(joinpath(@__DIR__,"../.."))
cd(ROOT)
module M
include(joinpath(pwd(),"scripts/production/reduce_electron_operating_points.jl"))
end
path="results/sipm_tpb_pilot/analysis.h5"
row=h5open(path,"r") do f
 g=f["tau12p5_n1024/ep2/pure_xe"]
 Dict(k=>read(g[k]) for k in ("manifest","event_id","sensor_total","asic_total","clip_loss","detected_corrected","sensor_corrected","asic_corrected"))
end
manifest=TOML.parsefile(String(row["manifest"]))
energies=h5open("events/0nubb_hitdf.h5","r") do f
 Dict(zip(read(f["events/event_id"]),Float64.(read(f["events/energy_MeV"]))))
end
noclip=Float64[];drifts=Float64[]
for shard in manifest["shards"]
 shard["gas"]=="pure_xe" || continue
 h5open(shard["sensor"],"r") do sf
  h5open(shard["asic"],"r") do af
   sensor=TOML.parse(String(read(sf["metadata/sensor_card_toml"])))
   dt=read(sf["metadata/bucket_width_s"])
   electron=h5open(String(read(sf["metadata/source_path"])),"r") do pf
    TOML.parse(String(read(pf["metadata/branch_card_toml"])))
   end
   for id in shard["event_ids"]
    s=sf["events/$id"];a=af["events/$id"]
    L=read(s["electron_drift_cm"]);push!(drifts,L)
    # Algebraic inversion of saved upper clipping: restores the signed noisy input sum.
    total=read(a["readout_pe_total"])+read(a["clip_loss_pe_equiv"])
    corrected,_=M.corrected_total(total,read(s["acquired_cells"]),dt,sensor,electron,L)
    push!(noclip,corrected/energies[id])
   end
  end
 end
end
width(x)=100std(x;corrected=false)/mean(x)
fraction=row["clip_loss"]./row["sensor_total"]
zero=iszero.(row["clip_loss"])
r=row["asic_corrected"];contrib=(r.-mean(r)).^2
worst=sortperm(fraction;rev=true)
output=Dict("no_clip_counterfactual_rms_percent"=>width(noclip),
 "final_rms_percent"=>width(r),"sensor_rms_percent"=>width(row["sensor_corrected"]),
 "detected_rms_percent"=>width(row["detected_corrected"]),
 "unclipped_subset_count"=>count(zero),"unclipped_subset_rms_percent"=>width(r[zero]),
 "unclipped_subset_gaussian_equiv_fwhm_percent"=>sqrt(8log(2))*width(r[zero]),
 "mean_event_clip_loss_percent"=>100mean(fraction),
 "std_event_clip_loss_percentage_points"=>100std(fraction;corrected=false),
 "sample_charge_loss_percent"=>100sum(row["clip_loss"])/sum(row["sensor_total"]),
 "worst_event_ids"=>Int64.(row["event_id"][worst[1:10]]),
 "worst_loss_percent"=>100fraction[worst[1:10]],"worst_drift_cm"=>drifts[worst[1:10]],
 "top2_variance_fraction_percent"=>100sum(sort(contrib;rev=true)[1:2])/sum(contrib),
 "top5_variance_fraction_percent"=>100sum(sort(contrib;rev=true)[1:5])/sum(contrib))
h5open("results/sipm_tpb_pilot/clipping_audit.h5","w") do f
 for (k,v) in output;f[k]=v;end
 f["source"]=read(@__FILE__,String);f["complete"]=true
end
for (k,v) in sort!(collect(output);by=first);println(k," = ",v);end
