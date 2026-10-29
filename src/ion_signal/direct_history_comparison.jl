module DirectHistoryComparison
using HDF5,TOML,Statistics,SHA,Test
sha(p)=open(io->bytes2hex(sha256(io)),p)
"""Compare historical and direct-collected input histories and completed CSA statistics.
Source counts remain exact integers; source hit coordinates and energy tables
are preserved. Pulse statistics refer to positive pad/frame cells.
"""
function run(root)
 cp=joinpath(root,"config/pad_history_direct_comparison.toml");c=TOML.parsefile(cp)
 old=TOML.parsefile(joinpath(root,c["historical_readout_card"]));new=TOML.parsefile(joinpath(root,c["direct_readout_card"]))
 dir=joinpath(root,c["output"]);mkpath(dir)
 h5open(joinpath(dir,"comparison.h5"),"w") do out
  attributes(out)["config_sha256"]=sha(cp);attributes(out)["analysis_sha256"]=sha(@__FILE__)
  @testset "Historical versus direct input comparison" begin
   for sample in new["samples"]
    name=sample["name"];before=only(filter(s->s["name"]==name,old["samples"]))
    g=create_group(out,name);attributes(g)["frame_s"]=sample["frame_s"];attributes(g)["pitch_m"]=sample["pitch_m"]
    h5open(joinpath(root,new["circuit_analysis"]),"r") do ref
     calibration=ref[new["circuit_name"]]["$(round(Int,1000sample["frame_s"]))ms"]
     index=read(attributes(calibration)["best_index"])
     attributes(g)["nominal_ENC_e"]=read(calibration["noise_rms_V"])[index]/read(calibration["gain_V_per_ion"])[index]
    end
    h5open(joinpath(root,before["path"]),"r") do a
     h5open(joinpath(root,sample["path"]),"r") do b
      @test read(a["events/event_id"])==read(b["events/event_id"])
      @test read(a["events/seed"])==read(b["events/seed"])
      @test read(a["events/drift_cm"])==read(b["events/drift_cm"])
      @test read(a["events/sigma_mm"])==read(b["events/sigma_mm"])
      @test read(a["events/produced_ions"])==read(b["events/produced_ions"])
      for (label,input,path) in (("historical",a,before["path"]),("direct",b,sample["path"]))
       dest=create_group(g,label);at=attributes(dest);n=read(input["cells/n_ions"])
       @test all(n.>0)
       collected=read(input["events/collected_ions"]);produced=read(input["events/produced_ions"]);outside=read(input["events/outside_ions"])
       @test sum(n)==sum(collected)
       @test produced==collected.+outside
       for (key,value) in ("source"=>path,"source_sha256"=>sha(joinpath(root,path)),"schema"=>read(attributes(input)["schema"]),"events"=>length(collected),"occupied_cells"=>length(n),"collected_ions"=>sum(n),"mean_ions_per_occupied_cell"=>sum(n)/length(n),"median_ions_per_occupied_cell"=>median(n),"maximum_bucket_ions"=>maximum(n),"mean_event_mean_ions"=>read(input["distribution_statistics/mean/mean_ions"]),"median_event_mean_ions"=>read(input["distribution_statistics/mean/median_ions"]))
        at[key]=value
       end
       dest["pulse_quantile_probability"]=[.1,.5,.9,.99,1.];dest["pulse_quantile_ions"]=quantile(n,[.1,.5,.9,.99,1.])
       si=findfirst(==(3790),read(input["events/event_id"]))
       sz=read(input["events/cells_with_charge"]);lo=sum(sz[1:si-1])+1;hi=lo+sz[si]-1
       edges=read(input["event_axes/3790/time_edges_ms"]);perframe=zeros(Int64,length(edges)-1)
       for (j,q) in zip(input["cells/time_slice"][lo:hi],input["cells/n_ions"][lo:hi]);perframe[j+1]+=q;end
       dest["event_3790/time_edges_ms"]=edges;dest["event_3790/ions_per_frame"]=perframe
       attributes(dest["event_3790"])["occupied_buckets"]=count(>(0),perframe)
       attributes(dest["event_3790"])["original_buckets"]=length(perframe)
       @test sum(perframe)==collected[si]
       screenpath=joinpath(root,(label=="historical" ? old : new)["output"],name*".h5")
       h5open(screenpath,"r") do screen
        for key in ("pads","pads_below_guard","events_with_guard_risk","minimum_linear_output_V")
         at[key]=read(attributes(screen)[key])
        end
       end
      end
      attributes(g)["same_retained_ions"]=read(a["events/retained_ions"])==read(b["events/retained_ions"])
      attributes(g)["same_outside_ions"]=read(a["events/outside_ions"])==read(b["events/outside_ions"])
     end
    end
    for (label,key) in (("historical","historical_statistics"),("direct","direct_statistics"))
     path=joinpath(root,c[key])
     h5open(path,"r") do stats
      dest=g[label];attributes(dest)["statistics_source"]=c[key];attributes(dest)["statistics_sha256"]=sha(path)
      sg=stats[name]
      @test read(attributes(sg)["completed"])==read(attributes(sg)["planned"])
      losses=read(stats["criterion_loss_fraction"]);noise=read(stats["criterion_noise_multiple"])
      dest["criterion_loss_fraction"]=losses;dest["criterion_noise_multiple"]=noise
      dest["compression_event_fraction"]=[read(attributes(sg["criterion_$i"])["event_fraction"]) for i in eachindex(losses)]
      dest["compression_event_count"]=[read(attributes(sg["criterion_$i"])["event_count"]) for i in eachindex(losses)]
      attributes(dest)["simulated_pads"]=read(attributes(sg)["completed"])
      attributes(dest)["control_pads"]=read(attributes(sg)["control_completed"])
     end
    end
   end
  end
 end
 println(joinpath(dir,"comparison.h5"))
end
end
