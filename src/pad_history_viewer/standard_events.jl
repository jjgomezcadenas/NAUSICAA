module PadHistoryNavigation
using HDF5,Random,SHA

"""Save a fixed-seed sample of events without a confirmed 10%-loss/one-ENC flag.
The initial pad is the event's lowest linear-voltage pad; all pads remain
available through the existing navigation. This is a browsing selection.
"""
function standard(c,name)
 source=joinpath(c.directory,"statistics.h5")
 isfile(source) || error("Complete the compression study to build the standard-event list")
 selection_hash=bytes2hex(sha256(read(@__FILE__)))
 stats_hash=bytes2hex(sha256(read(source)))
 destination=joinpath(c.directory,"navigation","standard_"*name*"_"*selection_hash[1:10]*"_"*stats_hash[1:10]*".h5")
 isfile(destination) && return destination
 h5open(source,"r") do stats
  g=stats[name];Bool(read(attributes(g)["full_candidates_complete"])) || error("The candidate study is incomplete")
  loss=read(stats["criterion_loss_fraction"]);noise=read(stats["criterion_noise_multiple"])
  criterion=findfirst((loss.==.1).&(noise.==1.))
  flags=read(g["criterion_$criterion/affected_pads"]);ids=read(g["event_id"])
  # Equal probability per eligible event, with a fixed seed for repeatable browsing.
  eligible=findall(==(0),flags);chosen=shuffle(MersenneTwister(c.card["random_seed"]),eligible)
  resize!(chosen,min(60,length(chosen)))
  screen=joinpath(c.root,c.readout["output"],name*".h5")
  h5open(screen,"r") do h
   npads=read(h["events/pad_count"]);ends=cumsum(npads)
   event=Int64[];px=Int64[];py=Int64[];peak=Int64[];voltage=Float64[]
   for j in chosen
    rows=(j==1 ? 1 : ends[j-1]+1):ends[j]
    v=h["pads/minimum_output_V"][rows];i=argmin(v);row=first(rows)+i-1
    push!(event,ids[j]);push!(voltage,v[i])
    push!(px,h["pads/pixel_x"][row]);push!(py,h["pads/pixel_y"][row]);push!(peak,h["pads/max_bucket_ions"][row])
   end
   mkpath(dirname(destination))
   h5open(destination,"w") do out
    attributes(out)["selection_sha256"]=selection_hash
    attributes(out)["statistics_sha256"]=stats_hash
    attributes(out)["seed"]=c.card["random_seed"]
    attributes(out)["eligible_events"]=length(eligible)
    attributes(out)["definition"]="Fixed-seed equal-probability sample of events without a confirmed 10% response loss and one-ENC deficit flag; initial pad has the lowest linear voltage in that event."
    for (key,value) in ("event_id"=>event,"pixel_x"=>px,"pixel_y"=>py,"max_bucket_ions"=>peak,"minimum_output_V"=>voltage);out[key]=value;end
   end
  end
 end
 destination
end
end
