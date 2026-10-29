module PadCompression
using HDF5,TOML,SHA,Random,Statistics,Test
include("pad_history_readout.jl")
const P=PadHistoryReadout
const B=P.B

"""Frozen source paths, calibrated kernels and cache identity for the stored pad histories."""
mutable struct Context
 root::String
 card::Dict{String,Any}
 readout::Dict{String,Any}
 directory::String
 identity::String
 kernels::Dict{String,Any}
 catalogues::Dict{String,Any}
end
function Context(root,cardpath)
 c=TOML.parsefile(cardpath);r=TOML.parsefile(joinpath(root,c["readout_card"]))
 hashes=Dict{String,String}()
 for path in (relpath(cardpath,root),c["readout_card"],"src/ion_signal/pad_compression.jl","src/ion_signal/pad_history_readout.jl","src/ion_signal/discharge_calibration.jl","src/ion_signal/bucket_study.jl","src/ion_signal/averaging_scan.jl",r["circuit_results"],r["circuit_analysis"],joinpath(r["source_decks"],"corrected.inc"),joinpath(r["source_decks"],"40ms_single_N400_r2.cir"))
  hashes[path]=P.sha(joinpath(root,path))
 end
 for sample in r["samples"]
  path=joinpath(r["output"],sample["name"]*".h5");hashes[path]=P.sha(joinpath(root,path))
 end
 identity=bytes2hex(sha256(join([k*":"*hashes[k] for k in sort(collect(keys(hashes)))],"\n")))
 dir=joinpath(root,c["output"],identity[1:12]);mkpath(dir)
 open(joinpath(dir,"provenance.toml"),"w") do io;TOML.print(io,Dict("identity"=>identity,"sha256"=>hashes));end
 Context(root,c,r,dir,identity,Dict{String,Any}(),Dict{String,Any}())
end
sample(c,name)=only(filter(s->s["name"]==name,c.readout["samples"]))
function kernel(c,name)
 get!(c.kernels,name) do
  s=sample(c,name);k=P.kernel(c.root,c.readout,s["frame_s"])
  noise=h5open(joinpath(c.root,c.readout["circuit_analysis"]),"r") do h
   g=h[c.readout["circuit_name"]]["$(round(Int,1000s["frame_s"]))ms"]
   j=read(attributes(g)["best_index"]);read(g["noise_rms_V"])[j]
  end
  baseline_error=h5open(joinpath(c.root,c.readout["circuit_results"]),"r") do h
   maximum(abs.(read(h["$(round(Int,1000s["frame_s"]))ms/single_N0_r2/V(out)"]).-k.baseline))
  end
  @assert baseline_error<1e-9
  (;k...,noise,enc=noise/k.gain,baseline_error)
 end
end
"""Event IDs and contiguous row offsets in source sparse cells and the screened pad table."""
function catalogue(c,name)
 get!(c.catalogues,name) do
  s=sample(c,name)
  h5open(joinpath(c.root,s["path"]),"r") do source
   h5open(joinpath(c.root,c.readout["output"],name*".h5"),"r") do screen
    @assert P.sha(joinpath(c.root,s["path"]))==read(attributes(screen)["source_sha256"])
    ids=read(source["events/event_id"]);sizes=read(source["events/cells_with_charge"])
    npads=read(screen["events/pad_count"])
    @assert ids==read(screen["events/event_id"])
    (;ids,cell_end=cumsum(sizes),pad_end=cumsum(npads),collected=read(source["events/collected_ions"]),positions=Dict(id=>i for (i,id) in enumerate(ids)))
   end
  end
 end
end
function eventrange(ends,i)
 (i==1 ? 1 : ends[i-1]+1):ends[i]
end
casekey(name,event,px,py)=name*"_"*string(event)*"_"*string(px)*"_"*string(py)
function casefolder(c,name,event,px,py)
 joinpath(c.directory,"cases",casekey(name,event,px,py))
end
"""Read one pad's counts without cuts, restoring every original zero bucket."""
function counts(c,name,event,px,py)
 cat=catalogue(c,name);haskey(cat.positions,event) || error("Event $event is absent")
 i=cat.positions[event];s=sample(c,name)
 h5open(joinpath(c.root,s["path"]),"r") do h
  rows=eventrange(cat.cell_end,i)
  x=h["cells/pixel_x"][rows];y=h["cells/pixel_y"][rows]
  sl=h["cells/time_slice"][rows];ions=h["cells/n_ions"][rows]
  edges=Float64.(read(h["event_axes/$event/time_edges_ms"])).*.001
  n=zeros(Int64,length(edges)-1);found=false
  for j in eachindex(x)
   if x[j]==px && y[j]==py;n[sl[j]+1]+=ions[j];found=true;end
  end
  found || error("Pad ($px, $py) has no collected charge in event $event")
  @assert all(n.>=0)
  n,edges
 end
end
function atomic_h5(writer,path)
 mkpath(dirname(path));temporary=path*".tmp"
 try
  h5open(writer,temporary,"w");mv(temporary,path;force=true)
 catch
  rm(temporary;force=true);rethrow()
 end
 path
end

"""Persist source counts and the calibrated linear prediction for one complete pad history."""
function history(c,name,event,px,py)
 folder=casefolder(c,name,event,px,py);path=joinpath(folder,"linear.h5")
 isfile(path) && return path
 n,edges=counts(c,name,event,px,py);k=kernel(c,name)
 extra=ceil(Int,c.card["recovery_tail_s"]/k.frame)
 truth=vcat(n,zeros(Int64,extra));times=vcat(edges,last(edges).+k.frame.*collect(1:extra))
 v,avg,d=P.response(truth,k)
 atomic_h5(path) do h
  for (key,value) in ("identity"=>c.identity,"sample"=>name,"event_id"=>event,"pixel_x"=>px,"pixel_y"=>py,"frame_s"=>k.frame,"baseline_V"=>k.baseline,"gain_V_per_ion"=>k.gain,"ENC_e"=>k.enc,"noise_rms_V"=>k.noise,"tau_s"=>k.tau,"alpha"=>k.alpha,"original_buckets"=>length(n),"averaging_s"=>c.readout["averaging_s"])
   attributes(h)[key]=value
  end
  h["counts"]=truth;h["edge_time_s"]=times
  h["linear_output_at_edges_V"]=v[1:length(times)]
  h["linear_average_V"]=avg[1:length(times)]
  h["linear_frame_signal_V"]=d[1:length(truth)]
  h["linear_estimated_ions"]=d[1:length(truth)]./k.gain
 end
end

"""Generate a finite-gap uniform-arrival current file for direct transistor simulation."""
function prepare(c,name,event,px,py;refine=false)
 path=history(c,name,event,px,py);folder=dirname(path);suffix=refine ? "refined" : "history"
 destination=joinpath(folder,suffix*".pwl")
 isfile(destination) && return path
 n,edges=counts(c,name,event,px,py);k=kernel(c,name)
 firstbin=findfirst(>(0),n);pre=c.card["prehistory_s"];offset=edges[firstbin]
 step=c.card["current_step_s"]/(refine ? 2 : 1)
 u=collect(range(-.005/3,0.,length=ceil(Int,.005/3/step)+1))
 ut=vcat(u,k.frame.+u);ui=B.current.(ut,k.frame,.00025,.005,3.)
 active=findall(>(0),n)
 stoptime=last(edges)-offset+pre+c.card["recovery_tail_s"]+c.readout["averaging_s"]
 t=sort!(unique(round.(vcat(0.,[x+pre+(j-firstbin)*k.frame for j in active for x in ut],stoptime);digits=12)))
 current=[sum(n[j]*((first(ut)<=x-pre-(j-firstbin)*k.frame<=last(ut)) ? B.interp(ut,ui,x-pre-(j-firstbin)*k.frame) : 0.) for j in active) for x in t]
 expected=sum(n)*B.E*(1-B.phi(.005,.00025));error=abs(B.integral(t,current)/expected-1)
 @assert error<1e-8
 open(destination*".tmp","w") do io
  for j in eachindex(t);println(io,t[j]," ",current[j]);end
 end
 mv(destination*".tmp",destination;force=true)
 open(joinpath(folder,suffix*"_input.toml"),"w") do io
  TOML.print(io,Dict("offset_s"=>offset,"prehistory_s"=>pre,"stoptime_s"=>stoptime,"charge_C"=>expected,"charge_error"=>error,"current_step_s"=>step))
 end
 path
end

"""Reduce a full SPICE transient to frame signals and measured compression criteria."""
function reduce(c,name,event,px,py;refine=false,rawpath=nothing)
 path=history(c,name,event,px,py);folder=dirname(path)
 suffix=refine ? "refined" : "history";destination=joinpath(folder,refine ? "refined.h5" : "spice.h5")
 rawpath===nothing && (rawpath=joinpath(folder,suffix*"_raw.h5"))
 input=TOML.parsefile(joinpath(folder,suffix*"_input.toml"));k=kernel(c,name)
 h5open(path,"r") do lin
  h5open(rawpath,"r") do raw
   t=read(raw["time"]);v=read(raw["V(out)"]);ii=read(raw["I(Is1)"])
   cumulative=vcat(0.,cumsum(diff(t).*((v[1:end-1].-k.baseline).+(v[2:end].-k.baseline))./2))
   function primitive(x)
    x<=first(t) && return 0.
    @assert x<=last(t)+1e-9
    j=clamp(searchsortedlast(t,x),1,length(t)-1);dx=x-t[j]
    cumulative[j]+(v[j]-k.baseline)*dx+(v[j+1]-v[j])/(t[j+1]-t[j])*dx^2/2
   end
   edges=read(lin["edge_time_s"]);n=read(lin["counts"]);times=edges.-input["offset_s"].+input["prehistory_s"]
   width=c.readout["averaging_s"]
   avg=[(primitive(x+width/2)-primitive(x-width/2))/width for x in times]
   point=[x<first(t) ? k.baseline : B.interp(t,v,x) for x in times]
   d=-(avg[2:end].-k.alpha.*avg[1:end-1]);linear=read(lin["linear_frame_signal_V"])
   deficit=max.(0.,(linear.-d)./k.gain)
   relative=[linear[j]>0 ? max(0.,1-d[j]/linear[j]) : 0. for j in eachindex(n)]
   cuts=[(r,s) for r in c.card["loss_fractions"] for s in c.card["deficit_noise_multiples"]]
   flags=[(n.>0).&(relative.>=r).&(deficit.>=s*k.enc) for (r,s) in cuts]
   charge=B.integral(t,ii);charge_error=abs(charge/input["charge_C"]-1);@assert charge_error<1e-7
   quiet=findall(==(0),n);positive=findall(>(0),n)
   # The plotted curve retains the exact voltage minimum and all frame boundaries.
   keep=sort!(unique(vcat(collect(1:10:length(t)),findmin(v)[2],length(t))))
   atomic_h5(destination) do h
    for key in keys(attributes(lin));attributes(h)[key]=read(attributes(lin)[key]);end
    for key in keys(attributes(raw));attributes(h)[key]=read(attributes(raw)[key]);end
    attributes(h)["baseline_validation_max_error_V"]=k.baseline_error
    for (key,value) in ("raw_sha256"=>P.sha(rawpath),"minimum_output_V"=>minimum(v),"maximum_output_V"=>maximum(v),"charge_error"=>charge_error,"maximum_relative_loss"=>maximum(relative[positive]),"maximum_deficit_e"=>maximum(deficit[positive]),"positive_deficit_sum_e"=>sum(deficit[positive]),"maximum_empty_residual_e"=>maximum(abs.((linear[quiet].-d[quiet])./k.gain)),"refined"=>refine)
     attributes(h)[key]=value
    end
    for (key,value) in ("edge_time_s"=>edges,"counts"=>n,"output_at_edges_V"=>point,"average_at_edges_V"=>avg.+k.baseline,"frame_signal_V"=>d,"estimated_ions"=>d./k.gain,"linear_frame_signal_V"=>linear,"linear_estimated_ions"=>linear./k.gain,"relative_loss"=>relative,"deficit_e"=>deficit,"criterion_loss_fraction"=>first.(cuts),"criterion_noise_multiple"=>last.(cuts),"affected_buckets"=>[count(f) for f in flags],"affected_ions"=>[sum(n[f]) for f in flags],"time_s"=>t[keep].-input["prehistory_s"].+input["offset_s"],"output_V"=>v[keep])
     h[key]=value
    end
   end
  end
 end
 destination
end

"""List all stored pads below the screen and fixed-seed controls in each higher-voltage stratum."""
function plan(c)
 path=joinpath(c.directory,"plan.h5");isfile(path) && return path
 rng=MersenneTwister(c.card["random_seed"])
 atomic_h5(path) do out
  attributes(out)["identity"]=c.identity
  for s in c.readout["samples"]
   name=s["name"];g=create_group(out,name)
   h5open(joinpath(c.root,c.readout["output"],name*".h5"),"r") do h
    v=read(h["pads/minimum_output_V"]);selected=findall(v.<c.card["screen_voltage_V"]);strata=fill(0,length(selected))
    edges=c.card["control_voltage_edges_V"]
    populations=Int64[];draws=Int64[]
    for j in 1:length(edges)-1
     possible=findall((v.>=edges[j]).&(v.<edges[j+1]));push!(populations,length(possible))
     chosen=shuffle(rng,possible)[1:min(length(possible),c.card["controls_per_stratum"])]
     append!(selected,chosen);append!(strata,fill(j,length(chosen)));push!(draws,length(chosen))
    end
    for key in ("event_id","pixel_x","pixel_y","max_bucket_ions","minimum_output_V")
     g[key]=read(h["pads/"*key])[selected]
    end
    g["screen_row"]=selected;g["stratum"]=strata
    g["control_population"]=populations;g["control_draws"]=draws
    attributes(g)["candidate_pads"]=count(==(0),strata);attributes(g)["events"]=read(attributes(h)["events"])
   end
  end
 end
end

"""Aggregate confirmed per-pad response loss to event fractions and affected event charge."""
function aggregate(c)
 planpath=plan(c);path=joinpath(c.directory,"statistics.h5")
 cuts=[(r,s) for r in c.card["loss_fractions"] for s in c.card["deficit_noise_multiples"]]
 atomic_h5(path) do out
  attributes(out)["identity"]=c.identity;attributes(out)["screen_voltage_V"]=c.card["screen_voltage_V"]
  attributes(out)["definition"]="Positive bucket: relative loss = max(0,1-D_spice/D_linear). Affected charge is collected ion count in flagged buckets / complete event collected count. Event fractions count confirmed cases in the tested candidates and controls. Quiet-frame signal spreading is present in both predictions."
  out["criterion_loss_fraction"]=first.(cuts);out["criterion_noise_multiple"]=last.(cuts)
  h5open(planpath,"r") do p
   for name in keys(p)
    g=create_group(out,name);cat=catalogue(c,name);pg=p[name]
    event=read(pg["event_id"]);px=read(pg["pixel_x"]);py=read(pg["pixel_y"]);stratum=read(pg["stratum"])
    buckets=zeros(Int64,length(cat.ids),length(cuts));ions=similar(buckets);fill!(ions,0);pads=similar(buckets);fill!(pads,0)
    complete=falses(length(event));maximumloss=fill(NaN,length(event));deficit=zeros(length(cat.ids));minvoltage=fill(NaN,length(event))
    controlsflag=zeros(Int64,length(cuts));refineerrors=Float64[]
    for i in eachindex(event)
     file=joinpath(casefolder(c,name,event[i],px[i],py[i]),"spice.h5")
     isfile(file) || continue
     h5open(file,"r") do h
      @assert read(attributes(h)["identity"])==c.identity
      complete[i]=true;maximumloss[i]=read(attributes(h)["maximum_relative_loss"]);minvoltage[i]=read(attributes(h)["minimum_output_V"])
      j=cat.positions[event[i]];b=read(h["affected_buckets"]);q=read(h["affected_ions"])
      buckets[j,:].+=b;ions[j,:].+=q;pads[j,:].+=(b.>0);deficit[j]+=read(attributes(h)["positive_deficit_sum_e"])
      if stratum[i]>0;controlsflag.+=(b.>0);end
     end
     refined=joinpath(dirname(file),"refined.h5")
     if isfile(refined)
      h5open(file,"r") do coarse
       h5open(refined,"r") do fine
        dd=read(coarse["frame_signal_V"]);df=read(fine["frame_signal_V"])
        push!(refineerrors,maximum(abs.(dd.-df))/max(maximum(abs.(df)),1e-12))
       end
      end
     end
    end
    g["event_id"]=cat.ids;g["event_collected_ions"]=cat.collected
    g["completed"]=Int8.(complete);g["plan_event_id"]=event;g["plan_pixel_x"]=px;g["plan_pixel_y"]=py;g["stratum"]=stratum
    g["pad_maximum_relative_loss"]=maximumloss;g["pad_minimum_output_V"]=minvoltage
    g["control_pads_with_flag"]=controlsflag;g["refinement_error"]=refineerrors
    g["event_positive_deficit_fraction"]=deficit./cat.collected
    for j in eachindex(cuts)
     cg=create_group(g,"criterion_"*string(j))
     cg["affected_pads"]=pads[:,j];cg["affected_buckets"]=buckets[:,j];cg["affected_ion_fraction"]=ions[:,j]./cat.collected
     n=count(>(0),pads[:,j]);attributes(cg)["event_count"]=n;attributes(cg)["event_fraction"]=n/length(cat.ids)
     # Wilson interval is conditional on the screened study coverage.
     z=1.959963984540054;ne=length(cat.ids);ph=n/ne;center=(ph+z^2/(2ne))/(1+z^2/ne);half=z*sqrt(ph*(1-ph)/ne+z^2/(4ne^2))/(1+z^2/ne)
     attributes(cg)["wilson_low"]=center-half;attributes(cg)["wilson_high"]=center+half
    end
    for (key,value) in ("planned"=>length(event),"completed"=>count(complete),"candidate_planned"=>count(==(0),stratum),"candidate_completed"=>count(complete.&(stratum.==0)),"control_planned"=>count(>(0),stratum),"control_completed"=>count(complete.&(stratum.>0)),"events"=>length(cat.ids),"full_candidates_complete"=>all(complete[stratum.==0]))
     attributes(g)[key]=value
    end
   end
  end
 end
end

"""Write UI catalogue or pad choices as HDF5 without changing event-analysis sources."""
function navigation(c,op,name,event,path)
 cat=catalogue(c,name)
 atomic_h5(path) do out
  if op=="events"
   out["event_id"]=cat.ids
  else
   i=get(cat.positions,event,0);i>0 || error("Event $event is absent")
   h5open(joinpath(c.root,c.readout["output"],name*".h5"),"r") do h
    rows=eventrange(cat.pad_end,i)
    for key in ("pixel_x","pixel_y","max_bucket_ions","total_ions","minimum_output_V")
     out[key]=h["pads/"*key][rows]
    end
   end
  end
 end
end
end
