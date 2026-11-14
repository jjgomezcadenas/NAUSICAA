#!/usr/bin/env julia
# Distribution of ions per pad over the accepted events of a detector run:
# what output range the VH1 readout needs. Reads the detector-response
# controls for the event selection and placements, resimulates the ion
# detector response with the recorded seeds (no electronics), and writes
# per-pad totals and per-event maxima to results/readout/<run>/vh1_pad_charge/.
#
#     JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia" julia --startup-file=no --project=. scripts/production/vh1_pad_charge_scan.jl config/readout/<detector card> [events]
include("run_framed_detector_response.jl")
using Statistics
function pad_charge_scan(cardpath;limit=100,selection=nothing)
 c=TOML.parsefile(cardpath);base=joinpath(ROOT,"results/readout",c["run_label"])
 source=selection===nothing ? joinpath(base,"detector_response/controls.h5") : abspath(selection);hp=joinpath(ROOT,c["hits"])
 ids,placements,frames=h5open(source,"r") do f
  ev=f["events"];n=min(limit,length(read(ev["event_number"])));ids=read(ev["event_number"])[1:n]
  ps=[EventPlacement(id,read(ev["ion_drift_cm"])[j],read(ev["electron_drift_cm"])[j]) for (j,id) in enumerate(ids)]
  fs=[EventFrame(mc_vertex_frame,ntuple(k->read(ev["frame_$("xyz"[k])_min"])[j],3),ntuple(k->read(ev["frame_$("xyz"[k])_max"])[j],3)) for j in eachindex(ids)]
  ids,ps,fs
 end
 totals=Int64[];maxima=Int64[];event_totals=Int64[];bucket_max=Int64[]
 h5open(hp,"r") do input
  rowmap=Dict(id=>j for (j,id) in enumerate(read(input["events/event_id"])))
  for (number,id) in enumerate(ids)
   h,info=read_mc_event(input,id;event_row=rowmap[id]);h=center_hits(h,info);p=placements[number];frame=frames[number]
   selected=hits_in_frame(h,frame);ds=event_settings(c,p)
   rng(tag)=Xoshiro(seedfor(c["seed"],id,tag))
   pairs=sample_fano_pairs(rng("pairs"),selected,Float32(c["W_eV"]),Float32(c["fano_factor"]))
   ions=sim_ion_detector_response(rng("ions"),selected,pairs,p,frame,ds;max_vector_bytes=c["max_vector_bytes"])
   t=vec(sum(ions.ions.counts;dims=1));append!(totals,t[t.>0]);push!(maxima,maximum(t));push!(event_totals,sum(t))
   push!(bucket_max,maximum(ions.ions.counts))
  end
 end
 out=joinpath(base,"vh1_pad_charge");mkpath(out)
 open(joinpath(out,"pad_totals.csv"),"w") do io;println(io,"ions_per_pad");foreach(x->println(io,x),totals);end
 open(joinpath(out,"event_maxima.csv"),"w") do io;println(io,"event,max_ions_per_pad,event_ions,max_ions_per_bucket");for j in eachindex(ids);println(io,ids[j],",",maxima[j],",",event_totals[j],",",bucket_max[j]);end;end
 s=sort(totals);q(x)=s[clamp(ceil(Int,x*length(s)),1,length(s))]
 println("events ",length(ids)," pads with ions ",length(totals)," ions ",sum(totals))
 println("ions per pad: median ",q(.5)," 90% ",q(.9)," 99% ",q(.99)," 99.9% ",q(.999)," max ",maximum(totals))
 for lim in (10370,20000,30000,50000,100000)
  println("limit ",lim,": pads above ",count(>(lim),totals)," (",round(100count(>(lim),totals)/length(totals);digits=2),"%), ions on them ",round(100sum(totals[totals.>lim])/sum(totals);digits=2),"%, events with any ",count(>(lim),maxima))
 end
 println("max ions in one bucket ",maximum(bucket_max))
end
pad_charge_scan(abspath(ARGS[1]);limit=length(ARGS)>1 ? parse(Int,ARGS[2]) : 100,selection=length(ARGS)>2 ? ARGS[3] : nothing)
