module HistoryTimeAudit
using HDF5,TOML,SHA,Test
const ROOT=normpath(joinpath(@__DIR__,"../.."))
sha(p)=open(io->bytes2hex(sha256(io)),p)
"""Audit the original event time axis against source-hit z and diffusion-grid resolution.
Preserve Float32 source coordinates/energies and integer counts in HDF5. Arrival
times are relative to the first retained event z edge, for increasing-z arrival.
"""
function run(event=3790)
 outdir=joinpath(ROOT,"results/ltspice/history_time_audit");mkpath(outdir)
 path=joinpath(outdir,"event_$event.h5")
 h5open(path,"w") do out
  attributes(out)["code_sha256"]=sha(@__FILE__)
  @testset "Event time-axis audit" begin
   for sample in ("bb0nu","electrons")
    file=joinpath(ROOT,"results/pixel_charge_1mm",sample*".h5")
    h5open(file,"r") do h
     ids=read(h["events/event_id"]);j=findfirst(==(event),ids);j===nothing && return
     sizes=read(h["events/cells_with_charge"]);lo=sum(sizes[1:j-1])+1;hi=lo+sizes[j]-1
     slices=h["cells/time_slice"][lo:hi];n=h["cells/n_ions"][lo:hi]
     zedges=read(h["event_axes/$event/z_edges_mm"]);tedges=read(h["event_axes/$event/time_edges_ms"])
     settings=TOML.parse(read(h["provenance/detector_toml"],String));speed=settings["gas"]["bulk_speed_m_s"]
     total=zeros(Int64,length(tedges)-1)
     for (k,q) in zip(slices,n);total[k+1]+=q;end
     g=create_group(out,sample)
     g["time_edges_ms"]=tedges;g["z_edges_mm"]=zedges;g["collected_ions_per_frame"]=total
     a=attributes(g);a["pixel_source_sha256"]=sha(file);a["bulk_speed_m_s"]=speed
     a["first_z_mm"]=first(zedges);a["last_z_mm"]=last(zedges);a["duration_s"]=Float64(last(tedges))/1000
     a["original_buckets"]=length(total);a["occupied_event_buckets"]=count(>(0),total)
     a["fraction_in_final_second"]=sum(total[tedges[1:end-1].>=last(tedges)-1000])/sum(total)
     @test sum(total)==read(h["events/collected_ions"])[j]
     @test maximum(abs.(Float64.(tedges).-Float64.(zedges.-first(zedges))./speed))<.01
     source=read(h["provenance/source_path"],String)
     h5open(source,"r") do hits
      i=findfirst(==(event),read(hits["events/event_id"]));start=read(hits["events/first_row"])[i];nh=read(hits["events/n_hits"])[i]
      rows=start:start+nh-1
      z=hits["hits/z_mm"][rows];energy=hits["hits/energy_MeV"][rows];pid=hits["hits/particle_id"][rows]
      for (key,value) in ("source_hit_z_mm"=>z,"source_hit_energy_MeV"=>energy,"source_hit_particle_id"=>pid);g[key]=value;end
      a["source_hit_sha256"]=sha(source);a["source_min_z_mm"]=minimum(z);a["source_max_z_mm"]=maximum(z)
      a["generation_vertex_z_mm"]=read(hits["events/vertex_z_mm"])[i]
     end
     h5open(joinpath(ROOT,"results/pixel_charge_1mm_4mm",sample*".h5"),"r") do other
      cache=read(other["provenance/diffusion_cache_path"],String)
      h5open(cache,"r") do d
       edges=read(d["events/$event/z_edges_mm"]);g["diffusion_z_edges_mm"]=edges
       a["diffusion_z_cell_width_mm"]=Float64(last(edges)-first(edges))/(length(edges)-1)
       a["diffusion_cache_sha256"]=sha(cache)
      end
     end
     println(sample," duration_s=",read(a["duration_s"])," buckets=",read(a["original_buckets"])," occupied=",read(a["occupied_event_buckets"])," final second fraction=",read(a["fraction_in_final_second"])," diffusion z cell mm=",read(a["diffusion_z_cell_width_mm"]))
    end
   end
  end
 end
 path
end
end
