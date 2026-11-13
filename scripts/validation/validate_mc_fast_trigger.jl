# Compare the reviewed library with the saved full-sample 5 keV prototype.

using HDF5,TOML,Test
include(joinpath(@__DIR__,"../../src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis
root=normpath(joinpath(@__DIR__,"../.."));c=TOML.parsefile(joinpath(root,"prototypes/mc_fast_trigger/fast_trigger_gamma_outside_frame.toml"))
s=FastTriggerSettings(1500f0,50f0,1450f0,(10f0,10f0,10f0),5f0)
h5open(joinpath(root,"work/fast_trigger/reference_5keV.h5"),"r") do ref
 @testset "Full-sample prototype decisions" begin
  for species in c["species"]
   g=ref[species];ids=read(g["all_events/event_number"]);vertex=read(g["all_events/vertex_z_mm"])
   contained=read(g["all_events/contained"]);fiducial=read(g["all_events/accepted"]);survivors=Set(read(g["after_gamma_cut_events/event_number"]))
   particles=read_mc_particles(expanduser(c[species]["particles"]),ids)
   stats=FastTriggerStatistics()
   h5open(joinpath(root,c[species]["hits"]),"r") do f
    for (i,id) in enumerate(ids)
     h,info=read_mc_event(f,id;event_row=i);v=center_hits(h,info)
     p=EventPlacement(id,vertex[i]/10f0,(1500f0-vertex[i])/10f0)
     o=mc_fast_trigger(v,info,particles[id],p,s);stats=update_statistics(stats,o)
     expected=contained[i]==0 ? outside_detector : fiducial[i]==0 ? outside_vertex_fiducial : id in survivors ? trigger_accepted : outside_frame_gamma
     @test o.decision==expected
     @test h.energy_keV==f["hits/energy_MeV"][f["events/first_row"][i]:f["events/first_row"][i]+f["events/n_hits"][i]-1].*1000f0
     if haskey(g["event_details"],string(id))
      d=g["event_details/$(id)"];lo=read(d["primary_low_mm"]);hi=read(d["primary_high_mm"]);lo[3]-=vertex[i];hi[3]-=vertex[i]
      @test all(isapprox.(collect(o.frame.lower_mm),lo.-10f0;atol=2f-4))
      @test all(isapprox.(collect(o.frame.upper_mm),hi.+10f0;atol=2f-4))
     end
    end
   end
   @test stats.accepted_events==read(g["selection/accepted_gamma_cut"])
   println(species,": ",stats.input_events," decisions matched; ",stats.accepted_events," accepted")
  end
 end
end
