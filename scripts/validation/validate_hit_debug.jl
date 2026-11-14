include("../production/run_hit_production.jl")
using Test
@testset "Control settings" begin
 @test hit_control_settings(Dict()).mode=="none"
 @test_throws ArgumentError hit_control_settings(Dict("controls"=>Dict("mode"=>"typo")))
 @test_throws ArgumentError hit_control_settings(Dict("controls"=>Dict("mode"=>"debug")))
 @test_throws ArgumentError hit_control_settings(Dict("controls"=>Dict("mode"=>"debug","event_limit"=>1,"diagnostics"=>["bad"])))
end
card=TOML.parsefile(ARGS[1]);card["accepted_events"]=3
card["output_subdir"]=joinpath(ROOT,"work/hit_output_debug_validation")
card["controls"]=Dict("mode"=>"debug","event_limit"=>1,"diagnostics"=>["noise","hit_cut","graph_scan","electronics_examples"])
path=joinpath(ROOT,"work/readout_timing/hit_debug.toml");open(io->TOML.print(io,card),path,"w")
result=run_hit_production(path)
@testset "Only selected debug events; identical production values" begin
 r=open_hit_input(result);reference=open_hit_input(abspath(ARGS[2]))
 for id in keys(r.event_rows)
  e=read_hit_event(r,id);ref=read_hit_event(reference,id)
  for arm in (:itrack,:etrack),key in fieldnames(SavedHits)
   @test isequal(getfield(getfield(e,arm),key),getfield(getfield(ref,arm),key))
  end
  for key in fieldnames(EventEnergy);@test isequal(getfield(e.energy,key),getfield(ref.energy,key));end
 end
 ids=sort(collect(keys(r.event_rows)));close_hit_input(r);close_hit_input(reference)
 h5open(joinpath(dirname(result),"debug/controls.h5"),"r") do f
  @test Set(keys(f))==Set(["card",string(first(ids))])
  @test Set(keys(f[string(first(ids))]))==Set(card["controls"]["diagnostics"])
 end
end
