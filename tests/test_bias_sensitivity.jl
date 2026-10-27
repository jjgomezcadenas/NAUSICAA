using Test,HDF5,TOML
root=normpath(joinpath(@__DIR__,".."));c=TOML.parsefile(joinpath(root,"config/bias_sensitivity.toml"))
@testset "Bias sensitivity: charge, refinement and selected working ranges" begin
 h5open(joinpath(root,c["output"],"analysis.h5"),"r") do h
  @test length(keys(h))==40
  for name in filter(!=("ranges"),collect(keys(h)))
   a=attributes(h[name]);family=read(a["family"]);kind=read(a["kind"])
   @test !haskey(a,"analysis_failure")
   @test read(a["ENC_e_rms"])>0
   @test abs(read(a["charge_relative_error"]))<1e-8
   @test 0<read(a["minimum_output_1000_V"])<read(a["maximum_output_1000_V"])<1.8
   if kind=="m4";@test read(a["m3_gate_V"])==c[family]["m3_gate_V"];end
   if kind=="m3";@test read(a["m4_gate_V"])==c[family]["m4_gate_V"];end
   if haskey(a,"signal_refinement_error")
    @test read(a["signal_refinement_error"])<1e-3
    @test read(a["noise_refinement_error"])<1e-3
    @test read(a["quadrature_error"])<1e-5
   end
  end
  for name in ("paper_rounded_4","paper_rounded_5","v2_rounded_3","v2_rounded_4",
               "paper_m3_5","paper_m3_7","paper_m4_3","paper_m4_9","v2_m4_3","v2_m4_9")
   a=attributes(h[name])
   @test read(a["ENC_relative_to_nominal"])<1.05
   @test abs(read(a["gain1000_relative_change"]))<.01
  end
 end
end
