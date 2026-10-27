using Test,HDF5
root=normpath(joinpath(@__DIR__,".."))
@testset "v3 actual pulse charge and amplitude" begin
 h5open(joinpath(root,"results/ltspice/rogelio_v3/analysis.h5"),"r") do h
  a=attributes(h);ex=attributes(h["exp_r2"]);un=attributes(h["uniform_1000_r2"])
  @test abs(read(a["exp_area_relative_error"]))<1e-5
  @test isapprox(read(ex["fitted_current_tau_s"]),1.7e-3;rtol=1e-8)
  @test .2<read(ex["first_positive_current_s"])<.200000001
  @test isapprox(read(ex["current_max_A"]),48e-15;rtol=1e-8)
  @test 500<read(ex["equivalent_electrons"])<520
  @test .013<read(ex["peak_voltage_V"])<.015
  @test .026<read(un["peak_voltage_V"])<.029
  for att in (ex,un)
   @test read(att["charge_refinement_error"])<1e-5
   @test read(att["voltage_refinement_error"])<1e-4
   @test isapprox(read(att["charge_C"]),read(att["behavioral_charge_C"]);rtol=1e-5)
  end
  uncomp=attributes(h["original_uncompressed"])
  @test isapprox(read(uncomp["charge_C"]),read(uncomp["behavioral_charge_C"]);rtol=1e-6)
  @test 1.9<read(a["uniform_to_exp_charge_ratio"])<2.1
  @test isapprox(read(a["uniform_to_exp_voltage_ratio"]),read(a["uniform_to_exp_charge_ratio"]);rtol=.01)
 end
end
