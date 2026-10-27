using Test,HDF5
include(joinpath(@__DIR__,"../src/ion_signal/discharge_calibration.jl"))
const D=DischargeCalibration
@testset "Decay fit and corrected readout" begin
 t=collect(0.:1e-5:.2);tau=2.7;v=-0.003 .* exp.(-t./tau)
 fitted,err,_,_,_=D.fit_tail(t,v,.021,.12,.0002)
 @test fitted≈tau rtol=1e-9
 @test err<1e-10
 @test D.average(t,3 .* t .+ 2,.05,.0003)≈2.15 atol=1e-12
 for center in (.04,.06,.1)
  now=D.average(t,v,center,.0003);previous=D.average(t,v,center-.02,.0003)
  @test abs(now-exp(-.02/tau)*previous)<1e-12
 end
 # Non-exponential decay must show a detectable model residual.
 w=-0.003 .* (0.7 .* exp.(-t./2.7).+0.3 .* exp.(-t./.05))
 _,mismatch,_,_,_=D.fit_tail(t,w,.021,.12,.0002)
 @test mismatch>.005
end
@testset "New current inputs conserve charge" begin
 root=normpath(joinpath(@__DIR__,".."));cfg,_=D.configuration(root)
 h5open(joinpath(root,cfg["output"],"input.h5"),"r") do h
  q=read(attributes(h)["single_ion_charge_C"])
  for r in (1,2)
   t=read(h["r$r/time_s"])
   @test D.B.integral(t,read(h["r$r/single_A_per_ion"]))≈q rtol=1e-10
   @test D.B.integral(t,read(h["r$r/train_A"]))≈sum(cfg["train_counts"])*q rtol=1e-10
  end
 end
end
@testset "Calibrated result and independent late tail" begin
 root=normpath(joinpath(@__DIR__,".."));cfg,_=D.configuration(root);dir=joinpath(root,cfg["output"])
 h5open(joinpath(dir,"analysis.h5"),"r") do h
  a=attributes(h);tau=read(a["tau_s"]);ts=read(h["tail_time_s"]);fit=read(h["tail_fit_V"])
  h5open(joinpath(dir,"raw.h5"),"r") do raw
   t,v=D.waveform(raw,cfg["calibration_ions"],2);t.-=cfg["prehistory_s"]
   times=collect(.14:.001:.199)
   actual=-[D.B.interp(t,v,x) for x in times]
   predicted=fit[end].*exp.(-(times.-ts[end])./tau)
   @test maximum(abs.(predicted./actual.-1))<1e-5
  end
  g=h["calibrated"];n=read(attributes(g)["verified_N_snr5"]);counts=read(g["counts"]);snrs=read(g["snr"])
  @test snrs[findfirst(==(n-1),counts)]<5<snrs[findfirst(==(n),counts)]
  @test abs(read(attributes(g)["previous_frame_residual_fraction"]))<.001
  @test maximum(abs.(read(g["train_error_ions"])))<.12*read(attributes(g)["ENC_e_rms"])
  @test maximum(read(g["train_refinement_ions"]))<.05
  @test maximum(abs.(read(h["local_tau_s"])./tau.-1))<1e-4
 end
end
