using Test
include("../production/run_framed_detector_response.jl")
include("../../src/electronics/ElectronicsResponse.jl")
using .ElectronicsResponse
include("../production/framed_electronics_settings.jl")
function validate_electronics(card)
 prepared_run=prepare_electronics_run(card)
 s=electronics_settings(prepared_run,IonReadoutTime(1000f0,20f0,Int32(25)))
 @testset "Run preparation and event lengths" begin
  shorter=electronics_settings(prepared_run,IonReadoutTime(1000f0,20f0,Int32(24)))
  @test shorter.ion===s.ion
  longer=electronics_settings(prepared_run,IonReadoutTime(1000f0,20f0,Int32(50)))
  @test longer.ion!==s.ion
  @test longer.ion.csa.current_pulse===s.ion.csa.current_pulse
  @test ER.exponential_parameters(s.ion.csa.linear_response)!==nothing
 end
 results=Dict{String,Any}()
 @testset "Existing CSA spectrum and 100 us samples" begin
  exact=ER.BucketNoiseCovariance.covariance(s.frequency,s.asd,0.,1e-4,0.,8)
  err=maximum(abs.(exact-s.ion.csa.noise.covariance_V2[1:9]))/exact[1]
  @test err<1e-5;results["noise_covariance_relative_error"]=err
  doubled=ER.voltage_noise_from_spectrum(s.frequency,s.asd,100f0,Int64(5200);oversampling=32)
  @test maximum(abs.(s.ion.csa.noise.covariance_V2[1:5000]-doubled.covariance_V2[1:5000]))/exact[1]<1e-4
  results["voltage_noise_rms_V"]=sqrt(exact[1])
  prepared=ER.prepare_csa_noise(Int64(5168),s.ion.csa.noise)
  eig=2 .*prepared.factors.^2 ./prepared.period;eig[1]/=2;eig[end]/=2
  reconstructed=ER.irfft(eig,prepared.period)
  error=maximum(abs.(reconstructed[1:5168]-s.ion.csa.noise.covariance_V2[1:5168]))/exact[1]
  @test error<1e-11
  results["prepared_noise_covariance_relative_error"]=error
  results["prepared_noise_period"]=prepared.period
  results["reference_noise_period"]=length(s.ion.csa.noise.covariance_V2)
 end
 @testset "Saved linear SPICE voltage" begin
  residuals=Float64[];coarse_errors=Float64[];fine_errors=Float64[]
  h5open(s.rawpath,"r") do h
   for n in [10,100],frame in [20,30,40]
    time=IonReadoutTime(0f0,Float32(frame),Int32(1));at=ion_adc_time(time,s.ion.csa.current_pulse,s.ion.acquisition,s.ion.adc)
    w=ER.pad_voltage(Int64[n],time,s.ion,at)
    coarse=ER.convolution_voltage_reference(w.current,w.time_s,s.ion.csa.linear_response;refinement=4)
    fine=ER.convolution_voltage_reference(w.current,w.time_s,s.ion.csa.linear_response;refinement=16)
    exact=linear_csa_voltage(w.current,w.time_s,s.ion.csa.linear_response)
    push!(coarse_errors,maximum(abs.(coarse-exact)));push!(fine_errors,maximum(abs.(fine-exact)))
    @test maximum(abs.(fine-exact))<maximum(abs.(coarse-exact))
    data=h["$(frame)ms/single_N$(n)_r2"];t=read(data["time"]);v=read(data["V(out)"])
    offset=Float64(time.start_ms-at.start_ms)*.001
    ref=[ER.interpolate(t,v,x-offset+.01) for x in w.time_s]
    err=maximum(abs.(w.voltage_V-ref));push!(residuals,err);@test err<1e-5
   end
  end
  results["linear_maximum_voltage_error_V"]=maximum(residuals)
  results["convolution_25us_error_V"]=maximum(coarse_errors)
  results["convolution_6p25us_error_V"]=maximum(fine_errors)
 end
 @testset "Saved nonlinear SPICE voltage at 100 us" begin
  paths,mins,frames,valid,hashes=h5open(joinpath(ROOT,card["nonlinear_calibration"]),"r") do f
   g=f["validation"];(read(g["raw_paths"]),read(g["minv"]),read(g["frame_ms"]),read(g["valid"]),read(g["raw_sha256"]))
  end
  selected=Int[]
  for frame in [20,30,40]
   indices=sort(findall(j->frames[j]==frame&&valid[j]==1&&isfile(paths[j]),eachindex(paths));by=j->mins[j])
   append!(selected,indices[unique(round.(Int,range(1,length(indices);length=min(8,length(indices)))))])
  end
  @test length(selected)>=12
  errs=Float64[]
  for j in selected
   @test bytes2hex(open(sha256,paths[j]))==hashes[j]
   folder=dirname(paths[j]);t=Float64[];i=Float64[]
   for line in eachline(joinpath(folder,"current.pwl"));a=split(line);push!(t,parse(Float64,a[1]));push!(i,parse(Float64,a[2]));end
   ts=collect(first(t):.0001:last(t));out=nonlinear_waveform(s.nonlinear.model,t,i;sample_time_s=ts)
   @test out.valid
   ref=h5open(paths[j],"r") do h
    rt=read(h["time"]);rv=read(h["V(out)"]);[ER.interpolate(rt,rv,x) for x in ts]
   end
   err=maximum(abs.(ref.-out.voltage_V));push!(errs,err)
   # Existing 0.5 ENC tolerance converted to voltage using the small-signal calibration.
   @test err<.5*50.179501626864536*s.gain*ER.E_CHARGE
  end
  results["nonlinear_cases"]=length(errs);results["nonlinear_maximum_voltage_error_V"]=maximum(errs)
 end
 results
end
if abspath(PROGRAM_FILE)==@__FILE__
 c=TOML.parsefile(joinpath(ROOT,"config/readout/bb0nu_electronics_pilot.toml"))
 out=validate_electronics(c);mkpath(joinpath(ROOT,"work/framed_electronics"))
 open(joinpath(ROOT,"work/framed_electronics/validation.toml"),"w") do io;TOML.print(io,out);end
end
