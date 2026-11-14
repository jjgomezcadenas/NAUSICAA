using Test
include("../production/run_framed_reconstruction.jl")
"""Validate the electron mean inverse against the existing statistical simulation.
Use steady and changing illumination over the pilot range, including full-device
dark activity and ASIC input noise. Report per-bucket bias and sampling error.
This validation does not change the SiPM simulation or its settings."""
function validate_electron_inverse()
 ec=TOML.parsefile(joinpath(ROOT,"config/readout/bb0nu_electronics_pilot_100.toml"))
 s=prepare_electronics_run(ec).sipm;sensor=s.sipm
 c=ElectronReconstructionSettings(sensor,s.illuminated_microcells,true)
 patterns=[("dark",zeros(Int64,8)),("weak",fill(Int64(100),8)),("steady",fill(Int64(30000),8)),
  ("rising_falling",Int64[0,0,3000,30000,100000,30000,3000,0]),("isolated",Int64[0,0,0,100000,0,0,0,0])]
 results=Dict{String,Any}();rng=Xoshiro(91425);arng=Xoshiro(42195);trials=500
 ch=ReadoutChannels(Int32[0],Int32[0],Float32[0],Float32[0]);time=ElectronReadoutTime(0f0,.25f0,Int32(8))
 q0=ER.SiPMResponse.charge_per_pe_C(sensor)
 @testset "Electron inverse: weak, steady and changing illumination" begin
  for (name,photons) in patterns
   samples=zeros(8,trials)
   for repeat in 1:trials
    simulated=ER.SiPMResponse.simulate_channel_illuminated(rng,photons,.25e-6,sensor,s.illuminated_microcells;options=s.statistical)
    # Unclipped charge isolates the sensor inverse; ASIC clipping is tested separately.
    charge=Float32.(simulated.charge_C./q0 .+s.asic.noise_rms_pe.*randn(arng,8))
    measured=SiPMASICResponse(1,ch,time,reshape(charge,8,1),falses(8,1))
    out=reco_electron_signal(measured,c);samples[:,repeat]=out.pe[:,1]
   end
   expected=sensor.pde.*photons;means=vec(mean(samples;dims=2));sigma=vec(std(samples;dims=2));error=sigma./sqrt(trials)
   tolerance=max.(.5,0.005 .*expected,3 .*error)
   @test all(abs.(means-expected).<=tolerance)
   results[name]=Dict("photons"=>photons,"expected_pe"=>expected,"mean_pe"=>means,"std_pe"=>sigma,"bias_pe"=>means-expected,"standard_error_pe"=>error,"tolerance_pe"=>tolerance)
  end
 end
 results["trials"]=trials
 results["criterion"]="absolute bias <= max(0.5 PE, 0.5% of expected PE, 3 standard errors)"
 open(joinpath(ROOT,"work/reconstruction/electron_validation.toml"),"w") do io;TOML.print(io,results);end
end
validate_electron_inverse()
