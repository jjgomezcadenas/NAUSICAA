include("runtests.jl")
@testset "Direct nonlinear inverse, noiseless voltages" begin
 ch,t,at,c,s=fixture();p=SR.prepare_ion_reconstruction(at,t,c,s)
 d=SR.prepare_direct_inverse!(p)
 @test d===SR.prepare_direct_inverse!(p)
 # F integrates reciprocal gain across table boundaries.
 for x in [.01,.2,.37,.5,.72]
  epsx=1e-7
  derivative=(SR.inverse_terms(d,x+epsx)[1]-SR.inverse_terms(d,x-epsx)[1])/(2epsx)
  gain=ER.interpolate(c.nonlinear_response.knots_V,c.nonlinear_response.gain_V_per_e,x)
  @test derivative≈1/gain rtol=1e-7
 end
 @test SR.reciprocal_gain_integral(.2,.001,0.)≈200.
 for q in [Int64[800,1200,400,200],Int64[0,0,2000,0],Int64[1800,0,0,1800]]
  current=ion_current(q,t,c.current_pulse,at)
  times=vcat(current.time_s,max(last(current.time_s),last(p.times))+1e-12)
  forward=nonlinear_waveform(c.nonlinear_response,times,vcat(current.current_A,0.);sample_time_s=p.times)
  @test forward.valid
  direct=SR.direct_ion_inverse(p,c.nonlinear_response.baseline_V.-forward.voltage_V)
  @test direct.valid && direct.iterations==0
  @test maximum(abs.(direct.ions-q))<.2
  @test abs(sum(direct.ions)-sum(q))<.2
 end
 @test !SR.direct_ion_inverse(p,fill(last(c.nonlinear_response.knots_V)+.01,length(p.times))).valid
 # Constant gain and linear recovery reduce to integrated linear deconvolution.
 gain=c.linear_response.voltage_per_charge_V_C[1]*-ER.E_CHARGE
 model=NonlinearCSA([0.,.8],[gain,gain],[0.,8.];baseline_V=.9)
 cal=IonReconstructionCalibration(c.adc,c.current_pulse,c.linear_response,model,c.voltage_noise)
 lp=SR.prepare_ion_reconstruction(at,t,cal,s);q=Int64[500,100,300,800]
 current=ion_current(q,t,c.current_pulse,at)
 current=IonCurrentPulse(vcat(current.time_s,max(last(current.time_s),last(lp.times))+1e-12),vcat(current.current_A,0.))
 voltage=linear_csa_voltage(current,lp.times,cal.linear_response)
 @test maximum(abs.(SR.direct_ion_inverse(lp,0.9 .-voltage).ions-q))<.2
end
