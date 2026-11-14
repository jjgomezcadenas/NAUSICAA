using Test,Random,Statistics,LinearAlgebra
include("../../src/event_analysis/EventAnalysis.jl")
include("../../src/detector_response/DetectorResponse.jl")
include("../../src/electronics/ElectronicsResponse.jl")
include("../../src/reconstruction/SignalReconstruction.jl")
using .EventAnalysis,.DetectorResponse,.ElectronicsResponse,.SignalReconstruction
const ER=ElectronicsResponse;const SR=SignalReconstruction
BLAS.set_num_threads(1)
@testset "Covariance of averaged and baseline-subtracted samples" begin
 c=[1e-4+1e-5*exp(-k/5)+(k==0 ? 2e-6 : 0.) for k in 0:99]
 groups=[1:10,11:20,21:30,31:37];quant=3e-8
 actual=SR.averaged_covariance(c,groups,quant)
 N=37;C=[c[abs(i-j)+1]+(i==j ? quant : 0.) for i in 1:N,j in 1:N]
 A=zeros(3,N)
 for i in 1:3;A[i,groups[1]].=1/10;A[i,groups[i+1]].=-1/length(groups[i+1]);end
 @test actual≈A*C*A' rtol=1e-11
 @test isposdef(Symmetric(actual))
end
function fixture()
 ch=ReadoutChannels(Int32[0],Int32[0],Float32[0],Float32[0])
 t=IonReadoutTime(1000f0,20f0,Int32(4));pulse=IonCurrentPulse([-.001,0.],[1e-16,1e-16])
 adc=ADCSettings(16,0f0,1.8f0,100f0,0f0);at=ion_adc_time(t,pulse,IonAcquisitionSettings(10f0,5f0),adc)
 times=collect(0.:1e-5:.2);gain=1e15;tau=.1
 linear=LinearCSAResponse(times,-gain.*exp.(-times/tau),.9)
 nonlinear=NonlinearCSA([0.,.2,.5,.8],[.00016,.00015,.00012,.0001],[0.,2.,5.,8.];baseline_V=.9)
 c=[1e-6*exp(-min(k,8192-k)/20)+(k==0 ? 1e-7 : 0.) for k in 0:8191]
 cal=IonReconstructionCalibration(adc,pulse,linear,nonlinear,CSAVoltageNoise(100f0,c))
 settings=IonReconstructionSettings(10f0,5f0,.63f0,20,1f-4)
 ch,t,at,cal,settings
end
@testset "Linear reconstruction and ADC baseline" begin
 ch,t,at,c,s=fixture();p=SR.prepare_ion_reconstruction(at,t,c,s)
 for q in ([1.,0.,0.,0.],[0.,3.,7.,0.],[0.,0.,0.,10.],[-2.,3.,-1.,5.])
  @test p.inverse*(p.response*q)≈q atol=1e-10
 end
 @test p.ion_covariance≈p.inverse*p.covariance*p.inverse' rtol=1e-12
 # Shift all ADC samples by a constant offset. Difference must cancel exactly.
 codes=fill(UInt16(32000),at.sample_count,1);r=IonADCResponse(1,ch,at,codes)
 result=reco_ion_signal(r,t,c,s)
 shifted=IonADCResponse(1,ch,at,codes.+UInt16(100))
 @test reco_ion_signal(shifted,t,c,s).ions==result.ions
 @test all(iszero,result.ions)
 @test all(result.valid)
 @test all(isfinite,adc_voltages(r,c.adc))
 bad=copy(codes);bad[3]=0
 @test !any(reco_ion_signal(IonADCResponse(1,ch,at,bad),t,c,s).valid)
 # Noise passed through the inverse has the predicted covariance.
 rng=Xoshiro(4);draws=p.inverse*(cholesky(Symmetric(p.covariance)).L*randn(rng,size(p.response,1),10000))
 @test var(vec(sum(draws;dims=1)))≈sum(p.ion_covariance) rtol=.04
end
@testset "Nonlinear inverse, derivatives and calibrated forward model" begin
 ch,t,at,c,s=fixture();p=SR.prepare_ion_reconstruction(at,t,c,s)
 q=[800.,1200.,400.,200.];prediction=SR.nonlinear_prediction(p,q)
 @test prediction.valid
 # Independent simulator integrates only physical currents, at its original knots.
 current=ion_current(Int64.(q),t,c.current_pulse,at)
 ts=p.times;tt=vcat(current.time_s,max(last(current.time_s),last(ts))+1e-12)
 forward=nonlinear_waveform(c.nonlinear_response,tt,vcat(current.current_A,0.);sample_time_s=ts)
 @test forward.valid
 expected=[.9-mean(view(forward.voltage_V,g)) for g in p.groups[2:end]]
 @test maximum(abs.(prediction.voltage-expected))<1e-8
 for j in 1:4
  perturbed=copy(q);perturbed[j]+=.01
  finite=(SR.nonlinear_prediction(p,perturbed).voltage-prediction.voltage)/.01
  @test finite≈prediction.jacobian[:,j] rtol=1e-3 atol=1e-8
 end
 fit=SR.fit_nonlinear(p,expected,p.inverse*expected)
 @test fit.valid
 @test maximum(abs.(fit.ions-q))<.01
 @test !SR.nonlinear_prediction(p,fill(1e8,4)).valid
end
@testset "Electron correction and charge-preserving rebin" begin
 ch,t,at,c,s=fixture();et=ElectronReadoutTime(123f0,.25f0,Int32(43))
 sensor=SiPM(.32,1000,3e5,12.5e-9,0.)
 settings=ElectronReconstructionSettings(sensor,Int64(500),true)
 input=Float32.(reshape(collect(-2.:40.),43,1));clipped=falses(43,1)
 response=SiPMASICResponse(1,ch,et,input,clipped)
 fine=reco_electron_signal(response,settings)
 @test fine.pe[1]<0
 @test all(fine.valid)
 coarse=rebin_electron_signal(fine,Int64(20))
 @test size(coarse.pe)==(3,1)
 @test sum(Float64,coarse.pe)≈sum(Float64,fine.pe) rtol=1e-7
 @test coarse.time.start_us==123f0
 @test coarse.time.bucket_us==5f0
 @test coarse.acquired_end_us==fine.acquired_end_us==133.75f0
 clipped[21]=true;invalid=rebin_electron_signal(reco_electron_signal(response,settings),Int64(20))
 @test invalid.asic_clipped[2]&&!invalid.valid[2]&&isnan(invalid.pe[2])
 @test invalid.valid[1]&&invalid.valid[3]
 @test_throws ArgumentError rebin_electron_signal(fine,Int64(0))
end
@testset "Linear charge recovered from digitised physical voltage" begin
 ch,t,at,c,s=fixture();q=Int64[50,100,200,50]
 current=ion_current(q,t,c.current_pulse,at)
 ts=collect(0:at.sample_count-1).*(Float64(at.sample_us)*1e-6)
 current=IonCurrentPulse(vcat(current.time_s,max(last(current.time_s),last(ts))+1e-12),vcat(current.current_A,0.))
 volts=linear_csa_voltage(current,ts,c.linear_response)
 codes=digitise_adc(Xoshiro(15),volts,c.adc).counts
 reconstructed=reco_ion_signal(IonADCResponse(1,ch,at,reshape(codes,:,1)),t,c,s)
 @test all(reconstructed.valid)
 @test maximum(abs.(reconstructed.ions[:,1]-q))<1.
 @test !any(reconstructed.checks.nonlinear)
end
include("../../scripts/production/reconstruction_controls.jl")
@testset "ZS control denominators and empty-channel labels" begin
 input=Float64[10 0 0;0 5 0;0 0 0]
 reco=Float32[3 -1 .2;-2 NaN .3;1 2 -1];valid=isfinite.(reco)
 result=zs_controls(reco,valid,input,[0.,5.],10)
 @test result["input_total"]==15
 @test result["input_invalid"]==5
 @test result["input_below"]==[0.,10.]
 @test result["noise_cells"]==7
 @test result["quiet_cells"]==3
 @test result["quiet_above"]==[2,0]
 @test sum(result["histograms"]["all"]["counts"])==8
 @test result["histograms"]["quiet_channels"]["entries"]==3
end
