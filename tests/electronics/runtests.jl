using Test,Random,Statistics,LinearAlgebra
include("../../src/event_analysis/EventAnalysis.jl")
include("../../src/detector_response/DetectorResponse.jl")
include("../../src/electronics/ElectronicsResponse.jl")
using .EventAnalysis,.DetectorResponse,.ElectronicsResponse
const ER=ElectronicsResponse
@testset "ADC transitions and noise" begin
 s=ADCSettings(12,0f0,1.8f0,100f0,0f0);lsb=Float64(s.maximum_V)/4096
 x=[-1.,0.,prevfloat(lsb),lsb,nextfloat(lsb),s.maximum_V,2.]
 a=digitise_adc(Xoshiro(1),x,s)
 @test a.counts==UInt16[0,0,0,1,1,4095,4095]
 @test (a.below_range,a.above_range)==(1,1)
 @test digitise_adc(Xoshiro(1),[1.],ADCSettings(16,0f0,1f0,100f0,0f0)).counts==[0xffff]
 @test_throws ArgumentError digitise_adc(Xoshiro(1),[1.],ADCSettings(17,0f0,1f0,100f0,0f0))
end
@testset "Current and linear voltage" begin
 pulse=IonCurrentPulse([-.001,0.],[1e-16,1e-16])
 t=IonReadoutTime(1000f0,20f0,2);a=IonAcquisitionSettings(10f0,5f0)
 adc=ADCSettings(12,0f0,1.8f0,100f0,0f0);at=ion_adc_time(t,pulse,a,adc)
 @test at.start_ms≈989f0
 c=ion_current(Int64[3,7],t,pulse,at)
 @test sum(diff(c.time_s).*(c.current_A[1:end-1]+c.current_A[2:end])/2)≈1e-18 rtol=1e-12
 ts=collect(0.:1e-4:.02);tau=.01;gain=1e12
 model=LinearCSAResponse(collect(0.:1e-5:.1),-gain.*exp.(-collect(0.:1e-5:.1)/tau),.9)
 current=IonCurrentPulse([0.,.02],[1e-12,1e-12])
 result=linear_csa_voltage(current,ts,model)
 expected=.9 .-gain*1e-12*tau.*(1 .-exp.(-ts/tau))
 @test maximum(abs.(result.-expected))<1e-8
end
@testset "CSA voltage noise spectrum and covariance" begin
 f=[.1,100.,1e4];asd=fill(1e-5,3)
 s=ER.voltage_noise_from_spectrum(f,asd,100f0,256)
 @test s.covariance_V2[1]≈(last(f)-first(f))*1e-10 rtol=1e-11
 exact=ER.BucketNoiseCovariance.covariance(f,asd,0.,1e-4,0.,5)
 @test maximum(abs.(s.covariance_V2[1:6]-exact))<1e-4*exact[1]
 draws=hcat([sample_csa_noise(Xoshiro(i),Int64(256),s) for i in 1:400]...)
 @test var(draws)≈exact[1] rtol=.03
 @test mean(draws[1:end-1,:].*draws[2:end,:])≈exact[2] atol=.02*exact[1]
 @test sample_csa_noise(Xoshiro(5),Int64(256),s)==sample_csa_noise(Xoshiro(5),Int64(256),s)
end
@testset "SiPM adapter matches existing calculation" begin
 channels=ReadoutChannels(Int32[0,1],Int32[0,0],Float32[0,4],Float32[0,0])
 time=ElectronReadoutTime(500f0,.25f0,Int32(20));counts=zeros(Int64,20,2);counts[3:5,1].=1000
 vectors=DetectorVectors(channels,time,counts)
 p=EventPlacement(1,75f0,75f0);frame=EventFrame(mc_vertex_frame,(-1f0,-1f0,-1f0),(1f0,1f0,1f0))
 losses=ReadoutLosses(3000,0,0,0,1,0)
 d=ElectronDetectorResponse(1,p,frame,vectors,1000,0,0,3000,0,losses,SparseCells(Int64[],Int64[],Int64[],Int64[]))
 s=SiPMReadoutSettings(SiPM(.32,100,3.6e5,12.5e-9,7e5),Int64(60),SiPMStatisticalOptions(128,16),SiPmAsic(.25e-6,.5,.001))
 response,control=sim_sipm_readout(Xoshiro(2),Xoshiro(3),d,s)
 rng=Xoshiro(2);arng=Xoshiro(3);q=ER.SiPMResponse.charge_per_pe_C(s.sipm)
 for j in 1:2
  a=ER.SiPMResponse.simulate_channel_illuminated(rng,counts[:,j],.25e-6,s.sipm,60;options=s.statistical)
  b=ER.SiPMElectronics.apply_electronics(arng,a.charge_C,q,s.asic)
  @test response.response_pe[:,j]==Float32.(b.current_A.*(.25e-6/q))
  @test response.clipped[:,j]==b.saturated
 end
 @test any(!iszero,response.response_pe[:,2])
 @test counts[3,1]==1000
end
@testset "Framed ion electronics and nonlinear selection" begin
 ch=ReadoutChannels(Int32[0,1],Int32[0,0],Float32[0,3],Float32[0,0])
 t=IonReadoutTime(1000f0,20f0,Int32(2));counts=Int64[100 0; 0 0]
 pulse=IonCurrentPulse([-.001,0.],[1e-16,1e-16]);a=IonAcquisitionSettings(10f0,5f0)
 adc=ADCSettings(12,0f0,1.8f0,100f0,0f0)
 at=ion_adc_time(t,pulse,a,adc)
 model=NonlinearCSA([0.,.8],[1e15*ER.E_CHARGE,1e15*ER.E_CHARGE],[0.,.8/.1];baseline_V=.9)
 times=collect(0.:1e-5:.1);linear=LinearCSAResponse(times,-1e15.*exp.(-times/.1),.9)
 noise=CSAVoltageNoise(100f0,zeros(nextpow(2,2at.sample_count)))
 s=IonElectronicsSettings(a,CSASettings(hybrid_nonlinear_csa,pulse,linear,model,.899f0,noise),adc)
 p=EventPlacement(1,75f0,75f0);f=EventFrame(mc_vertex_frame,(-1f0,-1f0,-1f0),(1f0,1f0,1f0))
 losses=ReadoutLosses(100,0,0,0,1,0)
 d=IonDetectorResponse(1,p,f,DetectorVectors(ch,t,counts),100,0,losses,SparseCells(Int64[],Int64[],Int64[],Int64[]))
 r,c=sim_ion_electronics_response(Xoshiro(4),Xoshiro(5),d,s)
 @test (c.linear_channels,c.nonlinear_channels)==(1,1)
 @test r.channels===ch
 @test size(r.counts)==(at.sample_count,2)
 @test length(unique(r.counts[:,2]))==1
 @test minimum(r.counts[:,1])<minimum(r.counts[:,2])
 @test sim_ion_electronics_response(Xoshiro(4),Xoshiro(5),d,s)[1].counts==r.counts
 w=ER.CSAWorkspace()
 @test ER.sim_ion_electronics_prepared(Xoshiro(4),Xoshiro(5),d,s,w)[1].counts==r.counts
 @test ER.sim_ion_electronics_prepared(Xoshiro(4),Xoshiro(5),d,s,w)[1].counts==r.counts
 @test r.counts!==ER.sim_ion_electronics_prepared(Xoshiro(4),Xoshiro(5),d,s,w)[1].counts
 @test counts==Int64[100 0;0 0]
 # Empty-pad baseline path must give exactly the same ADC codes as a zero-current CSA calculation.
 nc=fill(2e-6,length(noise.covariance_V2));nc[1]+=1e-6
 sn=IonElectronicsSettings(a,CSASettings(hybrid_nonlinear_csa,pulse,linear,model,.899f0,CSAVoltageNoise(100f0,nc)),adc)
 measured,_=sim_ion_electronics_response(Xoshiro(14),Xoshiro(15),d,sn)
 @test ER.sim_ion_electronics_prepared(Xoshiro(14),Xoshiro(15),d,sn,w)[1].counts==measured.counts
 nrng=Xoshiro(14);arng=Xoshiro(15)
 for j in 1:2
  full=ER.pad_voltage(view(counts,:,j),t,sn,at)
  nv=sample_csa_noise(nrng,at.sample_count,sn.csa.noise)
  @test measured.counts[:,j]==digitise_adc(arng,full.voltage_V.+nv,adc).counts
 end
 @test length(unique(measured.counts[:,2]))>1
 # The measured voltage accuracy follows the existing nonlinear integrator.
 wave=ER.pad_voltage(view(counts,:,1),t,s,at)
 @test maximum(abs.(wave.voltage_V-linear_csa_voltage(wave.current,wave.time_s,linear)))<1e-7
 # Selection immediately either side of the configured voltage boundary.
 minimum_linear=minimum(linear_csa_voltage(wave.current,wave.time_s,linear))
 for (guard,expected) in [(Float32(minimum_linear-1e-5),linear_csa),(Float32(minimum_linear+1e-5),hybrid_nonlinear_csa)]
  near=IonElectronicsSettings(a,CSASettings(hybrid_nonlinear_csa,pulse,linear,model,guard,noise),adc)
  @test ER.pad_voltage(view(counts,:,1),t,near,at).mode==expected
 end
 bad=NonlinearCSA([0.,.001],[1e15*ER.E_CHARGE,1e15*ER.E_CHARGE],[0.,.01];baseline_V=.9)
 sbad=IonElectronicsSettings(a,CSASettings(hybrid_nonlinear_csa,pulse,linear,bad,.899f0,noise),adc)
 @test_throws DomainError sim_ion_electronics_response(Xoshiro(1),Xoshiro(2),d,sbad)
end
@testset "Prepared noise preserves covariance and pad independence" begin
 n=Int64(128);m=4096
 # White plus correlated voltage noise, with a large slowly varying component.
 c=[1e-4+2e-5*exp(-min(k,m-k)/20)+(k==0 ? 1e-6 : 0.) for k in 0:m-1]
 s=CSAVoltageNoise(100f0,c);p=ER.prepare_csa_noise(n,s)
 @test p.period<m
 eig=2 .*p.factors.^2 ./p.period;eig[1]/=2;eig[end]/=2
 restored=ER.irfft(eig,p.period)
 @test restored[1:n]≈c[1:n] rtol=1e-12
 rng=Xoshiro(113);values=hcat([copy(ER.sample_csa_noise!(rng,p)) for _ in 1:4000]...)
 @test var(values[1,:])≈c[1] rtol=.06
 @test mean((values[2:end,:]-values[1:end-1,:]).^2)≈2(c[1]-c[2]) rtol=.015
 @test abs(cor(values[1,1:2:end],values[1,2:2:end]))<.08
 @test ER.sample_csa_noise!(Xoshiro(8),p)==sample_csa_noise(Xoshiro(8),n,s)
 # A narrow spectral line does not admit an arbitrary shorter embedding.
 tone=CSAVoltageNoise(100f0,[cospi(2*11*k/m) for k in 0:m-1]);pt=ER.prepare_csa_noise(n,tone)
 @test pt.period>nextpow(2,2(n-1))
 eig=2 .*pt.factors.^2 ./pt.period;eig[1]/=2;eig[end]/=2
 @test ER.irfft(eig,pt.period)[1:n]≈tone.covariance_V2[1:n] atol=1e-11
 @test_throws ArgumentError ER.prepare_csa_noise(n,CSAVoltageNoise(100f0,vcat(1.,fill(2.,m-1))))
end
@testset "Exact exponential integration and prepared convolution" begin
 tau=.01;gain=1e12;ts=collect(0.:1e-4:.02)
 table_t=collect(0.:1e-6:.03)
 s=LinearCSAResponse(table_t,-gain.*exp.(-table_t/tau),.9)
 # Ramp current: independent closed-form integral.
 slope=1e-10;c=IonCurrentPulse([0.,.02],[0.,slope*.02])
 expected=.9 .- gain*slope*tau.*(ts.-tau.*(1 .-exp.(-ts/tau)))
 exact=linear_csa_voltage(c,ts,s)
 @test maximum(abs.(exact.-expected))<2e-15
 coarse=ER.convolution_voltage_reference(c,ts,s;refinement=4)
 fine=ER.convolution_voltage_reference(c,ts,s;refinement=16)
 @test maximum(abs.(fine.-exact))<maximum(abs.(coarse.-exact))/10
 prep=ER.prepare_csa_convolution(ts,s)
 @test ER.convolution_voltage!(similar(ts),c,ts,s,prep)==coarse
 @test ER.convolution_voltage!(similar(ts),c,ts,s,prep)==coarse
 # A narrow triangular current pulse between ADC samples conserves its charge.
 narrow=IonCurrentPulse([0.,.000021,.000022,.000023,.02],[0.,0.,1e-9,0.,0.])
 v=linear_csa_voltage(narrow,ts,s)
 @test (.9-v[2])*exp((ts[2]-.000022)/tau)≈gain*1e-15 rtol=1e-8
 @test (.9-v[end])/(.9-v[2])≈exp(-(ts[end]-ts[2])/tau) rtol=1e-12
 @test ER.exponential_parameters(LinearCSAResponse([0.,.01,.02],[-1.,-.8,-.5],.9))===nothing
end
@testset "Prepared current matches direct sum" begin
 pulse=IonCurrentPulse([-.001,-.0003,0.],[1e-17,3e-16,8e-16])
 t=IonReadoutTime(2354f0,20f0,8)
 at=ion_adc_time(t,pulse,IonAcquisitionSettings(10f0,5f0),ADCSettings(12,0f0,1.8f0,100f0,0f0))
 counts=Int64[0,7,13,0,0,100,3,0]
 packet=ER.prepare_bucket_current(pulse,t.bucket_ms)
 current=ER.prepared_ion_current(counts,t,packet,at)
 offset=(Float64(t.start_ms)-at.start_ms)*.001;dt=Float64(t.bucket_ms)*.001
 direct=[sum((first(packet.time_s)<=x-offset-(j-1)*dt<=last(packet.time_s) ?
 counts[j]*ER.interpolate(packet.time_s,packet.current_A,x-offset-(j-1)*dt) : 0.) for j in eachindex(counts)) for x in current.time_s]
 @test current.current_A==direct
end
