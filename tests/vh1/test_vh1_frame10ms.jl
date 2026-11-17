# VH1 with 10 ms frames: the production 20 fF card with frame_ms = 10, two
# conversions of 5 ms. The circuit, gain and noise spectrum are those of the
# 20 ms card; only the frame integration changes. The existing VH1 tests are
# unchanged; this file checks the 10 ms card against them.
using Test,Random,Statistics,LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"))
include(joinpath(ROOT,"src/detector_response/DetectorResponse.jl"))
include(joinpath(ROOT,"src/electronics/ElectronicsResponse.jl"))
include(joinpath(ROOT,"src/reconstruction/SignalReconstruction.jl"))
include(joinpath(ROOT,"src/vh1/VH1Readout.jl"))
using .EventAnalysis,.DetectorResponse,.ElectronicsResponse,.SignalReconstruction,.VH1Readout
const ER=ElectronicsResponse;const V=VH1Readout

s10=load_vh1_settings(ROOT,joinpath(ROOT,"config/specs/response/ion_csa/vh1_device_changes_cf20f_10ms.toml");minimum_samples=4096)
s20=load_vh1_settings(ROOT,joinpath(ROOT,"config/specs/response/ion_csa/vh1_device_changes_cf20f.toml");minimum_samples=4096)
noiseless(s)=VH1Settings(s.gain_V_per_ion,s.charge_fraction,s.feedback_capacitance_F,s.v_ref_V,s.window_V,s.frame_ms,s.conversions_per_frame,s.window_ms,s.reference_frames,
    ADCSettings(16,0f0,1.8f0,s.adc.sample_us,0f0),s.noise_sample_us,CSAVoltageNoise(s.noise_sample_us,zeros(length(s.noise.covariance_V2))),s.current_pulse,s.provenance)
withwindow(s,w;k=1)=VH1Settings(s.gain_V_per_ion,s.charge_fraction,s.feedback_capacitance_F,s.v_ref_V,s.window_V,s.frame_ms,k,Float32(w),s.reference_frames,s.adc,s.noise_sample_us,s.noise,s.current_pulse,s.provenance)
channels(n)=ReadoutChannels(Int32.(1:n),ones(Int32,n),Float32.(1:n),zeros(Float32,n))
placement=EventPlacement(1,50f0,100f0)
frame=EventFrame(mc_vertex_frame,(0f0,0f0,0f0),(1f0,1f0,1f0))
losses=ReadoutLosses(0,0,0,0,0,0)
detector(counts::Matrix{Int64},t::IonReadoutTime)=IonDetectorResponse(1,placement,frame,DetectorVectors{IonReadoutTime}(channels(size(counts,2)),t,counts),sum(counts),0,losses,SparseCells(Int64[],Int64[],Int64[],Int64[]))

@testset "10 ms card: same circuit, shorter frame" begin
    @test s10.frame_ms==10f0 && s10.conversions_per_frame==2 && s10.window_ms==5f0 && s10.adc.bits==16
    @test V.samples_per_frame(s10)==100 && V.window_samples(s10)==50
    @test s10.gain_V_per_ion==s20.gain_V_per_ion && s10.charge_fraction==s20.charge_fraction && s10.window_V==s20.window_V
    @test s10.noise.covariance_V2==s20.noise.covariance_V2
    @test s10.provenance["noise_spectrum_sha256"]==s20.provenance["noise_spectrum_sha256"]
    @test s10.adc.sample_us==10000f0 && s20.adc.sample_us==20000f0
end

@testset "Exact recovery with 10 ms buckets" begin
    t=IonReadoutTime(1000f0,10f0,5);counts=reshape(Int64[0,120,3000,7,0],5,1)
    for (w,k) in ((5.0,2),(10.0,1),(2.5,4))
        sq=withwindow(noiseless(s10),w;k=k);ft=vh1_frame_time(t,sq)
        @test ft.frame_ms==10f0 && ft.start_ms+ft.frame_ms*ft.frame_count>=1000+50
        d=detector(counts,t);r=sim_vh1_electronics(Xoshiro(1),Xoshiro(2),d,sq)
        @test !r.saturated[1]
        plan=prepare_vh1_reconstruction(ft,t,sq)
        grid=V.internal_grid(ft,sq);packet=ER.prepare_bucket_current(sq.current_pulse,t.bucket_ms)
        current=ER.prepared_ion_current(vec(counts),t,packet,grid)
        reading=V.frame_readings!(zeros(V.readings(ft,sq)),current,ft,sq)
        i0=sq.conversions_per_frame
        q=plan.inverse*[reading[i0]-reading[i] for i in i0+1:length(reading)]
        @test maximum(abs.(q.-vec(counts)))<1e-7
        ri=reco_vh1_ions(r,plan);lsb_ions=(1.8/2^16)/sq.gain_V_per_ion
        @test maximum(abs.(vec(ri.ions).-vec(counts)))<2lsb_ions*sqrt(ft.frame_count) && all(ri.valid)
    end
end

@testset "Charge conservation: one ion in a 10 ms bucket" begin
    t=IonReadoutTime(1000f0,10f0,3);sq=withwindow(noiseless(s10),10.0;k=1);ft=vh1_frame_time(t,sq)
    A=prepare_vh1_reconstruction(ft,t,sq).response./sq.gain_V_per_ion
    inc=vcat(A[1:1,:],diff(A;dims=1))
    pulse=sq.current_pulse
    finite=sum(diff(pulse.time_s).*(pulse.current_A[1:end-1]+pulse.current_A[2:end])/2)/ER.E_CHARGE
    @test abs(sum(inc[:,1])-finite)<1e-9 && abs(sum(inc[:,2])-finite)<1e-9
    # the 1.67 ms pulse of an ion in bucket k reaches its own frame, the next and a trace in the one after
    for k in 1:3
        col=inc[:,k]
        @test col[k]+col[k+1]>.95*finite && all(abs(col[n])<1e-9 for n in eachindex(col) if n<k || n>k+2)
    end
    println("VH1 10 ms increments per frame for one ion in bucket 1: ",round.(inc[:,1];digits=4))
end

@testset "Response matrix and predicted noise: 10 ms against 20 ms" begin
    results=Dict{Float32,Any}()
    for (s,b) in ((s10,10f0),(s20,20f0))
        t=IonReadoutTime(1000f0,b,30);ft=vh1_frame_time(t,s);plan=prepare_vh1_reconstruction(ft,t,s)
        W=plan.lower\plan.response;sig=sqrt.(diag(plan.ion_covariance))
        results[b]=(cond_response=cond(plan.response),cond_whitened=cond(W),noise=sig,readings=size(plan.response,1))
        @test size(plan.response)==(2ft.frame_count-2,30) && isfinite(cond(W)) && all(isfinite,sig) && all(>(0),sig)
        println("frame $b ms: readings $(size(plan.response,1)), cond(response) $(round(cond(plan.response);digits=1)), cond(whitened) $(round(cond(W);digits=1)), ions per bucket median $(round(median(sig);digits=2)) max $(round(maximum(sig);digits=2))")
    end
    r10=results[10f0];r20=results[20f0]
    @test r10.readings==r20.readings   # 30 buckets, one reference frame: the same number of readings
    @test .5<median(r10.noise)/median(r20.noise)<2
    println("predicted ions per bucket, 10 ms / 20 ms: ",round(median(r10.noise)/median(r20.noise);digits=3))
end

@testset "Noise-only pads reproduce the 10 ms covariance" begin
    t=IonReadoutTime(1000f0,10f0,6);ft=vh1_frame_time(t,s10)
    npads=1500;counts=zeros(Int64,6,npads)
    r=sim_vh1_electronics(Xoshiro(11),Xoshiro(12),detector(counts,t),s10)
    plan=prepare_vh1_reconstruction(ft,t,s10);ri=reco_vh1_ions(r,plan)
    @test all(ri.valid)
    q=Float64.(ri.ions);predicted=sqrt.(diag(plan.ion_covariance));observed=vec(std(q;dims=2))
    println("VH1 10 ms noise-only ions per bucket: predicted ",round.(predicted;digits=2)," observed ",round.(observed;digits=2))
    @test all(abs.(observed./predicted.-1).<.08)
    @test abs(mean(q))<3*mean(predicted)/sqrt(6npads)
    Vr=frame_voltages(r,s10);d=diff(Vr[2:2:end,:];dims=1)./s10.gain_V_per_ion
    println("VH1 10 ms frame-difference ENC of the frame-end readings (5 ms windows): ",round(std(vec(d));digits=2)," ions")
end
println("VH1 10 ms tests passed")
