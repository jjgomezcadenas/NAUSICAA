using Test,Random,Statistics,LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"))
include(joinpath(ROOT,"src/detector_response/DetectorResponse.jl"))
include(joinpath(ROOT,"src/electronics/ElectronicsResponse.jl"))
include(joinpath(ROOT,"src/reconstruction/SignalReconstruction.jl"))
include(joinpath(ROOT,"src/vh1/VH1Readout.jl"))
using .EventAnalysis,.DetectorResponse,.ElectronicsResponse,.SignalReconstruction,.VH1Readout
const ER=ElectronicsResponse;const V=VH1Readout

card=joinpath(ROOT,"config/specs/response/ion_csa/vh1_device_changes.toml")
s=load_vh1_settings(ROOT,card;minimum_samples=4096)
noiseless(s)=VH1Settings(s.gain_V_per_ion,s.charge_fraction,s.feedback_capacitance_F,s.v_ref_V,s.window_V,s.frame_ms,s.conversions_per_frame,s.window_ms,s.reference_frames,
    ADCSettings(16,0f0,1.8f0,s.adc.sample_us,0f0),s.noise_sample_us,CSAVoltageNoise(s.noise_sample_us,zeros(length(s.noise.covariance_V2))),s.current_pulse,s.provenance)
withwindow(s,w;k=1)=VH1Settings(s.gain_V_per_ion,s.charge_fraction,s.feedback_capacitance_F,s.v_ref_V,s.window_V,s.frame_ms,k,Float32(w),s.reference_frames,s.adc,s.noise_sample_us,s.noise,s.current_pulse,s.provenance)
channels(n)=ReadoutChannels(Int32.(1:n),ones(Int32,n),Float32.(1:n),zeros(Float32,n))
placement=EventPlacement(1,50f0,100f0)
frame=EventFrame(mc_vertex_frame,(0f0,0f0,0f0),(1f0,1f0,1f0))
losses=ReadoutLosses(0,0,0,0,0,0)
function detector(counts::Matrix{Int64},t::IonReadoutTime)
    IonDetectorResponse(1,placement,frame,DetectorVectors{IonReadoutTime}(channels(size(counts,2)),t,counts),sum(counts),0,losses,SparseCells(Int64[],Int64[],Int64[],Int64[]))
end

@testset "Card and settings" begin
    @test s.gain_V_per_ion≈s.charge_fraction*ER.E_CHARGE/s.feedback_capacitance_F rtol=1e-3
    @test s.window_ms==10f0 && s.conversions_per_frame==2 && s.adc.bits==16 && s.frame_ms==20f0
    @test length(s.provenance["noise_spectrum_sha256"])==64
    @test_throws ArgumentError load_vh1_settings(ROOT,joinpath(ROOT,"config/readout/bb0nu_electronics_pilot.toml"))
end

@testset "Frame timing" begin
    t=IonReadoutTime(1000f0,20f0,4)
    ft=vh1_frame_time(t,s)
    onset=1000+1000first(s.current_pulse.time_s)
    @test ft.start_ms≈Float32(onset-20)
    @test ft.reference_frames==1
    # frames cover the end of the last bucket
    @test ft.start_ms+ft.frame_ms*ft.frame_count>=1000+80
    @test ft.start_ms+ft.frame_ms*(ft.frame_count-1)<1000+80
end

@testset "Exact recovery without noise and quantisation" begin
    t=IonReadoutTime(1000f0,20f0,5);counts=reshape(Int64[0,120,3000,7,0],5,1)
    for (w,k) in ((20.0,1),(5.0,1),(0.1,1),(5.0,4),(1.0,20))
        sq=withwindow(noiseless(s),w;k=k);ft=vh1_frame_time(t,sq)
        d=detector(counts,t);r=sim_vh1_electronics(Xoshiro(1),Xoshiro(2),d,sq)
        @test size(r.codes)==(ft.frame_count*sq.conversions_per_frame,1) && !r.saturated[1]
        # bypass the converter: reconstruct from the noiseless readings directly
        plan=prepare_vh1_reconstruction(ft,t,sq)
        grid=V.internal_grid(ft,sq);packet=ER.prepare_bucket_current(sq.current_pulse,t.bucket_ms)
        current=ER.prepared_ion_current(vec(counts),t,packet,grid)
        reading=V.frame_readings!(zeros(V.readings(ft,sq)),current,ft,sq)
        i0=sq.conversions_per_frame
        y=[reading[i0]-reading[i] for i in i0+1:length(reading)]
        q=plan.inverse*y
        @test maximum(abs.(q.-vec(counts)))<1e-7
        # with the 16-bit converter: within quantisation
        ri=reco_vh1_ions(r,plan)
        lsb_ions=(1.8/2^16)/sq.gain_V_per_ion
        @test maximum(abs.(vec(ri.ions).-vec(counts)))<2lsb_ions*sqrt(ft.frame_count)
        @test all(ri.valid)
    end
end

@testset "Charge per frame and spill into the next frame" begin
    t=IonReadoutTime(1000f0,20f0,3);sq=withwindow(noiseless(s),20.0;k=1);ft=vh1_frame_time(t,sq)
    plan=prepare_vh1_reconstruction(ft,t,sq)
    A=plan.response./sq.gain_V_per_ion   # one full-frame reading per frame
    inc=vcat(A[1:1,:],diff(A;dims=1))   # charge appearing in each frame
    pulse=sq.current_pulse
    finite=sum(diff(pulse.time_s).*(pulse.current_A[1:end-1]+pulse.current_A[2:end])/2)/ER.E_CHARGE
    @test abs(finite-(1-2/pi*atan((.00025)^2/(.005*sqrt(.005^2+2*(.00025)^2)))))<1e-3
    # the increments of one ion sum to its finite-path charge once the frames
    # have passed the end of its pulse (the last bucket completes inside the last frame)
    @test abs(sum(inc[:,1])-finite)<1e-9
    @test abs(sum(inc[:,2])-finite)<1e-9
    @test 0<finite-sum(inc[:,3])<.1
    # a bucket's charge appears in its own frame, the next one, and a trace in
    # the one after (integration weights late arrivals by the remaining frame)
    for k in 1:3
        col=inc[:,k]
        @test col[k]>.3*finite && col[k+1]>.3*finite && col[k+1]+col[k]>.95*finite
        @test all(abs(col[n])<1e-9 for n in eachindex(col) if n<k || n>k+2)
    end
    println("VH1 increments per frame for one ion in bucket 1: ",round.(inc[:,1];digits=4))
end

@testset "Window flag and converter limits" begin
    t=IonReadoutTime(1000f0,20f0,2);sq=noiseless(s)
    big=reshape(Int64[20000,0],2,1)   # 0.94 V drop leaves the 0.4 V window
    r=sim_vh1_electronics(Xoshiro(1),Xoshiro(2),detector(big,t),sq)
    @test r.saturated[1]
    ri=reco_vh1_ions(r,prepare_vh1_reconstruction(vh1_frame_time(t,sq),t,sq))
    @test !any(ri.valid) && all(isnan,ri.ions)
end

@testset "Noise-only pads reproduce the predicted covariance" begin
    t=IonReadoutTime(1000f0,20f0,6);ft=vh1_frame_time(t,s)
    npads=1500;counts=zeros(Int64,6,npads)
    r=sim_vh1_electronics(Xoshiro(11),Xoshiro(12),detector(counts,t),s)
    plan=prepare_vh1_reconstruction(ft,t,s)
    ri=reco_vh1_ions(r,plan)
    @test all(ri.valid)
    q=Float64.(ri.ions)
    predicted=sqrt.(diag(plan.ion_covariance));observed=vec(std(q;dims=2))
    println("VH1 noise-only ions per bucket: predicted ",round.(predicted;digits=2)," observed ",round.(observed;digits=2))
    @test all(abs.(observed./predicted.-1).<.08)
    @test abs(mean(q))<3*mean(predicted)/sqrt(6npads)
    # per-bucket noise against conversions per frame and window, 30 buckets
    t30=IonReadoutTime(1000f0,20f0,30)
    for (k,w) in ((1,20.0),(1,5.0),(2,10.0),(4,5.0),(10,2.0),(20,1.0),(4,2.0))
        sw=withwindow(s,w;k=k);ftw=vh1_frame_time(t30,sw);pw=prepare_vh1_reconstruction(ftw,t30,sw)
        sig=sqrt.(diag(pw.ion_covariance))
        println("conversions/frame ",k,", window ",w," ms: ions per bucket median ",round(median(sig);digits=1)," max ",round(maximum(sig);digits=1))
    end
    # frame-difference ENC of the raw readings, the quantity of the note
    Vr=frame_voltages(r,s);d=diff(Vr[s.conversions_per_frame:s.conversions_per_frame:end,:];dims=1)./s.gain_V_per_ion
    enc=std(vec(d))
    println("VH1 frame-difference ENC of the frame-end readings (",s.window_ms," ms windows): ",round(enc;digits=2)," ions")
    @test 16<enc<30
end
println("VH1 tests passed")
