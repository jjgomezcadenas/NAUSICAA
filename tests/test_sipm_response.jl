using Test, Random, Statistics, Distributions, HDF5, TOML
include("../scripts/production/sipm_pipeline.jl")
const SP=SiPMPipeline
const P=SP.P
const E=SP.E
const EA=SP.EA
const ROOT=normpath(joinpath(@__DIR__,".."))
const metrics=Dict{String,Float64}()

@testset "Sensor response and physical fluctuations" begin
    @test_throws ArgumentError P.Sensor(1.1,100,1.,1e-9,0.)
    @test_throws ArgumentError P.Sensor(.3,0,1.,1e-9,0.)
    @test P.recovered_fraction(0.,73.4e-9)==0
    @test P.recovered_fraction(73.4e-9,73.4e-9)≈1-exp(-1)
    @test P.recovered_fraction(Inf,73.4e-9)==1
    @test P.recovered_fraction(0.,0.)==1
    dt=250e-9
    ideal=P.Sensor(.3,100,2.55e6,0.,0.)
    counts=fill(Int64(100),20_000); original=copy(counts)
    x=P.simulate_channel(Xoshiro(12),counts,dt,ideal)
    @test counts==original
    @test abs(mean(x.detected_pe)-30)<.12
    @test abs(var(x.detected_pe)-21)<.6
    @test all(x.detected_pe .<= counts)
    @test x.charge_C≈x.detected_pe.*P.charge_per_pe_C(ideal)
    @test all(iszero,x.dark_pulses)
    metrics["binomial_mean"]=mean(x.detected_pe);metrics["binomial_variance"]=var(x.detected_pe)
    zero=P.simulate_channel(Xoshiro(1),counts,dt,P.Sensor(0.,1,1.,1e-7,0.))
    @test all(iszero,zero.charge_C)
    allpe=P.simulate_channel(Xoshiro(1),Int64[0,5,0],dt,P.Sensor(1.,10,1.,0.,0.))
    @test allpe.detected_pe==[0,5,0]
    @test isempty(P.simulate_channel(Xoshiro(1),Int64[],dt,ideal).charge_C)

    # Exact two-trigger result, one initially recovered cell: first charge=1,
    # second charge=1-exp(-D/tau), D has density 2(T-D)/T^2.
    # Independent midpoint quadrature supplies the expectation.
    tau=73.4e-9; single=P.Sensor(1.,1,1.,tau,0.); n=12_000
    draws=[sum(P.simulate_channel(Xoshiro(i),Int64[2],dt,single).charge_C)/P.charge_per_pe_C(single) for i in 1:n]
    steps=100_000; delta=dt/steps
    expected=1+sum((1-exp(-((j-.5)*delta)/tau))*2*(dt-(j-.5)*delta)/dt^2*delta for j in 1:steps)
    @test abs(mean(draws)-expected)<5std(draws)/sqrt(n)
    @test var(draws)>0.01 # random arrival separation fluctuates despite fixed charge law
    metrics["two_trigger_mean"]=mean(draws);metrics["two_trigger_quadrature"]=expected
    metrics["two_trigger_variance"]=var(draws)
    # Same repeated finite-bucket illumination: no artificial reset at a boundary.
    sustained=P.simulate_channel(Xoshiro(32),fill(Int64(20),20_000),dt,single)
    @test mean(sustained.charge_C)/P.charge_per_pe_C(single)<4
    @test all(sustained.charge_C .<= 20P.charge_per_pe_C(single))

    # Stationary Poisson arrivals independently generated at count level.
    stationary=P.Sensor(1.,40,1.,tau,0.)
    arrivals=rand(Xoshiro(28),Poisson(80.),30_000)
    response=P.simulate_channel(Xoshiro(39),arrivals,dt,stationary)
    measured=mean(response.charge_C[101:end])/P.charge_per_pe_C(stationary)
    prediction=P.axel_mean(80.,stationary,dt)
    @test abs(measured/prediction-1)<.005
    @test P.axel_inverse(prediction,stationary,dt).pe≈80
    @test !P.axel_inverse(prediction,stationary,dt;clipped=true).valid
    @test !P.axel_inverse(stationary.microcells*dt/tau*1.001,stationary,dt).valid
    @test P.axel_inverse(-.5,stationary,dt).valid
    metrics["stationary_mean_pe"]=measured;metrics["axel_mean_pe"]=prediction
    # Evaluate the proposed Float32 storage precision against Float64 charge sums.
    storage_error=abs(sum(Float64.(Float32.(response.charge_C)))/sum(response.charge_C)-1)
    metrics["float32_charge_storage_relative_error"]=storage_error
    @test storage_error<1e-6
    # A finite-bin inverse has a bias and dispersion; save both for review.
    inverse=[P.axel_inverse(q/P.charge_per_pe_C(stationary),stationary,dt).pe for q in response.charge_C[101:end]]
    metrics["stationary_inverse_mean_pe"]=mean(inverse)
    metrics["stationary_inverse_rms_pe"]=std(inverse)
    transient=[P.simulate_channel(Xoshiro(i),Int64[80],dt,stationary).charge_C[1]/P.charge_per_pe_C(stationary) for i in 1:4000]
    metrics["transient_mean_pe"]=mean(transient)
    metrics["transient_inverse_mean_pe"]=mean(P.axel_inverse(m,stationary,dt).pe for m in transient)
    @test mean(transient)>prediction # first bright bucket differs from stationary formula

    rate=P.dark_rate_Hz(1e6,9e-6,.005)
    @test P.dark_rate_Hz(1e6,9e-6,.010)≈4rate
    dark=P.simulate_channel(Xoshiro(89),zeros(Int64,150_000),dt,P.Sensor(0.,7854,1.,tau,rate))
    @test abs(mean(dark.dark_pulses)-rate*dt)<.008
    @test abs(var(dark.dark_pulses)-rate*dt)<.012
    @test all(iszero,dark.optical_charge_C)
    metrics["dark_mean"]=mean(dark.dark_pulses);metrics["dark_variance"]=var(dark.dark_pulses)
    @test P.simulate_channel(Xoshiro(42),Int64[10,0,10],dt,ideal)==
          P.simulate_channel(Xoshiro(42),Int64[10,0,10],dt,ideal)
end

@testset "ASIC units, noise and clipping" begin
    e=E.Electronics(250e-9,.5,1e-3);q0=P.elementary_charge_C*2.55e6
    x=E.apply_electronics(Xoshiro(2),zeros(100_000),q0,e)
    @test abs(mean(x.noise_pe))<.006
    @test abs(std(x.noise_pe)-.5)<.006
    @test abs(cor(x.noise_pe[1:end-1],x.noise_pe[2:end]))<.015
    @test any(x.current_A.<0) && !any(x.saturated)
    @test x.current_A.*e.bucket_width_s≈q0.*x.noise_pe
    metrics["electronics_noise_mean_pe"]=mean(x.noise_pe)
    metrics["electronics_noise_rms_pe"]=std(x.noise_pe)
    ideal=E.Electronics(e.bucket_width_s,0.,e.current_limit_A)
    q=[0.,.5,1.,2.].*e.current_limit_A.*e.bucket_width_s
    original=copy(q); trace=E.apply_electronics(Xoshiro(1),q,q0,ideal)
    @test q==original
    @test trace.saturated==[false,false,false,true]
    @test trace.current_A≈[0.,.5,1.,1.].*e.current_limit_A
    @test trace.current_unclipped_A.*e.bucket_width_s≈q
end

@testset "Signed buckets, empty acquisition and saved independent stages" begin
    s=SP.settings(ROOT,"config/specs/runs/sipm_response_pilot.toml")
    grid=EA.GridSpec(10f0,.25f0,(0f0,0f0,0f0),(0f0,0f0),Int8(-1))
    placement=EA.EventPlacement(Int64(0),50f0,100f0)
    cells=EA.SparseCells(Int64[-1,-1,2],Int64[3,3,1],Int64[-2,2,0],Int64[20,40,10])
    history=EA.PhotonSiPMHistory(placement,:pure_xe,grid,cells,Int64(100),Int64(90),Int64(1000))
    a=SP.acquisition(history,s.c["acquisition"])
    @test a.buckets==collect(-6:6) && sum(a.counts)==70
    @test a.counts[5,1]==20 && a.counts[9,1]==40 && a.counts[7,2]==10
    @test all(iszero,a.counts[6,:])
    emptyhistory=EA.PhotonSiPMHistory(placement,:pure_xe,grid,EA.SparseCells(Int64[],Int64[],Int64[],Int64[]),Int64(0),Int64(0),Int64(0))
    control=deepcopy(s.c["acquisition"]);control["additional_channels"]=[[0,0]]
    @test size(SP.acquisition(emptyhistory,control).counts)==(8,1)
    @test size(SP.acquisition(emptyhistory,s.c["acquisition"]).counts)==(8,0)
    mkpath(joinpath(ROOT,"work/sipm_tests"))
    mktempdir(joinpath(ROOT,"work/sipm_tests")) do dir
        input=joinpath(dir,"photons.h5")
        h5open(input,"w") do h
            m=create_group(h,"metadata");m["test_source"]=true
            g=create_group(h,"events/0");EA.write_photon_history(g,history);g["pair_checksum"]="test-pairs"
        end
        sourcehash=SP.sha(input);stream=SP.S.content_hash(s.sensor_model*sourcehash)
        src=(;path=input,sample="fixture",gas="pure_xe",sourcehash,dt=250e-9,
            card=Dict("drift_speed_mm_us"=>1.),stream,sensorid="fixture-sensor",eid="fixture-electronics",
            sensorpath=joinpath(dir,"sipm.h5"),response=joinpath(dir,"electronics.h5"))
        SP.sensor_file(s,src);SP.electronics_file(s,src)
        hashbefore=SP.sha(src.sensorpath)
        @test SP.sensor_file(s,src)==src.sensorpath
        @test SP.sha(src.sensorpath)==hashbefore
        # Persist a dark-only event with an explicitly selected acquisition channel.
        emptyinput=joinpath(dir,"empty_photons.h5")
        h5open(emptyinput,"w") do h
            create_group(h,"metadata")["test_source"]=true
            g=create_group(h,"events/0");EA.write_photon_history(g,emptyhistory);g["pair_checksum"]="empty-pairs"
        end
        emptycard=deepcopy(s.c);emptycard["acquisition"]=control
        emptys=merge(s,(;c=emptycard))
        emptysrc=merge(src,(;path=emptyinput,sourcehash=SP.sha(emptyinput),sensorid="empty-sensor",eid="empty-electronics",
            sensorpath=joinpath(dir,"empty_sensor.h5"),response=joinpath(dir,"empty_electronics.h5")))
        SP.sensor_file(emptys,emptysrc);SP.electronics_file(emptys,emptysrc)
        h5open(emptysrc.response,"r") do h
            @test all(iszero,read(h["events/0/incident_photons"]))
            @test length(read(h["events/0/current_A"]))==8
            @test all(isfinite,read(h["events/0/current_A"]))
        end
        h5open(src.sensorpath,"r") do sensor
            @test SP.schema(sensor)=="nausicaa.sipm_file.v1"
            @test read(sensor["metadata/complete"])
            @test read(sensor["events/0/pair_checksum"])=="test-pairs"
            @test sum(read(sensor["events/0/incident_photons"]))==70
            roundtrip=EA.read_photon_history(sensor["events/0/photon_input"])
            @test roundtrip.cells.bucket==cells.bucket
            @test roundtrip.cells.count==cells.count
            h5open(src.response,"r") do elec
                @test read(elec["events/0/charge_C"])==read(sensor["events/0/charge_C"])
                @test read(elec["events/0/bucket"])==read(sensor["events/0/bucket"])
                @test read(elec["metadata/source_sha256"])==hashbefore
            end
        end
        # Electronics changes reuse exactly the sensor file and its stochastic state.
        e2=merge(s,(;electronics=E.Electronics(250e-9,1.,1e-3),electronics_model="changed-noise"))
        src2=merge(src,(;eid="changed-electronics",response=joinpath(dir,"electronics2.h5")))
        SP.electronics_file(e2,src2)
        @test SP.sha(src.sensorpath)==hashbefore
        h5open(src.sensorpath,"r+") do h
            write(h["events/0/simulation_elapsed_s"],123.0)
        end
        # A regenerated sensor's runtime metadata may differ; stale electronics
        # must rebuild and bind to the new complete source file.
        SP.electronics_file(s,src)
        h5open(src.response,"r") do h
            @test read(h["metadata/source_sha256"])==SP.sha(src.sensorpath)
        end
        # Stable event/channel seeds are independent of the subset/order of events.
        seed(id)=SP.S.event_seed(stream,Int64(id),"sensor:-1:3:2026")
        forward=Dict(id=>P.simulate_channel(Xoshiro(seed(id)),Int64[10,0,20],250e-9,s.sensor) for id in [0,7])
        reverse=Dict(id=>P.simulate_channel(Xoshiro(seed(id)),Int64[10,0,20],250e-9,s.sensor) for id in [7,0])
        @test forward==reverse
    end
end

mkpath(joinpath(ROOT,"results/sipm_response_validation"))
write(joinpath(ROOT,"results/sipm_response_validation/.gitignore"),"*.h5\n")
h5open(joinpath(ROOT,"results/sipm_response_validation/validation.h5"),"w") do out
    attributes(out)["schema"]="nausicaa.sipm_validation.v1"
    for (k,v) in metrics;out[k]=v;end
    out["all_tests_passed"]=true
    m=create_group(out,"metadata");m["julia_version"]=string(VERSION)
    for (j,path) in enumerate([@__FILE__,joinpath(ROOT,"src/ion_signal/sipm_response.jl"),
        joinpath(ROOT,"src/ion_signal/sipm_electronics.jl"),joinpath(ROOT,"scripts/production/sipm_pipeline.jl")])
        g=create_group(m,string(j));g["path"]=relpath(path,ROOT)
        g["sha256"]=SP.sha(path);g["contents"]=read(path,String)
    end
end
