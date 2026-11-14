using Test,Random,Statistics,Distributions
include("../src/ion_signal/sipm_response.jl")
const P=SiPMResponse
@testset "Finite-count occupancy moments" begin
    r=[.1,.4,.7,1.];M=length(r)
    for k in 0:4
        charges=Float64[]
        for assignment in Iterators.product(ntuple(_->1:M,k)...)
            push!(charges,sum(r[i] for i in unique(assignment);init=0.))
        end
        moments=P.occupancy_moments(sum(r),sum(abs2,r),M,k)
        @test moments.mean≈mean(charges) atol=1e-12
        @test moments.variance≈var(charges;corrected=false) atol=1e-12
    end
    @test P.occupancy_moments(7854.,7854.,7854,1).variance==0
    # Independent arbitrary-precision expression checks cancellation at large M.
    M=7854;k=2
    v=setprecision(128) do
        a=(1-1/BigFloat(M))^k;b=(1-2/BigFloat(M))^k
        M*(a-a*a)+M*(M-1)*(b-a*a)
    end
    @test P.occupancy_moments(Float64(M),Float64(M),M,k).variance≈Float64(v) rtol=1e-7
end
@testset "Aggregate recovery and bounded sampling" begin
    r=[.03,.3,.8,1.];R=sum(r);S=sum(abs2,r)
    for duration in (0.,.01,1.,100.)
        recovered=1 .- (1 .-r).*exp(-duration)
        nextR,nextS=P.recover_totals(R,S,length(r),duration)
        @test nextR≈sum(recovered)
        @test nextS≈sum(abs2,recovered)
    end
    rng=Xoshiro(81)
    samples=[first(P.bounded_moment_draw(rng,.4,.03,0.,1.)) for _ in 1:30000]
    @test abs(mean(samples)-.4)<5sqrt(.03/length(samples))
    @test abs(var(samples)-.03)<.001
    @test P.bounded_moment_draw(rng,.4,0.,0.,1.)==(.4,false)
    @test last(P.bounded_moment_draw(rng,.4,1.,0.,1.))
    @test first(P.bounded_moment_draw(rng,.4,.24,0.,1.)) in (0.,1.)
end
@testset "Physical limiting cases and reproducibility" begin
    dt=250e-9
    @test_throws ArgumentError P.StatisticalOptions(0)
    ideal=P.Sensor(.3,7854,1.,0.,0.)
    input=fill(Int64(100),20000);copyinput=copy(input)
    t=P.simulate_channel_statistical(Xoshiro(2),input,dt,ideal)
    @test input==copyinput
    @test t.charge_C≈t.detected_pe.*P.charge_per_pe_C(ideal)
    @test abs(mean(t.detected_pe)-30)<.15
    @test abs(var(t.detected_pe)-21)<.7
    @test all(iszero,P.simulate_channel_statistical(Xoshiro(1),zeros(Int,10),dt,ideal).charge_C)
    @test isempty(P.simulate_channel_statistical(Xoshiro(1),Int[],dt,ideal).charge_C)
    s=P.Sensor(1.,1,1.,73.4e-9,0.)
    # Single-cell sparse mode carries exact recovery; compare independent quadrature.
    samples=[sum(P.simulate_channel_statistical(Xoshiro(i),[2],dt,s).charge_C)/P.charge_per_pe_C(s) for i in 1:10000]
    h=dt/100000
    expected=1+sum((1-exp(-(i-.5)*h/s.recovery_s))*2*(dt-(i-.5)*h)/dt^2*h for i in 1:100000)
    @test abs(mean(samples)-expected)<5std(samples)/sqrt(length(samples))
    @test var(samples)>.01
    second= sum((2-exp(-(i-.5)*h/s.recovery_s))^2*2*(dt-(i-.5)*h)/dt^2*h for i in 1:100000)
    @test abs(var(samples)-(second-expected^2))<.005
    sensor=P.Sensor(.3,7854,2.55e6,73.4e-9,2.18e6)
    input=[0,1,1000,100000,1000,0,0]
    a=P.simulate_channel_statistical(Xoshiro(10),input,dt,sensor)
    b=P.simulate_channel_statistical(Xoshiro(10),input,dt,sensor)
    @test a==b
    @test a.charge_C≈a.optical_charge_C+a.dark_charge_C
    @test all(0 .<= a.optical_charge_C .<= a.detected_pe.*P.charge_per_pe_C(sensor).+1e-20)
    @test all(0 .<= a.dark_charge_C .<= a.dark_pulses.*P.charge_per_pe_C(sensor).+1e-20)
    @test a.variance_limits==0
    @test a.state_repairs==0
end
@testset "Fully recovered initial state" begin
    sensor=P.Sensor(1.,39984,1.,12.5e-9,7e5)
    @test P.StatisticalOptions().substeps==64
    @test P.StatisticalOptions().sparse_trigger_limit==16
    @test P.StatisticalOptions().initial_state==:fully_recovered
    @test_throws ArgumentError P.StatisticalOptions(128,16;initial_state=:unknown)
    # With no buckets, fully recovered initialization consumes no random draws.
    rng=Xoshiro(151);control=Xoshiro(151)
    P.simulate_channel_statistical(rng,Int64[],250e-9,sensor)
    @test rand(rng)==rand(control)
    historical=Xoshiro(151)
    P.simulate_channel_statistical(historical,Int64[],250e-9,sensor;options=P.StatisticalOptions(;initial_state=:stationary_dark))
    @test rand(historical)!=rand(Xoshiro(151))
    one=P.simulate_channel_statistical(Xoshiro(7),Int64[1],1e-12,sensor)
    @test one.dark_pulses==[0]
    @test one.optical_charge_C[1]==P.charge_per_pe_C(sensor)
    # Dark activity continues throughout acquisition, including zero-input channels.
    dark=P.simulate_channel_statistical(Xoshiro(8),zeros(Int64,20000),250e-9,sensor)
    @test abs(mean(dark.dark_pulses)-.175)<.015
    @test sum(dark.dark_charge_C)>0
end
