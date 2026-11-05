using Test,Random,Statistics
include("../scripts/production/sipm_production.jl")
const Prod=SiPMProduction
const P=Prod.P
module Replay
include("../scripts/production/run_electron_250ns.jl")
end
@testset "ELM placement replay" begin
    e=Dict("placement"=>Dict("mode"=>"uniform_electron_reference","minimum_electron_drift_cm"=>10.,"maximum_electron_drift_cm"=>150.,"seed"=>2026))
    r=Dict("chamber_length_cm"=>150.)
    placements=[Replay.electron_placement(e,r,"fixed-primary-stream",Int64(i),nothing) for i in 1:10000]
    distances=[p.electron_drift_cm for p in placements]
    @test all(10 .<=distances.<150)
    @test all(p.ion_drift_cm+p.electron_drift_cm≈150 for p in placements)
    @test abs(mean(distances)-80)<2
    @test Replay.electron_placement(e,r,"fixed-primary-stream",Int64(1),nothing)==placements[1]
end
@testset "Illuminated disk and full-area dark response" begin
    settings=Prod.Pilot.settings(Prod.ROOT,"config/specs/runs/elm_sipm_response_tau12p5ns.toml")
    @test settings.sensor.microcells==39984
    @test settings.c["illuminated_microcell_count"]==round(Int,39984pi/4)
    @test settings.sensor.dark_rate_Hz≈700000
    n=100000;dt=250e-9
    sensor=P.Sensor(.32,39984,3.6e5,0.,700000.)
    input=fill(Int64(100),n)
    t=P.simulate_channel_illuminated(Xoshiro(12),input,dt,sensor,31403)
    @test all(input.==100)
    @test abs(mean(t.detected_pe)-32)<6sqrt(100*.32*.68/n)
    @test abs(mean(t.dark_pulses)-.175)<6sqrt(.175/n)
    @test all(isapprox.(t.charge_C./P.charge_per_pe_C(sensor),t.detected_pe+t.dark_pulses;atol=1e-10))
    @test all(t.expected_distinct_optical_cells.<=t.detected_pe)
    @test all(t.expected_distinct_optical_cells.<=31403)
    # Repeated optical triggers exceed the number of distinct cells; finite
    # recovery allows charge greater than the cell count over a long bucket.
    counts=fill(Int64(200000),10)
    s=P.Sensor(1.,39984,3.6e5,12.5e-9,0.)
    a=P.simulate_channel_illuminated(Xoshiro(22),counts,dt,s,31403)
    b=P.simulate_channel_illuminated(Xoshiro(22),counts,dt,s,31403;options=P.StatisticalOptions(256))
    qa=sum(a.charge_C)/P.charge_per_pe_C(s);qb=sum(b.charge_C)/P.charge_per_pe_C(s)
    @test all(a.detected_pe.==counts)
    @test all(a.charge_C./P.charge_per_pe_C(s).>31403)
    @test qa<sum(counts)
    @test abs(qa/qb-1)<.005
    @test a.state_repairs==a.variance_limits==0
    # Dark-only populations conserve the full rate and recover independently.
    @test_throws ArgumentError P.simulate_channel_illuminated(Xoshiro(1),[1],dt,s,0)
    @test_throws ArgumentError P.simulate_channel_illuminated(Xoshiro(1),[1],dt,s,40000)
end
