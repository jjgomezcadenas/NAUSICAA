using Test,Random,Statistics,Distributions
include("../scripts/production/run_elm_pilot.jl")
const R=ELMPilot.Replay
const EA=R.EventAnalysis
const PIPE=ELMPilot.SiPMProduction.Pilot
const P=PIPE.P
const E=PIPE.E
@testset "ELM analytical collection and linear response" begin
    args=(;tpb_efficiency=.65f0,forward_fraction=.5f0,guide_index=1.5f0,gap_index=1f0)
    a=EA.elm_light_acceptance(1.5f0,5f0;args...)
    n=100000;r=1.5;h=5.;theta=acos(1/1.5)
    axial=sum(.5*(1-z/hypot(z,r)) for z in ((j-.5)*h/n for j in 1:n))/n
    cone=sum(sin((j-.5)*theta/n) for j in 1:n)*theta/n
    @test a≈axial*.65*.5*cone rtol=3e-7
    @test EA.elm_light_acceptance(15f0,50f0;args...)≈a
    @test EA.elm_light_acceptance(1.5f0,5f0;args...,propagation_transmission=0f0)==0
    @test EA.elm_light_acceptance(1.5f0,5f0;merge(args,(gap_index=1.5f0,))...)==0
    @test_throws ArgumentError EA.elm_light_acceptance(1.5f0,5f0;merge(args,(tpb_efficiency=2f0,))...)
    sensor=P.LinearSensor(.32,3.6e5,7e5);q=P.charge_per_pe_C(sensor)
    incoming=fill(Int64(1000),50000);before=copy(incoming)
    trace=P.simulate_channel_linear(Xoshiro(11),incoming,250e-9,sensor)
    @test incoming==before
    @test abs(mean(trace.detected_pe)-320)<6sqrt(1000*.32*.68/length(incoming))
    @test var(trace.detected_pe)≈1000*.32*.68 rtol=.03
    @test mean(trace.dark_pulses)≈.175 atol=.01
    @test trace.charge_C≈q.*(trace.detected_pe+trace.dark_pulses)
    @test all(P.simulate_channel_linear(Xoshiro(1),[0,1000],250e-9,P.LinearSensor(1.,1.,0.)).detected_pe.==[0,1000])
    @test_throws ArgumentError P.LinearSensor(1.1,1.,0.)
    @test_throws ArgumentError P.simulate_channel_linear(Xoshiro(1),[-1],250e-9,sensor)
    # Analytic bucket clipping, including signed noise with no lower clipping.
    e=E.Electronics(250e-9,0.,1e-3)
    input=q.*[0.,4000.,5000.]
    t=E.apply_electronics(Xoshiro(2),input,q,e)
    @test t.saturated==[false,false,true]
    @test t.current_A≈min.(input./250e-9,1e-3)
    @test input==q.*[0.,4000.,5000.]
    e_noise=E.Electronics(250e-9,3.5416666667,1e-3)
    z=E.apply_electronics(Xoshiro(3),zeros(100000),q,e_noise)
    @test any(z.current_A.<0)
    @test std(z.noise_pe)≈e_noise.noise_rms_pe rtol=.01
    spec=R.load_specs(ELMPilot.ROOT,"config/specs/runs/elm_photons_ep2.toml")
    o=R.electron_optics(spec.electrons.data)
    @test o.band=="blue_tpb" && o.acceptance==a
    old=deepcopy(spec.electrons.data);delete!(old,"optical_model")
    @test R.electron_optics(old).acceptance==EA.axial_light_acceptance(1.5f0,5f0)
    s=PIPE.settings(ELMPilot.ROOT,"config/specs/runs/elm_sipm_response_ep2.toml")
    @test s.sensor isa P.LinearSensor
    @test s.sensor.pde==.32 && s.sensor.dark_rate_Hz≈7e5
end
