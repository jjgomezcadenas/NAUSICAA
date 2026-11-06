using Test,LinearAlgebra
include(joinpath(@__DIR__,"../src/ion_signal/ion_csa_nonlinear.jl"))
using .IonCSANonlinear
const E=IonCSANonlinear.ELEMENTARY_CHARGE_C
@testset "Nonlinear CSA independent analytic references" begin
    knots=[0.,.1,.2,.4,.8];g0=5e-5;tau=20.;base=.886
    # Linear gain rolloff makes the constant-current evolution exactly solvable.
    g=g0.*(1 .-knots);r=knots./tau
    m=NonlinearCSA(knots,g,r;baseline_V=base,step_s=.001)
    u=1000.;lambda=g0*u+1/tau;equilibrium=g0*u/lambda
    ts=collect(0.:.01:1.)
    exact=equilibrium.*(-expm1.(-lambda.*ts))
    wave=nonlinear_waveform(m,[0.,1.],[u*E,u*E];sample_time_s=ts)
    @test wave.valid
    @test maximum(abs.(wave.voltage_V.-(base.-exact)))<1e-11
    area=base.*ts.-equilibrium.*(ts.+expm1.(-lambda.*ts)./lambda)
    @test maximum(abs.(wave.voltage_integral_V_s.-area))<1e-11
    # Finite drive followed by zero input approaches baseline continuously.
    pulse=nonlinear_waveform(m,[0.,1.,1.01,100.],[u*E,u*E,0.,0.];sample_time_s=[1.01,20.,100.])
    @test pulse.valid
    @test pulse.voltage_V[1]<pulse.voltage_V[2]<pulse.voltage_V[3]<base
    @test (base-pulse.voltage_V[3])/(base-pulse.voltage_V[2])≈exp(-80/tau) rtol=1e-8
    # Fit independent analytic trajectories at several known constant currents.
    histories=map([500.,1000.,3000.,6000.,10000.]) do drive
        t=collect(0.:.002:4.);a=g0*drive+1/tau;q=g0*drive/a
        (;time_s=t,current_A=fill(drive*E,length(t)),voltage_V=base.-q.*(-expm1.(-a.*t)))
    end
    fit=fit_nonlinear_csa(histories,[0.,.05,.1,.2,.3];baseline_V=base,ridge=1e-12)
    @test maximum(abs.(fit.gain_V_per_e./(g0.*(1 .-fit.model.knots_V)).-1))<1e-4
    @test maximum(abs.(fit.recovery_V_s.-fit.model.knots_V./tau))<1e-6
    edges=collect(0.:.04:.8);width=.005;alpha=exp(-.04/tau);gain=4e-5
    sampled=sampled_response(m,[0.,1.],[u*E,u*E],edges;averaging_s=width,alpha,gain_V_per_ion=gain)
    primitive(t)=t<=0 ? base*t : base*t-equilibrium*(t+expm1(-lambda*t)/lambda)
    avg=[(primitive(t+width/2)-primitive(t-width/2))/width for t in edges]
    @test maximum(abs.(sampled.average_V.-avg))<1e-10
    @test maximum(abs.(sampled.estimated_ions.+((avg[2:end].-base).-alpha.*(avg[1:end-1].-base))./gain))<1e-6
end

include(joinpath(@__DIR__,"../src/ion_signal/ion_csa_nonlinear_io.jl"))
using HDF5,TOML
const IOModel=IonCSANonlinearIO
@testset "Validated calibration identity and production selection" begin
 mktempdir(joinpath(@__DIR__,"../work")) do dir
    p=joinpath(dir,"calibration.h5");circuit=joinpath(dir,"decks/fixture/corrected.inc")
    mkpath(dirname(circuit));write(circuit,"fixture circuit")
    h5open(p,"w") do h
        attributes(h)["schema"]="nausicaa.ion_csa_nonlinear.v1"
        for (key,value) in ("complete"=>true,"passed"=>true,
            "model_code_sha256"=>IOModel.sha(joinpath(@__DIR__,"../src/ion_signal/ion_csa_nonlinear.jl")),
            "knots_V"=>[0.,.2],"gain_V_per_e"=>[5e-5,4e-5],"recovery_V_s"=>[0.,.01],
            "baseline_V"=>.886,"step_s"=>1e-4,"maximum_current_A"=>1e-14,
            "validated_frame_ms"=>[20,30,40],"averaging_s"=>.005,"guard_V"=>.63,
            "maximum_error_ENC"=>.5,"circuit_sha256"=>IOModel.sha(circuit))
            h[key]=value
        end
    end
    c=Dict{String,Any}("nonlinear_calibration"=>p,"calibration_directory"=>dir,"variant"=>"fixture",
        "averaging_s"=>.005,"initial_spice_guard_V"=>.63,"linear_error_tolerance_ENC"=>.5,
        "input_timing"=>"uniform_mean_arrivals_with_finite_gap_current","pad_side_m"=>.0005,
        "gap_m"=>.005,"gap_speed_m_s"=>3.,"current_step_s"=>5e-6)
    @test IOModel.resolve_card(dir,c)==c # historical low-level default
    resolved=IOModel.resolve_card(dir,c;default_mode="hybrid_nonlinear")
    @test resolved["response_model"]=="hybrid_nonlinear"
    @test resolved["nonlinear_calibration_sha256"]==IOModel.sha(p)
    @test !haskey(c,"response_model")
    @test IOModel.resolve_card(dir,merge(resolved,Dict("mode"=>"production")))["mode"]=="production"
    @test IOModel.resolve_card(dir,merge(resolved,Dict("mode"=>"validation")))["mode"]=="validation"
    @test_throws ArgumentError IOModel.resolve_card(dir,merge(resolved,Dict("mode"=>"invalid")))
    @test_throws ArgumentError IOModel.resolve_card(dir,merge(c,Dict("mode"=>"production","response_model"=>"hybrid_spice")))
    @test IOModel.resolve_card(dir,merge(c,Dict("response_model"=>"hybrid_spice"));default_mode="hybrid_nonlinear")["response_model"]=="hybrid_spice"
    @test_throws ArgumentError IOModel.load_calibration(p;expected_sha256="different")
    @test_throws ArgumentError IOModel.resolve_card(dir,merge(resolved,Dict("gap_m"=>.01)))
    @test_throws ArgumentError IOModel.resolve_card(dir,merge(resolved,Dict("response_model"=>"unknown")))
    write(circuit,"changed circuit")
    @test_throws ArgumentError IOModel.resolve_card(dir,resolved)
    h5open(p,"r+") do h;write(h["passed"],false);end
    @test_throws ArgumentError IOModel.load_calibration(p)
 end
end
@testset "Nonlinear CSA domain and ownership" begin
    k=[0.,.1];g=[1e-4,1e-4];r=[0.,.005]
    m=NonlinearCSA(k,g,r;baseline_V=.9)
    k[2]=2.;g[1]=3.;r[2]=4.
    @test m.knots_V==[0.,.1] && m.gain_V_per_e==[1e-4,1e-4] && m.recovery_V_s==[0.,.005]
    t=[0.,1.];i=fill(1e5*E,2);before=copy(i)
    wave=nonlinear_waveform(m,t,i)
    @test !wave.valid && wave.reason=="excursion_outside_calibration"
    @test i==before && t==[0.,1.]
    @test_throws ArgumentError nonlinear_waveform(m,[1.,0.],i)
    @test_throws ArgumentError nonlinear_waveform(m,t,[-1.,0.])
    @test_throws ArgumentError nonlinear_waveform(m,t,i;sample_time_s=[2.])
    @test_throws ArgumentError NonlinearCSA([0.,0.],[1.,1.],[0.,1.];baseline_V=.9)
end
