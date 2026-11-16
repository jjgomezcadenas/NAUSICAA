using Test,Random
include("../../scripts/production/prepare_species_comparison.jl")
@testset "Separation by a cut" begin
    # Perfect separation and none
    c=cut_curve([10.,11.,12.],[1.,2.,3.]);@test c.area≈1
    c=cut_curve([1.,2.,3.,4.],[1.,2.,3.,4.]);@test c.area≈0.5
    # Hand-checked curve: signal {2,4}, background {1,3}
    c=cut_curve([2.,4.],[1.,3.])
    @test c.thresholds==[-Inf,1,2,3,4] && c.signal_efficiency==[1,1,.5,.5,0] && c.background_acceptance==[1,.5,.5,0,0]
    @test c.area≈0.75
    w=working_point(c,0.5);@test w.threshold==2 && w.signal_efficiency==0.5 && w.background_acceptance==0.5
    # Non-finite values fail the cut and stay in the denominator
    c=cut_curve([2.,NaN],[1.,3.]);@test c.signal_efficiency[1]==0.5
    # Two Gaussians one sigma apart: area = Phi(1/sqrt(2)) = 0.7602
    rng=Xoshiro(1);c=cut_curve(randn(rng,200_000).+1,randn(rng,200_000));@test isapprox(c.area,0.7602;atol=0.003)
    @test issorted(c.signal_efficiency;rev=true) && issorted(c.background_acceptance;rev=true)
end
@testset "ROI and track selection" begin
    lo,hi=roi_bounds(2457.83,5.0);@test hi-lo≈11.774 && (lo+hi)/2≈2457.83
    d=Dict("itrack_components"=>[1,1,2],"etrack_components"=>[1,2,1],"itrack_selected"=>Int8[0,1,1])
    @test track_selected(d,"itrack")==[false,true,true] && track_selected(d,"etrack")==[true,false,true]
    delete!(d,"itrack_selected");@test track_selected(d,"itrack")==[true,true,false]      # tables without the lattice selection
end
