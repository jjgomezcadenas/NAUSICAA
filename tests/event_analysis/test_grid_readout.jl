using Test, Random, DataFrames, Statistics
include(joinpath(@__DIR__,"../../src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis

@testset "EL yield and optical aperture" begin
    @test el_photon_yield(137f0,0.73f0,3f0,15f0,5f0) ≈ 2332.425f0 atol=0.002f0
    @test el_photon_yield(139f0,0.80f0,3f0,15f0,5f0) ≈ 2293.5f0 atol=0.002f0
    @test el_photon_yield(139f0,0.80f0,0.7f0,15f0,5f0) == 0f0
    radius,depth=2.5f0,5f0
    # Independent midpoint quadrature of the forward disk solid angle.
    n=100_000
    numerical=sum(0.5f0*(1f0-z/hypot(z,radius)) for z in
                  ((Float32(i)-0.5f0)*depth/Float32(n) for i in 1:n))/Float32(n)
    @test axial_light_acceptance(radius,depth) ≈ numerical atol=2f-5
    @test axial_light_acceptance(1f-5,depth) < 1f-5
    @test axial_light_acceptance(1f5,depth) > 0.49f0
end

@testset "Direct carrier collection" begin
    h=EventHits(Int64(3790),source_frame,Float32[0,0],Float32[0,0],
        Float32[0,1000],Float32[0.001,0.001],zeros(Float32,2),
        Int32[1,1],Int32[1,2])
    n=Int64[5000,5000]
    seed=UInt64(3790)
    kwargs=(;sigma_t_mm=1f0,sigma_l_mm=1.5f0,nsigma=3f0,
            reference_mm=(0f0,0f0,0f0),grid_offset_mm=(0.3f0,0.4f0),
            pitch_mm=2f0,bucket_width_mm=2f0,drift_speed_mm_ms=0.1f0,
            drift_direction=1,minimum_xy_span_mm=0f0)
    result=diffuse_to_grid(Xoshiro(seed),h,n;kwargs...)
    reference=diffused_carrier_reference(Xoshiro(seed),h,n;
        sigma_t_mm=1f0,sigma_l_mm=1.5f0,nsigma=3f0)
    @test reference.sampled_ions==result.sampled_ions
    @test reference.retained_ions==result.retained_ions
    @test reference.outside_ions==result.outside_ions
    @test result.sampled_ions==10000
    @test result.retained_ions+result.outside_ions==10000
    @test sum(result.counts.ions)==result.retained_ions
    @test maximum(result.counts.slice_index)-minimum(result.counts.slice_index)>490
    @test length(unique(result.counts.slice_index))>6
    @test result.time_edges_ms[2]-result.time_edges_ms[1] ≈ 20f0 atol=0.001f0
    opposite=diffuse_to_grid(Xoshiro(seed),h,n;merge(kwargs,(drift_direction=-1,))...)
    @test opposite.retained_ions==result.retained_ions
    @test opposite.time_edges_ms[2]-opposite.time_edges_ms[1] ≈ 20f0 atol=0.001f0
    # An independent replay assigns individual draws to the saved physical axes.
    rng=Xoshiro(seed)
    expected=Dict{NTuple{3,Int64},Int64}()
    sum_x=0.0;sum_y=0.0;min_z=Inf32;retained=0
    for i in eachindex(n),_ in 1:n[i]
        x=h.x_mm[i]+randn(rng,Float32)*1f0
        y=h.y_mm[i]+randn(rng,Float32)*1f0
        z=h.z_mm[i]+randn(rng,Float32)*1.5f0
        if -3f0<=x<=3f0 && -3f0<=y<=3f0 && -4.5f0<=z<=1004.5f0
            sum_x+=Float64(x);sum_y+=Float64(y);min_z=min(min_z,z);retained+=1
            key=ntuple(j->Int64(searchsortedlast((result.x_edges_mm,
                result.y_edges_mm,result.z_edges_mm)[j],(x,y,z)[j])),3)
            expected[key]=get(expected,key,0)+1
        end
    end
    actual=Dict((r.pad_ix,r.pad_iy,r.slice_index)=>r.ions for r in eachrow(result.counts))
    @test actual==expected
    @test reference.reference_mm==(Float32(sum_x/retained),Float32(sum_y/retained),min_z)
    @test_throws ArgumentError diffuse_to_grid(Xoshiro(seed),h,n;
        merge(kwargs,(bucket_width_mm=0f0,))...)
end

@testset "Direct grid boundaries" begin
    near=prevfloat(2f0)
    h=EventHits(Int64(9),source_frame,Float32[near,2f0],zeros(Float32,2),
        Float32[near,2f0],ones(Float32,2).*0.001f0,zeros(Float32,2),
        Int32[1,1],Int32[1,2])
    result=diffuse_to_grid(Xoshiro(1),h,Int64[1,1];
        sigma_t_mm=0f0,sigma_l_mm=0f0,nsigma=3f0,
        reference_mm=(0f0,0f0,0f0),grid_offset_mm=(0f0,0f0),
        pitch_mm=2f0,bucket_width_mm=2f0,drift_speed_mm_ms=0.1f0,
        drift_direction=1,minimum_xy_span_mm=0f0)
    @test result.counts.pad_ix==Int64[1,2]
    @test result.counts.slice_index==Int64[1,2]
    @test result.counts.ions==Int64[1,1]
end

@testset "PE sampling" begin
    rng=Xoshiro(2026)
    p=axial_light_acceptance(2.5f0,5f0)*0.3f0
    values=[sample_detected_pe(rng,Int64(1),2293.5f0,
                               axial_light_acceptance(2.5f0,5f0),0.3f0) for _ in 1:5000]
    trials=round(Int,2293.5)
    @test abs(mean(values)-trials*p) < 0.6
    @test abs(var(values)-trials*p*(1-p)) < 6
    @test sample_detected_pe(Xoshiro(1),0,2293.5f0,0.2f0,0.3f0)==0
    @test_throws ArgumentError sample_detected_pe(Xoshiro(1),1,2293.5f0,0.2f0,1.1f0)
end

@testset "Complete grid mapping and opposite arrival order" begin
    d=DataFrame(event_id=Int64[1,1,1],x=Float32[0,1,3],
                y=Float32[0,0,0],z=Float32[0,2,4],
                energy=Float32[1,2,3].*21.9f-6,ions=Int64[1,2,3])
    for (speed,direction,frame) in ((0.1f0,1,20f0),(1000f0,-1,0.002f0))
        result=collect_grid_history(d;reference_mm=(0f0,0f0,0f0),
            pitch_mm=2f0,z_bin_width_mm=2f0,grid_offset_mm=(0.3f0,0.4f0),
            drift_speed_mm_ms=speed,drift_direction=direction,
            minimum_xy_span_mm=0f0)
        @test sum(result.history.counts.ions)==6
        @test result.history.time_edges_ms[2]-result.history.time_edges_ms[1] ≈ frame atol=1f-5
        @test result.history.counts.slice_index[1] ==
              (direction==1 ? 1 : length(result.history.z_edges_mm)-1)
    end
end
