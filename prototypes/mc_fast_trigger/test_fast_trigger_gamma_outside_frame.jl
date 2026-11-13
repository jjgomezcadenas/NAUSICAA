using Test
include("fast_trigger_gamma_outside_frame_functions.jl")
using .MCFastTrigger
particle(id,parent,name;primary=0)=(particle_id=Int32(id),mother_id=Int32(parent),particle_name=name,primary=primary,initial_x=0f0,initial_y=0f0,initial_z=0f0,kin_energy=.01f0,creator_proc="test")
@testset "Ancestry, distances and deposited energy" begin
    p=[particle(1,0,"e-";primary=1),particle(2,1,"e-"),particle(3,1,"gamma"),particle(4,3,"e-"),particle(5,4,"gamma"),particle(6,5,"e-")]
    @test gamma_roots(p)==Dict(1=>0,2=>0,3=>3,4=>3,5=>3,6=>3)
    h=(x_mm=Float32[0,10,12,11,100],y_mm=zeros(Float32,5),z_mm=zeros(Float32,5),energy_MeV=Float32[1,1,.1,.002,.003],particle_id=Int32[1,1,2,4,6])
    r=mc_fast_trigger(h,p;expected_primaries=1)
    @test r.event.primary_dx_mm==10
    @test r.event.gamma_free_dx_mm==12
    @test length(r.gammas)==1
    @test isapprox(r.gammas[1].energy_keV,5;atol=1e-6)
    @test r.gammas[1].nearest_hit_min_mm==1
    @test r.gammas[1].nearest_hit_max_mm==90
    @test r.event.outside_frame_energy_keV≈3
    @test r.event.total_energy_keV≈2105
    @test_throws ErrorException gamma_roots([particle(2,9,"e-")])
    @test_throws ErrorException gamma_roots([particle(2,3,"e-"),particle(3,2,"e-")])
end
@testset "Placement and sequential fiducial cuts" begin
    s=detector_selection(Float32[-10,0,10],100f0)
    @test s.z==Float32[90,100,110]
    @test s.contained && s.accepted
    @test !detector_selection(Float32[-10,0,1500],100f0).contained
    @test detector_selection(Float32[-10,0,10],25f0).contained
    @test !detector_selection(Float32[-10,0,10],25f0).accepted
    @test detector_selection(Float32[-50,1450],50f0).accepted
end
@testset "Radial frame excursion and inside-frame gamma deposits" begin
    p=[particle(1,0,"e-";primary=1),particle(2,1,"gamma"),particle(3,2,"e-"),particle(4,1,"gamma"),particle(5,4,"e-")]
    h=(x_mm=Float32[0,10,20,10,11],y_mm=Float32[0,10,10,20,11],z_mm=zeros(Float32,5),energy_MeV=Float32[1,1,.002,.003,.004],particle_id=Int32[1,1,3,3,5])
    r=mc_fast_trigger(h,p;expected_primaries=1,margin_mm=2f0)
    # Excursions in X and Y occur at different hits: radial maximum is 8, not sqrt(128).
    @test r.gammas[1].outside_frame_r_mm==8
    @test r.gammas[1].inside_frame==0
    @test r.gammas[2].inside_frame==1
    @test r.gammas[2].outside_frame_r_mm==0
end
@testset "Components use the closest 3D pair" begin
    p=[particle(1,0,"e-";primary=1),particle(2,1,"gamma"),particle(3,2,"e-")]
    h=(x_mm=Float32[0,10,3,10],y_mm=zeros(Float32,4),z_mm=Float32[0,0,4,6],energy_MeV=Float32[1,1,.002,.003],particle_id=Int32[1,1,3,3])
    r=mc_fast_trigger(h,p;expected_primaries=1)
    d=only(r.gammas)
    @test d.nearest_hit_min_mm==5
    @test d.nearest_pair_z_mm==4
    @test d.nearest_pair_r_mm==3
    @test d.nearest_hit_min_mm^2≈d.nearest_pair_z_mm^2+d.nearest_pair_r_mm^2
end
@testset "Gamma energy cut rejects whole events" begin
    inside=(inside_frame=UInt8(1),energy_keV=500.0)
    low=(inside_frame=UInt8(0),energy_keV=19.0)
    boundary=(inside_frame=UInt8(0),energy_keV=20.0)
    high=(inside_frame=UInt8(0),energy_keV=20.01)
    @test passes_gamma_cut([inside,low,boundary],20)
    @test !passes_gamma_cut([inside,low,high],20)
    @test passes_gamma_cut(NamedTuple[],20)
    @test passes_gamma_cut([high],25)
end
@testset "Five keV pre-selection" begin
    @test passes_gamma_cut([(inside_frame=UInt8(0),energy_keV=5.0)],5.0)
    @test !passes_gamma_cut([(inside_frame=UInt8(0),energy_keV=5.001)],5.0)
    @test passes_gamma_cut([(inside_frame=UInt8(1),energy_keV=100.0)],5.0)
end
