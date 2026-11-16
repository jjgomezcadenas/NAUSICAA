using Test
include("../../scripts/production/prepare_batch_statistics.jl")
@testset "Main-track features" begin
    tc=Dict{String,Any}("ion_hit_cut_keV"=>5.0,"ion_voxel_width_mm"=>3.0,"ion_connection_mm"=>8.0,"ion_blob_radius_mm"=>7.0,
        "electron_hit_cut_keV"=>1.0,"electron_voxel_width_mm"=>4.0,"electron_connection_mm"=>8.0,"electron_blob_radius_mm"=>10.0)
    # A straight track of 20 voxels along X on the 3 mm lattice, 10 keV each, with 100 keV at one end
    x=Float32[1.5+3k for k in 0:19];e=fill(10f0,20);e[1]=100f0
    h=SavedHits(x,fill(1.5f0,20),fill(1.5f0,20),e)
    f=track_features(h,tc,"itrack")
    @test f.components==1 && f.hits==20 && f.main_voxels==20
    @test f.length_mm≈57 && f.extent_x_mm≈60 && f.extent_y_mm≈3 && f.extent_z_mm≈3
    @test f.main_keV≈290 && f.second_keV==0
    @test f.E1_keV≈120 && f.E2_keV≈30       # spheres of 7 mm: three voxels at each end
    # A detached deposit 30 mm away is a second component; hits below the cut are dropped
    h2=SavedHits(vcat(x,Float32[1.5,1.5,31.5]),vcat(fill(1.5f0,20),Float32[31.5,34.5,31.5]),fill(1.5f0,23),vcat(e,Float32[20,20,3]))
    f2=track_features(h2,tc,"itrack")
    @test f2.components==2 && f2.hits==22 && f2.second_keV≈40 && f2.length_mm≈57 && f2.main_keV≈290
    # No hit above the cut
    f0=track_features(SavedHits(Float32[0],Float32[0],Float32[0],Float32[2]),tc,"itrack")
    @test f0.components==0 && f0.hits==0 && isnan(f0.length_mm) && isnan(f0.E2_keV)
    # A bent track: 10 voxels along X then 10 along Y
    xb=vcat(Float32[1.5+3k for k in 0:9],fill(28.5f0,10));yb=vcat(fill(1.5f0,10),Float32[4.5+3k for k in 0:9])
    fb=track_features(SavedHits(xb,yb,fill(1.5f0,20),fill(10f0,20)),tc,"itrack")
    @test fb.length_mm≈57 && fb.extent_x_mm≈30 && fb.extent_y_mm≈33
end
@testset "Fits in slices" begin
    values=vcat(fill(10.,300).+repeat([-1.,0.,1.],100),fill(20.,50));by=vcat(fill(0.5,300),fill(1.5,50))
    s=sliced_fits(values,by,[0.,1.,2.];minimum_events=200)
    @test s["entries"]==[300,50] && isfinite(s["mu"][1]) && abs(s["mu"][1]-10)<0.2 && isnan(s["mu"][2])
end
@testset "Batch folders" begin
    @test batch_number("b00012")==12
    mktempdir() do d
        mkpath(joinpath(d,"vh1_hits"));mkpath(joinpath(d,"detector_response"));@test_throws ErrorException batch_output_subdir(d)
        touch(joinpath(d,"vh1_hits","features.h5"));@test batch_output_subdir(d)=="vh1_hits"
    end
end
