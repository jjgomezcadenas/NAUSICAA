using Test
include("../production/run_hit_production.jl")
@testset "Petit endpoint spheres" begin
    # A branched component has an independently known two-edge diameter.
    df=EventAnalysis.DataFrame(x=Float32[0,1,2,1],y=Float32[0,0,0,1],z=zeros(Float32,4),energy=Float32[1,2,3,4])
    tracks=make_tracks(df;method="kNN",k=1,max_distance_mm=1.1f0,energy_threshold_kev=0f0)
    t=only(tracks)
    u,v,path,_,_=EventAnalysis.find_extremes_mst_diameter(t,EventAnalysis.extract_coords(t))
    @test length(path)==3
    vox=t.voxels
    result=EventAnalysis.energy_in_spheres_around_extremes(t,(extremes=(vox[u,:],vox[v,:]),),1.0)
    for (idx,key) in ((u,:start_sphere_energy),(v,:end_sphere_energy))
        expected=sum(Float64(vox.energy[i]) for i in 1:4 if sqrt(sum((vox[i,a]-vox[idx,a])^2 for a in (:x,:y,:z)))<=1)
        @test result[key]==expected
    end
    @test result.blob1_energy>=result.blob2_energy
    source=read("/Users/jjgomezcadenas/Projects/Petit/src/track_blob_analysis.jl",String)
    localcode=read("src/event_analysis/petit/blobs.jl",String)
    a=findfirst("function energy_in_spheres_around_extremes",source).start
    b=findnext("\nend",source,a).stop
    @test occursin(source[a:b],localcode)
end
