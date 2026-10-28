@testset "Pixel-frame charge statistics" begin
    d=DataFrame(x=Float32[0,0,1,8],y=zeros(Float32,4),z=Float32[0,2,4,6],
        energy=Float32[10,11,30,100].*21.9f-6,ions=Int64[10,11,30,100])
    original=deepcopy(d)
    result=collect_ion_cloud(d;pitch_mm=2f0,z_bin_width_mm=2f0,grid_offset_mm=(.2f0,.4f0),
        drift_speed_mm_ms=.1f0,drift_direction=1,minimum_xy_span_mm=0f0)
    @test d==original
    @test result.history.total_ions==151
    @test sum(result.history.counts.ions)==151
    @test all(isapprox.(diff(result.history.time_edges_ms),20f0;atol=.001))
    s=pixel_frame_statistics(result.history;n_min=Int64(10))
    @test sort(s.counts)==Int64[10,11,30,100]
    @test s.minimum_ions==10 && s.maximum_ions==100
    s.counts[1]=999
    @test sum(result.history.counts.ions)==151
    @test !s.minimum_histogram_accepted
    @test pixel_frame_statistics(result.history;n_min=Int64(9)).minimum_histogram_accepted
    @test !pixel_frame_statistics(result.history;n_min=Int64(100)).minimum_histogram_accepted
    empty=pixel_frame_statistics((counts=DataFrame(ions=Int64[]),);n_min=Int64(10))
    @test isempty(empty.counts) && empty.minimum_ions===nothing && empty.maximum_ions===nothing
    @test_throws ArgumentError pixel_frame_statistics(result.history;n_min=Int64(-1))
end
