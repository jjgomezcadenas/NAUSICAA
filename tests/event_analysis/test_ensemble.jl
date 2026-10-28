using DataFrames,HDF5

@testset "Central-pad and first-arrival alignment" begin
    d = collection_table([-0.5,1.5],[0,4],[7,13])
    h = compute_occupancy_analysis(d;x_edges_mm=Float32[-2,-0,2,4],y_edges_mm=Float32[-1,1],
        z_edges_mm=Float32[0,2,4,6],drift_speed_mm_ms=0.1f0,drift_direction=1)
    a = align_pad_history(h,(0.25f0,0f0,0f0))
    @test a.pad_dx == Int64[-1,0] && a.arrival_bin == Int64[0,2]
    @test sum(a.ions)==20 && a.ions !== h.counts.ions
    # Translating the physical track and grid preserves the aligned pad labels.
    moved=copy(d);moved.x .+= 100f0;moved.y .-= 50f0;moved.z .+= 20f0
    hm = compute_occupancy_analysis(moved;x_edges_mm=h.x_edges_mm.+100f0,y_edges_mm=h.y_edges_mm.-50f0,
        z_edges_mm=h.z_edges_mm.+20f0,drift_speed_mm_ms=0.1f0,drift_direction=1)
    @test align_pad_history(hm,(100.25f0,-50f0,20f0)) == a
    @test_throws ArgumentError align_pad_history(h,(100f0,0f0,0f0))
end

@testset "Both ensemble denominators and full coverage" begin
    a = DataFrame(pad_dx=Int64[-1,0],pad_dy=Int64[0,0],arrival_bin=Int64[0,2],ions=Int64[7,13])
    b = DataFrame(pad_dx=Int64[0,3],pad_dy=Int64[0,-2],arrival_bin=Int64[1,0],ions=Int64[10,30])
    kw = (;processed_events=4,pitch_mm=2f0,z_bin_width_mm=2f0,drift_speed_mm_ms=0.1f0,minimum_half_cells=0)
    e = ensemble_pad_maps([a,b];kw...)
    @test e.accepted_events==2 && e.processed_events==4 && e.total_ions==60
    @test e.x_edges_mm==Float32[-3,-1,1,3,5,7] && e.y_edges_mm==Float32[-5,-3,-1,1]
    @test sum(e.xy_sum_ions)==sum(e.arrival_sum_ions)==60
    @test e.arrival_sum_ions==Int64[37,10,13]
    @test e.arrival_mean_accepted==Float32[18.5,5,6.5]
    @test e.arrival_mean_processed==Float32[9.25,2.5,3.25]
    @test e.time_edges_ms==Float32[0,20,40,60]
    @test all(e.xy_mean_accepted .== 2f0.*e.xy_mean_processed)
    reversed = ensemble_pad_maps([b,a];kw...)
    @test isequal(e,reversed)
    none = ensemble_pad_maps(DataFrame[];kw...)
    @test none.accepted_events==0 && none.total_ions==0
    @test all(isnan,none.xy_mean_accepted) && all(iszero,none.xy_mean_processed)
    @test_throws ArgumentError ensemble_pad_maps([a,b];kw...,processed_events=1)
    @test_throws ArgumentError ensemble_pad_maps([a];kw...,pitch_mm=0f0)
    bad = copy(a);bad.arrival_bin[1]=-1
    @test_throws ArgumentError ensemble_pad_maps([bad];kw...)
    bad = copy(a);bad.ions=[typemax(Int64),1]
    @test_throws OverflowError ensemble_pad_maps([bad];kw...)
    # Independent Float64 reference for a realistic-size integer numerator.
    large=copy(a);large.ions=Int64[123456789,987654321]
    result=ensemble_pad_maps([large];kw...,processed_events=8995)
    @test sum(Float64,result.xy_mean_processed) ≈ sum(large.ions)/8995 rtol=1e-7
end

@testset "Ensemble HDF5 round trip and corruption" begin
    mktempdir(joinpath(ROOT,"work")) do dir
        a=DataFrame(pad_dx=Int64[0],pad_dy=Int64[0],arrival_bin=Int64[0],ions=Int64[30])
        for tables in ([a],DataFrame[])
            e=ensemble_pad_maps(tables;processed_events=3,pitch_mm=2f0,z_bin_width_mm=2f0,drift_speed_mm_ms=0.1f0)
            path=joinpath(dir,"ensemble.h5")
            h5open(path,"w") do f;write_ensemble(f,e);end
            h5open(path,"r") do f
                r=read_ensemble(f)
                for name in keys(e);@test isequal(getproperty(r,name),getproperty(e,name));end
            end
        end
        e=ensemble_pad_maps([a];processed_events=3,pitch_mm=2f0,z_bin_width_mm=2f0,drift_speed_mm_ms=0.1f0)
        path=joinpath(dir,"bad.h5")
        h5open(path,"w") do f;write_ensemble(f,e);end
        h5open(path,"r+") do f
            f["arrival_sum_ions"][1:1]=Int64[31]
            @test_throws ErrorException read_ensemble(f)
        end
        h5open(path,"w") do f;write_ensemble(f,e);end
        h5open(path,"r+") do f
            f["arrival_mean_processed"][1:1]=Float32[11]
            @test_throws ErrorException read_ensemble(f)
        end
    end
end
