using DataFrames, Graphs, HDF5, Random

"""Small full-diffusion table with exact integer counts and common W."""
function collection_table(x,z,ions;y=zeros(Float32,length(x)),w=21.9f0)
    DataFrame(event_id=fill(Int64(42),length(x)),x=Float32.(x),y=Float32.(y),z=Float32.(z),
              energy=Float32.(ions).*w./1f6,ions=Int64.(ions))
end

@testset "Integer threshold after diffusion rebinning" begin
    # Two five-ion entries sum to the exact ten-ion threshold.
    d = collection_table([0.1,0.2,1.1,20.1],[0,0,0,0],[5,5,11,9])
    before = copy(d)
    r = select_diffused_event(d;bin_width_mm=1f0,max_distance_mm=2f0)
    @test r.accepted && length(r.tracks) == 1
    @test r.retained_ions == 30 && r.graph_retained_ions == 21 && r.graph_below_threshold_ions == 9
    @test sort(r.bins.ions) == [9,10,11]
    @test d == before
    @test sum(only(r.tracks).voxels.ions) == 21
    # A threshold-scale energy rounding change does not remove the ten-ion bin.
    d.energy[1] = prevfloat(d.energy[1])
    @test select_diffused_event(d;bin_width_mm=1f0,max_distance_mm=2f0).graph_retained_ions == 21
    d.ions[end] = 10
    @test !select_diffused_event(d;bin_width_mm=1f0,max_distance_mm=2f0).accepted
    @test !select_diffused_event(d;bin_width_mm=1f0,max_distance_mm=2f0,min_ions_per_bin=100).accepted
    @test isempty(select_diffused_event(d[1:0,:];bin_width_mm=1f0,max_distance_mm=2f0).tracks)
    @test_throws ArgumentError select_diffused_event(d;bin_width_mm=1f0,max_distance_mm=2f0,min_ions_per_bin=-1)
    for w in Float32[20,21.9,22,30], n in Int64[9,10,11]
        one = collection_table([0],[0],[n];w)
        @test select_diffused_event(one;bin_width_mm=1f0,max_distance_mm=2f0).accepted == (n>=10)
    end
    bad = collection_table([0,0],[0,0],[typemax(Int64),1])
    @test_throws OverflowError select_diffused_event(bad;bin_width_mm=1f0,max_distance_mm=2f0)
end

@testset "Full retained charge, coverage, phase and time direction" begin
    d = collection_table([0.1,1.1,100.1,-100.1],[0.1,1.1,8.1,-8.1],[14,14,6,3])
    r = select_diffused_event(d;bin_width_mm=1f0,max_distance_mm=3f0)
    @test r.accepted
    t = only(r.tracks)
    a = pad_axes(d,t;pitch_mm=2f0,z_bin_width_mm=2f0,grid_offset_mm=(0.25f0,0.5f0))
    @test first(a.x_edges_mm) <= minimum(d.x) <= maximum(d.x) <= last(a.x_edges_mm)
    @test first(a.z_edges_mm) <= minimum(d.z) <= maximum(d.z) <= last(a.z_edges_mm)
    @test all(diff(a.x_edges_mm) .≈ 2f0) && all(diff(a.z_edges_mm) .≈ 2f0)
    @test a.reference_mm[3] == minimum(t.voxels.z)
    h = compute_occupancy_analysis(d;x_edges_mm=a.x_edges_mm,y_edges_mm=a.y_edges_mm,z_edges_mm=a.z_edges_mm,
        drift_speed_mm_ms=0.1f0,drift_direction=1)
    p = pad_projections(h)
    @test h.total_ions == sum(h.counts.ions) == sum(p.xy_ions) == sum(p.slice_ions) == 37
    @test sum(t.voxels.ions) == 28 # below-threshold nine ions remain in pads
    @test all(diff(h.time_edges_ms) .≈ 20f0) && first(h.time_edges_ms) == 0f0
    @test h.x_edges_mm !== a.x_edges_mm
    reverse_time = compute_occupancy_analysis(d;x_edges_mm=a.x_edges_mm,y_edges_mm=a.y_edges_mm,z_edges_mm=a.z_edges_mm,
        drift_speed_mm_ms=0.1f0,drift_direction=-1)
    rp = pad_projections(reverse_time)
    @test rp.xy_ions == p.xy_ions && rp.slice_ions == reverse(p.slice_ions)
    @test reverse_time.time_origin_z_mm == last(a.z_edges_mm)
    # Integer-pitch phase shifts preserve physical occupied cells with relabelled local indices.
    b = pad_axes(d,t;pitch_mm=2f0,z_bin_width_mm=2f0,grid_offset_mm=(2.25f0,2.5f0))
    hb = compute_occupancy_analysis(d;x_edges_mm=b.x_edges_mm,y_edges_mm=b.y_edges_mm,z_edges_mm=b.z_edges_mm,
        drift_speed_mm_ms=0.1f0,drift_direction=1)
    physical(history) = sort([(history.x_edges_mm[row.pad_ix],history.y_edges_mm[row.pad_iy],row.slice_index,row.ions) for row in eachrow(history.counts)])
    @test physical(h) == physical(hb)
    narrow = Float32[-1,1]
    @test_throws ArgumentError compute_occupancy_analysis(d;x_edges_mm=narrow,y_edges_mm=narrow,z_edges_mm=narrow,
        drift_speed_mm_ms=0.1f0,drift_direction=1)
    @test_throws ArgumentError pad_axes(d,t;pitch_mm=0f0,z_bin_width_mm=2f0,grid_offset_mm=(0f0,0f0))
end

@testset "Sparse histories, exact edges, zero slices and ownership" begin
    d = collection_table([-1,0,1],[0,2,6],[4,5,6];y=Float32[-1,0,1])
    kw = (;x_edges_mm=Float32[-1,0,1],y_edges_mm=Float32[-1,0,1],z_edges_mm=Float32[0,2,4,6],drift_speed_mm_ms=0.1f0,drift_direction=1)
    h = compute_occupancy_analysis(d;kw...)
    @test h.counts.pad_ix == Int64[1,2,2] && h.counts.pad_iy == Int64[1,2,2]
    @test h.counts.slice_index == Int64[1,2,3]
    @test h.counts.ions == Int64[4,5,6]
    @test h.time_edges_ms == Float32[0,20,40,60]
    empty = compute_occupancy_analysis(d[1:0,:];kw...)
    @test isempty(empty.counts) && all(eltype(c)==Int64 for c in eachcol(empty.counts))
    @test pad_projections(empty).slice_ions == zeros(Int64,3)
    gaps = compute_occupancy_analysis(d[[1,3],:];kw...)
    @test pad_projections(gaps).slice_ions == Int64[4,0,6]
    h.x_edges_mm[1] = -99f0
    @test kw.x_edges_mm[1] == -1f0
    @test_throws ArgumentError compute_occupancy_analysis(d;kw...,drift_direction=0)
    @test_throws ArgumentError compute_occupancy_analysis(d;kw...,drift_speed_mm_ms=0f0)
    @test_throws ArgumentError compute_occupancy_analysis(d;kw...,z_edges_mm=Float32[0,0,2])
    @test_throws ArgumentError compute_occupancy_analysis(d;kw...,z_edges_mm=Float32[0,Inf])
    outside = copy(d);outside.x[end] = nextfloat(1f0)
    @test_throws ArgumentError compute_occupancy_analysis(outside;kw...)
end

@testset "Pinned Petit occupancy parity" begin
    h5open(joinpath(@__DIR__,"fixtures/petit_collection_reference.h5"),"r") do f
        d = collection_table(read(f["input/x"]),read(f["input/z"]),read(f["input/ions"]);y=read(f["input/y"]),w=22f0)
        # Original occupancy consumes track centres; use the same table for parity.
        t = Tracks(copy(d),SimpleGraph(nrow(d)),[collect(1:nrow(d))],DiffusionParams())
        a = pad_axes(d,t;pitch_mm=2f0,z_bin_width_mm=2f0,minimum_xy_span_mm=16f0,grid_offset_mm=(0f0,0f0))
        h = compute_occupancy_analysis(d;x_edges_mm=a.x_edges_mm,y_edges_mm=a.y_edges_mm,z_edges_mm=a.z_edges_mm,
            drift_speed_mm_ms=0.1f0,drift_direction=1)
        p = pad_projections(h)
        @test p.xy_ions ≈ read(f["total_counts"]) atol=2e-6
        @test p.xy_ions./length(p.slice_ions) ≈ read(f["avg_counts"]) atol=2e-6
        @test length(p.slice_ions) == read(f["n_time_slots"])
        occupancy = [count(==(i),h.counts.slice_index)/length(p.xy_ions) for i in eachindex(p.slice_ions)]
        @test occupancy ≈ read(f["occupancy_per_slot"])
    end
end

@testset "Hits through diffusion, selection and pads" begin
    function run_fixture(xs)
        h = EventHits(42,source_frame,Float32.(xs),zeros(Float32,length(xs)),zeros(Float32,length(xs)),
            fill(0.2f0,length(xs)),zeros(Float32,length(xs)),fill(Int32(7),length(xs)),Int32.(0:length(xs)-1))
        counts = ion_counts(h,21.9f0)
        d = diffuse_xyz_image_mc(Xoshiro(745),h,counts.sampled;sigma_t_mm=0.5f0,sigma_l_mm=0.5f0,
            nbins=100,nsigma=3f0,ionisation_energy_eV=21.9f0)
        b,_,distance = get_voxel_size_and_distance(0.5f0)
        s = select_diffused_event(d.bins;bin_width_mm=b,max_distance_mm=distance)
        return d,s
    end
    d,s = run_fixture([0])
    @test s.accepted
    a = pad_axes(d.bins,only(s.tracks);pitch_mm=2f0,z_bin_width_mm=2f0,grid_offset_mm=(0.1f0,0.2f0))
    h = compute_occupancy_analysis(d.bins;x_edges_mm=a.x_edges_mm,y_edges_mm=a.y_edges_mm,z_edges_mm=a.z_edges_mm,
        drift_speed_mm_ms=0.1f0,drift_direction=1)
    @test sum(h.counts.ions) == d.retained_ions
    @test sum(h.counts.ions)+d.outside_ions == d.sampled_ions
    @test s.graph_retained_ions+s.graph_below_threshold_ions == d.retained_ions
    _,background = run_fixture([0,100])
    @test !background.accepted && length(background.tracks) == 2
end
