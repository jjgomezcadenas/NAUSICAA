using DataFrames, Graphs, HDF5, LinearAlgebra
include("reconstruction_fixtures.jl")

@testset "Petit binning and adapter" begin
    h = EventHits(42,source_frame,Float32[-1,-0.01,0,0.9,1],zeros(Float32,5),
        zeros(Float32,5),Float32[1,2,4,8,16],Float32[4,3,2,1,0],fill(Int32(7),5),Int32.(0:4))
    d = to_petit_dataframe(h)
    @test d.time == h.time_ns
    @test d.x !== h.x_mm
    b = sort(voxelize_event(d,1f0),:x)
    @test b.x == Float32[-0.5,0.5,1.5]
    @test b.energy == Float32[3,12,16]
    @test sum(b.energy) == sum(h.energy_MeV)
    @test all(eltype(b[!,s]) == Float32 for s in (:x,:y,:z,:energy))
    @test all(b.y .== 0.5f0) && all(b.z .== 0.5f0)
    boundary = DataFrame(x=Float32[prevfloat(1f0),1f0,nextfloat(1f0)],y=zeros(Float32,3),z=zeros(Float32,3),energy=ones(Float32,3))
    @test sort(voxelize_event(boundary,1f0),:x).energy == Float32[1,2]
    @test nrow(voxelize_event(d[1:0,:],1f0)) == 0
    @test_throws ArgumentError voxelize_event(d,0f0)
    @test_throws ArgumentError voxelize_event(d,Inf32)
    @test_throws ArgumentError voxelize_event(transform(d,:x=>ByRow(Float64)=>:x),1f0)
    bad = copy(d); bad.energy[1] = -1f0
    @test_throws ArgumentError voxelize_event(bad,1f0)
    bad = copy(d); bad.event_id[1] = 99
    @test_throws ArgumentError voxelize_event(bad,1f0)
    @test get_voxel_size_and_distance(0.75f0) === (2.25f0,0.5f0,4.5f0)
    @test get_voxel_size_and_distance(1f0; bin_sigma_scale=2f0,edge_sigma_scale=4f0) === (2f0,0.5f0,4f0)
    @test_throws ArgumentError get_voxel_size_and_distance(0f0)
    @test_throws ArgumentError get_voxel_size_and_distance(floatmax(Float32))
end

@testset "Graph threshold, indexing and ownership" begin
    d = DataFrame(x=Float32[0,1,2,20],y=zeros(Float32,4),z=zeros(Float32,4),energy=Float32[1,0.002,0.001,3])
    for method in ("KDT","kNN","kNN_mutual")
        tracks = make_tracks(d; max_distance_mm=1f0,energy_threshold_kev=2f0,method,k=10)
        @test length(tracks) == 2
        @test [sum(t.voxels.energy) for t in tracks] == Float32[3,1.002]
        @test tracks[1].components == [[3]] # filtered parent table: x = 0,1,20
        @test tracks[2].components == [[1,2]]
        @test nv(tracks[2].graph) == 2 && ne(tracks[2].graph) == 1
        @test !is_single_track(tracks)
        tracks[2].voxels.x[1] = -99f0
        @test d.x[1] == 0f0
        @test isempty(make_tracks(d; max_distance_mm=1f0,energy_threshold_kev=4000f0,method))
        @test is_single_track(make_tracks(d[1:1,:];max_distance_mm=1f0,method))
        @test isempty(make_tracks(d[1:0,:];max_distance_mm=1f0,method))
    end
    @test !is_single_track(Tracks[])
    @test_throws ArgumentError make_tracks(d;max_distance_mm=Inf32)
    @test_throws ArgumentError make_tracks(d;max_distance_mm=1f0,k=0)
    @test_throws ArgumentError make_tracks(d;max_distance_mm=1f0,method="wrong")
    @test_throws ArgumentError make_tracks(d;max_distance_mm=1f0,energy_threshold_kev=-1f0)
    # Asymmetric nearest choices distinguish union and mutual kNN.
    d = DataFrame(x=Float32[0,1,3],y=zeros(Float32,3),z=zeros(Float32,3),energy=ones(Float32,3))
    @test length(make_tracks(d;max_distance_mm=4f0,k=1)) == 1
    @test length(make_tracks(d;max_distance_mm=4f0,k=1,method="kNN_mutual")) == 2
end

@testset "Pinned Petit graph/backbone parity" begin
    h5open(joinpath(@__DIR__,"fixtures/petit_reference.h5"),"r") do f
        for (name,case) in reconstruction_fixtures(), method in ("KDT","kNN","kNN_mutual")
            @testset "$name/$method" begin
                original = copy(case.data)
                tracks = make_tracks(case.data;max_distance_mm=case.cutoff,method,k=case.k)
                ref = f[name*"/"*method]
                @test [nv(t.graph) for t in tracks] == read(ref["component_sizes"])
                @test [sum(t.voxels.energy) for t in tracks] ≈ read(ref["component_energy_MeV"]) rtol=2f-6
                # Exact edge comparison is reserved for fixtures without distance ties.
                case.tied || @test(parent_edges(tracks) == read(ref["edges"]))
                # Independent brute-force nearest-neighbour criterion, allowing ties.
                coordinates = Matrix(case.data[:,[:x,:y,:z]])
                n = size(coordinates,1)
                distances = [norm(coordinates[i,:]-coordinates[j,:]) for i in 1:n,j in 1:n]
                kth = [sort([distances[i,j] for j in 1:n if i != j])[min(case.k,n-1)] for i in 1:n]
                actual = Set(Tuple(pair) for pair in eachcol(parent_edges(tracks)))
                for i in 1:n, j in i+1:n
                    near_i = distances[i,j] <= kth[i]
                    near_j = distances[i,j] <= kth[j]
                    allowed = distances[i,j] <= case.cutoff && (method == "KDT" ||
                        (method == "kNN" ? near_i || near_j : near_i && near_j))
                    (i,j) in actual && @test allowed
                    # Strictly nearer than both k-th neighbours must be present.
                    mandatory = distances[i,j] <= case.cutoff && (method == "KDT" ||
                        (method == "kNN" ? distances[i,j] < kth[i] || distances[i,j] < kth[j] :
                                              distances[i,j] < kth[i] && distances[i,j] < kth[j]))
                    mandatory && @test (i,j) in actual
                end
                for (i,t) in enumerate(tracks)
                    before = copy(t.voxels); graph_before = copy(t.graph)
                    r = walk_track_from_extremes(t)
                    @test r.total_length isa Float32 && r.confidence isa Float32
                    @test isfinite(r.confidence) && 0f0 <= r.confidence <= 1f0
                    @test allunique(r.path_indices)
                    @test all(has_edge(t.graph,a,b) for (a,b) in zip(r.path_indices[1:end-1],r.path_indices[2:end]))
                    @test all(norm([t.voxels.x[src(e)]-t.voxels.x[dst(e)],
                                    t.voxels.y[src(e)]-t.voxels.y[dst(e)],
                                    t.voxels.z[src(e)]-t.voxels.z[dst(e)]]) <= case.cutoff for e in edges(t.graph))
                    raw = get_raw_path(t,r.path_indices)
                    @test eltype(raw.s) == Float32
                    @test last(raw.s) ≈ r.total_length
                    @test last(get_raw_path(t,reverse(r.path_indices)).s) ≈ r.total_length
                    @test r.extremes[1] isa NamedTuple && r.extremes[2] isa NamedTuple
                    @test t.voxels == before && t.graph == graph_before
                    if !case.tied
                        rp = read(ref["$i/path"])
                        @test r.path_indices == rp || r.path_indices == reverse(rp)
                        @test r.total_length ≈ read(ref["$i/length_mm"]) rtol=3f-6
                        @test r.confidence ≈ read(ref["$i/confidence"]) atol=3f-6
                    end
                    mst = find_extremes_mst_diameter(t,EventAnalysis.extract_coords(t))
                    @test mst[5] isa Float32
                    @test all(has_edge(t.graph,a,b) for (a,b) in zip(mst[3][1:end-1],mst[3][2:end]))
                    @test t.graph == graph_before
                    r.path_voxels.x[1] = -999f0
                    @test t.voxels == before
                    @test r.extremes[1].x != -999f0
                end
                @test case.data == original
            end
        end
    end
    emptytable = DataFrame(x=Float32[],y=Float32[],z=Float32[],energy=Float32[])
    emptytrack = Tracks(emptytable,SimpleGraph(0),Vector{Int}[],DiffusionParams())
    @test isempty(walk_track_from_extremes(emptytrack).path_indices)
    @test isempty(get_raw_path(emptytrack,Int[]).s)
    one = only(make_tracks(DataFrame(x=Float32[1],y=Float32[2],z=Float32[3],energy=Float32[1]);max_distance_mm=1f0))
    @test walk_track_from_extremes(one).total_length === 0f0
    @test get_raw_path(one,[1]).s == Float32[0]
    @test_throws BoundsError get_raw_path(one,[2])
end

@testset "Time-ordered truth using actual primary IDs" begin
    # Primary 7 revisits its first cell after entering a different cell.
    h = EventHits(42,source_frame,Float32[0.1,1.2,0.3,0.7,-0.1,-1.2,99],zeros(Float32,7),zeros(Float32,7),
        Float32[1,2,3,1,0,0,10],Float32[0,2,3,1,0,1,0],Int32[7,7,7,7,19,19,44],Int32[0,1,2,3,0,1,0])
    info = EventInfo(42,(0f0,0f0,0f0),Int32[19,7])
    p = compute_mc_path(h,info,1f0)
    @test p.particle_id == Int32[7,7,7,19,19]
    @test p.x ≈ Float32[0.3,1.2,0.4,-0.1,-1.2]
    @test p.energy == Float32[3,2,2,0,0]
    @test p.s[1] === 0f0 && issorted(p.s)
    @test all(eltype(p[!,s]) == Float32 for s in (:x,:y,:z,:s,:energy))
    single = compute_mc_path(h,EventInfo(42,(0f0,0f0,0f0),Int32[7]),1f0)
    @test single.x == reverse(p.x[1:3])
    @test sum(single.energy) == sum(h.energy_MeV[h.particle_id .== 7])
    @test isempty(compute_mc_path(h,EventInfo(42,(0f0,0f0,0f0),Int32[90]),1f0))
    @test_throws ArgumentError compute_mc_path(h,EventInfo(43,(0f0,0f0,0f0),Int32[7]),1f0)
    @test_throws ArgumentError compute_mc_path(h,info,-1f0)
    @test h.x_mm[1] == 0.1f0
end
