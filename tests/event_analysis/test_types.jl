@testset "CoordinateFrame" begin
    @test instances(CoordinateFrame) == (source_frame, vertex_frame)
    @test source_frame isa CoordinateFrame
    @test vertex_frame isa CoordinateFrame
    @test source_frame != vertex_frame
end

@testset "EventInfo" begin
    ids = Int32[8, 3]
    info = EventInfo(Int64(42), (10f0, 20f0, 100f0), ids)
    @test ids == Int32[8, 3]
    @test info.primary_particle_ids == Int32[3, 8]
    @test info.primary_particle_ids !== ids
    @test info.event_id === Int64(42)
    @test info.vertex_mm === (10f0, 20f0, 100f0)
    ids[1] = 99
    @test info.primary_particle_ids == Int32[3, 8]
    for event_id in (Int64(0), Int64(1007))
        electron = EventInfo(event_id, (0f0, -5f0, 1500f0), Int32[7])
        @test electron.event_id == event_id
        @test electron.primary_particle_ids == Int32[7]
    end
    @test_throws ArgumentError EventInfo(Int64(0), (0f0, 0f0, 0f0), Int32[])
    @test_throws ArgumentError EventInfo(Int64(0), (0f0, 0f0, 0f0), Int32[7, 7])
    for bad in (NaN32, Inf32, -Inf32), axis in 1:3
        vertex = ntuple(i -> i == axis ? bad : 0f0, 3)
        @test_throws ArgumentError EventInfo(Int64(0), vertex, Int32[7])
    end
end

@testset "EventHits" begin
    hits = hit_fixture()
    @test hits.frame == source_frame
    @test hits.hit_id == Int32[0, 0]
    @test hits.particle_id == Int32[1, 2]
    @test hits.x_mm == Float32[10, 13]
    for f in HIT_VECTOR_FIELDS[1:5]
        @test getfield(hits, f) isa Vector{Float32}
    end
    @test hits.particle_id isa Vector{Int32}
    @test hits.hit_id isa Vector{Int32}
    @test hit_fixture(energy_MeV=Float32[0, 1]).energy_MeV[1] === 0f0
    @test hit_fixture(frame=vertex_frame).frame == vertex_frame
    @test hit_fixture(particle_id=Int32[9, 9], hit_id=Int32[4, 1]).hit_id == Int32[4, 1]
    @test_throws ArgumentError hit_fixture(particle_id=Int32[9, 9])

    # Identifiers can repeat in other events; validation is scoped to one event.
    @test hit_fixture(event_id=Int64(7)).hit_id == hit_fixture(event_id=Int64(100)).hit_id
    @test_throws ArgumentError hit_fixture(x_mm=Float32[], y_mm=Float32[], z_mm=Float32[],
        energy_MeV=Float32[], time_ns=Float32[], particle_id=Int32[], hit_id=Int32[])
    for f in HIT_VECTOR_FIELDS
        shorter = f in (:particle_id, :hit_id) ? Int32[1] : Float32[1]
        @test_throws DimensionMismatch hit_fixture(; (f => shorter,)...)
    end
    for f in HIT_VECTOR_FIELDS[1:5], bad in (NaN32, Inf32, -Inf32)
        @test_throws ArgumentError hit_fixture(; (f => Float32[0, bad],)...)
    end
    @test_throws ArgumentError hit_fixture(energy_MeV=Float32[-0.01, 1])

    # Each supplied array is copied, including when positions share a source array.
    shared = Float32[1, 2]
    energies = Float32[0.1, 0.2]
    times = Float32[0, 1]
    particles = Int32[1, 2]
    ids = Int32[0, 0]
    owned = hit_fixture(x_mm=shared, y_mm=shared, z_mm=shared,
                        energy_MeV=energies, time_ns=times, particle_id=particles, hit_id=ids)
    before = snapshot(owned)
    shared[1] = 200
    energies[1] = 50
    times[1] = 60
    particles[1] = 99
    ids[1] = 100
    check_snapshot(owned, before)
    owned.x_mm[1] = -5
    @test owned.y_mm == Float32[1, 2]
    @test owned.z_mm == Float32[1, 2]
end
