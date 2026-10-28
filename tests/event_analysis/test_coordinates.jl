@testset "center_hits" begin
    info = EventInfo(Int64(42), (10f0, 20f0, 100f0), Int32[1, 2])
    hits = hit_fixture()
    saved = snapshot(hits)
    local_hits = @inferred center_hits(hits, info)
    @test local_hits isa EventHits
    @test local_hits.frame == vertex_frame
    @test local_hits.event_id == hits.event_id
    @test local_hits.x_mm == Float32[0, 3]
    @test local_hits.y_mm == Float32[0, 4]
    @test local_hits.z_mm == Float32[0, 2]
    for f in HIT_VECTOR_FIELDS[1:5]
        @test getfield(local_hits, f) isa Vector{Float32}
    end
    for f in (:energy_MeV, :time_ns, :particle_id, :hit_id)
        @test isequal(getfield(local_hits, f), getfield(hits, f))
    end
    check_snapshot(hits, saved)
    @test info.vertex_mm === (10f0, 20f0, 100f0)
    @test info.primary_particle_ids == Int32[1, 2]
    for f in HIT_VECTOR_FIELDS
        @test getfield(local_hits, f) !== getfield(hits, f)
        getfield(local_hits, f)[1] += 1
    end
    check_snapshot(hits, saved)

    # Exact integer fixtures establish invariance independently of rounding tolerances.
    shift = (-256f0, 512f0, -1024f0)
    moved = hit_fixture(x_mm=hits.x_mm .+ shift[1], y_mm=hits.y_mm .+ shift[2],
                        z_mm=hits.z_mm .+ shift[3])
    moved_info = EventInfo(info.event_id, info.vertex_mm .+ shift, Int32[1, 2])
    local_moved = center_hits(moved, moved_info)
    @test local_moved.x_mm == Float32[0, 3]
    @test local_moved.y_mm == Float32[0, 4]
    @test local_moved.z_mm == Float32[0, 2]

    # Fractional source coordinates at metre scale: compare with an independent
    # Float64 subtraction of the same stored Float32 inputs, within one output ulp.
    fractional = hit_fixture(x_mm=Float32[1500.123, -10.37],
        y_mm=Float32[-1200.23, 51.7], z_mm=Float32[1499.93, 999.17])
    fractional_info = EventInfo(Int64(42), (1499.731f0, -1199.913f0, 1499.719f0), Int32[1])
    centered = center_hits(fractional, fractional_info)
    for (axis, f) in enumerate((:x_mm, :y_mm, :z_mm))
        reference = Float64.(getfield(fractional, f)) .- Float64(fractional_info.vertex_mm[axis])
        for (actual, expected) in zip(getfield(centered, f), reference)
            @test abs(Float64(actual) - expected) <= Float64(eps(actual))
        end
    end

    # A single-primary event includes secondary deposits and preserves their order.
    electron = hit_fixture(particle_id=Int32[7, 23], hit_id=Int32[5, 0])
    single = center_hits(electron, EventInfo(Int64(42), info.vertex_mm, Int32[7]))
    @test single.particle_id == Int32[7, 23]
    @test single.energy_MeV == electron.energy_MeV
    @test single.hit_id == Int32[5, 0]

    @test_throws ArgumentError center_hits(hit_fixture(event_id=Int64(43)), info)
    @test_throws ArgumentError center_hits(hit_fixture(frame=vertex_frame), info)
    for f in HIT_VECTOR_FIELDS
        invalid = hit_fixture()
        pop!(getfield(invalid, f))
        @test_throws DimensionMismatch center_hits(invalid, info)
    end
    for f in HIT_VECTOR_FIELDS[1:5]
        invalid = hit_fixture()
        getfield(invalid, f)[1] = NaN32
        @test_throws ArgumentError center_hits(invalid, info)
    end
    invalid = hit_fixture()
    invalid.energy_MeV[1] = -1
    @test_throws ArgumentError center_hits(invalid, info)
    invalid = hit_fixture()
    invalid.particle_id[2] = 1
    @test_throws ArgumentError center_hits(invalid, info)
    empty = hit_fixture()
    foreach(f -> empty!(getfield(empty, f)), HIT_VECTOR_FIELDS)
    @test_throws ArgumentError center_hits(empty, info)
    for ids in (Int32[], Int32[7, 7])
        invalid_info = EventInfo(Int64(42), info.vertex_mm, Int32[7])
        empty!(invalid_info.primary_particle_ids)
        append!(invalid_info.primary_particle_ids, ids)
        @test_throws ArgumentError center_hits(hits, invalid_info)
    end
end
