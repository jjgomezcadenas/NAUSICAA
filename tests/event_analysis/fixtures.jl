"""Two hits from different particles, with repeated hit ID zero and source order preserved."""
function hit_fixture(; event_id=Int64(42), frame=source_frame,
        x_mm=Float32[10, 13], y_mm=Float32[20, 24], z_mm=Float32[100, 102],
        energy_MeV=Float32[0.001, 0.002], time_ns=Float32[0, 0.5],
        particle_id=Int32[1, 2], hit_id=Int32[0, 0])
    return EventHits(event_id, frame, x_mm, y_mm, z_mm,
                     energy_MeV, time_ns, particle_id, hit_id)
end

const HIT_VECTOR_FIELDS = (:x_mm, :y_mm, :z_mm, :energy_MeV, :time_ns, :particle_id, :hit_id)

"""Independent copies of every hit field, used to check caller-visible preservation."""
snapshot(hits::EventHits) = (; (f => deepcopy(getfield(hits, f)) for f in fieldnames(EventHits))...)

"""Check value preservation across all scalar and vector fields."""
function check_snapshot(hits, saved)
    for f in fieldnames(EventHits)
        @test isequal(getfield(hits, f), getproperty(saved, f))
    end
end
