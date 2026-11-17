using Test
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "scripts/production/prepare_batch_features.jl"))   # loads EventAnalysis, TrackSelection and the saved-hit readers

const SETTINGS = (hit_cut_keV = 5.0, segment_connection_mm = 5.5, segment_min_hits = 6, deposit_min_keV = 10.0, max_segments = 3, connection_mm = 8.0)
"Cells of a 3 mm lattice along x, from cell `first` to cell `last`, at lattice row (j, k), 30 keV each."
line(first, last; j = 0, k = 0, keV = 30f0) = (Float32[3i for i in first:last], fill(Float32(3j), last - first + 1), fill(Float32(3k), last - first + 1), fill(keV, last - first + 1))
join4(parts...) = ntuple(c -> vcat((p[c] for p in parts)...), 4)
select(hits) = track_selection(hits...; SETTINGS...)

@testset "Components" begin
    x, y, z, _ = join4(line(0, 4), line(6, 9), line(0, 2; j = 10))
    near = point_components(x, y, z, 5.5); wide = point_components(x, y, z, 8.0)
    @test near == [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3]          # one empty cell separates at 5.5 mm
    @test wide == [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2]          # and is bridged at 8 mm
    @test point_components(Float32[0, 3], Float32[0, 3], Float32[0, 3], 5.5) == [1, 1]     # cells sharing a corner touch (5.2 mm)
    @test point_components(Float32[0, 6], Float32[0, 0], Float32[0, 0], 5.5) == [1, 2]
    @test isempty(point_components(Float32[], Float32[], Float32[], 5.5))
    order = [12, 3, 7, 1, 9, 5, 2, 11, 4, 8, 6, 10]             # labels do not depend on the order of the points
    shuffled = point_components(x[order], y[order], z[order], 5.5)
    @test all((near[order[i]] == near[order[j]]) == (shuffled[i] == shuffled[j]) for i in 1:12, j in 1:12)
end

@testset "Selection" begin
    track = line(0, 19)
    @test select(track) == (segments = 1, deposits = 0, deposit_keV = 0f0, segments_joined = true, selected = true)
    # a detached deposit of 12 keV, three cells away, rejects; below 10 keV it is ignored
    xray = (Float32[30], Float32[12], Float32[0], Float32[12])
    r = select(join4(track, xray))
    @test (r.segments, r.deposits, r.deposit_keV, r.selected) == (1, 1, 12f0, false)
    @test select(join4(track, (xray[1], xray[2], xray[3], Float32[9]))).selected
    # the same deposit touching the track belongs to it
    @test select(join4(track, (Float32[30], Float32[3], Float32[0], Float32[12]))) == select(track)
    # a deposit of two hits counts with its summed energy
    @test select(join4(track, (Float32[30, 33], Float32[12, 12], Float32[0, 0], Float32[6, 6]))).deposits == 1
    # hits at or below the hit cut are not used
    @test select(join4(track, (xray[1], xray[2], xray[3], Float32[5]))).selected
    @test select((Float32[0], Float32[0], Float32[0], Float32[4])) == (segments = 0, deposits = 0, deposit_keV = 0f0, segments_joined = false, selected = false)
    # two segments with one empty cell between them join at 8 mm; with three empty cells they stay apart
    gap1 = select(join4(line(0, 9), line(11, 19)))
    @test (gap1.segments, gap1.segments_joined, gap1.selected) == (2, true, true)
    gap3 = select(join4(line(0, 9), line(13, 22)))
    @test (gap3.segments, gap3.segments_joined, gap3.selected) == (2, false, false)
    # more segments than allowed reject even when they join
    four = join4(line(0, 5), line(7, 12), line(14, 19), line(21, 26))
    @test select(four).segments == 4 && select(four).segments_joined && !select(four).selected
    @test track_selection(four...; SETTINGS..., max_segments = 4).selected
    # a short component with enough energy is a deposit, not a segment
    @test select(line(0, 4)) == (segments = 0, deposits = 1, deposit_keV = 150f0, segments_joined = false, selected = false)
end

@testset "Features of a straight track" begin
    x, y, z, e = line(0, 19)
    f = TrackSelection.track_features(x, y, z, e; hit_cut_keV = 5.0, voxel_width_mm = 3.0, connection_mm = 8.0, blob_radius_mm = 7.0)
    @test (f.hits, f.components, f.main_voxels) == (20, 1, 20)
    @test f.main_keV ≈ 600 && f.second_keV == 0
    @test f.length_mm ≈ 57 && f.extent_x_mm ≈ 60 && f.extent_y_mm ≈ 3 && f.extent_z_mm ≈ 3
    @test f.E1_keV ≈ 90 && f.E2_keV ≈ 90                        # three cells within 7 mm of each endpoint
    none = TrackSelection.track_features(x, y, z, fill(1f0, 20); hit_cut_keV = 5.0, voxel_width_mm = 3.0, connection_mm = 8.0, blob_radius_mm = 7.0)
    @test none.hits == 0 && none.components == 0 && isnan(none.E2_keV)
end

# The feature tables written before the functions moved to the module are reproduced exactly.
@testset "Saved feature tables" begin
    found = 0
    for batch in ("b00000", "b00001")
        dir = joinpath(ROOT, "results/readout/bb0nu_itaca_pilot_10x100_ion15/batches", batch)
        isfile(joinpath(dir, "vh1_hits/features.h5")) && isfile(joinpath(dir, "vh1_hits/hits.h5")) && isfile(joinpath(dir, "detector_response/controls.h5")) || continue
        found += 1
        saved = read_batch_features(joinpath(dir, "vh1_hits/features.h5")); now = compute_batch_features(dir)
        @test Set(keys(now["events"])) == Set(keys(saved["events"]))
        for (k, v) in saved["events"]
            @test isequal(now["events"][k], v)
        end
    end
    found == 0 && @info "No pilot batch on disk: the comparison with saved feature tables was skipped"
end
