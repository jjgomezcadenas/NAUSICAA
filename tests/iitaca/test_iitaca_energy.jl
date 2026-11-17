using Test, Random, Statistics
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
include(joinpath(ROOT, "src/iitaca/IItaca.jl"))
include(joinpath(ROOT, "src/track_selection/TrackSelection.jl"))
using .EventAnalysis, .IItaca, .TrackSelection

INDEX = joinpath(ROOT, "config/iitaca/iitaca.toml")
detector = resolve_detector(INDEX; root = ROOT).detector
en = detector["energy"]; Q = en["qbb_keV"]

hits(x, y, z, e; frame = mc_vertex_frame) = MCEventHits(Int64(7), frame, Float32.(x), Float32.(y), Float32.(z), Float32.(e),
    Int32.(ones(length(x))), Int32.(1:length(x)))
box(lo, hi; frame = mc_vertex_frame) = EventFrame(frame, Float32.(lo), Float32.(hi))

@testset "Energy inside the frame" begin
    h = hits([0, 10, 79.9, 80, 80.1, -200], [0, 0, 0, 0, 0, 0], [0, 5, 0, 0, 0, 0], [1000, 400, 30, 20, 10, 157])
    cube = box((-80, -80, -80), (80, 80, 80))
    @test frame_energy_keV(h, cube) == 1450.0                       # edges included; 80.1 mm and the far deposit are outside
    @test frame_energy_keV(h, box((-300, -80, -80), (300, 80, 80))) == 1617.0
    @test frame_energy_keV(h, box((200, 200, 200), (300, 300, 300))) == 0.0
    @test_throws ArgumentError frame_energy_keV(h, box((-80, -80, -80), (80, 80, 80); frame = mc_source_frame))
    # the sum is carried in Float64: 10^5 hits of 24.5783 keV
    many = hits(zeros(100_000), zeros(100_000), zeros(100_000), fill(24.5783f0, 100_000))
    @test frame_energy_keV(many, cube) ≈ 100_000 * Float64(24.5783f0) rtol = 1e-12
end

@testset "Emulated energy" begin
    sigma = en["resolution_fwhm_at_qbb_keV"] / FWHM_PER_SIGMA
    @test sigma ≈ 4.9255 atol = 1e-4
    draw(E, n; seed = 1) = (rng = Xoshiro(seed); [emulated_energy_keV(rng, E, en) for _ in 1:n])
    n = 400_000; x = draw(Q, n)
    @test mean(x) ≈ Q atol = 5 * sigma / sqrt(n)                    # unbiased
    @test std(x) ≈ sigma rtol = 5 / sqrt(2n)                        # width of the card
    @test FWHM_PER_SIGMA * std(x) / Q ≈ 0.004719 rtol = 0.01
    @test std(draw(Q / 4, n)) ≈ sigma / 2 rtol = 5 / sqrt(2n)       # sqrt(E)
    @test std(draw(1592.5, n)) ≈ sigma * sqrt(1592.5 / Q) rtol = 5 / sqrt(2n)
    @test draw(Q, 5) == draw(Q, 5) && draw(Q, 5) != draw(Q, 5; seed = 2)     # reproducible with the seed
    @test emulated_energy_keV(Xoshiro(1), 0.0, en) == 0.0
    # one Gaussian deviate per call: the stream advances as randn does
    r1 = Xoshiro(3); emulated_energy_keV(r1, Q, en); r2 = Xoshiro(3); randn(r2)
    @test rand(r1) == rand(r2)
    # fraction inside the ROI of one FWHM for a line at Qbb: erf(sqrt(ln 2)) = 0.7610
    @test count(E -> in_roi(en, E), x) / n ≈ 0.7610 atol = 5 * sqrt(0.761 * 0.239 / n)
    # the 214Bi line at 2447.86 keV: 4.17 keV below the lower edge of the ROI
    bi = draw(2447.86, n)
    lower = (roi_keV(en).lower - 2447.86) / (sigma * sqrt(2447.86 / Q)); upper = (roi_keV(en).upper - 2447.86) / (sigma * sqrt(2447.86 / Q))
    expected = count(t -> lower < t < upper, randn(Xoshiro(9), n)) / n     # fraction of a unit Gaussian between the two bounds, about 0.20
    @test count(E -> in_roi(en, E), bi) / n ≈ expected atol = 7 * sqrt(expected * (1 - expected) / n)
    @test 0.15 < expected < 0.25
    @test !in_roi(en, roi_keV(en).lower) && !in_roi(en, roi_keV(en).upper) && in_roi(en, Q)      # strict bounds
    wide = resolve_detector(INDEX; overrides = Dict("energy.resolution_fwhm_at_qbb_keV" => 24.5783)).detector["energy"]
    @test std([emulated_energy_keV(r, Q, wide) for r in (Xoshiro(4),) for _ in 1:n]) ≈ 24.5783 / FWHM_PER_SIGMA rtol = 5 / sqrt(2n)
end

@testset "Track settings from the card" begin
    s = selection_settings(detector); f = feature_settings(detector)
    @test s == (hit_cut_keV = 5.0, segment_connection_mm = 5.5, segment_min_hits = 6, deposit_min_keV = 10.0, max_segments = 3, connection_mm = 8.0)
    @test f == (hit_cut_keV = 5.0, voxel_width_mm = 3.0, connection_mm = 8.0, blob_radius_mm = 7.0)
    x = Float32[3i for i in 0:19]; o = zeros(Float32, 20); e = fill(30f0, 20)
    @test track_selection(x, o, o, e; s...).selected
    @test track_features(x, o, o, e; f...).E2_keV ≈ 90
    loose = selection_settings(resolve_detector(INDEX; overrides = Dict("tracks.deposit_min_keV" => 20.0)).detector)
    xray = (vcat(x, 30f0), vcat(o, 12f0), vcat(o, 0f0), vcat(e, 12f0))
    @test !track_selection(xray...; s...).selected && track_selection(xray...; loose...).selected
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("tracks.saved_hit_cut_keV" => 6.0))
    @test_throws ArgumentError resolve_detector(INDEX; overrides = Dict("tracks.segment_connection_mm" => 9.0))
end
