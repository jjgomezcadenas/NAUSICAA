# Features, statistics and comparison of iItaca runs.
using Test
include(joinpath(@__DIR__, "../../scripts/iitaca/run_iitaca_batch.jl"))
include(joinpath(ROOT, "scripts/iitaca/prepare_iitaca_statistics.jl"))
include(joinpath(ROOT, "scripts/iitaca/prepare_iitaca_comparison.jl"))

@testset "Cut on E2" begin
    signal = Float64.(1:100); background = Float64.(1:2:99) ./ 2            # 50 values from 0.5 to 49.5
    t = threshold_for(signal, 0.8)
    @test accepted(signal, t).fraction == 0.8 && t == 20.5
    @test accepted(signal, threshold_for(signal, 1.0)).fraction == 1.0
    a = accepted(background, t)
    @test (a.events, a.fraction) == (29, 0.58) && a.error ≈ sqrt(0.58 * 0.42 / 50)
    @test isnan(accepted(Float64[], 1.0).fraction)
    c = curve(signal, background)
    @test first(c.efficiency) == 1 && first(c.acceptance) == 1 && last(c.efficiency) == 0 && last(c.acceptance) == 0
    @test issorted(c.efficiency; rev = true) && issorted(c.acceptance; rev = true)
    @test 0.5 < c.area < 1
    @test curve(signal, signal).area ≈ 0.5 atol = 0.01                      # no separation
    @test curve(signal .+ 1000, signal).area ≈ 1                            # complete separation
end

# The pilot read with the ITACA ion settings: the selection and the track features of the
# iItaca tables equal the itrack columns of the ITACA feature tables, event by event.
const CARD = joinpath(ROOT, "config/iitaca/runs/hybrid_equivalence_bb0nu_pilot.toml")
const REFERENCE = joinpath(ROOT, "results/readout/bb0nu_itaca_pilot_10x100_ion15/batches")
data = get(ENV, "NAUSICAA_DATA", joinpath(homedir(), "NAUSICAA"))
batches = [b for b in (0, 1) if isfile(joinpath(REFERENCE, batch_tag(b), "vh1_hits/features.h5")) &&
           isfile(joinpath(data, "pilot_10x100/bb0nu/bb0nu_$(batch_tag(b))_hitdf.h5"))]
if isempty(batches)
    @info "No pilot files on this machine: the tests on real events were skipped"
else
    ENV["NAUSICAA_DATA"] = data
    dir = mktempdir(joinpath(ROOT, "work"); prefix = "iitaca_analysis_"); run_dir = joinpath(dir, "hybrid_equivalence_bb0nu_pilot")
    for b in batches; run_iitaca_batch(CARD; batch = b, output_root = dir); end
    @testset "Feature tables" begin
        for b in batches
            folder = joinpath(run_dir, "batches", batch_tag(b)); path = write_features(folder)
            @test_throws ErrorException write_features(folder)                      # an existing table is preserved
            t = read_features(path); e = t["events"]
            old = h5open(f -> Dict(k => read(f["events"][k]) for k in keys(f["events"])), joinpath(REFERENCE, batch_tag(b), "vh1_hits/features.h5"), "r")
            @test e["event_number"] == old["event_number"]
            for (new, reference) in (("selected", "itrack_selected"), ("segments", "itrack_segments"), ("deposits", "itrack_deposits"), ("deposit_keV", "itrack_deposit_keV"),
                    ("segments_joined", "itrack_segments_joined"), ("hits", "itrack_hits"), ("components", "itrack_components"), ("main_keV", "itrack_main_keV"),
                    ("length_mm", "itrack_length_mm"), ("E1_keV", "itrack_E1_keV"), ("E2_keV", "itrack_E2_keV"), ("hit_count", "itrack_saved_hits"))
                @test isequal(e[new], old[reference])
            end
            @test e["in_roi"] == Int8.(t["roi_lower_keV"] .< e["energy_keV"] .< t["roi_upper_keV"])
            @test (t["roi_lower_keV"], t["roi_upper_keV"]) == (2457.83 - 11.5988 / 2, 2457.83 + 11.5988 / 2)
            @test t["sample"] == "bb0nu" && t["generated_events"] == 100 && all(==(b), e["batch"])
        end
    end
    @testset "Run statistics" begin
        s = run_statistics(run_dir); flow = Dict(s.flow)
        @test s.run_label == "hybrid_equivalence_bb0nu_pilot" && length(s.batches) == length(batches)
        @test flow["generated"] == 100 * length(batches)
        @test flow["generated"] >= flow["read"] >= flow["trigger"] >= flow["roi"] >= flow["track"] > 0
        @test flow["trigger"] == length(s.events["event_number"]) == s.counts["accepted_events"]
        @test flow["read"] == sum(s.counts[k] for k in keys(s.counts) if k != "input_events")       # every event read has one decision
        @test flow["track"] == count((s.events["in_roi"] .== 1) .& (s.events["selected"] .== 1))
        path = write_statistics(run_dir)
        @test h5open(f -> read(f["cut_flow/events"]), path, "r") == [n for (_, n) in s.flow]
        target = write_comparison(joinpath(dir, "comparison"), run_dir, [run_dir])                 # the sample against itself
        h5open(target, "r") do f
            @test read(f["backgrounds/01/acceptance"]) ≈ read(f["signal_efficiency"])
            @test read(f["backgrounds/01/curve/area"]) ≈ 0.5 atol = 0.02
            @test all(abs.(read(f["signal_efficiency"]) .- read(f["targets"])) .< 1 / flow["track"])
        end
    end
    # Track settings of a second variant from the same hit file (prepare_iitaca_features.jl --tracks-variant)
    @testset "Features of a tracks variant" begin
        card(name, overrides) = (p = joinpath(dir, "$name.toml"); write(p, "schema = \"nausicaa.iitaca_variant.v1\"\nname = \"$name\"\n[overrides]\n$overrides"); p)
        base = card("T0", ""); other = card("T1", "\"tracks.segment_connection_mm\" = 4.0\n\"tracks.connection_mm\" = 6.0\n")
        run = joinpath(dir, "T0_bb0nu.toml")
        write(run, replace(read(CARD, String), "run_label = \"hybrid_equivalence_bb0nu_pilot\"" => "run_label = \"T0_bb0nu\"\nvariant = \"$(relpath(base, ROOT))\""))
        run_iitaca_batch(run; batch = batches[1], output_root = dir)
        folder = joinpath(dir, "T0_bb0nu", "batches", batch_tag(batches[1])); write_features(folder)
        path = write_features(folder; tracks_variant = other)
        @test path == joinpath(dir, "T1_bb0nu", "batches", batch_tag(batches[1]), "features.h5") == features_target(folder; tracks_variant = other)
        @test islink(joinpath(dirname(path), "hits.h5")) && read(joinpath(dirname(path), "hits.h5")) == read(joinpath(folder, "hits.h5"))
        @test_throws ErrorException write_features(folder; tracks_variant = other)
        a = read_features(joinpath(folder, "features.h5")); b = read_features(path)
        @test (b["variant"], b["run_label"], b["hits_run_label"]) == ("T1", "T1_bb0nu/batches/$(batch_tag(batches[1]))", a["run_label"])
        @test b["hit_file_sha256"] == a["hit_file_sha256"] && b["events"]["event_number"] == a["events"]["event_number"]
        @test TOML.parse(b["settings"])["tracks"]["segment_connection_mm"] == 4.0 && TOML.parse(a["settings"])["tracks"]["segment_connection_mm"] == 5.5
        d = TOML.parse(a["resolved"])["detector"]; s = merge(selection_settings(d), (segment_connection_mm = 4.0, connection_mm = 6.0))
        x, y, z, e, off, cnt = h5open(f -> (read(f["itrack/x_mm"]), read(f["itrack/y_mm"]), read(f["itrack/z_mm"]), read(f["itrack/energy_keV"]),
            read(f["events/hit_offset"]), read(f["events/hit_count"])), joinpath(folder, "hits.h5"), "r")
        @test b["events"]["selected"] == [Int8(track_selection(x[k], y[k], z[k], e[k]; s...).selected) for k in (off[j]+1:off[j]+cnt[j] for j in eachindex(off))]
        @test run_statistics(joinpath(dir, "T1_bb0nu")).variant == "T1"
        # a variant that changes the readout or the saved hits is refused
        @test_throws ErrorException compute_features(joinpath(folder, "hits.h5"); tracks_variant = card("T2", "\"tracks.saved_hit_cut_keV\" = 2.0\n"))
        @test_throws ErrorException compute_features(joinpath(folder, "hits.h5"); tracks_variant = card("T3", "\"ions.pad_side_mm\" = 0.4\n"))
    end
    rm(dir; recursive = true)
end
