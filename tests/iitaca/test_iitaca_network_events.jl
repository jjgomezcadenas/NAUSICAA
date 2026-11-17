# Training files of iItaca against those of the ITACA builder: on the pilot read with
# the ITACA ion settings, the itrack cells of every event common to both are the same.
# The two chains define modules with common names, so the iItaca batches are read in a separate process.
using Test
include(joinpath(@__DIR__, "../../scripts/iitaca/build_iitaca_network_events.jl"))

const CARD = joinpath(ROOT, "config/iitaca/runs/hybrid_equivalence_bb0nu_pilot.toml")
const REFERENCE = joinpath(ROOT, "results/readout/bb0nu_itaca_pilot_10x100_ion15")
data = get(ENV, "NAUSICAA_DATA", joinpath(homedir(), "NAUSICAA"))
if !(isfile(joinpath(REFERENCE, "batches/b00000/vh1_hits/features.h5")) && isfile(joinpath(data, "pilot_10x100/bb0nu/bb0nu_b00000_hitdf.h5")))
    @info "No pilot files on this machine: the comparison of training files was skipped"
else
    ENV["NAUSICAA_DATA"] = data
    dir = mktempdir(joinpath(ROOT, "work"); prefix = "iitaca_network_")
    run_dir = joinpath(dir, "hybrid_equivalence_bb0nu_pilot")
    code = "include(\"scripts/iitaca/run_iitaca_batch.jl\"); include(\"scripts/iitaca/prepare_iitaca_features.jl\"); " *
           "for b in (0, 1); write_features(dirname(run_iitaca_batch(\"$CARD\"; batch = b, output_root = \"$dir\"))); end"
    run(pipeline(Cmd(`$(Base.julia_cmd()) --startup-file=no --project=$ROOT -e $code`; dir = ROOT); stdout = devnull))
    new = build_iitaca_network_events(run_dir, joinpath(dir, "iitaca.h5"))
    old = build_network_events(REFERENCE, joinpath(dir, "itaca.h5"))           # no ROI cut: every event passing the track selection
    cells(path) = h5open(path, "r") do f
        off = read(f["itrack/offsets"]); ids = read(f["events/event_number"]); c = [read(f["itrack/$k"]) for k in ("ix", "iy", "iz", "energy_keV")]
        (ids = ids, pitch = read(f["itrack/pitch_mm"]), cells = read(f["itrack/cells"]), side = read(f["frame_side_mm"]), has_etrack = haskey(f, "etrack"),
         event = Dict(ids[j] => sort(collect(zip((v[off[j]+1:off[j+1]] for v in c)...))) for j in eachindex(ids)),
         columns = keys(f["events"]), selection = read(f["selection"]), species = read(f["species"]))
    end
    a = cells(new); b = cells(old)
    @testset "Training files" begin
        @test a.pitch ≈ b.pitch && a.cells == b.cells && a.side == b.side
        @test !a.has_etrack && b.has_etrack && a.species == "bb0nu"
        common = intersect(a.ids, b.ids)
        @test issubset(a.ids, b.ids) && length(common) == length(a.ids) > 50     # the ROI cut only removes events
        @test all(a.event[id] == b.event[id] for id in common)
        @test all(in(a.columns), ("itrack_E2_keV", "itrack_E1_keV", "itrack_length_mm", "ion_drift_cm", "energy_keV", "batch", "event_number"))
        @test occursin("ROI cut on the event energy", a.selection)
        @test_throws ErrorException build_iitaca_network_events(run_dir, joinpath(dir, "iitaca.h5"))   # an existing file is preserved
        w = build_iitaca_network_events(run_dir, joinpath(dir, "window.h5"); window = (2420.0, 2490.0))    # box-energy window instead of the ROI
        h5open(w, "r") do f
            @test occursin("window cut on the box energy (energy_true_keV), 2420.0 to 2490.0 keV", read(f["selection"]))
            E = read(f["events/energy_true_keV"]); @test all(2420 .<= E .<= 2490) && all(==(1), read(f["events/selected"]))
            @test read(f["roi/selected_events"]) == length(a.ids)
        end
    end
    rm(dir; recursive = true)
end
