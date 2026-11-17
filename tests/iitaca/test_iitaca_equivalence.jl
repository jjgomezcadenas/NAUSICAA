# The iItaca chain against the ITACA production on the same Monte Carlo events.
# With the ion settings of the ITACA pilot as overrides, the iItaca driver must
# select the same events and write the same itrack hits, bit for bit: the two
# chains use the same modules and the same random streams per event.
# Needs the pilot NEXUS files (NAUSICAA_DATA) and the ITACA hit files of
# results/readout/bb0nu_itaca_pilot_10x100_ion15; skipped when they are absent.
using Test
include(joinpath(@__DIR__, "../../scripts/iitaca/run_iitaca_batch.jl"))

const CARD = joinpath(ROOT, "config/iitaca/runs/hybrid_equivalence_bb0nu_pilot.toml")
const REFERENCE = joinpath(ROOT, "results/readout/bb0nu_itaca_pilot_10x100_ion15/batches")
data = get(ENV, "NAUSICAA_DATA", joinpath(homedir(), "NAUSICAA"))
batches = [b for b in (0, 1) if isfile(joinpath(REFERENCE, batch_tag(b), "vh1_hits/hits.h5")) &&
           isfile(joinpath(data, "pilot_10x100/bb0nu/bb0nu_$(batch_tag(b))_hitdf.h5"))]

if isempty(batches)
    @info "No pilot files on this machine: the comparison with the ITACA production was skipped"
else
    ENV["NAUSICAA_DATA"] = data
    for b in batches
        @testset "Batch $(batch_tag(b))" begin
            dir = mktempdir(joinpath(ROOT, "work"); prefix = "iitaca_equivalence_")
            path = run_iitaca_batch(CARD; batch = b, output_root = dir)
            new = h5open(path, "r") do f
                (ids = read(f["events/event_number"]), offset = read(f["events/hit_offset"]), count = read(f["events/hit_count"]),
                 hits = Dict(k => read(f["itrack/$k"]) for k in ("x_mm", "y_mm", "z_mm", "energy_keV")),
                 E_true = read(f["events/energy_true_keV"]), E = read(f["events/energy_keV"]), ion = read(f["events/ion_energy_keV"]),
                 drift = read(f["events/ion_drift_cm"]), decisions = read(f["selection/events/decision"]), all_ids = read(f["selection/events/event_number"]),
                 accepted = read(f["selection/counts/accepted_events"]), schema = read(f["schema"]), sim = keys(f["configuration/simulation_inputs"]))
            end
            @test new.schema == "iitaca.hits.v1"
            @test !any(occursin("documentation", k) for k in new.sim)             # no documentation value reached the simulation
            @test new.accepted == length(new.ids) == count(==(0), new.decisions)
            old = open_hit_input(joinpath(REFERENCE, batch_tag(b), "vh1_hits/hits.h5"))
            try
                old_ids = read(old.file["events/event_number"])
                @test new.ids == old_ids                                          # same accepted events, same order
                same_hits = 0; total_hits = 0
                for (j, id) in enumerate(new.ids)
                    ev = read_hit_event(old, id); range = new.offset[j]+1:new.offset[j]+new.count[j]
                    @test new.hits["x_mm"][range] == ev.itrack.x_mm
                    @test new.hits["y_mm"][range] == ev.itrack.y_mm
                    @test new.hits["z_mm"][range] == ev.itrack.z_mm
                    @test new.hits["energy_keV"][range] == ev.itrack.energy_keV
                    @test Float32(new.ion[j]) === ev.energy.ion_energy_keV || (isnan(new.ion[j]) && isnan(ev.energy.ion_energy_keV))
                    @test new.drift[j] == ev.info.ion_drift_cm
                    total_hits += length(ev.itrack.x_mm)
                end
                @info "Batch $(batch_tag(b)): $(length(new.ids)) events, $total_hits itrack hits compared"
            finally
                close_hit_input(old)
            end
            # the energy inside the box is the one recorded by the ITACA detector-response stage
            h5open(joinpath(REFERENCE, batch_tag(b), "detector_response/controls.h5"), "r") do f
                @test read(f["events/event_number"]) == new.ids
                @test new.E_true ≈ Float64.(read(f["events/energy_inframe_keV"])) rtol = 1e-6
            end
            @test all(isfinite, new.E) && new.E != new.E_true
            rm(dir; recursive = true)
        end
    end
end
