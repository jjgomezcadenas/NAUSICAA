#!/usr/bin/env julia
#=
eITACA: Monte Carlo truth of the selected events of a run, against the network score.

  julia --project=. scripts/eitaca/prepare_eitaca_truth.jl <run card> <scores> <output.h5> [--batches=A:B]

As scripts/iitaca/prepare_iitaca_truth.jl (its functions event_truth,
read_truth_particles and scored_events), on the statistics of the eITACA run,
results/eitaca/<run_label>/statistics/events.h5. NAUSICAA_DATA must point to
the folder holding IE_L1500_D400/prod/ (NEXUS files and reduced hit tables).

The electron file does not hold the simulated box, so step 3 of the iItaca
script takes it from the production code, for every scored event:

 1. Fast trigger and frame (EventAnalysis.mc_fast_trigger) with the resolved
    trigger and geometry of the run (statistics/events.h5, resolved), on the
    hits centred on the start of the main track: for trigger.frame_track =
    "leading_gamma_electron" the leading electron of the gamma, as
    run_eitaca_batch.jl. The event must be accepted.
 2. Box: the main track and the gamma deposits inside the frame, with the
    margin of the card (EventAnalysis.simulation_box), relative to that start.

The energy in the box recomputed from the hits must equal energy_true_keV of
the production for every event. Output: schema iitaca.truth.v1.
=#
include(joinpath(@__DIR__, "../iitaca/prepare_iitaca_truth.jl"))

const EITACA_COLUMNS = ("energy_keV", "energy_true_keV", "electron_energy_keV", "E1_keV", "E2_keV", "length_mm", "electron_drift_cm", "segments",
    "segments_joined", "deposits", "components", "main_keV", "second_keV", "deposit_keV", "hit_count")

"""Feature columns of the eITACA run, indexed by event number."""
function eitaca_columns(statistics::AbstractString, ids)
    h5open(statistics, "r") do f
        row = Dict(id => i for (i, id) in enumerate(read(f["events/event_number"]))); r = [row[id] for id in ids]
        Dict(n => read(f["events/$n"])[r] for n in EITACA_COLUMNS)
    end
end

"""Fast-trigger settings from the resolved detector, as run_eitaca_batch.jl."""
trigger_settings(d) = (g = d["geometry"]; t = d["trigger"]; FixedFrameTriggerSettings(Float32(10g["fiducial_length_cm"]), Float32(5g["fiducial_diameter_cm"]),
    Float32(t["frame_side_mm"]), Float32(t["gamma_energy_cut_keV"]), t["frame_placement"] == "inside_detector", Float32(t["el_exclusion_mm"])))

function prepare_eitaca_truth(card::AbstractString, scores::AbstractString, output::AbstractString; batches = nothing)
    c = TOML.parsefile(card); label = c["run_label"]
    statistics = joinpath(ROOT, "results/eitaca", label, "statistics/events.h5")
    d, gamma_event = h5open(statistics, "r") do f
        TOML.parse(read(f["resolved"]))["detector"], read(f["frame_track"]) == "leading_gamma_electron"
    end
    trigger = trigger_settings(d); L = trigger.chamber_length_mm; margin = Float32(d["trigger"]["simulation_margin_mm"])
    s = scored_events(scores, statistics)
    keep = batches === nothing ? trues(length(s.event_number)) : in.(s.batch, Ref(batches))
    ids = s.event_number[keep]; bs = s.batch[keep]; score = s.score[keep]
    feats = eitaca_columns(statistics, ids); order = Dict(id => i for (i, id) in enumerate(ids))
    rows = Vector{Any}(undef, length(ids)); deps = Vector{Any}(undef, length(ids)); identity = Dict{String,String}()
    started = time()
    for b in sort(unique(bs))
        tag = batch_tag(b); bid = ids[bs .== b]
        hp = input_path(ROOT, replace(c["hits"], "{batch}" => tag)); pp = input_path(ROOT, replace(c["particles"], "{batch}" => tag))
        truth = read_truth_particles(pp, bid); parts = read_mc_particles(pp, bid)
        leading = gamma_event ? read_leading_gamma_electrons(pp, bid) : nothing
        h5open(hp, "r") do f
            identity[tag] = read(f["metadata/source_sha256"])
            row = Dict(id => i for (i, id) in enumerate(read(f["events/event_id"])))
            for id in bid
                h, info = read_mc_event(f, id; event_row = row[id]); i = order[id]
                start = leading === nothing ? info : leading_electron_info(info, leading[id])
                hc = center_hits(h, start); o = mc_fast_trigger(hc, start, parts[id], placement_from_vertex(start, L), trigger)
                o.decision == trigger_accepted || error("Event $id of $tag is not accepted by the fast trigger")
                box = simulation_box(hc, start, parts[id], o.frame, margin)
                rows[i], deps[i] = event_truth(h, info, truth[id], (Tuple(box.lower_mm), Tuple(box.upper_mm)), gamma_event)
            end
        end
        println("Batch $tag: $(length(bid)) events; $(round(time() - started; digits = 1)) s"); flush(stdout)
    end
    dE = maximum(abs(rows[i]["box_keV_truth"] - feats["energy_true_keV"][i]) for i in eachindex(ids))
    dE < 0.05 || error("Energy in the box differs from the production by $dE keV")
    mkpath(dirname(output)); isfile(output) && error("Preserve existing output: $output"); tmp = output * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = TRUTH_SCHEMA; f["run_label"] = label; f["sample"] = c["sample"]
            f["scores_file"] = basename(scores); f["scores_sha256"] = bytes2hex(open(sha256, scores))
            f["statistics_sha256"] = bytes2hex(open(sha256, statistics)); f["threshold_80"] = s.threshold
            f["selection"] = "selected events of $(label) (ROI and track selection) with a network score; passed: score above the 80 % signal threshold"
            f["code_sha256"] = bytes2hex(open(sha256, @__FILE__)); f["code"] = "scripts/eitaca/prepare_eitaca_truth.jl"
            f["iitaca_code_sha256"] = bytes2hex(open(sha256, joinpath(ROOT, "scripts/iitaca/prepare_iitaca_truth.jl")))
            f["batches"] = sort(collect(keys(identity))); f["particle_file_sha256"] = [identity[k] for k in sort(collect(keys(identity)))]
            f["deposit_types"] = DEPOSIT_TYPES; f["histories"] = HISTORIES; f["leading_processes"] = LEADING_PROCESSES
            f["blob_radius_mm"] = BLOB_RADIUS_MM; f["max_box_energy_difference_keV"] = dE
            f["events/event_number"] = ids; f["events/batch"] = bs; f["events/score"] = score; f["events/passed"] = Int8.(score .> s.threshold)
            for (k, v) in feats; f["events/$k"] = v; end
            for k in sort(collect(keys(rows[1]))); f["events/$k"] = [r[k] for r in rows]; end
            n = [length(x) for x in deps]; f["events/deposit_offset"] = cumsum(n) .- n; f["events/deposit_count"] = n
            flat = reduce(vcat, deps); f["deposits/event_number"] = reduce(vcat, [fill(ids[i], n[i]) for i in eachindex(ids)])
            for k in keys(first(flat)); f["deposits/$k"] = k == :in_box ? Int8[x[k] for x in flat] : [x[k] for x in flat]; end
        end
        Base.rename(tmp, output)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("OUTPUT $output events=$(length(ids)) passed=$(count(score .> s.threshold)) max_box_dE_keV=$dE time_s=$(round(time() - started; digits = 1))")
    output
end

if abspath(PROGRAM_FILE) == @__FILE__
    k = findall(startswith("--batches="), ARGS); args = ARGS[setdiff(eachindex(ARGS), k)]
    batches = isempty(k) ? nothing : (r = parse.(Int, split(split(ARGS[k[1]], "=")[2], ":")); r[1]:r[end])
    length(args) == 3 || error("Usage: prepare_eitaca_truth.jl <run card> <scores> <output.h5> [--batches=A:B]")
    prepare_eitaca_truth(abspath(args[1]), abspath(args[2]), abspath(args[3]); batches)
end
