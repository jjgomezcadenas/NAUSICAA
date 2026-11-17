#!/usr/bin/env julia
#=
iItaca: the same Monte Carlo events read on 3 mm and on 2 mm pads.

  julia --project=. scripts/iitaca/prepare_iitaca_pitch_comparison.jl <output folder> <settings> ...

Each <settings> is a run prefix: II_L1000_D500 (3 mm), II_L1000_D500_P2 (2 mm,
radii scaled to the lattice), II_L1000_D500_P2_R3 (2 mm, radii of the 3 mm
selection). For a prefix P it reads

  results/iitaca/P_comparison/comparison.h5             cut flow, E₂ thresholds and acceptances
  results/network/P/trainings/itrack_sparse_e100/       training.h5 (single electrons of the test sample),
                                                        scores_bi214.h5, scores_tl208.h5
  results/iitaca/P_bi214_truth/                         truth_bi214.h5 (interaction history),
                                                        summary.h5 (second electron)

and writes <output folder>/comparison.h5 (schema iitaca.pitch_comparison.v1):

 1. Cut flow of every sample (generated, read, trigger, ROI, track selection).
 2. At 90, 80, 70 and 60 % signal: E₂ threshold and acceptance of every
    background; network threshold, signal efficiency on the test sample and
    acceptance of single electrons (test sample), ²¹⁴Bi and ²⁰⁸Tl.
 3. ²¹⁴Bi by truth group, with the E₂ and network acceptances at each target:
    interaction history (photoelectric absorption only; Compton scattering,
    then absorption or escape; other) and distance of the second electron to
    the main track (below 4 mm, 4 to 8 mm, beyond 8 mm, no second electron of
    at least 10 keV in the box). Definitions: summarize_iitaca_truth.jl.
 4. Events selected in both of two settings (²¹⁴Bi): count, for the page.

Errors are binomial, sqrt(max(k, 1) (1 - k/n)) / n, as in the truth summary.
=#
using HDF5, Statistics
const ROOT = normpath(joinpath(@__DIR__, "../.."))
const PITCH_SCHEMA = "iitaca.pitch_comparison.v1"
const TARGETS = [0.9, 0.8, 0.7, 0.6]
const SAMPLES = ("bb0nu", "single_electrons", "bi214", "tl208")
const HISTORY_GROUPS = ["photoelectric", "compton", "other"]            # from truth history 1; 2 and 3; the rest
const DISTANCE_GROUPS = ["lt4mm", "4to8mm", "gt8mm", "no_second_electron"]
const TRAINING = "trainings/itrack_sparse_e100"

fraction(k, n) = n > 0 ? k / n : NaN
binomial_error(k, n) = n > 0 ? sqrt(max(k, 1) * (1 - k / n)) / n : NaN
index_of(targets, t) = (i = findfirst(x -> isapprox(x, t), targets); i === nothing ? error("No target $t") : i)

"""Cut flow, E₂ thresholds and background acceptances of the comparison of one setting."""
function read_comparison(path)
    h5open(path, "r") do f
        flow = Dict(read(f["signal/sample"]) => (read(f["signal/flow_step"]), read(f["signal/flow_events"])))
        e2 = Dict{String,Any}(); targets = read(f["targets"]); thr = read(f["thresholds_keV"])
        for k in keys(f["backgrounds"])
            g = f["backgrounds/$k"]; s = read(g["sample"]); flow[s] = (read(g["flow_step"]), read(g["flow_events"]))
            e2[s] = (events = read(g["accepted_events"]), acceptance = read(g["acceptance"]), error = read(g["acceptance_error"]),
                     E2 = read(g["E2_keV"]), event_number = read(g["event_number"]))
        end
        i = [index_of(targets, t) for t in TARGETS]
        (flow = flow, threshold = thr[i], efficiency = read(f["signal_efficiency"])[i],
         e2 = Dict(s => (events = v.events[i], acceptance = v.acceptance[i], error = v.error[i], E2 = v.E2, event_number = v.event_number) for (s, v) in e2),
         settings = read(f["settings"]), variant = read(f["variant"]))
    end
end

"""Network: thresholds and acceptances at TARGETS; scores of ²¹⁴Bi per event."""
function read_network(dir)
    t = h5open(joinpath(dir, "training.h5"), "r") do f
        s = f["test/separation"]; i = [index_of(read(s["target_signal_efficiency"]), x) for x in TARGETS]
        (threshold = read(s["threshold"])[i], efficiency = read(s["signal_efficiency"])[i], acceptance = read(s["background_acceptance"])[i],
         error = read(s["background_acceptance_error"])[i], events = read(s["background_events"]), signal_events = read(s["signal_events"]),
         area = read(s["area"]), best_epoch = read(f["settings/best_epoch"]))
    end
    scores = Dict{String,Any}()
    for s in ("bi214", "tl208")
        scores[s] = h5open(joinpath(dir, "scores_$s.h5"), "r") do f
            i = [index_of(read(f["target_signal_efficiency"]), x) for x in TARGETS]
            read(f["threshold"])[i] ≈ t.threshold || error("Thresholds of scores_$s.h5 and training.h5 differ in $dir")
            (accepted = read(f["accepted_events"])[i], acceptance = read(f["acceptance"])[i], error = read(f["acceptance_error"])[i],
             event_number = read(f["events/event_number"]), score = read(f["events/score"]))
        end
    end
    (training = t, scores = scores)
end

"""History group and distance group (indices into HISTORY_GROUPS, DISTANCE_GROUPS) of the ²¹⁴Bi events of a truth folder."""
function read_truth(folder)
    ids, history = h5open(f -> (read(f["events/event_number"]), read(f["events/history"])), joinpath(folder, "truth_bi214.h5"), "r")
    ids2, distance = h5open(f -> (read(f["bi214/event_number"]), read(f["bi214/second_distance_mm"])), joinpath(folder, "summary.h5"), "r")
    ids == ids2 || error("Events of truth_bi214.h5 and summary.h5 differ in $folder")
    hg = [h == 1 ? 1 : h in (2, 3) ? 2 : 3 for h in history]
    dg = [isnan(d) ? 4 : d < 4 ? 1 : d < 8 ? 2 : 3 for d in distance]
    Dict(ids[i] => (hg[i], dg[i]) for i in eachindex(ids))
end

"""Events, E₂ passes and network passes at each target, per group (rows) and target (columns)."""
function grouped(groups, m, e2pass, netpass)
    n = [count(==(g), groups) for g in 1:m]
    e2 = [count(i -> groups[i] == g && e2pass[i, j], eachindex(groups)) for g in 1:m, j in eachindex(TARGETS)]
    net = [count(i -> groups[i] == g && netpass[i, j], eachindex(groups)) for g in 1:m, j in eachindex(TARGETS)]
    (events = n, e2_passed = e2, net_passed = net)
end

function prepare(out::AbstractString, prefixes::Vector{String})
    mkpath(out); target = joinpath(out, "comparison.h5"); isfile(target) && error("Preserve existing comparison: $target")
    tmp = target * ".tmp"; selected = Dict{String,Set{Int}}()
    try
        h5open(tmp, "w") do f
            f["schema"] = PITCH_SCHEMA; f["settings"] = prefixes; f["targets"] = TARGETS
            f["samples"] = collect(SAMPLES); f["history_groups"] = HISTORY_GROUPS; f["distance_groups"] = DISTANCE_GROUPS
            for p in prefixes
                c = read_comparison(joinpath(ROOT, "results/iitaca", "$(p)_comparison", "comparison.h5"))
                c.variant == p || error("Comparison of $p holds variant $(c.variant)")
                n = read_network(joinpath(ROOT, "results/network", p, TRAINING)); truth = read_truth(joinpath(ROOT, "results/iitaca", "$(p)_bi214_truth"))
                g = create_group(f, p); g["track_settings"] = c.settings
                for s in SAMPLES; g["flow/$s/step"] = c.flow[s][1]; g["flow/$s/events"] = c.flow[s][2]; end
                g["e2/threshold_keV"] = c.threshold; g["e2/signal_efficiency"] = c.efficiency
                for (s, v) in c.e2; g["e2/$s/accepted"] = v.events; g["e2/$s/acceptance"] = v.acceptance; g["e2/$s/error"] = v.error; end
                t = n.training
                g["network/threshold"] = t.threshold; g["network/signal_efficiency"] = t.efficiency; g["network/area"] = t.area
                g["network/signal_test_events"] = t.signal_events; g["network/best_epoch"] = t.best_epoch
                g["network/single_electrons/events"] = t.events; g["network/single_electrons/acceptance"] = t.acceptance; g["network/single_electrons/error"] = t.error
                g["network/single_electrons/accepted"] = round.(Int, t.acceptance .* t.events)
                for (s, v) in n.scores
                    g["network/$s/events"] = length(v.event_number); g["network/$s/accepted"] = v.accepted; g["network/$s/acceptance"] = v.acceptance; g["network/$s/error"] = v.error
                end
                # ²¹⁴Bi by truth group: the events of the comparison, their E₂ and their scores
                b = c.e2["bi214"]; sc = n.scores["bi214"]
                sort(b.event_number) == sort(sc.event_number) || error("Selected ²¹⁴Bi events of the comparison and of the scores differ for $p")
                Set(keys(truth)) == Set(b.event_number) || error("Truth table of $p does not hold the selected ²¹⁴Bi events")
                score = Dict(zip(sc.event_number, sc.score)); ids = b.event_number; selected[p] = Set(ids)
                e2pass = [b.E2[i] > c.threshold[j] for i in eachindex(ids), j in eachindex(TARGETS)]
                netpass = [score[ids[i]] > t.threshold[j] for i in eachindex(ids), j in eachindex(TARGETS)]
                for (name, groups, m) in (("history", [truth[id][1] for id in ids], length(HISTORY_GROUPS)), ("distance", [truth[id][2] for id in ids], length(DISTANCE_GROUPS)))
                    r = grouped(groups, m, e2pass, netpass)
                    g["bi214_$name/events"] = r.events; g["bi214_$name/e2_passed"] = r.e2_passed; g["bi214_$name/net_passed"] = r.net_passed
                    g["bi214_$name/e2_fraction"] = fraction.(r.e2_passed, r.events); g["bi214_$name/e2_error"] = binomial_error.(r.e2_passed, r.events)
                    g["bi214_$name/net_fraction"] = fraction.(r.net_passed, r.events); g["bi214_$name/net_error"] = binomial_error.(r.net_passed, r.events)
                end
                println("$p: bi214 selected $(length(ids)), network at 80 % $(round(100 * n.scores["bi214"].acceptance[2]; digits = 2)) %")
            end
            for (i, p) in enumerate(prefixes), q in prefixes[i+1:end]
                f["bi214_common/$(p)__$(q)"] = length(intersect(selected[p], selected[q]))
            end
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("PITCH COMPARISON $target")
    target
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("Usage: prepare_iitaca_pitch_comparison.jl <output folder> <run prefix> ...")
    prepare(abspath(ARGS[1]), String.(ARGS[2:end]))
end
