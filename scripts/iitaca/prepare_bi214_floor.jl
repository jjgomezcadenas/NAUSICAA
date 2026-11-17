#!/usr/bin/env julia
#=
iItaca II_L1000_D500: numerical results of the note on the 214Bi events that
pass the selection (notes/physics/bi214_compton_floor.tex).

  julia --project=. scripts/iitaca/prepare_bi214_floor.jl <extract folder> <output.h5>

<extract folder>: the files copied from Hyperion (work/bi214_note/out):
  roi_mc_hits.h5         Geant4 hits (detector coordinates, mm, MeV) of the
                         31,641 214Bi events in the ROI;
  display_particles.h5   NEXUS MC/particles of the ROI events of batch 301;
  b00301_hits.h5         hit file of batch 301 (iitaca.hits.v1).

Steps:
 1. Cut flow of the production (statistics/events.h5), acceptance of the E2
    cut at 80 % signal (comparison.h5) and of the network (scores_bi214.h5,
    training.h5).
 2. Monte Carlo truth classes and acceptances (truth summary.h5), and the
    distance from the start of the leading electron to the energy centre of
    the second electron (truth_bi214.h5, deposits).
 3. Network acceptance of single electrons in thirds of the drift length, for
    II_L1000_D500 and for the Xe/He ITACA run (bb0nu_itaca_100k_ion15), with
    the threshold of 80 % signal of all test events; thermal ion spread
    sigma = sqrt(2 k T L / (e E)) at the median drift of each third.
 4. Attenuation length of photons in the gas: NIST XCOM mu/rho of xenon and
    the gas density of the detector cards; energy of the scattered photon of
    a Compton scattering of the 2447.86 keV line.
 5. Voxelised Monte Carlo: the Geant4 hits of every ROI event inside its
    simulated box, summed in 3 mm cells of a lattice fixed in detector
    coordinates, cells above the saved-hit cut (3 keV); the track selection
    and track features of the production applied to the cells.
 6. Event displays: Geant4 hits of selected events labelled by particle
    (main track, second electron, fluorescence photon, other) and their
    reconstructed itrack hits.
=#
using HDF5, TOML, SHA, Statistics, Printf
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl")); using .EventAnalysis
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
include(joinpath(ROOT, "src/iitaca/IItaca.jl")); using .IItaca
include(joinpath(ROOT, "src/track_selection/TrackSelection.jl")); using .TrackSelection

const SCHEMA = "iitaca.bi214_floor.v1"
const STATS = joinpath(ROOT, "results/iitaca/II_L1000_D500_bi214/statistics/events.h5")
const COMPARISON = joinpath(ROOT, "results/iitaca/II_L1000_D500_comparison/comparison.h5")
const SCORES = joinpath(ROOT, "results/network/II_L1000_D500/trainings/itrack_sparse_e100/scores_bi214.h5")
const TRAINING = joinpath(ROOT, "results/network/II_L1000_D500/trainings/itrack_sparse_e100/training.h5")
const TRAINING_XEHE = joinpath(ROOT, "results/network/II_L1500_D400/trainings/itrack_sparse_e100/training.h5")
const TRUTH = joinpath(ROOT, "results/iitaca/II_L1000_D500_bi214_truth/truth_bi214.h5")
const SUMMARY = joinpath(ROOT, "results/iitaca/II_L1000_D500_bi214_truth/summary.h5")
const BOLTZMANN_OVER_E_V_K = 8.617333262e-5          # k/e (V/K), CODATA 2018
const ELECTRON_MASS_KEV = 510.99895                   # CODATA 2018
# Xenon (Z = 54), total attenuation with coherent scattering, mu/rho (cm^2/g):
# J. H. Hubbell and S. M. Seltzer, NIST Standard Reference Database 126,
# physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z54.html (read 2026-11-17).
const XCOM_MEV = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0]
const XCOM_MU_RHO = [2.011, 0.7202, 0.3760, 0.1797, 0.1223, 0.09699, 0.08281, 0.06696, 0.05785, 0.05054, 0.04594, 0.04078, 0.03681]
# Events shown (batch 301): passing Compton events with the second electron at
# 2.4 mm and 6.5 mm, a passing photoelectric event; a Compton event rejected
# by the track selection is chosen below.
const DISPLAY = [(30101089, "compton_lt4mm"), (30100721, "compton_4to8mm"), (30100952, "photoelectric")]

sha(path) = bytes2hex(open(sha256, path))
binomial(k, n) = (f = k / n; (f, sqrt(f * (1 - f) / n)))
loglog(x, xs, ys) = (i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1);
    exp(log(ys[i]) + (log(x) - log(xs[i])) * (log(ys[i+1]) - log(ys[i])) / (log(xs[i+1]) - log(xs[i]))))

function thirds(path, field_V_m, T_K)
    s, l, d = h5open(f -> (read(f["test/score"]), read(f["test/label"]), Float64.(read(f["test/ion_drift_cm"]))), path, "r")
    thr = quantile(s[l .== 1], 0.2); edges = quantile(d, [0, 1 / 3, 2 / 3, 1])
    third = [clamp(searchsortedlast(edges, x), 1, 3) for x in d]
    sigma(Lcm) = 1000 * sqrt(2 * BOLTZMANN_OVER_E_V_K * T_K * Lcm / 100 / field_V_m)   # mm
    k = [count((third .== t) .& (l .== 0) .& (s .> thr)) for t in 1:3]; n = [count((third .== t) .& (l .== 0)) for t in 1:3]
    ks = [count((third .== t) .& (l .== 1) .& (s .> thr)) for t in 1:3]; ns = [count((third .== t) .& (l .== 1)) for t in 1:3]
    all_k = count((l .== 0) .& (s .> thr)); all_n = count(l .== 0)
    Dict("threshold" => thr, "edges_cm" => edges, "passed" => k, "events" => n, "fraction" => k ./ n,
         "error" => [binomial(k[t], n[t])[2] for t in 1:3], "signal_efficiency" => ks ./ ns,
         "median_drift_cm" => [median(d[third .== t]) for t in 1:3], "sigma_median_mm" => [sigma(median(d[third .== t])) for t in 1:3],
         "sigma_lower_mm" => [sigma(edges[t]) for t in 1:3], "sigma_upper_mm" => [sigma(edges[t+1]) for t in 1:3],
         "overall_fraction" => all_k / all_n, "overall_error" => binomial(all_k, all_n)[2], "overall_passed" => all_k, "overall_events" => all_n,
         "field_V_m" => field_V_m, "temperature_K" => T_K, "sha256" => sha(path))
end

"""Cells of a 3 mm lattice fixed in detector coordinates: centres relative to the event position and energies (keV)."""
function voxelise(x, y, z, e_keV, vertex, box, pitch)
    cells = Dict{NTuple{3,Int},Float64}()
    for j in eachindex(x)
        r = (x[j] - vertex[1], y[j] - vertex[2], z[j] - vertex[3])
        all(box[1][a] <= r[a] <= box[2][a] for a in 1:3) || continue
        key = (floor(Int, x[j] / pitch), floor(Int, y[j] / pitch), floor(Int, z[j] / pitch))
        cells[key] = get(cells, key, 0.0) + e_keV[j]
    end
    ks = collect(keys(cells))
    c(a, i) = Float32((i + 0.5) * pitch - vertex[a])
    (Float32[c(1, k[1]) for k in ks], Float32[c(2, k[2]) for k in ks], Float32[c(3, k[3]) for k in ks], Float32[cells[k] for k in ks])
end

function voxelised_mc(extract, detector)
    sel = selection_settings(detector); feat = feature_settings(detector); saved = detector["tracks"]["saved_hit_cut_keV"]
    pitch = detector["ions"]["pad_pitch_mm"]
    st = h5open(STATS, "r") do f
        Dict(k => read(f["events/$k"]) for k in ("event_number", "in_roi", "selected", "vertex_x_mm", "vertex_y_mm", "vertex_z_mm",
            "box_x_min_mm", "box_y_min_mm", "box_z_min_mm", "box_x_max_mm", "box_y_max_mm", "box_z_max_mm", "E2_keV", "segments", "deposits", "energy_true_keV"))
    end
    row = Dict(e => i for (i, e) in enumerate(st["event_number"]))
    path = joinpath(extract, "roi_mc_hits.h5")
    h5open(path, "r") do f
        ids = read(f["events/event_number"]); first = read(f["events/first_row"]); n = read(f["events/n_hits"])
        X = read(f["hits/x_mm"]); Y = read(f["hits/y_mm"]); Z = read(f["hits/z_mm"]); E = read(f["hits/energy_MeV"])
        N = length(ids); out = Dict{String,Any}(
            "event_number" => ids, "selected" => zeros(Int8, N), "segments" => zeros(Int, N), "deposits" => zeros(Int, N),
            "cells" => zeros(Int, N), "E2_keV" => fill(NaN32, N), "box_keV" => zeros(N),
            "readout_selected" => zeros(Int8, N), "readout_segments" => zeros(Int, N), "readout_deposits" => zeros(Int, N), "readout_E2_keV" => fill(NaN32, N))
        for i in 1:N
            r = row[ids[i]]; st["in_roi"][r] == 1 || error("Event $(ids[i]) is not in the ROI")
            k = first[i]+1:first[i]+n[i]
            v = (Float64(st["vertex_x_mm"][r]), Float64(st["vertex_y_mm"][r]), Float64(st["vertex_z_mm"][r]))
            box = ((st["box_x_min_mm"][r], st["box_y_min_mm"][r], st["box_z_min_mm"][r]), (st["box_x_max_mm"][r], st["box_y_max_mm"][r], st["box_z_max_mm"][r]))
            x, y, z, e = voxelise(Float64.(X[k]), Float64.(Y[k]), Float64.(Z[k]), 1000 .* Float64.(E[k]), v, box, pitch)
            out["box_keV"][i] = sum(e; init = 0f0)
            keep = e .> saved; x = x[keep]; y = y[keep]; z = z[keep]; e = e[keep]
            s = track_selection(x, y, z, e; sel...)
            out["selected"][i] = s.selected; out["segments"][i] = s.segments; out["deposits"][i] = s.deposits; out["cells"][i] = length(e)
            s.selected && (out["E2_keV"][i] = track_features(x, y, z, e; feat...).E2_keV)
            out["readout_selected"][i] = st["selected"][r]; out["readout_segments"][i] = st["segments"][r]; out["readout_deposits"][i] = st["deposits"][r]
            st["selected"][r] == 1 && (out["readout_E2_keV"][i] = st["E2_keV"][r])
        end
        # the box energy of the cells equals the production energy in the box
        dE = maximum(abs(out["box_keV"][i] - st["energy_true_keV"][row[ids[i]]]) for i in 1:N)
        dE < 0.05 || error("Energy in the box differs from the production by $dE keV")
        out["max_box_energy_difference_keV"] = dE; out["sha256"] = sha(path); out["pitch_mm"] = pitch; out["saved_hit_cut_keV"] = saved
        out
    end
end

"""Label of every particle: 1 main track, 2 second electron (Compton electron or photoelectron produced by the gamma), 3 fluorescence photon of the gamma, 4 other."""
function particle_labels(P, ev)
    m = findall(==(ev), P["event_id"]); id = P["particle_id"][m]; mo = P["mother_id"][m]
    name = strip.(String.(P["particle_name"][m])); cre = strip.(String.(P["creator_proc"][m])); ke = 1000 .* P["kin_energy"][m]; prim = P["primary"][m]
    byid = Dict(id[i] => i for i in eachindex(id))
    g = only(i for i in eachindex(id) if prim[i] == 1 && name[i] == "gamma")
    ch = [i for i in eachindex(id) if mo[i] == id[g]]
    el = [i for i in ch if name[i] == "e-"]; lead = el[argmax([(ke[i], -id[i]) for i in el])]
    label = Dict{Int32,Int8}()
    function root_label(i)
        haskey(label, id[i]) && return label[id[i]]
        l = if i == lead; Int8(1)
        elseif mo[i] == id[g]; name[i] == "e-" && ke[i] >= 40 ? Int8(2) : name[i] == "gamma" && ke[i] < 40 ? Int8(3) : Int8(4)
        elseif i == g; Int8(4)
        else
            p = root_label(byid[mo[i]]); p == 1 && name[i] == "gamma" ? Int8(4) : p
        end
        label[id[i]] = l
    end
    for i in eachindex(id); root_label(i); end
    info = (lead_creator = cre[lead], lead_keV = ke[lead], lead_start = Tuple(Float64(P["initial_$a"][m][lead]) for a in ("x", "y", "z")),
            second = [(ke[i], Tuple(Float64(P["initial_$a"][m][i]) for a in ("x", "y", "z")), cre[i]) for i in el if i != lead && ke[i] >= 40],
            compton = count(i -> name[i] == "e-" && cre[i] == "compt", ch))
    label, info
end

function displays(extract)
    P = h5open(f -> Dict(k => read(f["particles"][k]) for k in keys(f["particles"])), joinpath(extract, "display_particles.h5"), "r")
    st = h5open(f -> Dict(k => read(f["events/$k"]) for k in ("event_number", "in_roi", "selected", "batch", "vertex_x_mm", "vertex_y_mm", "vertex_z_mm", "segments", "deposits")), STATS, "r")
    row = Dict(e => i for (i, e) in enumerate(st["event_number"]))
    mc = h5open(f -> Dict(k => read(f[k]) for k in ("events/event_number", "events/first_row", "events/n_hits", "hits/x_mm", "hits/y_mm", "hits/z_mm", "hits/energy_MeV", "hits/particle_id")),
                joinpath(extract, "roi_mc_hits.h5"), "r")
    mrow = Dict(e => i for (i, e) in enumerate(mc["events/event_number"]))
    hp = joinpath(extract, "b00301_hits.h5")
    rec = h5open(f -> Dict(k => read(f[k]) for k in ("events/event_number", "events/hit_offset", "events/hit_count", "itrack/x_mm", "itrack/y_mm", "itrack/z_mm", "itrack/energy_keV")), hp, "r")
    rrow = Dict(e => i for (i, e) in enumerate(rec["events/event_number"]))
    # a Compton event of the ROI rejected by the track selection, with the second electron farthest from the start
    cand = [(e, info) for e in unique(P["event_id"]) if haskey(row, e) && st["in_roi"][row[e]] == 1 && st["selected"][row[e]] == 0 && haskey(rrow, e)
            for (_, info) in (particle_labels(P, e),) if info.lead_creator == "compt" && length(info.second) == 1]
    dist(info) = sqrt(sum((info.second[1][2][a] - info.lead_start[a])^2 for a in 1:3))
    ok = [c for c in cand if 30 <= dist(c[2]) <= 150]
    rej = ok[argmin([dist(c[2]) for c in ok])][1]
    out = Dict{String,Any}()
    for (ev, tag) in (DISPLAY..., (rej, "compton_rejected"))
        label, info = particle_labels(P, ev); i = mrow[ev]; k = mc["events/first_row"][i]+1:mc["events/first_row"][i]+mc["events/n_hits"][i]
        v = info.lead_start; r = row[ev]; j = rrow[ev]; q = rec["events/hit_offset"][j]+1:rec["events/hit_offset"][j]+rec["events/hit_count"][j]
        d = Dict{String,Any}("event_number" => ev, "lead_creator" => info.lead_creator, "lead_keV" => info.lead_keV,
            "mc_x_mm" => mc["hits/x_mm"][k] .- Float32(v[1]), "mc_y_mm" => mc["hits/y_mm"][k] .- Float32(v[2]), "mc_z_mm" => mc["hits/z_mm"][k] .- Float32(v[3]),
            "mc_energy_keV" => 1000f0 .* mc["hits/energy_MeV"][k], "mc_label" => Int8[label[p] for p in mc["hits/particle_id"][k]],
            "reco_x_mm" => rec["itrack/x_mm"][q], "reco_y_mm" => rec["itrack/y_mm"][q], "reco_z_mm" => rec["itrack/z_mm"][q], "reco_energy_keV" => rec["itrack/energy_keV"][q],
            "readout_selected" => st["selected"][r], "readout_segments" => st["segments"][r], "readout_deposits" => st["deposits"][r],
            "vertex_offset_mm" => collect(Float64.((st["vertex_x_mm"][r], st["vertex_y_mm"][r], st["vertex_z_mm"][r])) .- collect(v)))
        if !isempty(info.second)
            s = info.second[argmax([x[1] for x in info.second])]
            d["second_keV"] = s[1]; d["second_start_mm"] = collect(s[2] .- v); d["second_start_distance_mm"] = sqrt(sum((s[2] .- v) .^ 2))
        end
        out[tag] = d
    end
    out["sha256_particles"] = sha(joinpath(extract, "display_particles.h5")); out["sha256_hits"] = sha(hp)
    out
end

function second_centre(truth)
    h5open(truth, "r") do f
        ev = Dict(k => read(f["events/$k"]) for k in ("event_number", "deposit_offset", "deposit_count", "passed", "leading_process"))
        dep = Dict(k => read(f["deposits/$k"]) for k in ("type", "deposit_in_box_keV", "centre_start_mm", "min_distance_mm"))
        n = length(ev["event_number"]); c = fill(NaN, n); e = fill(NaN, n)
        for i in 1:n, k in ev["deposit_offset"][i] .+ (1:ev["deposit_count"][i])
            dep["type"][k] in (1, 2, 6) && dep["deposit_in_box_keV"][k] >= 10 || continue
            (isnan(e[i]) || dep["deposit_in_box_keV"][k] > e[i]) && (e[i] = dep["deposit_in_box_keV"][k]; c[i] = dep["centre_start_mm"][k])
        end
        (centre_start_mm = c, leading_process = ev["leading_process"])
    end
end

"""214Bi in thirds of the drift length (tertiles of the ROI events): ROI events, selected, passing the network; Compton class; voxelised selection."""
function bi214_thirds(vox)
    st = h5open(f -> (read(f["events/ion_drift_cm"]), read(f["events/in_roi"]), read(f["events/event_number"])), STATS, "r")
    roi = st[2] .== 1; d = Float64.(st[1][roi]); edges = quantile(d, [0, 1 / 3, 2 / 3, 1])
    third(x) = clamp(searchsortedlast(edges, x), 1, 3)
    sc = h5open(f -> (Float64.(read(f["events/ion_drift_cm"])), read(f["events/score"]), read(f["threshold"])[findfirst(≈(0.8), read(f["target_signal_efficiency"]))],
                      read(f["events/event_number"])), SCORES, "r")
    cls = h5open(f -> (read(f["bi214/class"]), read(f["bi214/event_number"])), SUMMARY, "r")
    cls[2] == sc[4] || error("Event order of the scores and of the truth summary differ")
    t = third.(sc[1]); p = sc[2] .> sc[3]; c = cls[1] .== 3
    vt = third.(Float64.(st[1][[findfirst(==(e), st[3]) for e in vox["event_number"]]]))
    Dict("edges_cm" => edges, "roi" => [count(third.(d) .== i) for i in 1:3], "selected" => [count(t .== i) for i in 1:3],
         "passed" => [count((t .== i) .& p) for i in 1:3], "compton_selected" => [count((t .== i) .& c) for i in 1:3],
         "compton_passed" => [count((t .== i) .& c .& p) for i in 1:3],
         "voxelised_selected" => [count((vt .== i) .& (vox["selected"] .== 1)) for i in 1:3], "voxelised_roi" => [count(vt .== i) for i in 1:3],
         "readout_only" => [count((vt .== i) .& (vox["readout_selected"] .== 1) .& (vox["selected"] .== 0)) for i in 1:3],
         "voxelised_only" => [count((vt .== i) .& (vox["readout_selected"] .== 0) .& (vox["selected"] .== 1)) for i in 1:3])
end

"""Scaling estimate for a pad pitch p: the network acceptance of the Compton class taken as a function of
distance / pitch, with the acceptances measured at 3 mm in the bins 0-4, 4-8, > 8 mm; the selected events are kept."""
function pitch_estimate(p, p0 = 3.0)
    s = h5open(f -> Dict(k => read(f["bi214/$k"]) for k in ("class", "second_distance_mm", "passed")), SUMMARY, "r")
    acc = h5open(f -> read(f["bi214/by_second_distance/fraction"]), SUMMARY, "r")
    c = s["class"] .== 3; d = s["second_distance_mm"][c]; e = [0.0, 4, 8, Inf] .* (p / p0)
    n = [count(e[i] .<= d .< e[i+1]) for i in 1:3]
    compton = sum(n .* acc); others = count((s["passed"] .== 1) .& .!(s["class"] .== 3))
    Dict("pitch_mm" => p, "edges_mm" => e[1:3], "events" => n, "acceptance_3mm" => acc, "compton_passed" => compton,
         "compton_fraction" => compton / count(c), "total_fraction" => (compton + others) / length(c), "others_passed" => others)
end

function main(extract, output)
    resolved = h5open(f -> read(f["resolved"]), STATS, "r"); detector = TOML.parse(resolved)["detector"]
    flow = h5open(f -> (read(f["cut_flow/step"]), read(f["cut_flow/events"])), STATS, "r")
    comp = h5open(COMPARISON, "r") do f
        i = findfirst(≈(0.8), read(f["targets"]))
        (threshold = read(f["thresholds_keV"])[i], bi = (read(f["backgrounds/02/accepted_events"])[i], read(f["backgrounds/02/flow_events"])[end]),
         se = (read(f["backgrounds/01/accepted_events"])[i], read(f["backgrounds/01/flow_events"])[end]), bi_label = read(f["backgrounds/02/run_label"]), se_label = read(f["backgrounds/01/run_label"]))
    end
    sc = h5open(f -> (i = findfirst(≈(0.8), read(f["target_signal_efficiency"])); (read(f["accepted_events"])[i], length(read(f["events/score"])))), SCORES, "r")
    tr = h5open(f -> (i = findfirst(≈(0.8), read(f["test/separation/target_signal_efficiency"]));
        (round(Int, read(f["test/separation/background_acceptance"])[i] * read(f["test/separation/background_events"])), read(f["test/separation/background_events"]))), TRAINING, "r")
    gas = detector["gas"]; T = gas["temperature_K"]
    cond = IItaca.derived_conditions(detector); rho = only(c[3] for c in cond if c[1] == "gas" && c[2] == "density_kg_m3")
    E_line = detector["generation"]["bi214_gamma_keV"]
    θ = collect(range(0, π; length = 181)); Eθ = E_line ./ (1 .+ (E_line / ELECTRON_MASS_KEV) .* (1 .- cos.(θ)))
    egrid = collect(range(200.0, 2447.86; length = 200))
    λ(EkeV) = 10 / (loglog(EkeV / 1000, XCOM_MEV, XCOM_MU_RHO) * rho / 1000)        # mm
    sc2 = second_centre(TRUTH)
    vox = voxelised_mc(extract, detector)
    disp = displays(extract)
    tmp = output * ".tmp"; mkpath(dirname(output))
    try
        h5open(tmp, "w") do f
            f["schema"] = SCHEMA; f["code"] = "scripts/iitaca/prepare_bi214_floor.jl"; f["code_sha256"] = sha(@__FILE__)
            for (k, p) in (("statistics", STATS), ("comparison", COMPARISON), ("scores", SCORES), ("training", TRAINING), ("training_xehe", TRAINING_XEHE), ("truth", TRUTH), ("summary", SUMMARY))
                f["inputs/$k/path"] = relpath(p, ROOT); f["inputs/$k/sha256"] = sha(p)
            end
            f["flow/step"] = flow[1]; f["flow/events"] = flow[2]
            h5open(STATS, "r") do s
                for k in keys(s["selection/counts"]); f["flow/trigger_counts/$k"] = read(s["selection/counts/$k"]); end
            end
            f["e2/threshold_keV"] = comp.threshold; f["e2/bi214_passed"] = comp.bi[1]; f["e2/bi214_events"] = comp.bi[2]
            f["e2/single_electrons_passed"] = comp.se[1]; f["e2/single_electrons_events"] = comp.se[2]
            f["network/bi214_passed"] = sc[1]; f["network/bi214_events"] = sc[2]; f["network/single_electrons_passed"] = tr[1]; f["network/single_electrons_events"] = tr[2]
            h5open(SUMMARY, "r") do s
                for g in ("all", "by_history", "by_second_distance", "by_xray", "by_class", "by_second_energy"), k in ("events", "passed", "fraction", "error", "e2_passed")
                    f["truth/bi214/$g/$k"] = read(s["bi214/$g/$k"])
                end
                for k in ("events", "passed", "fraction", "error"); f["truth/single_electrons/all/$k"] = read(s["single_electrons/all/$k"]); end
                for k in ("expected_bi214_fraction", "expected_bi214_fraction_error", "second_distance_edges_mm", "second_energy_edges_keV"); f["truth/$k"] = read(s[k]); end
                f["truth/second_distance_mm"] = read(s["bi214/second_distance_mm"]); f["truth/second_keV"] = read(s["bi214/second_keV"])
                f["truth/class"] = read(s["bi214/class"]); f["truth/passed"] = read(s["bi214/passed"]); f["truth/event_number"] = read(s["bi214/event_number"])
            end
            f["truth/second_centre_start_mm"] = sc2.centre_start_mm; f["truth/leading_process"] = sc2.leading_process
            for (tag, p, E, TT) in (("II_L1000_D500", TRAINING, detector["drift"]["field_V_cm"] * 100, T), ("xehe_itaca", TRAINING_XEHE, 350.0 * 100, 290.0))
                for (k, v) in thirds(p, E, TT); f["thirds/$tag/$k"] = v; end
            end
            f["thirds/xehe_itaca/note"] = "bb0nu_itaca_100k_ion15 / single_electrons_itaca_100k_ion15: Xe-15% He, 16 bar, 290 K, 350 V/cm, chamber 1300 mm (config/detector at commit 4e09ccf80)"
            f["thirds/xehe_itaca/full_drift_sigma_mm"] = 1000 * sqrt(2 * BOLTZMANN_OVER_E_V_K * 290.0 * 1.3 / 35000.0)
            f["thirds/II_L1000_D500/full_drift_sigma_mm"] = 1000 * sqrt(2 * BOLTZMANN_OVER_E_V_K * T * detector["geometry"]["fiducial_length_cm"] / 100 / (detector["drift"]["field_V_cm"] * 100))
            f["attenuation/xcom_energy_MeV"] = XCOM_MEV; f["attenuation/xcom_mu_rho_cm2_g"] = XCOM_MU_RHO
            f["attenuation/source"] = "Hubbell and Seltzer, NIST Standard Reference Database 126, xenon Z = 54, mu/rho total with coherent scattering"
            f["attenuation/density_kg_m3"] = rho; f["attenuation/energy_keV"] = egrid; f["attenuation/length_mm"] = λ.(egrid)
            f["attenuation/line_keV"] = E_line; f["attenuation/line_length_mm"] = λ(E_line)
            f["attenuation/backscatter_keV"] = minimum(Eθ); f["attenuation/backscatter_length_mm"] = λ(minimum(Eθ))
            f["attenuation/compton_angle_rad"] = θ; f["attenuation/compton_scattered_keV"] = Eθ
            for EkeV in (231.0, 300.0, 400.0, 500.0)
                f["attenuation/at_$(Int(EkeV))keV_length_mm"] = λ(EkeV); f["attenuation/at_$(Int(EkeV))keV_p10mm"] = 1 - exp(-10 / λ(EkeV))
            end
            for (k, v) in vox; f["voxelised/$k"] = v; end
            for (k, v) in bi214_thirds(vox); f["bi214_thirds/$k"] = v; end
            for (k, v) in pitch_estimate(2.0); f["pitch_estimate/$k"] = v; end
            for (tag, d) in disp
                d isa AbstractDict || (f["displays/$tag"] = d; continue)
                for (k, v) in d; f["displays/$tag/$k"] = v; end
            end
        end
        Base.rename(tmp, output)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("OUTPUT $output voxelised_events=$(length(vox["event_number"])) voxelised_selected=$(sum(vox["selected"])) readout_selected=$(sum(vox["readout_selected"]))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error("Usage: prepare_bi214_floor.jl <extract folder> <output.h5>")
    main(abspath(ARGS[1]), abspath(ARGS[2]))
end
