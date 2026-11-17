#!/usr/bin/env julia
#=
iItaca: Monte Carlo truth of the selected events of a run, against the network score.

  julia --project=. scripts/iitaca/prepare_iitaca_truth.jl <run card> <scores> <output.h5> [--batches=A:B]

<run card>: config/iitaca/runs/II_L1000_D500_bi214.toml (or _single_electrons).
<scores>: the scores of a trained network, either scores_<sample>.h5 (every
selected event of the sample) or training.h5 (its test background events).
NAUSICAA_DATA must point to the folder holding <variant>/prod/ (NEXUS files and
reduced hit tables). --batches limits the batches read, for tests.

For every scored event:

 1. Read its hits (reduced hit table, mm, keV) and its particles (NEXUS
    MC/particles: names, mothers, creator processes, kinetic energies, creation
    and end points).
 2. Main track as in the production: the leading electron of the gamma (the
    most energetic electron produced directly by it), or the generated
    electron, with its descendants not passing through a gamma
    (EventAnalysis.gamma_map). Every other hit belongs to a secondary deposit:
    a particle produced directly by the generated gamma, or a photon emitted
    by the main track, with its descendants.
 3. Distance of every secondary hit to the nearest main-track hit; energy of
    the secondary deposits in the simulated box of the event (vertex-relative
    box of the feature table) by that distance: below 4 mm, 4 to 8 mm, beyond.
 4. Per deposit: type, energy of the particle that starts it, deposited energy,
    smallest distance and energy-weighted centre distance to the main track,
    distances of the centre to the start and to the end of the leading electron.
 5. Interaction history of the gamma: Compton electrons, photoelectric
    absorption, pair production, fluorescence photons.
 6. Energy within 7 mm (blob radius) of the start and of the end of the
    leading electron: main-track hits only and every hit in the box.
 7. Join the feature columns of the event (statistics/events.h5 of the run)
    and its score; passed = score above the threshold of 80 % signal.

Output (schema iitaca.truth.v1): events/<column> one row per event,
deposits/<column> one row per secondary deposit, and the identity of inputs
and code. The energy in the box recomputed from the hits must equal the
energy_true_keV of the production for every event.
=#
using HDF5, TOML, SHA, NearestNeighbors, Statistics
const ROOT = get(ENV, "NAUSICAA_ROOT", normpath(joinpath(@__DIR__, "../..")))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl")); using .EventAnalysis
include(joinpath(ROOT, "scripts/production/batch_inputs.jl"))

const TRUTH_SCHEMA = "iitaca.truth.v1"
const BLOB_RADIUS_MM = 7.0
const XRAY_MAX_KEV = 40.0          # photons of the gamma below this energy: fluorescence
const AUGER_MAX_KEV = 40.0         # electrons of the gamma below this energy, other than Compton: atomic relaxation
const DISTANCE_EDGES_MM = (4.0, 8.0)

# Deposit types (deposits/type)
const DEPOSIT_TYPES = ["compton_electron", "photoelectron", "relaxation_electron", "fluorescence_photon",
                       "other_photon_of_gamma", "positron", "photon_of_main_track", "other", "generated_gamma"]
# Interaction history of the generated gamma (events/history)
const HISTORIES = ["generated_electron", "photoelectric", "compton_photoelectric", "compton_only", "pair", "other"]
const LEADING_PROCESSES = ["primary", "phot", "compt", "conv", "other"]

"""One row of NEXUS MC/particles: positions in mm, kinetic energy in keV, track length in mm."""
struct TruthParticle
    particle_id::Int32
    mother_id::Int32
    name::String
    primary::Bool
    initial_mm::NTuple{3,Float32}
    final_mm::NTuple{3,Float32}
    initial_t_ns::Float32
    kinetic_keV::Float32
    length_mm::Float32
    creator::String
end

"""Particles of the requested events of a NEXUS file, read in blocks."""
function read_truth_particles(path::AbstractString, ids)
    wanted = Set(ids); result = Dict(id => TruthParticle[] for id in ids)
    h5open(path, "r") do f
        ds = f["MC/particles"]
        for first in 1:65536:length(ds), p in ds[first:min(first + 65535, length(ds))]
            p.event_id in wanted || continue
            push!(result[p.event_id], TruthParticle(p.particle_id, p.mother_id, String(p.particle_name), p.primary == 1,
                (p.initial_x, p.initial_y, p.initial_z), (p.final_x, p.final_y, p.final_z), p.initial_t, 1000f0 * p.kin_energy, p.length,
                String(p.creator_proc)))
        end
    end
    for id in ids; isempty(result[id]) && throw(KeyError(id)); end
    result
end

"""Scored events of the sample: event number, batch, score, and the 80 % signal threshold."""
function scored_events(scores::AbstractString, statistics::AbstractString)
    batch_of = h5open(statistics, "r") do f
        Dict(zip(read(f["events/event_number"]), read(f["events/batch"])))
    end
    h5open(scores, "r") do f
        if haskey(f, "test")
            k = read(f["test/label"]) .== 0
            ids = read(f["test/event_number"])[k]; s = read(f["test/score"])[k]
            eff = read(f["test/separation/target_signal_efficiency"]); thr = read(f["test/separation/threshold"])
            b = read(f["test/batch"])[k]; all(b[i] == batch_of[ids[i]] for i in eachindex(ids)) || error("Batches of the scores and of the statistics differ")
        else
            ids = read(f["events/event_number"]); s = read(f["events/score"])
            eff = read(f["target_signal_efficiency"]); thr = read(f["threshold"])
        end
        i = findfirst(e -> isapprox(e, 0.8), eff); i === nothing && error("No 80 % threshold in $scores")
        (event_number = ids, batch = [batch_of[id] for id in ids], score = s, threshold = thr[i])
    end
end

"""Feature columns of the run, indexed by event number."""
function feature_columns(statistics::AbstractString, ids)
    names = ("energy_keV", "energy_true_keV", "E1_keV", "E2_keV", "length_mm", "ion_drift_cm", "segments", "segments_joined", "deposits",
             "components", "main_keV", "second_keV", "deposit_keV", "hit_count",
             "box_x_min_mm", "box_x_max_mm", "box_y_min_mm", "box_y_max_mm", "box_z_min_mm", "box_z_max_mm")
    h5open(statistics, "r") do f
        row = Dict(id => i for (i, id) in enumerate(read(f["events/event_number"])))
        r = [row[id] for id in ids]
        Dict(n => read(f["events/$n"])[r] for n in names)
    end
end

distance(a, b) = sqrt(sum((Float64(a[i]) - Float64(b[i]))^2 for i in 1:3))

"""Truth of one event. `box` is (lower, upper) in mm relative to the event vertex (start of the main track)."""
function event_truth(h::MCEventHits, info::MCEventInfo, parts::Vector{TruthParticle}, box, gamma_event::Bool)
    byid = Dict(p.particle_id => p for p in parts)
    if gamma_event
        g = only(p for p in parts if p.primary && p.name == "gamma")
        children = [p for p in parts if p.mother_id == g.particle_id]
        el = [p for p in children if p.name == "e-"]
        lead = el[argmax([(p.kinetic_keV, -p.particle_id) for p in el])]
        info = MCEventInfo(info.event_id, lead.initial_mm, [lead.particle_id])
    else
        lead = byid[only(info.primary_particle_ids)]; children = TruthParticle[]; g = nothing
    end
    roots = gamma_map([MCParticle(p.particle_id, p.mother_id, p.name) for p in parts], info)
    v = info.generation_vertex_mm
    X = [h.x_mm[j] - v[1] for j in eachindex(h.x_mm)]; Y = [h.y_mm[j] - v[2] for j in eachindex(h.x_mm)]; Z = [h.z_mm[j] - v[3] for j in eachindex(h.x_mm)]
    inbox = [box[1][1] <= X[j] <= box[2][1] && box[1][2] <= Y[j] <= box[2][2] && box[1][3] <= Z[j] <= box[2][3] for j in eachindex(X)]
    E = Float64.(h.energy_keV); root = [roots[p] for p in h.particle_id]; main = root .== 0
    any(main) || error("No main-track hits in event $(h.event_id)")
    tree = KDTree(vcat(X[main]', Y[main]', Z[main]') .|> Float64)
    other = findall(.!main); d = fill(0.0, length(X))
    isempty(other) || (d[other] = nn(tree, vcat(X[other]', Y[other]', Z[other]') .|> Float64)[2])
    start = (0f0, 0f0, 0f0); stop = lead.final_mm .- v
    pos(j) = (X[j], Y[j], Z[j])
    near(p, R) = [distance(pos(j), p) <= R for j in eachindex(X)]
    ns = near(start, BLOB_RADIUS_MM); ne = near(stop, BLOB_RADIUS_MM)
    sec = .!main .& inbox
    ev = Dict{String,Any}(
        "leading_keV" => Float64(lead.kinetic_keV), "leading_length_mm" => Float64(lead.length_mm),
        "leading_extent_mm" => distance(start, stop),
        "main_keV_truth" => sum(E[main]), "box_keV_truth" => sum(E[inbox]), "outside_box_keV_truth" => sum(E[.!inbox]),
        "secondary_keV" => sum(E[sec]; init = 0.0),
        "secondary_lt4_keV" => sum(E[sec .& (d .< DISTANCE_EDGES_MM[1])]; init = 0.0),
        "secondary_4to8_keV" => sum(E[sec .& (DISTANCE_EDGES_MM[1] .<= d .< DISTANCE_EDGES_MM[2])]; init = 0.0),
        "secondary_gt8_keV" => sum(E[sec .& (d .>= DISTANCE_EDGES_MM[2])]; init = 0.0),
        "start_main_keV" => sum(E[ns .& main]; init = 0.0), "start_all_keV" => sum(E[ns .& inbox]; init = 0.0),
        "end_main_keV" => sum(E[ne .& main]; init = 0.0), "end_all_keV" => sum(E[ne .& inbox]; init = 0.0))
    # interaction history of the gamma
    proc = [p.creator for p in children]
    ncompt = count(p -> p.name == "e-" && p.creator == "compt", children); phot = any(==("phot"), proc); conv = any(==("conv"), proc)
    ev["compton_electrons"] = ncompt
    ev["history"] = !gamma_event ? 0 : conv ? 4 : (phot && ncompt == 0) ? 1 : (phot && ncompt > 0) ? 2 : ncompt > 0 ? 3 : 5
    ev["leading_process"] = !gamma_event ? 0 : something(findfirst(==(lead.creator), LEADING_PROCESSES[2:4]), 4)
    xr = [p for p in children if p.name == "gamma" && p.kinetic_keV < XRAY_MAX_KEV]
    ev["fluorescence_photons"] = length(xr)
    # deposits
    groups = Dict{Int32,Vector{Int}}()
    for j in other; push!(get!(groups, root[j], Int[]), j); end
    dep = NamedTuple[]
    for r in sort(collect(keys(groups)))
        ix = groups[r]; p = byid[r]; e = sum(E[ix]); e > 0 || continue
        t = if g !== nothing && p.mother_id == g.particle_id
            p.name == "e-" ? (p.creator == "compt" ? 1 : p.kinetic_keV < AUGER_MAX_KEV ? 3 : 2) :
            p.name == "gamma" ? (p.kinetic_keV < XRAY_MAX_KEV ? 4 : 5) : p.name == "e+" ? 6 : 8
        elseif g !== nothing && r == g.particle_id
            9
        elseif p.name == "gamma" && haskey(roots, p.mother_id) && roots[p.mother_id] == 0
            7
        else
            8
        end
        w = E[ix] ./ e; c = (sum(w .* X[ix]), sum(w .* Y[ix]), sum(w .* Z[ix]))
        dc = nn(tree, Float64[c[1], c[2], c[3]])[2]
        push!(dep, (root_id = Int64(r), type = t, particle_keV = Float64(p.kinetic_keV), deposit_keV = e, deposit_in_box_keV = sum(E[ix[inbox[ix]]]; init = 0.0),
            hits = length(ix), min_distance_mm = minimum(d[ix]), centre_distance_mm = dc,
            centre_start_mm = distance(c, start), centre_end_mm = distance(c, stop),
            centre_x_mm = Float64(c[1]), centre_y_mm = Float64(c[2]), centre_z_mm = Float64(c[3]), in_box = all(inbox[ix])))
    end
    # largest deposit in the box, and the largest fluorescence photon
    inb = [x for x in dep if x.deposit_in_box_keV > 0]
    big = isempty(inb) ? nothing : inb[argmax([x.deposit_in_box_keV for x in inb])]
    ev["largest_deposit_keV"] = big === nothing ? 0.0 : big.deposit_in_box_keV
    ev["largest_deposit_type"] = big === nothing ? 0 : big.type
    ev["largest_deposit_distance_mm"] = big === nothing ? NaN : big.min_distance_mm
    ev["largest_deposit_start_mm"] = big === nothing ? NaN : big.centre_start_mm
    for (cut, name) in ((3.0, "3"), (10.0, "10"))
        ev["deposits_ge$(name)_lt4"] = count(x -> x.deposit_in_box_keV >= cut && x.min_distance_mm < 4, dep)
        ev["deposits_ge$(name)_4to8"] = count(x -> x.deposit_in_box_keV >= cut && 4 <= x.min_distance_mm < 8, dep)
        ev["deposits_ge$(name)_gt8"] = count(x -> x.deposit_in_box_keV >= cut && x.min_distance_mm >= 8, dep)
    end
    xd = [x for x in dep if x.type == 4]
    xb = isempty(xd) ? nothing : xd[argmax([x.deposit_keV for x in xd])]
    ev["xray_keV"] = isempty(xr) ? NaN : Float64(maximum(p.kinetic_keV for p in xr))
    ev["xray_deposit_keV"] = xb === nothing ? 0.0 : xb.deposit_in_box_keV
    ev["xray_distance_mm"] = xb === nothing ? NaN : xb.min_distance_mm
    ev["xray_centre_distance_mm"] = xb === nothing ? NaN : xb.centre_distance_mm
    ev["xray_start_mm"] = xb === nothing ? NaN : xb.centre_start_mm
    ev["xray_end_mm"] = xb === nothing ? NaN : xb.centre_end_mm
    ev, dep
end

function prepare_truth(card::AbstractString, scores::AbstractString, output::AbstractString; batches = nothing)
    c = TOML.parsefile(card); label = c["run_label"]; gamma_event = c["frame_track"] == "leading_gamma_electron"
    statistics = joinpath(ROOT, "results/iitaca", label, "statistics/events.h5")
    s = scored_events(scores, statistics)
    keep = batches === nothing ? trues(length(s.event_number)) : in.(s.batch, Ref(batches))
    ids = s.event_number[keep]; bs = s.batch[keep]; score = s.score[keep]
    feats = feature_columns(statistics, ids)
    order = Dict(id => i for (i, id) in enumerate(ids))
    rows = Vector{Any}(undef, length(ids)); deps = Vector{Any}(undef, length(ids)); identity = Dict{String,String}()
    started = time()
    for b in sort(unique(bs))
        tag = batch_tag(b); bid = ids[bs .== b]
        hp = input_path(ROOT, replace(c["hits"], "{batch}" => tag)); pp = input_path(ROOT, replace(c["particles"], "{batch}" => tag))
        parts = read_truth_particles(pp, bid)
        h5open(hp, "r") do f
            identity[tag] = read(f["metadata/source_sha256"])
            row = Dict(id => i for (i, id) in enumerate(read(f["events/event_id"])))
            for id in bid
                h, info = read_mc_event(f, id; event_row = row[id]); i = order[id]
                box = ((feats["box_x_min_mm"][i], feats["box_y_min_mm"][i], feats["box_z_min_mm"][i]),
                       (feats["box_x_max_mm"][i], feats["box_y_max_mm"][i], feats["box_z_max_mm"][i]))
                rows[i], deps[i] = event_truth(h, info, parts[id], box, gamma_event)
            end
        end
        println("Batch $tag: $(length(bid)) events; $(round(time() - started; digits = 1)) s"); flush(stdout)
    end
    # the box energy recomputed from the hits must match the production
    dE = maximum(abs(rows[i]["box_keV_truth"] - feats["energy_true_keV"][i]) for i in eachindex(ids))
    dE < 0.05 || error("Energy in the box differs from the production by $dE keV")
    mkpath(dirname(output)); tmp = output * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = TRUTH_SCHEMA; f["run_label"] = label; f["sample"] = c["sample"]
            f["scores_file"] = basename(scores); f["scores_sha256"] = bytes2hex(open(sha256, scores))
            f["statistics_sha256"] = bytes2hex(open(sha256, statistics)); f["threshold_80"] = s.threshold
            f["selection"] = "selected events of $(label) (ROI and track selection) with a network score; passed: score above the 80 % signal threshold"
            f["code_sha256"] = bytes2hex(open(sha256, @__FILE__)); f["code"] = "scripts/iitaca/prepare_iitaca_truth.jl"
            f["batches"] = sort(collect(keys(identity))); f["particle_file_sha256"] = [identity[k] for k in sort(collect(keys(identity)))]
            f["deposit_types"] = DEPOSIT_TYPES; f["histories"] = HISTORIES; f["leading_processes"] = LEADING_PROCESSES
            f["blob_radius_mm"] = BLOB_RADIUS_MM; f["max_box_energy_difference_keV"] = dE
            f["events/event_number"] = ids; f["events/batch"] = bs; f["events/score"] = score; f["events/passed"] = Int8.(score .> s.threshold)
            for (k, v) in feats; f["events/$k"] = v; end
            for k in sort(collect(keys(rows[1]))); f["events/$k"] = [r[k] for r in rows]; end
            n = [length(d) for d in deps]; f["events/deposit_offset"] = cumsum(n) .- n; f["events/deposit_count"] = n
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
    length(args) == 3 || error("Usage: prepare_iitaca_truth.jl <run card> <scores> <output.h5> [--batches=A:B]")
    prepare_truth(abspath(args[1]), abspath(args[2]), abspath(args[3]); batches)
end
