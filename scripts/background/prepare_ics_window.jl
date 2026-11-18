#!/usr/bin/env julia
#=
Radiogenic background: fraction of the copper-shield decays that reach the saved
energy window, from the Geant4 files of one isotope.

  julia --project=. scripts/background/prepare_ics_window.jl <shielding card> <isotope> <folder of .next.h5 files> [OUTPUT.h5]

The folder holds the files <isotope>_bNNNNN.next.h5 written by
background_nexus_array.slurm (region ICS, full decays, events saved when the
energy deposited in the fiducial gas lies in the window). For each saved event:

 1. The decay vertex is the initial position of the primary ion (MC/particles,
    primary = 1). Its piece of the copper (barrel, cathode disk, anode disk) and
    its depth from the inner face give the weight that restores the uniform
    distribution of the decays (Background.weight); the piece is checked
    against the volume name Geant4 recorded.
 2. The energy is the sum of the hits (MC/hits, MeV to keV).
 3. The particles that entered the gas: for every particle with hits, the first
    ancestor created outside the fiducial gas (initial_volume not ITACA_GAS).
    Their number and kinds classify the event: one gamma, one electron (e- or
    e+), one other particle, or several entering particles (two separate
    deposits, rejected by the track selection).
 4. The number of decays generated is the sum of num_events over the files
    (MC/configuration).

Writes OUTPUT.h5 (default results/background/<name>/<isotope>/window.h5, schema
nausicaa.background_window.v1): the event table, the sums of the weights per
class, the fraction of the decays in the window per class with its error
sqrt(sum w^2) / N, the copper mass, its activity, the decays per year and the
events per year in the window per class. An existing output is preserved.
=#
include(joinpath(@__DIR__, "resolve_background.jl"))
using HDF5, SHA, Statistics

const WINDOW_SCHEMA = "nausicaa.background_window.v1"
const GAS_VOLUME = "ITACA_GAS"
const PIECE_CODE = Dict(:barrel => Int8(1), :cathode => Int8(2), :anode => Int8(3))
const PIECE_VOLUME = Dict(:barrel => "ITACA_ICS_BARREL", :cathode => "ITACA_ICS_CATHODE", :anode => "ITACA_ICS_ANODE")
const CLASSES = ("gamma", "electron", "other", "multiple")      # codes 1 to 4
const ISOTOPE_KEY = Dict("bi214" => "bi214_Bq_kg", "tl208" => "tl208_Bq_kg")

"""Configuration of a NEXUS file as a dictionary of its parameter keys and values."""
function nexus_configuration(f)
    c = read(f["MC/configuration"])
    Dict(String(strip(String(r.param_key))) => String(strip(String(r.param_value))) for r in c)
end

"""The entering particle of a hit particle: the first ancestor created outside the fiducial gas."""
function entering_particle(particles::Dict{Int32,<:NamedTuple}, id::Int32)
    seen = 0
    while true
        p = particles[id]
        p.initial_volume == GAS_VOLUME || return p
        p.mother_id == 0 && error("Particle $id of event $(p.event_id) was created in the gas without a mother")
        id = p.mother_id; (seen += 1) < 10_000 || error("Mother chain of particle $id does not end")
    end
end

kind(name) = name == "gamma" ? "gamma" : name in ("e-", "e+") ? "electron" : "other"

"""One file: per-event rows and the number of decays generated."""
function analyse_file(path::AbstractString, g::ShieldGeometry, batch::Int)
    h5open(path, "r") do f
        conf = nexus_configuration(f)
        haskey(conf, "num_events") || error("No num_events in the configuration of $path")
        generated = parse(Int, conf["num_events"])
        particles = read(f["MC/particles"]); hits = read(f["MC/hits"])
        by_event = Dict{Int64,Dict{Int32,eltype(particles)}}()
        for p in particles
            push!(get!(by_event, Int64(p.event_id), Dict{Int32,eltype(particles)}()), p.particle_id => p)
        end
        energy = Dict{Int64,Float64}(); hit_particles = Dict{Int64,Set{Int32}}()
        for h in hits
            e = Int64(h.event_id); energy[e] = get(energy, e, 0.0) + Float64(h.energy) * 1000
            push!(get!(hit_particles, e, Set{Int32}()), h.particle_id)
        end
        rows = NamedTuple[]
        for (e, ps) in sort!(collect(by_event); by = first)
            prim = [p for p in values(ps) if p.primary == 1]
            length(prim) == 1 || error("Event $e of $path has $(length(prim)) primaries")
            ion = prim[1]
            x, y, z = Float64(ion.initial_x) / 10, Float64(ion.initial_y) / 10, Float64(ion.initial_z) / 10      # mm to cm
            pc = piece(g, x, y, z)
            String(ion.initial_volume) == PIECE_VOLUME[pc] || error("Event $e of $path: vertex in $(ion.initial_volume), geometry says $(PIECE_VOLUME[pc])")
            entering = Set{Int32}(); kinds = String[]
            for id in get(hit_particles, e, Set{Int32}())
                q = entering_particle(ps, id)
                q.particle_id in entering && continue
                push!(entering, q.particle_id); push!(kinds, kind(String(q.particle_name)))
            end
            class = length(entering) == 0 ? Int8(0) : length(entering) > 1 ? Int8(4) : Int8(findfirst(==(kinds[1]), CLASSES))
            push!(rows, (event_number = e, batch = batch, x_cm = x, y_cm = y, z_cm = z, piece = PIECE_CODE[pc],
                         depth_cm = depth_cm(g, x, y, z), weight = weight(g, x, y, z), energy_keV = get(energy, e, 0.0),
                         entering = Int32(length(entering)), class = class, ion = String(ion.particle_name)))
        end
        rows, generated, conf
    end
end

function prepare_ics_window(card_path::AbstractString, isotope::AbstractString, folder::AbstractString, target::AbstractString = "")
    haskey(ISOTOPE_KEY, isotope) || error("Unknown isotope $isotope (bi214 or tl208)")
    r = resolved_shield(card_path); g = r.geometry; name = r.card.card["name"]
    files = sort(filter(p -> occursin(Regex("^$(isotope)_b\\d{5}\\.next\\.h5\$"), basename(p)), readdir(abspath(folder); join = true)))
    isempty(files) && error("No $(isotope)_bNNNNN.next.h5 under $folder")
    started = time(); rows = NamedTuple[]; generated = 0; shas = String[]; batches = Int[]; window = nothing
    for path in files
        batch = parse(Int, match(r"_b(\d{5})\.next\.h5$", basename(path)).captures[1])
        fr, n, conf = analyse_file(path, g, batch)
        append!(rows, fr); generated += n; push!(shas, bytes2hex(open(sha256, path))); push!(batches, batch)
        w = map(("min_energy", "max_energy")) do k                    # "2420. keV" as stored by NEXUS
            v = split(conf["/Actions/DefaultEventAction/$k"]); length(v) == 2 && v[2] == "keV" || error("Energy $k of $path is not in keV: $(conf["/Actions/DefaultEventAction/$k"])")
            parse(Float64, v[1])
        end
        window === nothing ? (window = w) : (w == window || error("Saved window differs in $path"))
        println("$(basename(path)): $(length(fr)) events in the window of $n decays, $(round(time() - started; digits = 1)) s"); flush(stdout)
    end
    allunique(getfield.(rows, :event_number)) || error("Duplicate event numbers across files")
    W = [x.weight for x in rows]; C = [x.class for x in rows]
    sumw = [sum(W[C .== c]; init = 0.0) for c in 1:4]; sumw2 = [sum(W[C .== c] .^ 2; init = 0.0) for c in 1:4]
    fraction = sumw ./ generated; error_ = sqrt.(sumw2) ./ generated
    Bq_kg = r.card.card["ics"][ISOTOPE_KEY[isotope]]; mass = copper_mass_kg(g); per_year = decays_per_year(g, Bq_kg)
    isempty(target) && (target = joinpath(ROOT, "results/background", name, isotope, "window.h5")); target = abspath(target)
    isfile(target) && error("Preserve existing results: $target"); mkpath(dirname(target))
    tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = WINDOW_SCHEMA; f["name"] = name; f["isotope"] = isotope; f["family"] = r.family; f["variant"] = r.variant
            f["card_path"] = relpath(r.card.path, ROOT); f["card_text"] = r.card.text; f["card_sha256"] = r.card.sha256
            f["files"] = relpath.(files, ROOT); f["file_sha256"] = shas; f["batches"] = batches
            f["generated_decays"] = generated; f["saved_events"] = length(rows)
            f["window/lower_keV"] = window[1]; f["window/upper_keV"] = window[2]; f["window/source"] = "MC/configuration of the files (DefaultEventAction min_energy, max_energy, keV)"
            for k in fieldnames(ShieldGeometry); f["geometry/$k"] = getfield(g, k); end
            f["geometry/outer_cm"] = outer_cm(g); f["geometry/volume_cm3"] = volume_cm3(g); f["geometry/volume_barrel_cm3"] = volume_barrel_cm3(g); f["geometry/volume_disk_cm3"] = volume_disk_cm3(g)
            f["copper/mass_kg"] = mass; f["copper/specific_activity_Bq_kg"] = Bq_kg; f["copper/activity_Bq"] = activity_Bq(g, Bq_kg); f["copper/decays_per_year"] = per_year
            f["classes"] = collect(CLASSES); f["class_codes"] = Int8.(1:4)
            f["window_fraction"] = fraction; f["window_fraction_error"] = error_; f["sum_weights"] = sumw; f["sum_weights_squared"] = sumw2
            f["events_per_year"] = fraction .* per_year; f["events_per_year_error"] = error_ .* per_year
            f["counts"] = [count(==(c), C) for c in 1:4]                       # saved events per class, unweighted
            for k in (:event_number, :batch, :x_cm, :y_cm, :z_cm, :piece, :depth_cm, :weight, :energy_keV, :entering, :class)
                f["events/$k"] = [getfield(x, k) for x in rows]
            end
            f["events/ion"] = [x.ion for x in rows]
            f["time_s"] = time() - started
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("WINDOW $target: $(length(rows)) events of $generated decays; window fraction per class " *
            join(("$(CLASSES[c]) $(round(fraction[c]; sigdigits = 3)) ± $(round(error_[c]; sigdigits = 2)) ($(count(==(c), C)))" for c in 1:4), ", ") *
            "; copper $(round(mass; digits = 0)) kg, $(round(per_year; sigdigits = 3)) decays per year; events per year in the window " *
            join(("$(CLASSES[c]) $(round(fraction[c] * per_year; sigdigits = 3))" for c in 1:4), ", "))
    target
end

if abspath(PROGRAM_FILE) == @__FILE__
    3 <= length(ARGS) <= 4 || error("Usage: prepare_ics_window.jl <shielding card> <isotope> <folder> [OUTPUT.h5]")
    prepare_ics_window(ARGS...)
end
