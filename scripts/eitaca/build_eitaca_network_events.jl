#!/usr/bin/env julia
#=
eITACA: training events of one run for the classification network, from the etrack.

  julia --project=. scripts/eitaca/build_eitaca_network_events.jl <run folder> [OUTPUT.h5] [--window=LO,HI]

For each batch it reads features.h5 and electrons.h5 (prepare_eitaca_features.jl,
run_eitaca_batch.jl), keeps the events whose event energy is in the ROI and that
pass the track selection (in_roi and selected), and places their etrack hits on
the readout lattice (SiPM pitch in x and y, rebinned bucket in z) inside a cube
of the frame side of the trigger card centred on the hits of the event, one hit
per cell. The electron file does not hold the frame; the network sees the hits
relative to each other, so the centring changes no input. Writes
results/network/<variant>/events/<sample>/events.h5 by default, schema
itaca.network_events.v1 with the etrack only; network_trainer.py reads it with
--input=etrack. The track features also appear with the etrack_ prefix.

With --window=LO,HI (keV) the energy cut is the box energy energy_true_keV in
[LO, HI] instead of in_roi; the default output is then
results/network/<variant>/events/<run label>_window/events.h5.
=#
include(joinpath(@__DIR__, "../network/network_events_builder.jl"))     # lattice_cells, lattice_indices, one_hit_per_cell; GuardedFiles

const EITACA_FEATURES = "eitaca.features.v1"
const TRACK_COLUMNS = ("hits", "components", "main_keV", "second_keV", "main_voxels", "length_mm", "extent_x_mm", "extent_y_mm", "extent_z_mm", "E1_keV", "E2_keV")

"""Pitches (x, y, z in mm) of the etrack lattice and the frame side (mm), from the resolved configuration of an electron file."""
function eitaca_lattice(resolved::AbstractString)
    d = TOML.parse(resolved)["detector"]; pitch = Float64(d["el"]["hole_pitch_mm"])
    dz = d["sipm"]["electron_bucket_us"] * d["reconstruction"]["rebin_factor"] * d["drift"]["electron_drift_velocity_mm_us"]
    (pitch = (pitch, pitch, Float64(dz)), side = Float64(d["trigger"]["frame_side_mm"]))
end

function build_eitaca_network_events(run_dir::AbstractString, target::AbstractString = ""; window = nothing)
    window === nothing || (length(window) == 2 && window[1] < window[2]) || error("The window needs LO < HI, found $window")
    run_dir = abspath(run_dir); root = joinpath(run_dir, "batches")
    tags = sort(filter(t -> occursin(r"^b\d{5}$", t) && isfile(joinpath(root, t, "features.h5")), readdir(root))); isempty(tags) && error("No feature tables under $root")
    started = time(); columns = Dict{String,Vector}(); accepted = 0; in_roi_events = 0; roi_selected = 0; shas = String[]
    ref = nothing; lattice = nothing; n = (0, 0, 0)
    a = (offsets = Int64[0], ix = UInt8[], iy = UInt8[], iz = UInt8[], energy = Float32[]); dropped = 0
    for t in tags
        dir = joinpath(root, t); electrons = joinpath(dir, "electrons.h5")
        h5open(joinpath(dir, "features.h5"), "r") do ft
            read(ft["schema"]) == EITACA_FEATURES || error("Unexpected schema in $dir/features.h5")
            ident = Dict(k => read(ft[k]) for k in ("sample", "variant", "frame_track", "settings", "roi_lower_keV", "roi_upper_keV", "run_label"))
            ident["run_label"] = dirname(dirname(ident["run_label"]))
            ref === nothing ? (ref = ident) : for k in keys(ident); ident[k] == ref[k] || error("Setting $k differs in batch $t"); end
            bytes2hex(open(sha256, electrons)) == read(ft["electron_file_sha256"]) || error("Feature table and electron file of $t do not match")
            push!(shas, read(ft["electron_file_sha256"]))
            e = Dict(k => read(ft["events"][k]) for k in keys(ft["events"]))
            accepted += length(e["event_number"]); in_roi_events += count(==(1), e["in_roi"])
            roi_selected += count((e["in_roi"] .== 1) .& (e["selected"] .== 1))
            energy_cut = window === nothing ? e["in_roi"] .== 1 : window[1] .<= e["energy_true_keV"] .<= window[2]
            keep = findall(energy_cut .& (e["selected"] .== 1))
            h5open(electrons, "r") do f
                read(f["events/event_number"]) == e["event_number"] || error("Event order differs between features and electrons in $t")
                l = eitaca_lattice(read(f["configuration/resolved"]))
                lattice === nothing ? (lattice = l; n = ntuple(i -> lattice_cells(l.side, l.pitch[i]), 3)) : (l == lattice || error("Lattice or frame differs in $t"))
                all(<=(256), n) || error("Grid too large for 8-bit indices: $n")
                off = read(f["events/hit_offset"]); cnt = read(f["events/hit_count"]); xyz = [read(f["etrack/$(c)_mm"]) for c in "xyz"]; energy = read(f["etrack/energy_keV"])
                for row in keep
                    r = off[row]+1:off[row]+cnt[row]
                    lower = [(minimum(view(xyz[i], r)) + maximum(view(xyz[i], r))) / 2 - lattice.side / 2 for i in 1:3]
                    idx = [lattice_indices(view(xyz[i], r), lower[i], lattice.pitch[i], n[i]; tolerance = i == 3 ? 0.5 : 0.02) for i in 1:3]   # z: the last time group of an event can be shorter, as network_events_builder.jl
                    en = view(energy, r); kept = one_hit_per_cell(idx[1], idx[2], idx[3], en); dropped += length(en) - length(kept)
                    append!(a.ix, idx[1][kept]); append!(a.iy, idx[2][kept]); append!(a.iz, idx[3][kept]); append!(a.energy, en[kept]); push!(a.offsets, length(a.energy))
                end
            end
            for (k, v) in e; append!(get!(columns, k, similar(v, 0)), v[keep]); end
        end
        println("batch $t: $(length(columns["event_number"])) of $accepted events kept, $(round(time() - started; digits = 1)) s"); flush(stdout)
    end
    allunique(columns["event_number"]) || error("Duplicate event numbers across batches")
    for k in TRACK_COLUMNS; columns["etrack_$k"] = columns[k]; end
    isempty(target) && (target = joinpath(ROOT, "results/network", ref["variant"], "events", window === nothing ? ref["sample"] : ref["run_label"] * "_window", "events.h5"))
    target = abspath(target)
    isfile(target) && error("Preserve existing events: $target"); mkpath(dirname(target))
    budget = 16 * 1024^2 + 2 * (Base.summarysize(columns) + Base.summarysize(a)); GuardedFiles.space_check(dirname(target), budget, 1024^3); tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = NETWORK_EVENTS_SCHEMA; f["species"] = ref["sample"]; f["run_label"] = ref["run_label"]; f["batches"] = tags; f["electron_file_sha256"] = shas
            f["selection"] = window === nothing ? "eITACA $(ref["variant"]): ROI cut on the event energy, $(ref["roi_lower_keV"]) to $(ref["roi_upper_keV"]) keV; track selection of the etrack on the readout lattice (selected)" :
                "eITACA $(ref["variant"]): window cut on the box energy (energy_true_keV), $(Float64(window[1])) to $(Float64(window[2])) keV; track selection of the etrack on the readout lattice (selected)"
            window === nothing || (f["window/lower_keV"] = Float64(window[1]); f["window/upper_keV"] = Float64(window[2]); f["window/source"] = "features.h5 of the batches (energy_true_keV)")
            f["accepted_events"] = accepted; f["frame_side_mm"] = lattice.side
            f["roi/lower_keV"] = ref["roi_lower_keV"]; f["roi/upper_keV"] = ref["roi_upper_keV"]; f["roi/events"] = in_roi_events; f["roi/source"] = "features.h5 of the batches (in_roi)"; f["roi/selected_events"] = roi_selected
            f["settings"] = ref["settings"]; f["variant"] = ref["variant"]
            for (k, v) in columns; f["events/$k"] = v; end
            f["etrack/pitch_mm"] = collect(lattice.pitch); f["etrack/cells"] = collect(n); f["etrack/offsets"] = a.offsets; f["etrack/hits_dropped_from_shared_cells"] = dropped
            f["etrack/placement"] = "cube of frame_side_mm centred on the hits of each event"
            for (k, v) in (("ix", a.ix), ("iy", a.iy), ("iz", a.iz), ("energy_keV", a.energy)); f["etrack/$k", chunk = (max(1, min(length(v), 65536)),), shuffle = true, deflate = 4] = v; end
            f["time_s"] = time() - started; flush(f); GuardedFiles.check_size(tmp, budget)
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("NETWORK EVENTS $target events=$(length(columns["event_number"])) of $accepted ($in_roi_events in the ROI) hits dropped from shared cells: $dropped bytes=$(filesize(target)) time_s=$(round(time() - started; digits = 1))")
    target
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = filter(a -> startswith(a, "--"), ARGS); args = filter(a -> !startswith(a, "--"), ARGS)
    all(o -> startswith(o, "--window="), opts) && length(opts) <= 1 && 1 <= length(args) <= 2 ||
        error("Usage: build_eitaca_network_events.jl <run folder> [OUTPUT.h5] [--window=LO,HI]")
    window = isempty(opts) ? nothing : parse.(Float64, split(opts[1][length("--window=")+1:end], ","))
    build_eitaca_network_events(args...; window)
end
