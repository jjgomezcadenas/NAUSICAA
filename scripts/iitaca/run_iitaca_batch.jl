#!/usr/bin/env julia
#=
iItaca: read one batch of Monte Carlo events and write its hit file.

  julia --project=. scripts/iitaca/run_iitaca_batch.jl <run card> --batch=N [--limit=M]

For every event of the batch:

 1. Fast trigger and fixed frame from the Monte Carlo (EventAnalysis.mc_fast_trigger).
 2. Box of the ion response: main track and the gamma deposits inside the frame, with the margin.
 3. Event energy: energy deposited in that box, smeared with the resolution of the energy card.
 4. Ions of each deposit (W, Fano), drift with the thermal spread to the pads, frames.
 5. VH1 electronics and reconstruction of the ions per pad and frame.
 6. Hits above the saved hit cut, and the noise hits of the rest of the frame.

Output: results/iitaca/<run_label>/batches/<tag>/hits.h5 (schema iitaca.hits.v1)
with the decision of every event, the table of accepted events, their itrack
hits, the configuration and the hashes of inputs and code. The electrons are
not simulated: the driver reads `IItaca.simulation_inputs` and nothing else of
the detector. --limit reads only the first M events and writes under
<run_label>_first<M>.
=#
using HDF5, TOML, Random, SHA, Statistics, LinearAlgebra
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl")); using .EventAnalysis
include(joinpath(ROOT, "src/detector_response/DetectorResponse.jl")); using .DetectorResponse
include(joinpath(ROOT, "src/ion_signal/guarded_files.jl"))
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
include(joinpath(ROOT, "src/iitaca/IItaca.jl")); using .IItaca
include(joinpath(ROOT, "src/electronics/ElectronicsResponse.jl")); using .ElectronicsResponse
include(joinpath(ROOT, "src/reconstruction/SignalReconstruction.jl")); using .SignalReconstruction
include(joinpath(ROOT, "src/vh1/VH1Readout.jl")); using .VH1Readout
include(joinpath(ROOT, "scripts/production/batch_inputs.jl"))
include(joinpath(ROOT, "scripts/production/frame_noise_hits.jl"))

const HIT_SCHEMA = "iitaca.hits.v1"
const CODE_DIRECTORIES = ("src/event_analysis", "src/detector_response", "src/electronics", "src/reconstruction", "src/ion_signal", "src/vh1", "src/iitaca", "src/itaca")

"""Fast-trigger settings from the simulation inputs."""
trigger_settings(sim) = FixedFrameTriggerSettings(Float32(10sim["geometry.fiducial_length_cm"]), Float32(5sim["geometry.fiducial_diameter_cm"]),
    Float32(sim["trigger.frame_side_mm"]), Float32(sim["trigger.gamma_energy_cut_keV"]), sim["trigger.frame_placement"] == "inside_detector",
    Float32(sim["trigger.el_exclusion_mm"]))

"""Ion-response settings of one event from the simulation inputs.

The pad lattice has a random offset per event (stream `ion_grid`). The ion
spread is the thermal one for the drift length of the event. The electron
entries of `DetectorResponseSettings` are not used: they are NaN, except the
electron spreads, zero, which the ion response requires to be finite.
"""
function ion_event_settings(sim, seed, p)
    F = Float32; pitch = F(sim["ions.pad_pitch_mm"])
    grid = ReadoutGrid(pitch, sample_grid_offset_mm(Xoshiro(event_seed(seed, p.event_id, "ion_grid")), pitch))
    sigma = sigma_t_ion_mm(F(sim["gas.temperature_K"]), p.ion_drift_cm, F(sim["drift.field_V_cm"]))
    gas = sim["gas.helium_mole_fraction"] == 0 ? Xe_100 : Xe_85_He_15
    DetectorResponseSettings(F(10sim["geometry.fiducial_length_cm"]), F(sim["energy.W_eV"]), F(sim["energy.fano_factor"]), gas, grid, grid,
        F(sim["ions.frame_ms"]), NaN32, IonTransportSettings(F(sim["drift.ion_drift_velocity_cm_s"] / 100), sigma, sigma),
        ElectronTransportSettings(NaN32, 0f0, 0f0), NaN32, NaN32, NaN32, NaN32, :not_simulated)
end

"""SHA-256 of the driver and of every Julia file of the modules it uses, keyed by path."""
function code_hashes()
    files = [@__FILE__, joinpath(ROOT, "scripts/production/batch_inputs.jl"), joinpath(ROOT, "scripts/production/frame_noise_hits.jl")]
    for dir in CODE_DIRECTORIES, (folder, _, names) in walkdir(joinpath(ROOT, dir)), name in names
        endswith(name, ".jl") && push!(files, joinpath(folder, name))
    end
    Dict(relpath(f, ROOT) => bytes2hex(open(sha256, f)) for f in files)
end

columns(g, rows) = (isempty(rows) || for k in keys(first(rows)); g[string(k)] = [r[k] for r in rows]; end)

"""
    run_iitaca_batch(card; batch, limit=0, output_root=joinpath(ROOT, "results/iitaca"))::String

Process one batch of the run card and write its hit file; return the path.
An existing hit file is preserved: the function stops.
"""
function run_iitaca_batch(card::AbstractString; batch::Integer, limit::Integer = 0, output_root::AbstractString = joinpath(ROOT, "results/iitaca"))
    BLAS.set_num_threads(1)
    started = time(); r = resolve_run(card; root = ROOT); c = r.run; select_batch!(c, batch)
    limit > 0 && (c["max_events"] = limit; c["run_label"] = replace(c["run_label"], "/batches/" => "_first$(limit)/batches/"))
    sim = simulation_inputs(r.resolved.detector); en = r.resolved.detector["energy"]
    out = joinpath(output_root, c["run_label"]); target = joinpath(out, "hits.h5")
    isfile(target) && error("Preserve existing hits: $target")
    hp = input_path(ROOT, c["hits"]); pp = input_path(ROOT, c["particles"])
    source_sha = bytes2hex(open(sha256, hp)); particle_sha = bytes2hex(open(sha256, pp))
    ids = h5open(hp, "r") do f
        read(f["metadata/source_sha256"]) == particle_sha || error("Particle source mismatch")
        read(f["events/event_id"])
    end
    get(c, "max_events", 0) > 0 && (ids = ids[1:min(length(ids), c["max_events"])])
    parts = read_mc_particles(pp, ids)
    leading = c["frame_track"] == "leading_gamma_electron" ? read_leading_gamma_electrons(pp, ids) : nothing
    trigger = trigger_settings(sim); L = trigger.chamber_length_mm; margin = Float32(sim["trigger.simulation_margin_mm"])
    W = Float32(sim["energy.W_eV"]); fano = Float32(sim["energy.fano_factor"]); scale = Float64(W) / 1000      # keV per ion
    cut = Float64(sim["tracks.saved_hit_cut_keV"]); seed = c["seed"]; eseed = c["electronics_seed"]
    max_bytes = get(c, "max_vector_bytes", 268435456)
    vh = load_vh1_settings(ROOT, joinpath(ROOT, sim["ions.vh1_card"]); minimum_samples = get(c, "vh1_minimum_samples", 32768))
    stats = FixedFrameStatistics(); decisions = NamedTuple[]; events = NamedTuple[]
    x = Float32[]; y = Float32[]; z = Float32[]; e = Float32[]; timing = zeros(5)
    h5open(hp, "r") do f
        for (row, id) in enumerate(ids)
            local h, info, p, o
            timing[1] += @elapsed begin
                h, info = read_mc_event(f, id; event_row = row)
                leading === nothing || (info = leading_electron_info(info, leading[id]))
                h = center_hits(h, info); p = placement_from_vertex(info, L); o = mc_fast_trigger(h, info, parts[id], p, trigger)
                stats = update_statistics(stats, o); push!(decisions, (event_number = id, decision = Int(o.decision), ion_drift_cm = p.ion_drift_cm))
            end
            o.decision == trigger_accepted || continue
            rng(tag) = Xoshiro(event_seed(seed, id, tag)); erng(tag) = Xoshiro(event_seed(eseed, id, tag))
            local box, selected, ds, ions, adc, ri, hits, E_true, E, sigma_ions
            timing[2] += @elapsed begin
                box = simulation_box(h, info, parts[id], o.frame, margin); selected = hits_in_frame(h, box)
                gm = gamma_map(parts[id], info)
                iszero(sum(h.energy_keV[j] for j in eachindex(h.x_mm) if gm[h.particle_id[j]] == 0 && !contains(box, (h.x_mm[j], h.y_mm[j], h.z_mm[j])); init = 0f0)) ||
                    error("Main-track energy outside the box in event $id")
                E_true = frame_energy_keV(h, box); E = emulated_energy_keV(rng("emulated_energy"), E_true, en)
                pairs = sample_fano_pairs(rng("pairs"), selected, W, fano); ds = ion_event_settings(sim, seed, p)
                ions = sim_ion_detector_response(rng("ions"), selected, pairs, p, box, ds; max_vector_bytes = max_bytes)
            end
            timing[3] += @elapsed adc = sim_vh1_electronics(erng("vh1_noise"), erng("vh1_adc"), ions, vh)
            timing[4] += @elapsed begin
                plan = prepare_vh1_reconstruction(adc.time, ions.ions.time, vh); ri = reco_vh1_ions(adc, plan)
                sigma_ions = sqrt(sum(plan.ion_covariance[k, k] for k in axes(plan.ion_covariance, 1)) / size(plan.ion_covariance, 1))
            end
            timing[5] += @elapsed begin
                hits = reconstruct_hits(ri, p, ds; keV_per_count = scale, minimum_energy_keV = cut)
                hits = join_hits(hits, frame_noise_hits(erng("frame_noise"), o.frame, ri, p, ds, sigma_ions, scale, cut))
            end
            valid = ri.valid .& isfinite.(ri.ions); total = sum(Float64, ri.ions[valid]; init = 0.0)
            v = info.generation_vertex_mm
            push!(events, (event_number = id, energy_keV = E, energy_true_keV = E_true, energy_outside_box_keV = Float64(sum(Float64, h.energy_keV)) - E_true,
                ion_energy_keV = all(valid) ? total * scale : NaN, invalid_ion_buckets = count(.!valid), saturated_pads = count(adc.saturated),
                ion_drift_cm = p.ion_drift_cm, ion_sigma_mm = ds.ion_transport.sigma_t_mm, ion_noise_ions = sigma_ions,
                vertex_x_mm = v[1], vertex_y_mm = v[2], vertex_z_mm = v[3],
                frame_x_min_mm = o.frame.lower_mm[1], frame_x_max_mm = o.frame.upper_mm[1], frame_y_min_mm = o.frame.lower_mm[2],
                frame_y_max_mm = o.frame.upper_mm[2], frame_z_min_mm = o.frame.lower_mm[3], frame_z_max_mm = o.frame.upper_mm[3],
                box_x_min_mm = box.lower_mm[1], box_x_max_mm = box.upper_mm[1], box_y_min_mm = box.lower_mm[2],
                box_y_max_mm = box.upper_mm[2], box_z_min_mm = box.lower_mm[3], box_z_max_mm = box.upper_mm[3],
                pads = length(ri.channels.ix), frames = Int(ri.time.bucket_count), hit_offset = length(x), hit_count = length(hits.x_mm)))
            append!(x, hits.x_mm); append!(y, hits.y_mm); append!(z, hits.z_mm); append!(e, hits.energy_keV)
            if length(events) % 100 == 0
                println("Accepted $(length(events)); input $row/$(length(ids)); $(round(time() - started; digits = 1)) s"); flush(stdout); GC.gc(false)
            end
        end
    end
    mkpath(out); budget = 16 * 1024^2 + 16 * length(x) + 512 * length(events) + 64 * length(decisions)
    GuardedFiles.space_check(out, budget, 1024^3); tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = HIT_SCHEMA
            g = create_group(f, "configuration")
            g["run_label"] = c["run_label"]; g["sample"] = c["sample"]; g["variant"] = r.resolved.variant; g["batch"] = batch_tag(batch)
            g["seed"] = seed; g["electronics_seed"] = eseed; g["frame_track"] = c["frame_track"]
            g["resolved"] = resolved_text(r, ROOT); g["run_card"] = read(card, String); g["vh1_card"] = read(joinpath(ROOT, sim["ions.vh1_card"]), String)
            for (name, text) in r.resolved.set.texts; g["cards/$name"] = text; end
            r.resolved.variant_text == "" || (g["cards/variant"] = r.resolved.variant_text)
            sims = create_group(g, "simulation_inputs")
            for (k, value) in sim; sims[k] = value isa Real && !(value isa Bool) ? Float64(value) : value isa AbstractDict ? sprint(io -> TOML.print(io, value)) : string(value); end
            g["sha256/hits"] = source_sha; g["sha256/particles"] = particle_sha
            g["sha256/vh1_card"] = bytes2hex(open(sha256, joinpath(ROOT, sim["ions.vh1_card"])))
            code = code_hashes(); names = sort(collect(keys(code))); g["code_files"] = names; g["code_sha256"] = [code[n] for n in names]
            columns(create_group(f, "selection/events"), decisions)
            for k in fieldnames(typeof(stats)); f["selection/counts/" * string(k)] = getfield(stats, k); end
            columns(create_group(f, "events"), events)
            t = create_group(f, "itrack")
            for (name, values) in (("x_mm", x), ("y_mm", y), ("z_mm", z), ("energy_keV", e))
                isempty(values) ? (t[name] = values) : (t[name, chunk = (min(length(values), 65536),), deflate = 4] = values)
            end
            attrs(t)["coordinates"] = "mm relative to the event vertex; z along the drift"
            for (name, value) in zip(("selection", "transport", "vh1", "reconstruction", "hits"), timing); f["timing/$(name)_s"] = value; end
            f["timing/total_s"] = time() - started; f["timing/input_events"] = length(ids); f["timing/accepted_events"] = length(events)
            flush(f); GuardedFiles.check_size(tmp, budget)
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    n = max(length(events), 1)
    println("TIMING input=$(length(ids)) accepted=$(length(events)) total_s=$(round(time() - started; digits = 1))")
    for (name, value) in zip(("selection", "transport", "vh1", "reconstruction", "hits"), timing)
        println("TIMING $name total_s=$(round(value; digits = 2)) per_accepted_event_s=$(round(value / n; digits = 4))")
    end
    println("OUTPUT $target bytes=$(filesize(target))")
    target
end

if abspath(PROGRAM_FILE) == @__FILE__
    batch, args = batch_argument(ARGS)
    k = findall(startswith("--limit="), args); limit = isempty(k) ? 0 : parse(Int, split(args[k[1]], "=")[2]); args = args[setdiff(eachindex(args), k)]
    (batch !== nothing && length(args) == 1) || error("Usage: run_iitaca_batch.jl run_card.toml --batch=N [--limit=M]")
    run_iitaca_batch(abspath(args[1]); batch, limit)
end
