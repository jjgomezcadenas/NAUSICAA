#!/usr/bin/env julia
# Replay the completed primary-pair realisation into 250 ns photon histories.
include(joinpath(@__DIR__, "run_detector_pilot.jl"))

"""Hash the code that determines this electron-only replay and its output."""
function electron_replay_code_hash()::String
    files = ("src/event_analysis/EventAnalysis.jl", "src/event_analysis/types.jl",
             "src/event_analysis/io.jl", "src/event_analysis/production.jl",
             "src/event_analysis/production_readout.jl", "src/event_analysis/grid_readout.jl",
             "src/ion_signal/event_io.jl", "scripts/production/specs.jl",
             "scripts/production/run_detector_pilot.jl",
             "scripts/production/run_electron_250ns.jl")
    return content_hash(join((content_hash(joinpath(ROOT, file); file=true)
                              for file in files), "|"))
end

"""Resolve the saved realisation and new photon files without touching ion files."""
function electron_paths(spec, sample::String, code_hash::String, n::Int)
    r = spec.realisation.data
    source = spec.sample_paths[sample]
    source_sha = content_hash(source; file=true)
    realisation_values = Dict(k => v for (k, v) in r if k != "samples")
    stream_values = Dict(k => v for (k, v) in realisation_values
                         if k != "max_events_per_sample")
    stream_identity = content_hash(join((sample, source_sha,
                                         canonical(stream_values)), "|"))
    realisation_id = content_hash(join((stream_identity,
                                        canonical(realisation_values)), "|"))
    identity_tail = n == Int(r["max_events_per_sample"]) ? "" : "|limit=$n"
    electron_id = content_hash(join((realisation_id,
                                     canonical(spec.electrons.data), code_hash), "|") *
                               identity_tail)
    base = joinpath(spec.results_dir, sample, realisation_id)
    output_root = n == Int(r["max_events_per_sample"]) ? base :
                  joinpath(ROOT, "work/electron_250ns_pilot", sample, realisation_id)
    photons_dir = joinpath(output_root, "detector", "photons", electron_id)
    return (; source, source_sha, stream_identity, realisation_id, electron_id,
            realisation=joinpath(base, "realisation.h5"),
            pure_xe=joinpath(photons_dir, "photons_pure_xe.h5"),
            xe_he15=joinpath(photons_dir, "photons_xe_he15.h5"))
end

"""Check the original 8,000-event primary realisation before replay."""
function check_common_realisation(file::HDF5.File, paths, expected_events::Int,
                                  seed::Int)::Nothing
    String(read(attributes(file)["schema"])) ==
        "nausicaa.event_realisation.pilot.v1" || error("Realisation schema mismatch")
    m = file["metadata"]
    Bool(read(m["complete"])) &&
        Int(read(m["event_count"])) == expected_events ||
        error("Primary realisation is incomplete")
    String(read(m["source_sha256"])) == paths.source_sha ||
        error("Source content differs from the saved realisation")
    String(read(m["card_hash"])) == paths.realisation_id ||
        error("Realisation card differs")
    String(read(m["stream_identity"])) == paths.stream_identity ||
        error("Primary RNG stream differs")
    Int(read(m["seed"])) == seed || error("Primary seed differs")
    return nothing
end

"""Check an existing completed photon file and recover its aggregate ledgers."""
function completed_photon_totals(path::String, paths, gas::String, n::Int)
    return h5open(path, "r") do f
        String(read(attributes(f)["schema"])) ==
            "nausicaa.photon_sipm_file.pilot.v1" || error("Photon schema mismatch: $path")
        m = f["metadata"]
        Bool(read(m["complete"])) && Int(read(m["event_count"])) == n &&
            String(read(m["card_hash"])) == paths.electron_id &&
            String(read(m["branch"])) == gas &&
            String(read(m["source_sha256"])) == paths.source_sha ||
            error("Existing photon output does not match this run: $path")
        return (; pairs=Int64(read(m["initial_pairs"])),
                survivors=Int64(read(m["surviving_electrons"])),
                emitted=Int64(read(m["emitted_photons"])),
                incident=Int64(read(m["incident_photons"])),
                elapsed_s=Float64(read(m["elapsed_s"])))
    end
end

"""Replay one sample at the configured electron bucket width for both gases."""
function run_electron_sample(spec, sample::String, code_hash::String, n::Int)
    r = spec.realisation.data
    e = spec.electrons.data
    paths = electron_paths(spec, sample, code_hash, n)
    isfile(paths.realisation) || error("Missing common realisation: $(paths.realisation)")
    existing = (isfile(paths.pure_xe), isfile(paths.xe_he15))
    if all(existing)
        xe = completed_photon_totals(paths.pure_xe, paths, "pure_xe", n)
        he = completed_photon_totals(paths.xe_he15, paths, "xe_he15", n)
        xe.pairs == he.pairs && xe.survivors == he.survivors ||
            error("Existing gas branches disagree on common primary totals")
        println("$sample: reused complete 250 ns photon files")
        flush(stdout)
        return (; sample, paths, events=n, pairs=xe.pairs,
                survivors=xe.survivors, xe_emitted=xe.emitted,
                xe_incident=xe.incident, he_emitted=he.emitted,
                he_incident=he.incident, elapsed_s=xe.elapsed_s)
    end
    any(existing) && error("Only one gas output exists for $sample; inspect the partial run")
    mkpath(dirname(paths.pure_xe))
    temporary = (pure_xe=paths.pure_xe * ".tmp", xe_he15=paths.xe_he15 * ".tmp")
    for path in values(temporary)
        isfile(path) && rm(path)
    end
    seed = Int(r["seed"])
    started = time_ns()
    pairs_total = Int64(0)
    survivors_total = Int64(0)
    xe_emitted = Int64(0)
    xe_incident = Int64(0)
    he_emitted = Int64(0)
    he_incident = Int64(0)
    h5open(paths.source, "r") do source
        ids = read(source["events/event_id"])
        length(ids) >= n || error("Source has fewer than $n events")
        h5open(paths.realisation, "r") do common
            check_common_realisation(common, paths, Int(r["max_events_per_sample"]), seed)
            h5open(temporary.pure_xe, "w") do fx
                h5open(temporary.xe_he15, "w") do fh
                    for (file, gas) in ((fx, "pure_xe"), (fh, "xe_he15"))
                        write_header(file, "nausicaa.photon_sipm_file.pilot.v1",
                                     sample, paths.source, paths.source_sha,
                                     paths.electron_id, code_hash, seed, gas,
                                     read(spec.realisation.path, String),
                                     read(spec.electrons.path, String))
                        file["metadata/stream_identity"] = paths.stream_identity
                    end
                    for (j, raw_id) in enumerate(ids[1:n])
                        id = Int64(raw_id)
                        hits, _ = read_event(source, id)
                        saved = common["events/$(id)"]
                        pairs = sample_fano_pairs(
                            Xoshiro(event_seed(paths.stream_identity, id, "pairs")),
                            hits, Float32(r["w_eV"]), Float32(r["fano_factor"]))
                        pair_count = foldl(Base.Checked.checked_add, pairs; init=Int64(0))
                        checksum = pair_checksum(pairs)
                        pair_count == Int64(read(saved["initial_pairs"])) &&
                            checksum == String(read(saved["pair_checksum"])) ||
                            error("Primary-pair replay failed at $sample event $id")
                        placement = EventPlacement(id,
                            Float32(read(saved["ion_drift_cm"])),
                            Float32(read(saved["electron_drift_cm"])))
                        isapprox(placement.ion_drift_cm + placement.electron_drift_cm,
                                 Float32(r["chamber_length_cm"]); atol=1f-4) ||
                            error("Drift lengths fail to sum to chamber length")
                        rx, ry = transverse_reference(hits)
                        phase = sample_grid_offset_mm(
                            Xoshiro(event_seed(paths.stream_identity, id, "electron_phase")),
                            Float32(e["hole_pitch_mm"]))
                        bucket_mm = Float32(e["drift_speed_mm_us"] * e["bucket_width_us"])
                        grid = GridSpec(Float32(e["hole_pitch_mm"]), bucket_mm,
                                        (rx, ry, maximum(hits.z_mm)), phase, Int8(-1))
                        drift_ms = Float32(10 * placement.electron_drift_cm /
                                           (1000 * e["drift_speed_mm_us"]))
                        survivors = sample_attached_electrons(
                            Xoshiro(event_seed(paths.stream_identity, id, "attachment")),
                            pairs, drift_ms, Float32(e["lifetime_ms"]))
                        survivor_count = foldl(Base.Checked.checked_add, survivors;
                                               init=Int64(0))
                        pairs_total = Base.Checked.checked_add(pairs_total, pair_count)
                        survivors_total = Base.Checked.checked_add(survivors_total,
                                                                  survivor_count)
                        for (file, gas) in ((fx, "pure_xe"), (fh, "xe_he15"))
                            model = e[gas]
                            scale = sqrt(placement.electron_drift_cm /
                                         Float32(e["pressure_bar"]))
                            sigma_t = Float32(model["d_T_mm_sqrtbar_per_sqrtcm"]) * scale
                            sigma_l = Float32(model["d_L_mm_sqrtbar_per_sqrtcm"]) * scale
                            yield_per_e = el_photon_yield(
                                Float32(model["yield_gain_photons_per_electron_kV"]),
                                Float32(model["yield_threshold_kV_cm_bar"]),
                                Float32(e["reduced_field_kV_cm_bar"]),
                                Float32(e["pressure_bar"]), Float32(e["hole_depth_mm"]))
                            acceptance = axial_light_acceptance(
                                Float32(e["hole_diameter_mm"]) / 2f0,
                                Float32(e["hole_depth_mm"]))
                            history = collect_photon_sipm_history(
                                Xoshiro(event_seed(paths.stream_identity, id,
                                                  "$(gas)_diffusion_el")),
                                hits, survivors, placement, grid;
                                sigma_t_mm=sigma_t, sigma_l_mm=sigma_l,
                                yield_per_e=yield_per_e,
                                J=Float32(e["j_variance_ratio"]),
                                acceptance=acceptance, gas_model=Symbol(gas),
                                initial_pairs=pair_count)
                            group = create_group(file["events"], string(id))
                            write_photon_history(group, history)
                            group["pair_checksum"] = checksum
                            incident = foldl(Base.Checked.checked_add,
                                             history.cells.count; init=Int64(0))
                            if gas == "pure_xe"
                                xe_emitted = Base.Checked.checked_add(
                                    xe_emitted, history.emitted_photons)
                                xe_incident = Base.Checked.checked_add(xe_incident, incident)
                            else
                                he_emitted = Base.Checked.checked_add(
                                    he_emitted, history.emitted_photons)
                                he_incident = Base.Checked.checked_add(he_incident, incident)
                            end
                        end
                        if j % 250 == 0 || j == n
                            println("$sample photons 250 ns: $j/$n events")
                            flush(stdout)
                        end
                    end
                    elapsed_s = (time_ns() - started) / 1e9
                    for (file, gas, emitted, incident) in
                        ((fx, "pure_xe", xe_emitted, xe_incident),
                         (fh, "xe_he15", he_emitted, he_incident))
                        m = file["metadata"]
                        m["initial_pairs"] = pairs_total
                        m["surviving_electrons"] = survivors_total
                        m["emitted_photons"] = emitted
                        m["incident_photons"] = incident
                        m["event_count"] = n
                        m["elapsed_s"] = elapsed_s
                        write(m["complete"], true)
                    end
                end
            end
        end
    end
    mv(temporary.pure_xe, paths.pure_xe)
    mv(temporary.xe_he15, paths.xe_he15)
    elapsed_s = (time_ns() - started) / 1e9
    println("$sample photons 250 ns: $n events in $(round(elapsed_s; digits=2)) s")
    flush(stdout)
    return (; sample, paths, events=n, pairs=pairs_total,
            survivors=survivors_total, xe_emitted, xe_incident,
            he_emitted, he_incident, elapsed_s)
end

"""Run both gas branches for the selected samples and save a numerical index."""
function main_electrons(card::String, limit::Union{Nothing,Int}=nothing)
    spec = load_specs(ROOT, card)
    Float64(spec.electrons.data["bucket_width_us"]) == 0.25 ||
        error("This runner requires 0.25 µs electron buckets")
    maximum_events = Int(spec.realisation.data["max_events_per_sample"])
    n = isnothing(limit) ? maximum_events : limit
    1 <= n <= maximum_events || error("Limit must be in 1:$maximum_events")
    hash = electron_replay_code_hash()
    run_id = content_hash(join((canonical(spec.run_spec),
                                canonical(spec.configuration),
                                canonical(spec.realisation.data),
                                canonical(spec.electrons.data), hash, string(n)), "|"))[1:16]
    output_dir = n == maximum_events ? spec.results_dir :
                 joinpath(ROOT, "work/electron_250ns_pilot")
    summary_path = joinpath(output_dir, "electron_250ns_summary_$(run_id).h5")
    results = [run_electron_sample(spec, sample, hash, n) for sample in
               get(spec.run_spec, "selected_samples", ["bb0nu", "single_electrons"])]
    mkpath(dirname(summary_path))
    temporary = summary_path * ".tmp"
    isfile(temporary) && rm(temporary)
    h5open(temporary, "w") do f
        attributes(f)["schema"] = "nausicaa.electron_rebin_summary.v1"
        m = create_group(f, "metadata")
        m["complete"] = false
        m["run_id"] = run_id
        m["run_card"] = read(spec.run_path, String)
        m["realisation_card"] = read(spec.realisation.path, String)
        m["electron_card"] = read(spec.electrons.path, String)
        m["code_hash"] = hash
        m["events_per_sample"] = n
        m["bucket_width_us"] = Float64(spec.electrons.data["bucket_width_us"])
        m["bucket_width_mm"] = Float64(spec.electrons.data["drift_speed_mm_us"] *
                                        spec.electrons.data["bucket_width_us"])
        for result in results
            group = create_group(f, result.sample)
            group["events"] = result.events
            group["initial_pairs"] = result.pairs
            group["surviving_electrons"] = result.survivors
            group["realisation_path"] = result.paths.realisation
            group["realisation_id"] = result.paths.realisation_id
            group["source_sha256"] = result.paths.source_sha
            group["stream_identity"] = result.paths.stream_identity
            group["electron_id"] = result.paths.electron_id
            group["pure_xe_path"] = result.paths.pure_xe
            group["xe_he15_path"] = result.paths.xe_he15
            group["pure_xe_emitted"] = result.xe_emitted
            group["pure_xe_incident"] = result.xe_incident
            group["xe_he15_emitted"] = result.he_emitted
            group["xe_he15_incident"] = result.he_incident
            group["elapsed_s"] = result.elapsed_s
            result.xe_incident <= result.xe_emitted &&
                result.he_incident <= result.he_emitted &&
                result.survivors <= result.pairs || error("Electron ledger failed")
        end
        write(m["complete"], true)
    end
    mv(temporary, summary_path; force=true)
    println("Summary: $summary_path")
    return summary_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (1, 2) ||
        error("Usage: julia --project=. run_electron_250ns.jl RUN.toml [event_limit]")
    main_electrons(ARGS[1], length(ARGS) == 2 ? parse(Int, ARGS[2]) : nothing)
end
