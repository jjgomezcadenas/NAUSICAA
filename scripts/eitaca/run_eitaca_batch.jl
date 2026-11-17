#!/usr/bin/env julia
#=
eITACA: read one batch of Monte Carlo events, simulate the electron readout
and write the electron file of the batch.

  julia --project=. scripts/eitaca/run_eitaca_batch.jl <run card> --batch=N [--limit=M] [--output=DIR]

For every event of the batch:

 1. Fast trigger and fixed frame from the Monte Carlo (EventAnalysis.mc_fast_trigger),
    with the exclusion of tracks near the EL entry. With trigger.frame_track =
    "leading_gamma_electron" (gamma samples) the main track is the electron of
    highest energy produced by the generated gamma, with its descendants, as
    in scripts/iitaca/run_iitaca_batch.jl.
 2. Box of the response: main track and the gamma deposits inside the frame, with the margin.
 3. Electron-ion pairs of each deposit (W, Fano); attachment, diffusion to the EL
    holes and 250 ns buckets, EL photons and blue photons at the SiPM
    (DetectorResponse.sim_electron_detector_response), with the eITACA transport,
    EL yield and optics (EItaca).
 4. SiPM microcells and ASIC (ElectronicsResponse.sim_sipm_readout), with the
    SiPM and ASIC cards of the hybrid production. For an event with clipped
    buckets the readout is repeated with the same random streams and no ASIC
    limit, which gives the charge the clipping removed.
 5. Reconstruction: recovery inverse (reco_electron_signal), event energy,
    rebinning and hits (rebin_electron_signal, reconstruct_hits).

Saturation measures of each event: the largest blue-photon count and expected
photoelectrons (photons x PDE) in one SiPM and one bucket, the ASIC clipped
buckets and the charge they lose, the recovery loss of the microcells, and the
event energy with and without clipping. The cells with photons are also
histogrammed in expected photoelectrons.

Output: results/eitaca/<run_label>/batches/<tag>/electrons.h5 (schema
eitaca.electrons.v1). No ion is simulated. --limit reads only the first M
events and writes under <run_label>_first<M>.
=#
using HDF5, TOML, Random, SHA, Statistics, LinearAlgebra
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl")); using .EventAnalysis
include(joinpath(ROOT, "src/detector_response/DetectorResponse.jl")); using .DetectorResponse
include(joinpath(ROOT, "src/ion_signal/guarded_files.jl"))
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
include(joinpath(ROOT, "src/eitaca/EItaca.jl")); using .EItaca
include(joinpath(ROOT, "src/electronics/ElectronicsResponse.jl")); using .ElectronicsResponse
include(joinpath(ROOT, "src/reconstruction/SignalReconstruction.jl")); using .SignalReconstruction
include(joinpath(ROOT, "scripts/production/batch_inputs.jl"))

const SCHEMA = "eitaca.electrons.v1"
const CODE_DIRECTORIES = ("src/event_analysis", "src/detector_response", "src/electronics", "src/reconstruction", "src/eitaca", "src/itaca")
"Edges of the histogram of expected photoelectrons per SiPM and bucket: 0.01 to 10^5 PE, 20 bins per decade."
const PE_EDGES = 10 .^ range(-2, 5; length = 141)

"""Fast-trigger settings from the resolved detector."""
trigger_settings(d) = (g = d["geometry"]; t = d["trigger"]; FixedFrameTriggerSettings(Float32(10g["fiducial_length_cm"]), Float32(5g["fiducial_diameter_cm"]),
    Float32(t["frame_side_mm"]), Float32(t["gamma_energy_cut_keV"]), t["frame_placement"] == "inside_detector", Float32(t["el_exclusion_mm"])))

"""Electron-response settings of one event.

The SiPM lattice has a random offset per event (stream `electron_grid`, as in
the hybrid production). Spreads D sqrt(L) for the drift length of the event
(`electron_spreads_mm`); EL yield, optics and lifetime from the resolved
detector. The ion entries are not used: the ion lattice is the SiPM lattice
and the ion velocity and bucket are NaN, the ion spreads zero, which the
response requires to be finite. The gas label Xe_85_He_15 is the closest
entry of `GasMixture`; it is a label and enters no calculation.
"""
function electron_event_settings(d, seed, p)
    F = Float32; pitch = F(d["el"]["hole_pitch_mm"])
    grid = ReadoutGrid(pitch, sample_grid_offset_mm(Xoshiro(event_seed(seed, p.event_id, "electron_grid")), pitch))
    s = electron_spreads_mm(d, p.electron_drift_cm)
    DetectorResponseSettings(F(10d["geometry"]["fiducial_length_cm"]), F(d["energy"]["W_eV"]), F(d["energy"]["fano_factor"]), Xe_85_He_15, grid, grid,
        NaN32, F(d["sipm"]["electron_bucket_us"]), IonTransportSettings(NaN32, 0f0, 0f0),
        ElectronTransportSettings(F(d["drift"]["electron_drift_velocity_mm_us"]), F(s.transverse), F(s.longitudinal)),
        F(d["drift"]["electron_lifetime_ms"]), el_photons_per_electron(d), F(d["el"]["variance_ratio_J"]), light_acceptance(d), Symbol(d["el"]["photon_band"]))
end

"ASIC current limit of the repeated readout without clipping: 10^6 times the card limit."
const UNLIMITED_FACTOR = 1e6

"""SiPM and ASIC settings from the response cards, as `prepare_electronics_run`
in scripts/production/framed_electronics_settings.jl (its SiPM part); `limit_A`
replaces the ASIC current limit."""
function sipm_settings(sc, ac; limit_A = ac["current_limit_A"])
    SiPMReadoutSettings(SiPM(sc["pde"], sc["microcell_count"], sc["gain"], sc["recovery_s"], sc["dark_reference_rate_Hz"]),
        Int64(sc["illuminated_microcell_count"]), SiPMStatisticalOptions(sc["statistical_substeps"], sc["sparse_trigger_limit"];
        initial_state = Symbol(get(sc, "initial_state", "stationary_dark"))), SiPmAsic(ac["bucket_width_s"], ac["noise_rms_pe"], limit_A))
end

"""Event energy in keV from reconstructed PE: attachment correction exp(t/tau)
x valid PE sum x W / (PDE x photons per electron x acceptance), as
`event_charge_estimates` (src/reconstruction/hits.jl). Returns the energy of
the valid buckets and the count of invalid ones."""
function electron_energy(fine, p, ds, pde)
    correction = exp(10Float64(p.electron_drift_cm) / ds.electron_transport.electron_drift_velocity_mm_us / 1000 / ds.electron_lifetime_ms)
    scale = Float64(ds.W_eV) / 1000 / (pde * ds.el_photons_per_electron * ds.optical_acceptance)
    ok = fine.valid .& isfinite.(fine.pe)
    (energy_keV = correction * sum(Float64, fine.pe[ok]; init = 0.0) * scale, invalid = count(.!ok), correction, keV_per_pe = scale)
end

"""SHA-256 of the driver and of every Julia file of the modules it uses, keyed by path."""
function code_hashes()
    files = [@__FILE__, joinpath(ROOT, "scripts/production/batch_inputs.jl")]
    for dir in CODE_DIRECTORIES, (folder, _, names) in walkdir(joinpath(ROOT, dir)), name in names
        endswith(name, ".jl") && push!(files, joinpath(folder, name))
    end
    Dict(relpath(f, ROOT) => bytes2hex(open(sha256, f)) for f in files)
end

columns(g, rows) = (isempty(rows) || for k in keys(first(rows)); g[string(k)] = [r[k] for r in rows]; end)

"""
    run_eitaca_batch(card; batch, limit=0, output_root=joinpath(ROOT, "results/eitaca"))::String

Process one batch of the run card and write its electron file; return the path.
An existing file is preserved: the function stops.
"""
function run_eitaca_batch(card::AbstractString; batch::Integer, limit::Integer = 0, output_root::AbstractString = joinpath(ROOT, "results/eitaca"))
    BLAS.set_num_threads(1)
    started = time(); r = resolve_run(card; root = ROOT); c = r.run; select_batch!(c, batch); d = r.resolved.detector
    limit > 0 && (c["max_events"] = limit; c["run_label"] = replace(c["run_label"], "/batches/" => "_first$(limit)/batches/"))
    out = joinpath(output_root, c["run_label"]); target = joinpath(out, "electrons.h5")
    isfile(target) && error("Preserve existing output: $target")
    hp = input_path(ROOT, c["hits"]); pp = input_path(ROOT, c["particles"])
    source_sha = bytes2hex(open(sha256, hp)); particle_sha = bytes2hex(open(sha256, pp))
    ids = h5open(hp, "r") do f
        read(f["metadata/source_sha256"]) == particle_sha || error("Particle source mismatch")
        read(f["events/event_id"])
    end
    get(c, "max_events", 0) > 0 && (ids = ids[1:min(length(ids), c["max_events"])])
    parts = read_mc_particles(pp, ids)
    leading = d["trigger"]["frame_track"] == "leading_gamma_electron" ? read_leading_gamma_electrons(pp, ids) : nothing
    trigger = trigger_settings(d); L = trigger.chamber_length_mm; margin = Float32(d["trigger"]["simulation_margin_mm"])
    W = Float32(d["energy"]["W_eV"]); fano = Float32(d["energy"]["fano_factor"]); seed = c["seed"]; eseed = c["electronics_seed"]
    max_bytes = get(c, "max_vector_bytes", 268435456); rc = d["reconstruction"]
    sc, ac = response_cards(d, ROOT); settings = sipm_settings(sc, ac); unlimited = sipm_settings(sc, ac; limit_A = UNLIMITED_FACTOR * ac["current_limit_A"])
    es = ElectronReconstructionSettings(settings.sipm, settings.illuminated_microcells, rc["subtract_mean_dark"])
    pde = settings.sipm.pde; limit_pe = asic_limit_pe(sc, ac)
    stats = FixedFrameStatistics(); decisions = NamedTuple[]; events = NamedTuple[]
    histogram = zeros(Int64, length(PE_EDGES) - 1); below = 0
    x = Float32[]; y = Float32[]; z = Float32[]; e = Float32[]; timing = zeros(6)
    h5open(hp, "r") do f
        for (row, id) in enumerate(ids)
            local h, info, p, o
            event_started = time()
            timing[1] += @elapsed begin
                h, info = read_mc_event(f, id; event_row = row)
                leading === nothing || (info = leading_electron_info(info, leading[id]))
                h = center_hits(h, info); p = placement_from_vertex(info, L); o = mc_fast_trigger(h, info, parts[id], p, trigger)
                stats = update_statistics(stats, o)
                push!(decisions, (event_number = id, decision = Int(o.decision), electron_drift_cm = p.electron_drift_cm))
            end
            o.decision == trigger_accepted || continue
            rng(tag) = Xoshiro(event_seed(seed, id, tag)); erng(tag) = Xoshiro(event_seed(eseed, id, tag))
            local box, selected, ds, photons, asic, ctl, free, fine, pe, hits, E_true, main
            timing[2] += @elapsed begin
                box = simulation_box(h, info, parts[id], o.frame, margin); selected = hits_in_frame(h, box)
                gm = gamma_map(parts[id], info)
                iszero(sum(h.energy_keV[j] for j in eachindex(h.x_mm) if gm[h.particle_id[j]] == 0 && !contains(box, (h.x_mm[j], h.y_mm[j], h.z_mm[j])); init = 0f0)) ||
                    error("Main-track energy outside the box in event $id")
                E_true = sum(Float64, selected.energy_keV; init = 0.0); main = main_track_box(h, info, parts[id])
                pairs = sample_fano_pairs(rng("pairs"), selected, W, fano); ds = electron_event_settings(d, seed, p)
                photons = sim_electron_detector_response(rng("electrons"), selected, pairs, p, box, ds; max_vector_bytes = max_bytes)
            end
            timing[3] += @elapsed asic, ctl = sim_sipm_readout(erng("sipm"), erng("asic"), photons, settings)
            # The same streams without the ASIC limit: the response the clipping removed
            timing[4] += @elapsed free = ctl.asic_clipped_buckets > 0 ? first(sim_sipm_readout(erng("sipm"), erng("asic"), photons, unlimited)) : asic
            local E, E_kept, E_free
            timing[5] += @elapsed begin
                fine = reco_electron_signal(asic, es); E = electron_energy(fine, p, ds, pde)
                E_kept = ctl.asic_clipped_buckets > 0 ? electron_energy(reco_electron_signal(SiPMASICResponse(asic.event_id, asic.channels, asic.time, asic.response_pe,
                    falses(size(asic.clipped))), es), p, ds, pde) : E
                E_free = ctl.asic_clipped_buckets > 0 ? electron_energy(reco_electron_signal(free, es), p, ds, pde) : E
                pe = rebin_electron_signal(fine, Int64(rc["rebin_factor"]))
                hits = reconstruct_hits(pe, p, ds; charge_correction = E.correction, keV_per_count = E.keV_per_pe, minimum_energy_keV = rc["electron_hit_cut_keV"])
            end
            v = photons.photons.counts; positive = 0; clipped_with_photons = 0
            @inbounds for k in eachindex(v)
                v[k] > 0 || continue
                positive += 1; clipped_with_photons += asic.clipped[k]
                q = pde * v[k]; b = searchsortedlast(PE_EDGES, q)
                1 <= b < length(PE_EDGES) ? (histogram[b] += 1) : (b < 1 && (below += 1))
            end
            peak = isempty(v) ? 0 : maximum(v); jpeak = isempty(v) ? 0 : argmax(v)
            free_peak = isempty(free.response_pe) ? 0f0 : maximum(free.response_pe)
            vz = 10p.ion_drift_cm                                   # vertex z, mm from the cathode
            push!(events, (event_number = id, energy_true_keV = E_true, electron_energy_keV = E.invalid == 0 ? E.energy_keV : NaN,
                electron_energy_valid_keV = E.energy_keV, electron_energy_clipped_kept_keV = E_kept.energy_keV, electron_energy_unclipped_keV = E_free.energy_keV,
                invalid_electron_buckets = E.invalid, electron_drift_cm = p.electron_drift_cm,
                track_distance_to_el_mm = Float64(L) - (main.upper_mm[3] + vz), track_z_length_mm = main.upper_mm[3] - main.lower_mm[3],
                sigma_t_mm = ds.electron_transport.sigma_t_mm, sigma_l_mm = ds.electron_transport.sigma_l_mm,
                e_ion_pairs = photons.e_ion_pairs, e_lost_attachment = photons.e_lost_attachment, EL_photons = photons.EL_photons,
                blue_photons = sum(v), photons_outside_acquisition = photons.losses.excluded_total,
                channels = size(v, 2), buckets = size(v, 1), positive_cells = positive,
                peak_photons = peak, peak_expected_pe = pde * peak, peak_bucket = jpeak == 0 ? 0 : jpeak[1],
                peak_asic_pe = isempty(asic.response_pe) ? 0f0 : maximum(asic.response_pe), peak_unclipped_asic_pe = free_peak,
                detected_pe = ctl.detected_pe, dark_pulses = ctl.dark_pulses, sipm_charge_pe = ctl.sipm_charge_pe,
                recovery_loss_fraction = 1 - ctl.sipm_charge_pe / max(ctl.detected_pe + ctl.dark_pulses, 1),
                clipped_buckets = ctl.asic_clipped_buckets, clipped_cells_with_photons = clipped_with_photons,
                asic_charge_lost_pe = ctl.asic_charge_lost_pe, charge_retained_fraction = 1 - ctl.asic_charge_lost_pe / max(ctl.sipm_charge_pe, 1e-300),
                hit_offset = length(x), hit_count = length(hits.x_mm), event_s = time() - event_started))
            append!(x, hits.x_mm); append!(y, hits.y_mm); append!(z, hits.z_mm); append!(e, hits.energy_keV)
            timing[6] += time() - event_started
            if length(events) % 10 == 0
                println("Accepted $(length(events)); input $row/$(length(ids)); $(round(time() - started; digits = 1)) s"); flush(stdout); GC.gc(false)
            end
        end
    end
    mkpath(out); budget = 16 * 1024^2 + 16 * length(x) + 1024 * length(events) + 64 * length(decisions)
    GuardedFiles.space_check(out, budget, 1024^3); tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = SCHEMA
            g = create_group(f, "configuration")
            g["run_label"] = c["run_label"]; g["sample"] = c["sample"]; g["variant"] = r.resolved.variant; g["batch"] = batch_tag(batch)
            g["seed"] = seed; g["electronics_seed"] = eseed; g["generated_events"] = Int(get(c, "generated_events_per_batch", 0))
            g["resolved"] = resolved_text(r, ROOT); g["run_card"] = read(card, String)
            g["sipm_card"] = read(joinpath(ROOT, d["sipm"]["response_card"]), String); g["asic_card"] = read(joinpath(ROOT, d["sipm"]["asic_card"]), String)
            for (name, text) in r.resolved.set.texts; g["cards/$name"] = text; end
            r.resolved.variant_text == "" || (g["cards/variant"] = r.resolved.variant_text)
            g["optics"] = d["el"]["optics"]; g["light_acceptance"] = Float64(light_acceptance(d)); g["photons_per_electron"] = Float64(el_photons_per_electron(d))
            g["pde"] = pde; g["asic_limit_pe"] = limit_pe; g["illuminated_microcells"] = settings.illuminated_microcells
            g["recovery_s"] = settings.sipm.recovery_s; g["bucket_s"] = settings.asic.bucket_width_s
            g["sha256/hits"] = source_sha; g["sha256/particles"] = particle_sha
            code = code_hashes(); names = sort(collect(keys(code))); g["code_files"] = names; g["code_sha256"] = [code[n] for n in names]
            columns(create_group(f, "selection/events"), decisions)
            for k in fieldnames(typeof(stats)); f["selection/counts/" * string(k)] = getfield(stats, k); end
            columns(create_group(f, "events"), events)
            hg = create_group(f, "cell_pe_histogram"); hg["edges_pe"] = collect(PE_EDGES); hg["counts"] = histogram; hg["below_first_edge"] = below
            attrs(hg)["definition"] = "SiPM-bucket cells with at least one blue photon, in expected photoelectrons = photons x PDE"
            t = create_group(f, "etrack")
            for (name, values) in (("x_mm", x), ("y_mm", y), ("z_mm", z), ("energy_keV", e))
                isempty(values) ? (t[name] = values) : (t[name, chunk = (min(length(values), 65536),), deflate = 4] = values)
            end
            attrs(t)["coordinates"] = "mm relative to the event vertex; z along the drift"
            for (name, value) in zip(("selection", "transport", "sipm_asic", "unclipped_repeat", "reconstruction", "accepted_events"), timing); f["timing/$(name)_s"] = value; end
            f["timing/total_s"] = time() - started; f["timing/input_events"] = length(ids); f["timing/accepted_events"] = length(events)
            flush(f); GuardedFiles.check_size(tmp, budget)
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    n = max(length(events), 1)
    println("TIMING input=$(length(ids)) accepted=$(length(events)) total_s=$(round(time() - started; digits = 1))")
    for (name, value) in zip(("selection", "transport", "sipm_asic", "unclipped_repeat", "reconstruction", "accepted_events"), timing)
        println("TIMING $name total_s=$(round(value; digits = 2)) per_accepted_event_s=$(round(value / n; digits = 4))")
    end
    println("OUTPUT $target bytes=$(filesize(target))")
    target
end

if abspath(PROGRAM_FILE) == @__FILE__
    batch, args = batch_argument(ARGS)
    k = findall(startswith("--limit="), args); limit = isempty(k) ? 0 : parse(Int, split(args[k[1]], "=")[2]); args = args[setdiff(eachindex(args), k)]
    k = findall(startswith("--output="), args); output_root = isempty(k) ? joinpath(ROOT, "results/eitaca") : abspath(split(args[k[1]], "=")[2])
    args = args[setdiff(eachindex(args), k)]
    (batch !== nothing && length(args) == 1) || error("Usage: run_eitaca_batch.jl run_card.toml --batch=N [--limit=M] [--output=DIR]")
    run_eitaca_batch(abspath(args[1]); batch, limit, output_root)
end
