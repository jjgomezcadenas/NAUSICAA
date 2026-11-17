#!/usr/bin/env julia
#=
eITACA pilot: merge the batches of each pilot run and compute the saturation
summary that results/eitaca/pilot/index.html shows.

  julia --project=. scripts/eitaca/prepare_eitaca_pilot.jl [results/eitaca/pilot]

Steps:

 1. Read every batches/*/electrons.h5 of each run folder (IE_pilot_<sample>[_guides][_no_exclusion]).
 2. Concatenate the event tables; add the selection counts, the histograms of
    expected photoelectrons per SiPM and bucket, and the timing.
 3. Compute, per run: the distribution of the event peak (largest expected
    photoelectrons in one SiPM and one 250 ns bucket), the fraction of cells
    with photons beyond the ASIC limit, the events with clipped buckets, the
    charge they lose, the recovery loss and the reconstructed energy.
 4. Pair each eITACA run with its guide reference by event number: ratio of
    the event peaks and of the detected photoelectrons.
 5. Read the hybrid production of the same events (results/readout/
    <sample>_itaca_pilot_10x100) for its count of invalid electron buckets.

Output: <pilot>/pilot_summary.h5, read by scripts/eitaca/present_eitaca_pilot.py.
=#
using HDF5, Statistics, TOML, SHA
const ROOT = normpath(joinpath(@__DIR__, "../.."))
const SAMPLES = ("bb0nu", "single_electrons")
const ARMS = ("", "_guides", "_no_exclusion", "_guides_no_exclusion")
"Edges of the histogram of event peaks, expected photoelectrons: 10 to 10^5, 20 bins per decade."
const PEAK_EDGES = 10 .^ range(1, 5; length = 81)

"""Read and concatenate the batches of one run folder."""
function read_run(dir)
    files = sort([joinpath(dir, "batches", b, "electrons.h5") for b in readdir(joinpath(dir, "batches")) if isfile(joinpath(dir, "batches", b, "electrons.h5"))])
    isempty(files) && error("No batches in $dir")
    ev = Dict{String,Vector}(); counts = Dict{String,Int}(); hist = nothing; edges = nothing; below = 0
    timing = Dict{String,Float64}(); config = Dict{String,Any}(); shas = String[]
    for p in files
        h5open(p, "r") do f
            read(f["schema"]) == "eitaca.electrons.v1" || error("Unknown schema in $p")
            for k in keys(f["events"]); append!(get!(ev, k, Any[]), read(f["events"][k])); end
            for k in keys(f["selection/counts"]); counts[k] = get(counts, k, 0) + read(f["selection/counts"][k]); end
            h = read(f["cell_pe_histogram/counts"]); hist = hist === nothing ? h : hist .+ h
            edges = read(f["cell_pe_histogram/edges_pe"]); below += read(f["cell_pe_histogram/below_first_edge"])
            for k in keys(f["timing"]); timing[k] = get(timing, k, 0.0) + Float64(read(f["timing"][k])); end
            g = f["configuration"]
            for k in ("optics", "light_acceptance", "photons_per_electron", "pde", "asic_limit_pe", "illuminated_microcells", "recovery_s", "bucket_s", "sample")
                config[k] = read(g[k])
            end
            config["el_exclusion_mm"] = TOML.parse(read(g["resolved"]))["detector"]["trigger"]["el_exclusion_mm"]
            push!(shas, bytes2hex(open(sha256, p)))
        end
    end
    ev = Dict(k => [x for x in v] for (k, v) in ev)
    (; events = ev, counts, hist, edges, below, timing, config, files, shas)
end

quantile_or_nan(x, q) = isempty(x) ? NaN : quantile(x, q)

"""Saturation summary of one run."""
function summary(run)
    e = run.events; c = run.config; limit = c["asic_limit_pe"]; n = length(e["event_number"])
    peak = Float64.(e["peak_expected_pe"]); clipped = e["clipped_buckets"] .> 0
    beyond_cells = sum(run.hist[k] for k in eachindex(run.hist) if run.edges[k] >= limit; init = 0)    # bins entirely above the limit
    fe = e["electron_energy_unclipped_keV"]; kept = e["electron_energy_clipped_kept_keV"]
    occupancy = peak .* c["recovery_s"] ./ (c["illuminated_microcells"] * c["bucket_s"])
    ok = .!clipped
    ratio = e["electron_energy_valid_keV"][ok] ./ e["energy_true_keV"][ok]
    loss = 1 .- kept[clipped] ./ fe[clipped]
    Dict{String,Any}(
        "input_events" => run.counts["input_events"], "accepted_events" => n,
        "main_track_near_el_events" => run.counts["main_track_near_el_events"],
        "seconds_per_accepted_event" => run.timing["accepted_events_s"] / max(n, 1),
        "seconds_per_input_event" => run.timing["total_s"] / max(run.counts["input_events"], 1),
        "asic_limit_pe" => limit, "light_acceptance" => c["light_acceptance"], "photons_per_electron" => c["photons_per_electron"],
        "peak_pe_mean" => mean(peak), "peak_pe_median" => median(peak), "peak_pe_q99" => quantile_or_nan(peak, 0.99), "peak_pe_max" => maximum(peak),
        "peak_pe_over_limit_max" => maximum(peak) / limit,
        "cells_with_photons" => sum(e["positive_cells"]), "cells_clipped_with_photons" => sum(e["clipped_cells_with_photons"]),
        "cells_clipped_fraction" => sum(e["clipped_cells_with_photons"]) / sum(e["positive_cells"]),
        "cells_expected_beyond_limit" => beyond_cells,
        "events_with_clipping" => count(clipped), "events_with_clipping_fraction" => count(clipped) / n,
        "events_with_invalid_energy" => count(isnan, e["electron_energy_keV"]),
        "sample_charge_retained_fraction" => 1 - sum(e["asic_charge_lost_pe"]) / sum(e["sipm_charge_pe"]),
        "clipped_event_charge_retained_min" => any(clipped) ? minimum(e["charge_retained_fraction"][clipped]) : 1.0,
        "clipped_event_energy_loss_mean" => any(clipped) ? mean(loss) : 0.0, "clipped_event_energy_loss_max" => any(clipped) ? maximum(loss) : 0.0,
        "recovery_loss_mean" => mean(e["recovery_loss_fraction"]), "recovery_loss_max" => maximum(e["recovery_loss_fraction"]),
        "peak_occupancy_max" => maximum(occupancy),
        "clipped_event_distance_to_el_max_mm" => any(clipped) ? maximum(e["track_distance_to_el_mm"][clipped]) : NaN,
        "clipped_event_drift_max_cm" => any(clipped) ? maximum(e["electron_drift_cm"][clipped]) : NaN,
        # accepted events as close to the EL as the farthest clipped one: those a larger exclusion would remove
        "events_within_clipped_distance" => any(clipped) ? count(<=(maximum(e["track_distance_to_el_mm"][clipped])), e["track_distance_to_el_mm"]) : 0,
        "detected_pe_mean" => mean(Float64.(e["detected_pe"])),
        "energy_ratio_unclipped_mean" => mean(ratio), "energy_ratio_unclipped_std" => std(ratio; corrected = false),
        "energy_ratio_unclipped_entries" => length(ratio))
end

"""Per-event ratios of a run to its reference, paired by event number."""
function paired(a, b)
    ia = Dict(id => j for (j, id) in enumerate(a.events["event_number"])); ids = [id for id in b.events["event_number"] if haskey(ia, id)]
    ib = Dict(id => j for (j, id) in enumerate(b.events["event_number"]))
    pa = [a.events["peak_expected_pe"][ia[i]] for i in ids]; pb = [b.events["peak_expected_pe"][ib[i]] for i in ids]
    da = [a.events["detected_pe"][ia[i]] for i in ids]; db = [b.events["detected_pe"][ib[i]] for i in ids]
    (ids = ids, peak_ratio = Float64.(pa ./ pb), detected_ratio = Float64.(da ./ db))
end

"""Invalid electron buckets of the hybrid production of the same events."""
function hybrid_reference(sample)
    n = 0; invalid = 0; nan = 0; files = String[]
    for b in 0:9
        p = joinpath(ROOT, "results/readout/$(sample)_itaca_pilot_10x100/batches/b$(lpad(b, 5, '0'))/vh1_hits/hits.h5")
        isfile(p) || continue
        push!(files, relpath(p, ROOT))
        h5open(p, "r") do f
            v = read(f["events/invalid_electron_buckets"]); n += length(v); invalid += count(>(0), v); nan += count(isnan, read(f["events/electron_energy_keV"]))
        end
    end
    (events = n, events_with_invalid_electron_buckets = invalid, events_with_nan_energy = nan, files = files)
end

histogram(x, edges) = (h = zeros(Int64, length(edges) - 1); for v in x; k = searchsortedlast(edges, v); 1 <= k < length(edges) && (h[k] += 1); end; h)

function main(pilot = joinpath(ROOT, "results/eitaca/pilot"))
    target = joinpath(pilot, "pilot_summary.h5"); isfile(target) && error("Preserve existing summary: $target")
    tmp = target * ".tmp"
    h5open(tmp, "w") do f
        f["schema"] = "eitaca.pilot_summary.v1"; f["peak_edges_pe"] = collect(PEAK_EDGES)
        f["source"] = read(@__FILE__, String)
        for s in SAMPLES
            runs = Dict(arm => read_run(joinpath(pilot, "IE_pilot_$(s)$(arm)")) for arm in ARMS)
            for arm in ARMS
                run = runs[arm]; g = create_group(f, "$(s)$(arm)")
                for (k, v) in summary(run); g["summary/$k"] = v; end
                for (k, v) in run.events; g["events/$k"] = v; end
                g["cell_histogram/edges_pe"] = run.edges; g["cell_histogram/counts"] = run.hist
                g["peak_histogram/counts"] = histogram(run.events["peak_expected_pe"], PEAK_EDGES)
                g["optics"] = run.config["optics"]; g["el_exclusion_mm"] = Float64(run.config["el_exclusion_mm"])
                g["files"] = relpath.(run.files, ROOT); g["sha256"] = run.shas
            end
            for (a, b) in (("", "_guides"), ("_no_exclusion", "_guides_no_exclusion"))
                p = paired(runs[a], runs[b]); g = create_group(f, "$(s)$(a)/paired")
                g["event_number"] = p.ids; g["peak_ratio"] = p.peak_ratio; g["detected_ratio"] = p.detected_ratio
                g["peak_ratio_mean"] = mean(p.peak_ratio); g["peak_ratio_median"] = median(p.peak_ratio)
                g["detected_ratio_mean"] = mean(p.detected_ratio)
            end
            h = hybrid_reference(s); g = create_group(f, "hybrid_production/$s")
            g["events"] = h.events; g["events_with_invalid_electron_buckets"] = h.events_with_invalid_electron_buckets
            g["events_with_nan_energy"] = h.events_with_nan_energy; g["files"] = isempty(h.files) ? [""] : h.files
        end
    end
    Base.rename(tmp, target); println("OUTPUT $target")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(isempty(ARGS) ? joinpath(ROOT, "results/eitaca/pilot") : abspath(ARGS[1]))
