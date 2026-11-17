#!/usr/bin/env julia
#=
iItaca: signal against backgrounds after the selection, with the cut on the endpoint energy.

  julia --project=. scripts/iitaca/prepare_iitaca_comparison.jl <output folder> <signal run folder> <background run folder>...

Reads statistics/events.h5 of each run (prepare_iitaca_statistics.jl) and
writes comparison.h5 in the output folder (schema iitaca.comparison.v1). The
runs must share the variant, the ROI and the track settings.

For each run: the cut flow and, for the events inside the ROI that pass the
track selection, E1 and E2, the endpoint energies of the main track (E1 >= E2).

Classical cut: an event is accepted when E2 exceeds a threshold. For each
signal efficiency of TARGETS the threshold is the corresponding quantile of the
signal E2; each background is given its accepted fraction with its binomial
error. The full curve (signal efficiency against background acceptance for
every threshold) and its area are saved too. Python makes the page
(scripts/iitaca/present_iitaca_comparison.py).
=#
using HDF5, Statistics
const ROOT = normpath(joinpath(@__DIR__, "../.."))
isdefined(@__MODULE__, :GuardedFiles) || include(joinpath(ROOT, "src/ion_signal/guarded_files.jl"))

const COMPARISON_SCHEMA = "iitaca.comparison.v1"
const TARGETS = (0.9, 0.8, 0.7, 0.6)

"""Selected events of a run and its identity, from statistics/events.h5."""
function read_run(run_dir::AbstractString)
    h5open(joinpath(run_dir, "statistics/events.h5"), "r") do f
        read(f["schema"]) == "iitaca.statistics.v1" || error("Unexpected schema in $run_dir")
        keep = (read(f["events/in_roi"]) .== 1) .& (read(f["events/selected"]) .== 1)
        column(k) = Float64.(read(f["events/$k"])[keep])
        (run_label = read(f["run_label"]), sample = read(f["sample"]), variant = read(f["variant"]), settings = read(f["settings"]),
         roi = (read(f["roi_lower_keV"]), read(f["roi_upper_keV"])), flow_step = read(f["cut_flow/step"]), flow_events = read(f["cut_flow/events"]),
         E1 = column("E1_keV"), E2 = column("E2_keV"), drift = column("ion_drift_cm"), event_number = read(f["events/event_number"])[keep])
    end
end

"""Fraction of `values` above `threshold` and its binomial error."""
function accepted(values::AbstractVector, threshold::Real)
    n = length(values); n == 0 && return (fraction = NaN, error = NaN, events = 0)
    k = count(>(threshold), values); p = k / n
    (fraction = p, error = sqrt(p * (1 - p) / n), events = k)
end

"""Threshold on E2 that keeps the fraction `target` of the signal: the event of rank ceil(target n) from the top lies just above it."""
function threshold_for(signal::AbstractVector, target::Real)
    s = sort(signal; rev = true); k = clamp(ceil(Int, target * length(s)), 1, length(s))
    k == length(s) ? prevfloat(s[k]) : (s[k] + s[k+1]) / 2
end

"""Signal efficiency and background acceptance for every threshold between consecutive values, and the area under the curve."""
function curve(signal::AbstractVector, background::AbstractVector)
    t = sort(unique(vcat(signal, background))); t = vcat(prevfloat(first(t)), (t[1:end-1] .+ t[2:end]) ./ 2, nextfloat(last(t)))
    eff = [count(>(x), signal) / length(signal) for x in t]; acc = [count(>(x), background) / length(background) for x in t]
    area = sum((acc[i] - acc[i+1]) * (eff[i] + eff[i+1]) / 2 for i in 1:length(t)-1)      # efficiency against acceptance, trapezoids
    (threshold = t, efficiency = eff, acceptance = acc, area = area)
end

function write_comparison(out_dir::AbstractString, signal_dir::AbstractString, background_dirs::Vector{<:AbstractString})
    signal = read_run(signal_dir); backgrounds = [read_run(d) for d in background_dirs]
    for b in backgrounds
        b.variant == signal.variant && b.roi == signal.roi && b.settings == signal.settings || error("Runs $(signal.run_label) and $(b.run_label) differ in variant, ROI or settings")
    end
    isempty(signal.E2) && error("No selected signal events")
    mkpath(out_dir); target = joinpath(out_dir, "comparison.h5")
    budget = 16 * 1024^2 + 64 * (length(signal.E2) + sum(length(b.E2) for b in backgrounds; init = 0)) * 4
    GuardedFiles.space_check(out_dir, budget, 1024^3); tmp = target * ".tmp"
    thresholds = [threshold_for(signal.E2, t) for t in TARGETS]
    try
        h5open(tmp, "w") do f
            f["schema"] = COMPARISON_SCHEMA; f["variant"] = signal.variant; f["settings"] = signal.settings
            f["roi_lower_keV"] = signal.roi[1]; f["roi_upper_keV"] = signal.roi[2]
            f["targets"] = collect(TARGETS); f["thresholds_keV"] = thresholds
            f["signal_efficiency"] = [accepted(signal.E2, t).fraction for t in thresholds]
            for (role, runs) in (("signal", [signal]), ("background", backgrounds)), (j, r) in enumerate(runs)
                g = create_group(f, role == "signal" ? "signal" : "backgrounds/$(lpad(j, 2, '0'))")
                g["run_label"] = r.run_label; g["sample"] = r.sample; g["flow_step"] = r.flow_step; g["flow_events"] = r.flow_events
                g["E1_keV"] = r.E1; g["E2_keV"] = r.E2; g["ion_drift_cm"] = r.drift; g["event_number"] = r.event_number
                role == "signal" && continue
                rows = [accepted(r.E2, t) for t in thresholds]
                g["acceptance"] = [x.fraction for x in rows]; g["acceptance_error"] = [x.error for x in rows]; g["accepted_events"] = [x.events for x in rows]
                if !isempty(r.E2)
                    c = curve(signal.E2, r.E2)
                    g["curve/threshold_keV"] = c.threshold; g["curve/signal_efficiency"] = c.efficiency; g["curve/background_acceptance"] = c.acceptance; g["curve/area"] = c.area
                end
            end
            flush(f); GuardedFiles.check_size(tmp, budget)
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("COMPARISON $target signal=$(length(signal.E2)) backgrounds=" * join((string(b.sample, ":", length(b.E2)) for b in backgrounds), ","))
    target
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 3 || error("Usage: prepare_iitaca_comparison.jl <output folder> <signal run folder> <background run folder>...")
    write_comparison(abspath(ARGS[1]), abspath(ARGS[2]), abspath.(ARGS[3:end]))
end
