#!/usr/bin/env julia
#=
iItaca: classes of the Monte Carlo truth of the selected events and their network acceptance.

  julia --project=. scripts/iitaca/summarize_iitaca_truth.jl <folder> [comparison.h5]

Reads <folder>/truth_bi214.h5 and <folder>/truth_single_electrons.h5
(prepare_iitaca_truth.jl) and the E₂ thresholds of the comparison of the
variant (default results/iitaca/II_L1000_D500_comparison/comparison.h5), and writes <folder>/summary.h5 (schema iitaca.truth_summary.v1).
Python makes the page (scripts/iitaca/present_iitaca_truth.py).

 1. Secondary deposit: a deposit of the truth table with at least 10 keV (the
    deposit threshold of the track selection) inside the simulated box. The
    deposits of the atomic relaxation at the absorption point (relaxation
    electrons and the local energy of the gamma) are counted apart: they lie
    at the start of the photoelectron.
 2. Class of every event, the first that applies:
      second_electron   a Compton electron, a photoelectron that is not the
                        leading electron, or a positron, as secondary deposit;
      fluorescence      a fluorescence photon of the gamma as secondary deposit;
      photon            a photon of the main track or another photon;
      electron_only     no secondary deposit.
 3. Distance class of the fluorescence photon of the gamma: no photon, photon
    with less than 10 keV in the box, its smallest distance to the main track
    below 4 mm, 4 to 8 mm, beyond 8 mm. Its acceptance table counts only the
    events whose gamma was absorbed without Compton scattering.
 4. Second electron: the most energetic second-electron deposit of the event,
    its smallest distance to the main track and its energy in the box.
 5. Acceptance of the network (score above the 80 % signal threshold) and of
    the E₂ cut at 80 % signal, per class and per bin of the energy within
    7 mm of the start of the leading electron (every hit in the box), with
    binomial errors.
 6. Expected ²¹⁴Bi acceptance if each bin of that start energy were accepted
    as single electrons are: sum over bins of the ²¹⁴Bi fraction times the
    single-electron acceptance.
=#
using HDF5, Statistics
const ROOT = normpath(joinpath(@__DIR__, "../.."))

const SUMMARY_SCHEMA = "iitaca.truth_summary.v1"
const DEPOSIT_MIN_KEV = 10.0
const CLASSES = ["electron_only", "fluorescence", "second_electron", "photon"]
const XRAY_CLASSES = ["no_photon", "below_10keV_in_box", "lt4mm", "4to8mm", "gt8mm"]
const HISTORIES = ["generated_electron", "photoelectric", "compton_photoelectric", "compton_only", "pair", "other"]
const SECOND_DISTANCE_EDGES_MM = [0.0, 4, 8, 1e9]
const SECOND_ENERGY_EDGES_KEV = [10.0, 200, 400, 800, 2500]
const START_EDGES_KEV = [0.0, 50, 75, 100, 150, 200, 300, 400, 600, 2500]

readgroup(g) = Dict(k => read(g[k]) for k in keys(g))

"""Class of every event (index into CLASSES) from its secondary deposits."""
function event_classes(ev, dep)
    n = length(ev["event_number"]); cls = ones(Int8, n)
    rank = Dict(1 => 3, 2 => 3, 6 => 3, 4 => 2, 5 => 4, 7 => 4, 8 => 4)       # deposit type => class
    priority = Dict(1 => 0, 2 => 2, 3 => 3, 4 => 1)                           # class => order of precedence
    for i in 1:n, k in ev["deposit_offset"][i] .+ (1:ev["deposit_count"][i])
        dep["deposit_in_box_keV"][k] >= DEPOSIT_MIN_KEV || continue
        c = get(rank, dep["type"][k], 0); c == 0 && continue
        priority[c] > priority[cls[i]] && (cls[i] = c)
    end
    cls
end

"""Smallest distance (mm) and box energy (keV) of the most energetic second electron of every event; NaN without one."""
function second_electrons(ev, dep)
    n = length(ev["event_number"]); d = fill(NaN, n); e = fill(NaN, n)
    for i in 1:n, k in ev["deposit_offset"][i] .+ (1:ev["deposit_count"][i])
        dep["type"][k] in (1, 2, 6) && dep["deposit_in_box_keV"][k] >= DEPOSIT_MIN_KEV || continue
        (isnan(e[i]) || dep["deposit_in_box_keV"][k] > e[i]) && (e[i] = dep["deposit_in_box_keV"][k]; d[i] = dep["min_distance_mm"][k])
    end
    d, e
end

binof(x, edges) = isnan(x) ? 0 : clamp(searchsortedlast(edges, x), 1, length(edges) - 1)

function xray_classes(ev)
    map(eachindex(ev["event_number"])) do i
        isnan(ev["xray_keV"][i]) ? Int8(1) : ev["xray_deposit_keV"][i] < DEPOSIT_MIN_KEV ? Int8(2) :
        ev["xray_distance_mm"][i] < 4 ? Int8(3) : ev["xray_distance_mm"][i] < 8 ? Int8(4) : Int8(5)
    end
end

"""Events, network passes, E₂ passes and fractions with binomial errors, per value of `label` in 1:m."""
function acceptance(label, m, passed, e2)
    n = [count(==(c), label) for c in 1:m]; k = [count(i -> label[i] == c && passed[i], eachindex(label)) for c in 1:m]
    k2 = [count(i -> label[i] == c && e2[i], eachindex(label)) for c in 1:m]
    f(k, n) = n > 0 ? k / n : NaN; err(k, n) = n > 0 ? sqrt(max(k, 1) * (1 - k / n)) / n : NaN
    (events = n, passed = k, fraction = f.(k, n), error = err.(k, n), e2_passed = k2, e2_fraction = f.(k2, n), e2_error = err.(k2, n))
end

function summarize(folder::AbstractString, comparison::AbstractString = joinpath(ROOT, "results/iitaca/II_L1000_D500_comparison/comparison.h5"))
    e2_threshold = h5open(comparison, "r") do f
        read(f["thresholds_keV"])[findfirst(t -> isapprox(t, 0.8), read(f["targets"]))]
    end
    out = joinpath(folder, "summary.h5"); tmp = out * ".tmp"
    samples = Dict{String,Any}()
    for sample in ("bi214", "single_electrons")
        ev, dep, thr = h5open(joinpath(folder, "truth_$sample.h5"), "r") do f
            readgroup(f["events"]), readgroup(f["deposits"]), read(f["threshold_80"])
        end
        passed = ev["passed"] .== 1; e2 = ev["E2_keV"] .> e2_threshold
        cls = event_classes(ev, dep); xc = xray_classes(ev); sd, se = second_electrons(ev, dep)
        sbin = [clamp(searchsortedlast(START_EDGES_KEV, x), 1, length(START_EDGES_KEV) - 1) for x in ev["start_all_keV"]]
        samples[sample] = (ev = ev, passed = passed, e2 = e2, threshold = thr, classes = cls, xray = xc, start_bin = sbin,
            by_class = acceptance(cls, length(CLASSES), passed, e2), by_xray = acceptance(ifelse.(ev["history"] .== 1, xc, Int8(0)), length(XRAY_CLASSES), passed, e2),
            by_history = acceptance(ev["history"] .+ 1, length(HISTORIES), passed, e2),
            second_distance_mm = sd, second_keV = se,
            by_second_distance = acceptance(binof.(sd, Ref(SECOND_DISTANCE_EDGES_MM)), length(SECOND_DISTANCE_EDGES_MM) - 1, passed, e2),
            by_second_energy = acceptance(binof.(se, Ref(SECOND_ENERGY_EDGES_KEV)), length(SECOND_ENERGY_EDGES_KEV) - 1, passed, e2),
            by_start = acceptance(sbin, length(START_EDGES_KEV) - 1, passed, e2), all = acceptance(ones(Int, length(passed)), 1, passed, e2))
    end
    b = samples["bi214"]; e = samples["single_electrons"]
    w = b.by_start.events ./ sum(b.by_start.events)
    fe = replace(e.by_start.fraction, NaN => 0.0); fe2 = replace(e.by_start.e2_fraction, NaN => 0.0)
    expected = sum(w .* fe); expected_e2 = sum(w .* fe2)
    expected_error = sqrt(sum((w .* replace(e.by_start.error, NaN => 0.0)) .^ 2))
    try
        h5open(tmp, "w") do f
            f["schema"] = SUMMARY_SCHEMA; f["classes"] = CLASSES; f["xray_classes"] = XRAY_CLASSES; f["histories"] = HISTORIES
            f["start_edges_keV"] = START_EDGES_KEV; f["second_distance_edges_mm"] = SECOND_DISTANCE_EDGES_MM; f["second_energy_edges_keV"] = SECOND_ENERGY_EDGES_KEV; f["deposit_min_keV"] = DEPOSIT_MIN_KEV; f["e2_threshold_keV"] = e2_threshold
            f["expected_bi214_fraction"] = expected; f["expected_bi214_fraction_error"] = expected_error; f["expected_bi214_e2_fraction"] = expected_e2
            for (sample, s) in samples
                g = create_group(f, sample)
                g["event_number"] = s.ev["event_number"]; g["class"] = s.classes; g["xray_class"] = s.xray; g["start_bin"] = s.start_bin; g["second_distance_mm"] = s.second_distance_mm; g["second_keV"] = s.second_keV
                g["passed"] = Int8.(s.passed); g["e2_passed"] = Int8.(s.e2); g["threshold_80"] = s.threshold
                for name in (:by_class, :by_xray, :by_history, :by_second_distance, :by_second_energy, :by_start, :all), (k, v) in pairs(getfield(s, name))
                    g["$name/$k"] = v
                end
            end
        end
        Base.rename(tmp, out)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("OUTPUT $out expected_bi214=$(round(100expected; digits = 2)) % (measured $(round(100b.all.fraction[1]; digits = 2)) %)")
    out
end

if abspath(PROGRAM_FILE) == @__FILE__
    1 <= length(ARGS) <= 2 || error("Usage: summarize_iitaca_truth.jl <folder> [comparison.h5]")
    summarize(abspath.(ARGS)...)
end
