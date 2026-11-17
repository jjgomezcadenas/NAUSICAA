#!/usr/bin/env julia
#=
iItaca II_L1000_D500_P2: network acceptance of the gamma backgrounds measured on a
box-energy window and applied to the ROI.

  julia --project=. scripts/iitaca/prepare_iitaca_gamma_window.jl

Within 2420-2490 keV of box energy the event topology does not depend on energy,
so the network acceptance measured there applies to the events in the ROI.

 1. Read the scores of the window samples (network_scorer.py --other-selection):
    events passing the track selection with energy_true_keV in [2420, 2490] keV.
    208Tl pools II_L1000_D500_P2_tl208 (40 M gammas) and II_L1000_D500_P2_tl208_w
    (100 M); 214Bi is II_L1000_D500_P2_bi214 (40 M).
 2. Acceptance at each signal efficiency of the training: events above the
    threshold over the events of the window, with its binomial error.
 3. Read the cut flow of the ROI runs (statistics/events.h5): 208Tl from
    II_L1000_D500_P2_tl208_w only, whose saved window (2420-2490 keV) covers every
    event that can reach the ROI; 214Bi from II_L1000_D500_P2_bi214. ROI events
    passing the track selection per generated gamma, with its binomial error.
 4. Background per generated gamma after selection and network: the product of
    2 and 3, errors added in quadrature. Rejection factors relative to the events
    saved by Geant4 and to the generated gammas.

Writes results/iitaca/II_L1000_D500_P2_gamma_window/gamma_window.h5; an existing
file is preserved.
=#
using HDF5, SHA

const ROOT = normpath(joinpath(@__DIR__, "../.."))
const TRAINING = "results/network/II_L1000_D500_P2/trainings/itrack_sparse_e100"
const WINDOW_KEV = (2420.0, 2490.0)
const TARGETS = (0.9, 0.8, 0.7, 0.6)
const ISOTOPES = (("tl208", ["II_L1000_D500_P2_tl208", "II_L1000_D500_P2_tl208_w"], "II_L1000_D500_P2_tl208_w"),
                  ("bi214", ["II_L1000_D500_P2_bi214"], "II_L1000_D500_P2_bi214"))
const OUTPUT = joinpath(ROOT, "results/iitaca/II_L1000_D500_P2_gamma_window/gamma_window.h5")

sha(path) = bytes2hex(open(sha256, path))

"""Scores of one window sample: events, events above each threshold, threshold, signal efficiency."""
function window_scores(species, label)
    path = joinpath(ROOT, TRAINING, "scores_$(species)_window_$label.h5")
    h5open(path, "r") do f
        occursin("window cut on the box energy (energy_true_keV), $(WINDOW_KEV[1]) to $(WINDOW_KEV[2]) keV", read(f["selection"])) || error("Unexpected selection in $path")
        target = read(f["target_signal_efficiency"]); rows = [findfirst(t -> isapprox(t, x; atol = 1e-6), target) for x in TARGETS]
        any(isnothing, rows) && error("Missing signal efficiencies in $path")
        (n = length(read(f["events/score"])), accepted = Int.(read(f["accepted_events"])[rows]), threshold = read(f["threshold"])[rows],
         efficiency = read(f["signal_efficiency"])[rows], path = path, sha = sha(path), model = read(f["model_sha256"]))
    end
end

"""Cut flow of a run: step names and events."""
cut_flow(label) = (p = joinpath(ROOT, "results/iitaca", label, "statistics/events.h5");
    h5open(f -> (step = read(f["cut_flow/step"]), events = Int.(read(f["cut_flow/events"])), roi = (read(f["roi_lower_keV"]), read(f["roi_upper_keV"])), path = p, sha = sha(p)), p, "r"))

binomial(k, n) = (p = k / n; (p, sqrt(p * (1 - p) / n)))

function prepare()
    isfile(OUTPUT) && error("Preserve the existing result: $OUTPUT"); mkpath(dirname(OUTPUT)); tmp = OUTPUT * ".tmp"
    h5open(tmp, "w") do f
        f["schema"] = "nausicaa.iitaca_gamma_window.v1"; f["variant"] = "II_L1000_D500_P2"; f["training"] = TRAINING
        f["window_lower_keV"] = WINDOW_KEV[1]; f["window_upper_keV"] = WINDOW_KEV[2]; f["target_signal_efficiency"] = collect(TARGETS)
        for (iso, window_runs, roi_run) in ISOTOPES
            g = create_group(f, iso); s = [window_scores(iso, l) for l in window_runs]
            all(x -> x.model == s[1].model && x.threshold == s[1].threshold, s) || error("The window samples of $iso were scored with different networks")
            n = sum(x.n for x in s); k = sum(x.accepted for x in s); A = binomial.(k, n)
            c = cut_flow(roi_run); steps = c.step; ev(name) = c.events[findfirst(==(name), steps)]
            N = ev("generated"); saved = ev("read"); track = ev("track"); r, σr = binomial(track, N)
            b = [r * a for (a, _) in A]; σb = [a == 0 ? r * σa : r * a * sqrt((σr / r)^2 + (σa / a)^2) for (a, σa) in A]
            g["window_runs"] = window_runs; g["window_events_per_run"] = [x.n for x in s]; g["window_events"] = n
            g["accepted_events"] = k; g["acceptance"] = first.(A); g["acceptance_error"] = last.(A)
            g["threshold"] = s[1].threshold; g["signal_efficiency"] = s[1].efficiency; g["model_sha256"] = s[1].model
            g["score_files"] = [relpath(x.path, ROOT) for x in s]; g["score_sha256"] = [x.sha for x in s]
            g["roi_run"] = roi_run; g["roi_lower_keV"] = c.roi[1]; g["roi_upper_keV"] = c.roi[2]
            g["cut_flow/step"] = steps; g["cut_flow/events"] = c.events; g["statistics_file"] = relpath(c.path, ROOT); g["statistics_sha256"] = c.sha
            g["generated_gammas"] = N; g["saved_events"] = saved; g["roi_selected_events"] = track
            g["roi_selected_per_gamma"] = r; g["roi_selected_per_gamma_error"] = σr
            g["background_per_gamma"] = b; g["background_per_gamma_error"] = σb
            g["rejection_saved"] = [saved / (track * a) for (a, _) in A]; g["rejection_generated"] = 1 ./ b
            println("$iso: window $n events; ", join(["$(round(Int, 100 * e)) %: $(k[i]) ($(round(100A[i][1]; sigdigits = 3)) ± $(round(100A[i][2]; sigdigits = 2)) %), b = $(b[i]) ± $(σb[i]) per gamma" for (i, e) in enumerate(TARGETS)], "; "))
        end
    end
    Base.rename(tmp, OUTPUT); println(OUTPUT)
end

abspath(PROGRAM_FILE) == @__FILE__() && prepare()
