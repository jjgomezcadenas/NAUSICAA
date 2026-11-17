#!/usr/bin/env julia
#=
eITACA IE_L1500_D400: numbers of the results page (notes/html/eitaca/results.html).

  julia --project=. scripts/eitaca/prepare_eitaca_results.jl

 1. Cut flow of each sample (statistics/events.h5 of prepare_eitaca_features.jl):
    generated, saved, fast trigger, track selection, track selection in the ROI.
 2. Network of the etrack (results/network/IE_L1500_D400/trainings/etrack_sparse_e100):
    the test events of the training at each signal efficiency (bb0nu, single
    electrons), and the scores of the 214Bi and 208Tl events passing the track
    selection with box energy in 2420-2470 keV (the gamma files stop at
    2470 keV). Within this window the topology does not depend on the energy;
    the two isotopes are pooled into one acceptance, with its binomial error.
 3. Noise of one hit, from the saved columns of each bb0nu event: one SiPM over
    rebin_factor buckets, sigma_pe^2 = rebin x (noise_rms_pe^2 + mean dark
    pulses per SiPM and bucket), the dark pulses Poisson; in keV with the
    scale of the event energy, W / (PDE x EL photons per electron x
    acceptance), and the attachment correction exp(t / tau).
 4. Hits per mm of the selected bb0nu tracks in the ROI (hits over length_mm),
    the basis of segment_min_hits.

Writes results/eitaca/IE_L1500_D400_results/results.h5; an existing file is preserved.
=#
using HDF5, TOML, SHA, Statistics

const ROOT = normpath(joinpath(@__DIR__, "../.."))
const RUN = "IE_L1500_D400"
const TRAINING = "results/network/$RUN/trainings/etrack_sparse_e100"
const WINDOW_KEV = (2420.0, 2470.0)
const TARGETS = (0.9, 0.8, 0.7, 0.6)
const SAMPLES = ("bb0nu", "single_electrons", "bi214", "tl208")
const OUTPUT = joinpath(ROOT, "results/eitaca/$(RUN)_results/results.h5")

sha(path) = bytes2hex(open(sha256, path))
binomial(k, n) = (p = k / n; (p, sqrt(p * (1 - p) / n)))
rows_of(target) = (r = [findfirst(t -> isapprox(t, x; atol = 1e-6), target) for x in TARGETS]; any(isnothing, r) && error("Missing signal efficiencies"); r)

function prepare()
    isfile(OUTPUT) && error("Preserve the existing result: $OUTPUT"); mkpath(dirname(OUTPUT)); tmp = OUTPUT * ".tmp"
    inputs = String[]
    try
    h5open(tmp, "w") do out
        out["schema"] = "nausicaa.eitaca_results.v1"; out["variant"] = RUN; out["training"] = TRAINING
        out["window_lower_keV"] = WINDOW_KEV[1]; out["window_upper_keV"] = WINDOW_KEV[2]; out["target_signal_efficiency"] = collect(TARGETS)
        resolved = nothing
        for s in SAMPLES
            p = joinpath(ROOT, "results/eitaca/$(RUN)_$s/statistics/events.h5"); push!(inputs, p)
            h5open(p, "r") do f
                flow = Dict(zip(read(f["cut_flow/step"]), read(f["cut_flow/events"])))
                sel = read(f["events/selected"]) .== 1; roi = read(f["events/in_roi"]) .== 1
                g = create_group(out, s)
                for k in ("generated", "read", "trigger", "roi", "track"); g["flow/$k"] = Int(flow[k]); end
                g["selected"] = count(sel); g["selected_in_roi"] = count(sel .& roi)
                g["window_selected"] = count(sel .& (WINDOW_KEV[1] .<= read(f["events/energy_true_keV"]) .<= WINDOW_KEV[2]))
                s == "bb0nu" || return
                resolved = TOML.parse(read(f["resolved"]))["detector"]; settings = TOML.parse(read(f["settings"]))
                hits = read(f["events/hits"]); len = read(f["events/length_mm"]); k = sel .& roi
                g["hits_per_mm_median"] = median(hits[k] ./ len[k])
                # 3. noise of one hit
                asic = TOML.parsefile(joinpath(ROOT, resolved["sipm"]["asic_card"])); sipm = TOML.parsefile(joinpath(ROOT, resolved["sipm"]["response_card"]))
                ep = h5open(joinpath(ROOT, "results/eitaca/$(RUN)_bb0nu/batches/b00000/electrons.h5"), "r") do e
                    (pde = read(e["configuration/pde"]), photons = read(e["configuration/photons_per_electron"]), acceptance = read(e["configuration/light_acceptance"]))
                end
                rebin = resolved["reconstruction"]["rebin_factor"]; noise = asic["noise_rms_pe"]
                dark = read(f["events/dark_pulses"]) ./ (read(f["events/channels"]) .* read(f["events/buckets"]))
                t_ms = 10 .* read(f["events/electron_drift_cm"]) ./ resolved["drift"]["electron_drift_velocity_mm_us"] ./ 1000
                keV_per_pe = resolved["energy"]["W_eV"] / 1000 / (ep.pde * ep.photons * ep.acceptance)
                sigma = keV_per_pe .* exp.(t_ms ./ resolved["drift"]["electron_lifetime_ms"]) .* sqrt.(rebin .* (noise^2 .+ dark))
                g["noise/sigma_hit_keV"] = sigma; g["noise/median_keV"] = median(sigma); g["noise/max_keV"] = maximum(sigma)
                g["noise/dark_per_bucket_median"] = median(dark); g["noise/asic_rms_pe"] = noise; g["noise/rebin_factor"] = rebin
                g["noise/keV_per_pe"] = keV_per_pe; g["noise/hit_cut_keV"] = settings["tracks"]["hit_cut_keV"]
                g["noise/hit_cut_sigma"] = settings["tracks"]["hit_cut_keV"] / maximum(sigma)
            end
        end
        # 2. network
        tr = joinpath(ROOT, TRAINING, "training.h5"); push!(inputs, tr)
        h5open(tr, "r") do f
            sep = f["test/separation"]; r = rows_of(read(sep["target_signal_efficiency"]))
            g = create_group(out, "network")
            g["area"] = read(sep["area"]); g["threshold"] = read(sep["threshold"])[r]
            g["signal_events"] = read(sep["signal_events"]); g["signal_efficiency"] = read(sep["signal_efficiency"])[r]
            g["background_events"] = read(sep["background_events"]); g["background_acceptance"] = read(sep["background_acceptance"])[r]
            g["background_accepted"] = round.(Int, read(sep["background_acceptance"])[r] .* read(sep["background_events"]))
        end
        pooled_n = 0; pooled_k = zeros(Int, length(TARGETS)); thresholds = nothing
        for iso in ("bi214", "tl208")
            p = joinpath(ROOT, TRAINING, "scores_$(iso)_window_$(RUN)_$(iso).h5"); push!(inputs, p)
            h5open(p, "r") do f
                occursin("window cut on the box energy (energy_true_keV), $(WINDOW_KEV[1]) to $(WINDOW_KEV[2]) keV", read(f["selection"])) || error("Unexpected selection in $p")
                r = rows_of(read(f["target_signal_efficiency"])); n = length(read(f["events/score"])); k = Int.(read(f["accepted_events"])[r])
                thr = read(f["threshold"])[r]; thresholds === nothing ? (thresholds = thr) : (thr == thresholds || error("Different thresholds in $p"))
                n == read(out["$iso/window_selected"]) || error("Window events of $iso differ from the statistics")
                g = out[iso]; g["window/events"] = n; g["window/accepted"] = k
                A = binomial.(k, n); g["window/acceptance"] = first.(A); g["window/acceptance_error"] = last.(A)
                pooled_n += n; pooled_k .+= k
            end
        end
        A = binomial.(pooled_k, pooled_n)
        out["pooled/events"] = pooled_n; out["pooled/accepted"] = pooled_k; out["pooled/acceptance"] = first.(A); out["pooled/acceptance_error"] = last.(A)
        out["inputs"] = [relpath(p, ROOT) for p in inputs]; out["input_sha256"] = sha.(inputs)
        out["code_sha256"] = sha(@__FILE__)
    end
    Base.rename(tmp, OUTPUT)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println(OUTPUT)
end

abspath(PROGRAM_FILE) == @__FILE__() && prepare()
