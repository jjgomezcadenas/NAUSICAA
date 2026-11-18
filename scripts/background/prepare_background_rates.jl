#!/usr/bin/env julia
#=
Radiogenic background of the copper shield in the ROI, per year and per tonne-year,
for iItaca (2 mm pads) and eITACA.

  julia --project=. scripts/background/prepare_background_rates.jl [EITACA_TL208_COUNT.toml] [OUTPUT.h5]

For each detector and isotope, in the order of the ITACA revisited budget:

 1. Copper: mass of the shield (Background.copper_mass_kg of the card), its
    activity from the specific activity of the card, decays per year.
 2. Window: the fraction of decays with 2420-2490 keV in the fiducial gas and
    one entering gamma, from the pilot tables results/background/<name>/pilot/<isotope>/window.h5
    (prepare_ics_window.jl); events per year in the window.
 3. After the window, per saved window event: the track selection in the ROI
    (eps_track) and the network at 80 % signal (eps_net), from the gamma
    productions of the same detector: iItaca from results/iitaca/II_L1000_D500_P2_gamma_window/gamma_window.h5;
    eITACA from results/eitaca/IE_L1500_D400_results/results.h5, with the 208Tl
    saved-window count of count_window_events.jl (its production was saved on a
    wider window) and the 214Bi + 208Tl pooled network acceptance.
 4. Events per year: decays x window fraction x eps_track x eps_net; per
    tonne-year with the xenon mass of the detector (derived conditions of the
    variant). Errors: the window fraction (weights), eps_track (binomial) and
    eps_net (binomial), in quadrature; the specific activities have none.

Writes OUTPUT.h5 (default results/background/rates.h5, schema
nausicaa.background_rates.v1). An existing output is preserved.
=#
include(joinpath(@__DIR__, "resolve_background.jl"))
using HDF5, TOML, SHA, Statistics

const RATES_SCHEMA = "nausicaa.background_rates.v1"
const ISOTOPES = ("bi214", "tl208")
const DETECTORS = (
    (key = "iitaca", label = "iItaca, 2 mm pads", card = "config/background/ics_iitaca.toml", pilot = "results/background/II_L1500_ICS20/pilot",
     selection = "II_L1000_D500_P2", efficiencies = "results/iitaca/II_L1000_D500_P2_gamma_window/gamma_window.h5"),
    (key = "eitaca", label = "eITACA", card = "config/background/ics_eitaca.toml", pilot = "results/background/IE_L1500_ICS20/pilot",
     selection = "IE_L1500_D400", efficiencies = "results/eitaca/IE_L1500_D400_results/results.h5"))
"ITACA revisited (srcx/sensi.tex), for comparison: ion / electron track, all components, 15 cm ICS, 2400-2500 keV."
const IR = Dict("bi214_single_track_per_year" => (1.42, 4.42), "tl208_single_track_per_year" => (4.04, 24.0),
                "total_after_topology_per_tonne_year" => 0.05, "total_after_roi_per_tonne_year" => (0.02, 0.7), "resolution_fwhm_percent" => 0.7)

"""Xenon mass of a detector from its derived conditions."""
function xenon_mass_kg(family, variant)
    if family == "iitaca"
        d = IItaca.resolve_detector(joinpath(ROOT, "config/iitaca/iitaca.toml"); variant, root = ROOT).detector
        rows = IItaca.derived_conditions(d)
    else
        d = EItaca.resolve_detector(joinpath(ROOT, "config/eitaca/eitaca.toml"); variant, root = ROOT).detector
        rows = EItaca.derived_conditions(d; root = ROOT)
    end
    only(r[3] for r in rows if r[1] == "geometry" && r[2] == "xenon_mass_kg")
end

"""Window table of a pilot: fractions per class, the one-gamma class first."""
function window_table(path)
    h5open(path, "r") do f
        read(f["schema"]) == "nausicaa.background_window.v1" || error("Unexpected schema in $path")
        (fraction = read(f["window_fraction"]), error = read(f["window_fraction_error"]), counts = read(f["counts"]), classes = read(f["classes"]),
         generated = read(f["generated_decays"]), saved = read(f["saved_events"]), window = (read(f["window/lower_keV"]), read(f["window/upper_keV"])),
         mass = read(f["copper/mass_kg"]), Bq_kg = read(f["copper/specific_activity_Bq_kg"]), sha = bytes2hex(open(sha256, path)))
    end
end

index80(eff) = (i = findfirst(e -> isapprox(e, 0.8), eff); i === nothing && error("No 80 % signal point"); i)

"""Efficiencies after the window of iItaca: saved window events, ROI track-selected events, network acceptance at 80 %."""
function iitaca_efficiencies(path)
    h5open(path, "r") do f
        i = index80(read(f["target_signal_efficiency"]))
        Dict(iso => (saved = read(f["$iso/saved_events"]), track = read(f["$iso/roi_selected_events"]), net = read(f["$iso/acceptance"])[i],
                     net_error = read(f["$iso/acceptance_error"])[i], net_events = read(f["$iso/accepted_events"])[i], net_sample = read(f["$iso/window_events"]),
                     note = "gamma_window.h5: saved on 2420-2490 keV; network acceptance of the track-selected events of that window")
             for iso in ISOTOPES)
    end
end

"""Efficiencies after the window of eITACA: the 214Bi production was saved on the window, the 208Tl one on a wider window (count given)."""
function eitaca_efficiencies(path, count_toml)
    c = TOML.parsefile(count_toml); c["window_keV"] == [2420.0, 2490.0] || error("The 208Tl count is not on 2420-2490 keV")
    h5open(path, "r") do f
        eff = haskey(f, "target_signal_efficiency") ? read(f["target_signal_efficiency"]) : [0.9, 0.8, 0.7, 0.6]; i = index80(eff)
        acc = sum(read(f["$iso/window/accepted"])[i] for iso in ISOTOPES); n = sum(read(f["$iso/window/events"]) for iso in ISOTOPES)
        net = acc / n; net_error = sqrt(net * (1 - net) / n)
        Dict("bi214" => (saved = read(f["bi214/flow/read"]), track = read(f["bi214/flow/track"]), net, net_error, net_events = acc, net_sample = n,
                         note = "results.h5: saved on 2420-2490 keV; network acceptance pooled over 214Bi and 208Tl window events"),
             "tl208" => (saved = c["in_window"], track = read(f["tl208/flow/track"]), net, net_error, net_events = acc, net_sample = n,
                         note = "results.h5 with the 2420-2490 keV count of count_window_events.jl ($(c["in_window"]) of $(c["read_events"]) read events); pooled network acceptance"))
    end
end

function prepare_background_rates(count_toml::AbstractString = joinpath(ROOT, "results/background/eitaca_tl208_window_count.toml"),
                                  target::AbstractString = joinpath(ROOT, "results/background/rates.h5"))
    started = time(); target = abspath(target)
    isfile(target) && error("Preserve existing results: $target"); mkpath(dirname(target))
    tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = RATES_SCHEMA; f["roi"] = "one FWHM around Qbb, 2452.0-2463.6 keV; network at 80 % signal"
            for (k, v) in IR; f["reference/itaca_revisited/$k"] = collect(v); end
            f["reference/itaca_revisited/note"] = "ITACA revisited, srcx/sensi.tex: events per year for 1 t after the single-track cut (ion, electron), 2400-2500 keV; after topology; after the ROI at 0.7 % FWHM (ion track; conventional pure-xenon HPXeEL)"
            for det in DETECTORS
                r = resolved_shield(joinpath(ROOT, det.card)); g = r.geometry; p = "detectors/$(det.key)"
                mass_xe = xenon_mass_kg(r.family, joinpath(ROOT, r.card.card["variant"]))
                eff = det.key == "iitaca" ? iitaca_efficiencies(joinpath(ROOT, det.efficiencies)) : eitaca_efficiencies(joinpath(ROOT, det.efficiencies), count_toml)
                f["$p/label"] = det.label; f["$p/shield"] = r.card.card["name"]; f["$p/variant"] = r.variant; f["$p/selection_run"] = det.selection
                f["$p/card_sha256"] = r.card.sha256; f["$p/efficiency_file"] = det.efficiencies; f["$p/efficiency_sha256"] = bytes2hex(open(sha256, joinpath(ROOT, det.efficiencies)))
                f["$p/xenon_mass_kg"] = mass_xe; f["$p/copper_mass_kg"] = copper_mass_kg(g); f["$p/copper_thickness_cm"] = g.thickness_cm
                f["$p/copper_inner_radius_cm"] = g.inner_cm; f["$p/copper_cavity_length_cm"] = cavity_length_cm(g)
                total = 0.0; total_var = 0.0
                for iso in ISOTOPES
                    w = window_table(joinpath(ROOT, det.pilot, iso, "window.h5")); e = eff[iso]; q = "$p/$iso"
                    w.mass ≈ copper_mass_kg(g) || error("Copper mass of the pilot differs from the card for $(det.key) $iso")
                    Bq_kg = w.Bq_kg; activity = activity_Bq(g, Bq_kg); decays = decays_per_year(g, Bq_kg)
                    fw = w.fraction[1]; fw_err = w.error[1]                                       # one entering gamma
                    eps_track = e.track / e.saved; eps_track_err = sqrt(eps_track * (1 - eps_track) / e.saved)
                    eps = eps_track * e.net; rel = sqrt((fw_err / fw)^2 + (eps_track_err / eps_track)^2 + (e.net_error / e.net)^2)
                    rate = decays * fw * eps; rate_err = rate * rel
                    total += rate; total_var += rate_err^2
                    f["$q/specific_activity_Bq_kg"] = Bq_kg; f["$q/activity_Bq"] = activity; f["$q/decays_per_year"] = decays
                    f["$q/window/classes"] = w.classes; f["$q/window/fraction"] = w.fraction; f["$q/window/fraction_error"] = w.error; f["$q/window/counts"] = w.counts
                    f["$q/window/generated_decays"] = w.generated; f["$q/window/saved_events"] = w.saved; f["$q/window/keV"] = collect(w.window); f["$q/window/sha256"] = w.sha
                    f["$q/window/events_per_year"] = decays * fw; f["$q/window/events_per_year_error"] = decays * fw_err
                    f["$q/selection/saved_window_events"] = e.saved; f["$q/selection/track_selected_roi"] = e.track
                    f["$q/selection/eps_track"] = eps_track; f["$q/selection/eps_track_error"] = eps_track_err
                    f["$q/selection/eps_net"] = e.net; f["$q/selection/eps_net_error"] = e.net_error; f["$q/selection/net_events"] = e.net_events; f["$q/selection/net_sample"] = e.net_sample
                    f["$q/selection/eps_after_window"] = eps; f["$q/selection/eps_after_window_error"] = eps * sqrt((eps_track_err / eps_track)^2 + (e.net_error / e.net)^2)
                    f["$q/selection/note"] = e.note
                    f["$q/events_per_year"] = rate; f["$q/events_per_year_error"] = rate_err
                    f["$q/events_per_tonne_year"] = rate / (mass_xe / 1000); f["$q/events_per_tonne_year_error"] = rate_err / (mass_xe / 1000)
                    println("$(det.label) $iso: $(round(decays; sigdigits = 3)) decays/yr, window $(round(fw; sigdigits = 3)), eps_track $(round(eps_track; sigdigits = 3)), eps_net $(round(e.net; sigdigits = 3)): $(round(rate; sigdigits = 3)) ± $(round(rate_err; sigdigits = 2)) /yr, $(round(rate / (mass_xe / 1000); sigdigits = 3)) /(t yr)")
                end
                f["$p/total_events_per_year"] = total; f["$p/total_events_per_year_error"] = sqrt(total_var)
                f["$p/total_events_per_tonne_year"] = total / (mass_xe / 1000); f["$p/total_events_per_tonne_year_error"] = sqrt(total_var) / (mass_xe / 1000)
            end
            f["tl208_count_toml"] = relpath(abspath(count_toml), ROOT); f["tl208_count_sha256"] = bytes2hex(open(sha256, count_toml))
            f["code_sha256"] = bytes2hex(open(sha256, @__FILE__)); f["time_s"] = time() - started
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("RATES $target"); target
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 2 || error("Usage: prepare_background_rates.jl [EITACA_TL208_COUNT.toml] [OUTPUT.h5]")
    prepare_background_rates(ARGS...)
end
