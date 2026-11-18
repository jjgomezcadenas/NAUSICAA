#!/usr/bin/env julia
#=
iItaca: bound on the gain of a diffusion-limited ion readout, from the Monte
Carlo truth of the events passing the 2 mm selection and network.

  julia --project=. scripts/iitaca/prepare_diffusion_limit.jl <truth folder> [OUTPUT.h5]

Reads truth_bb0nu.h5 (test signal events of the training), truth_single_electrons.h5
(its test background events) and truth_bi214.h5 (selected events) written by
prepare_iitaca_truth.jl with the secondary energy profile (schema iitaca.truth.v2),
and keeps the events that pass the network at 80 % signal. For each of them it
asks whether a readout blurred only by ion diffusion would show a detached
deposit: secondary energy at a distance of at least k sigma from the main track
of at least E_thr (the profile at the first grid point not below k sigma, so the
energy is never overestimated). sigma is the thermal ion spread
sigma_t_ion_mm(T, L, E), with the drift L of the event:

 a. the production: ion_drift_cm of the event, uniform over the 1 m of
    II_L1000_D500 at 500 V/cm (sigma up to 1.0 mm);
 b. II_L1500_D350: a drift uniform over 1.5 m at 350 V/cm, drawn per event
    (sigma up to 1.5 mm);
 c. fixed sigma, a scan.

Writes OUTPUT.h5 (default <truth folder>/../II_L1000_D500_P2_diffusion_limit/diffusion_limit.h5,
schema iitaca.diffusion_limit.v1): for each sample the number of passing events
and the flagged fraction with its binomial error on the grids of k and E_thr
for a and b, and on the grid of sigma, k and E_thr for c; for 214Bi also by the
distance of its largest secondary deposit (< 4 mm, 4-8 mm, > 8 mm, none). The
flagged fraction of bb0nu is the signal loss, that of the backgrounds the
further rejection: 1 / (1 - fraction).
=#
using HDF5, SHA, Random, Statistics
const ROOT = normpath(joinpath(@__DIR__, "../.."))
if !isdefined(@__MODULE__, :IItaca)
    include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
    include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
    include(joinpath(ROOT, "src/iitaca/IItaca.jl")); using .IItaca
end
using .EventAnalysis: sigma_t_ion_mm

const LIMIT_SCHEMA = "iitaca.diffusion_limit.v1"
const SAMPLES = ("bb0nu", "single_electrons", "bi214")
const K_GRID = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
const E_GRID_KEV = [5.0, 10.0, 20.0]
const SIGMA_GRID_MM = collect(0.25:0.25:2.0)
const CLASSES = ("lt4mm", "4to8mm", "gt8mm", "none")
const CLASS_MIN_KEV = 3.0          # a largest deposit below this is "none"

"""Drift conditions of a variant: temperature, field and fiducial length."""
function drift_conditions(variant)
    r = resolve_detector(joinpath(ROOT, "config/iitaca/iitaca.toml"); variant = joinpath(ROOT, "config/iitaca/variants/$variant.toml"), root = ROOT)
    d = r.detector
    (T_K = Float32(d["gas"]["temperature_K"]), E_V_cm = Float32(d["drift"]["field_V_cm"]), L_cm = Float32(d["geometry"]["fiducial_length_cm"]), name = variant)
end

"""Passing events of a truth file: profiles (edges x events), drift, classes."""
function load_passing(path)
    h5open(path, "r") do f
        read(f["schema"]) == "iitaca.truth.v2" || error("$path needs the truth schema v2 (with the secondary energy profile)")
        passed = read(f["events/passed"]) .== 1
        prof = read(f["events/secondary_profile_keV"])[:, passed]
        dist = read(f["events/largest_deposit_distance_mm"])[passed]; dep = read(f["events/largest_deposit_keV"])[passed]
        class = [dep[i] < CLASS_MIN_KEV || !isfinite(dist[i]) ? 4 : dist[i] < 4 ? 1 : dist[i] < 8 ? 2 : 3 for i in eachindex(dist)]
        (profile = prof, edges = read(f["profile_edges_mm"]), drift_cm = Float32.(read(f["events/ion_drift_cm"])[passed]),
         event_number = read(f["events/event_number"])[passed], class = class, selected = length(passed), passed = count(passed),
         sha = bytes2hex(open(sha256, path)), selection = read(f["selection"]))
    end
end

"""Secondary energy at a distance of at least d: the profile at the first grid point not below d."""
function energy_beyond(profile::AbstractVector, edges, d)
    j = searchsortedfirst(edges, d); j > length(edges) ? 0.0 : profile[j]
end

"""Flagged fraction with its binomial error of events whose energy beyond k sigma is at least E, per k and E (k x E)."""
function flagged(s, sigma_mm::AbstractVector)
    n = size(s.profile, 2); frac = zeros(length(K_GRID), length(E_GRID_KEV)); err = similar(frac); flags = falses(length(K_GRID), length(E_GRID_KEV), n)
    for (a, k) in enumerate(K_GRID), (b, E) in enumerate(E_GRID_KEV)
        c = 0
        for i in 1:n
            f = energy_beyond(view(s.profile, :, i), s.edges, k * sigma_mm[i]) >= E; flags[a, b, i] = f; c += f
        end
        p = n == 0 ? 0.0 : c / n; frac[a, b] = p; err[a, b] = n == 0 ? 0.0 : sqrt(p * (1 - p) / n)
    end
    frac, err, flags
end

function prepare_diffusion_limit(folder::AbstractString, target::AbstractString = "")
    folder = abspath(folder); started = time()
    prod = drift_conditions("II_L1000_D500"); long = drift_conditions("II_L1500_D350")
    isempty(target) && (target = joinpath(dirname(folder), "II_L1000_D500_P2_diffusion_limit", "diffusion_limit.h5")); target = abspath(target)
    isfile(target) && error("Preserve existing results: $target"); mkpath(dirname(target))
    tmp = target * ".tmp"
    try
        h5open(tmp, "w") do f
            f["schema"] = LIMIT_SCHEMA; f["truth_folder"] = relpath(folder, ROOT); f["k_grid"] = K_GRID; f["E_grid_keV"] = E_GRID_KEV; f["sigma_grid_mm"] = SIGMA_GRID_MM
            f["classes"] = collect(CLASSES); f["class_min_keV"] = CLASS_MIN_KEV
            for (name, c) in (("production", prod), ("long", long))
                f["drift/$name/variant"] = c.name; f["drift/$name/T_K"] = c.T_K; f["drift/$name/E_V_cm"] = c.E_V_cm; f["drift/$name/L_cm"] = c.L_cm
                f["drift/$name/sigma_max_mm"] = sigma_t_ion_mm(c.T_K, c.L_cm, c.E_V_cm)
            end
            for sample in SAMPLES
                path = joinpath(folder, "truth_$sample.h5"); isfile(path) || error("Missing $path")
                s = load_passing(path); n = s.passed
                f["$sample/truth_sha256"] = s.sha; f["$sample/selection"] = s.selection; f["$sample/selected"] = s.selected; f["$sample/passed"] = n
                f["$sample/event_number"] = s.event_number; f["$sample/class"] = Int8.(s.class)
                # a: the production drift; b: a drift uniform over the long detector, drawn per event
                sig_a = [sigma_t_ion_mm(prod.T_K, L, prod.E_V_cm) for L in s.drift_cm]
                rng = MersenneTwister(20261118 + sum(Int64, s.event_number) % 1_000_000)
                drift_b = Float32.(long.L_cm .* rand(rng, n)); sig_b = [sigma_t_ion_mm(long.T_K, L, long.E_V_cm) for L in drift_b]
                for (name, sig) in (("production", sig_a), ("long", sig_b))
                    frac, err, flags = flagged(s, sig)
                    f["$sample/$name/sigma_mm"] = sig; f["$sample/$name/flagged"] = frac; f["$sample/$name/flagged_error"] = err
                    f["$sample/$name/flags"] = Int8.(flags)
                    if sample == "bi214"
                        for (ci, cname) in enumerate(CLASSES)
                            m = s.class .== ci; nc = count(m)
                            cnt = [count(flags[a, b, m]) for a in eachindex(K_GRID), b in eachindex(E_GRID_KEV)]
                            f["$sample/$name/class/$cname/events"] = nc; f["$sample/$name/class/$cname/flagged_events"] = cnt
                            f["$sample/$name/class/$cname/flagged"] = nc == 0 ? zeros(size(cnt)) : cnt ./ nc
                        end
                    end
                end
                # c: fixed sigma
                frac_c = zeros(length(SIGMA_GRID_MM), length(K_GRID), length(E_GRID_KEV)); err_c = similar(frac_c)
                for (si, sg) in enumerate(SIGMA_GRID_MM)
                    frac, err, _ = flagged(s, fill(sg, n)); frac_c[si, :, :] = frac; err_c[si, :, :] = err
                end
                f["$sample/fixed/flagged"] = frac_c; f["$sample/fixed/flagged_error"] = err_c
                println("$sample: $n passing of $(s.selected); flagged at k = 3, 10 keV: production $(round(100 * f["$sample/production/flagged"][4, 2]; digits = 1)) %, long $(round(100 * f["$sample/long/flagged"][4, 2]; digits = 1)) %")
            end
            f["code_sha256"] = bytes2hex(open(sha256, @__FILE__)); f["time_s"] = time() - started
        end
        Base.rename(tmp, target)
    catch
        isfile(tmp) && rm(tmp); rethrow()
    end
    println("DIFFUSION LIMIT $target time_s=$(round(time() - started; digits = 1))")
    target
end

if abspath(PROGRAM_FILE) == @__FILE__
    1 <= length(ARGS) <= 2 || error("Usage: prepare_diffusion_limit.jl <truth folder> [OUTPUT.h5]")
    prepare_diffusion_limit(ARGS...)
end
