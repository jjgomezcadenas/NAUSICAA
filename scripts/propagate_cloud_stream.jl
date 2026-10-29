#!/usr/bin/env julia
# Prototype one independently writable sample/carrier diffusion stream.
# Usage: julia --project=. scripts/propagate_cloud_stream.jl SAMPLE CARRIER LIMIT OUTPUT
# SAMPLE is bb0nu or electrons; CARRIER is ion, pure_xe or xe_he.
using HDF5, Random, SHA, TOML
const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis

function event_seed(path::String, id::Int64, seed::Int, label::String)::UInt64
    digest = bytes2hex(sha256(join((realpath(path), string(id), string(seed), label), "|")))
    return parse(UInt64, digest[1:16]; base=16)
end

"""Diffuse one sample and carrier, writing independent per-event HDF5 clouds.

Stored coordinates are Float32 mm, occupied-cell counts Int64, and explicit
Float32 histogram edges. The common sampled event location sets complementary
ion/electron drift lengths. Ion thermal width and Petit's electron coefficients
are evaluated from the repository cards. The write is complete only after
the stream-level `complete` dataset is written.
"""
function propagate_stream(sample::String, carrier::String, limit::Int, output::String)
    sample in ("bb0nu", "electrons") || error("Unknown sample: $sample")
    carrier in ("ion", "pure_xe", "xe_he") || error("Unknown carrier: $carrier")
    limit > 0 || error("LIMIT must be positive")
    card = TOML.parsefile(joinpath(ROOT, "config/pixel_charge_1mm.toml"))
    ec = TOML.parsefile(joinpath(ROOT, "config/electron_transport.toml"))
    det = TOML.parsefile(joinpath(ROOT, "config/detector.toml"))
    tc = card["transport"]
    input = realpath(joinpath(ROOT, card["samples"][sample]))
    w = Float32(tc["ionisation_energy_eV"])
    pressure = Float32(ec["gas"]["pressure_bar"])
    maximum_cm = Float32(ec["placement"]["l_max_cm"])
    minimum_cm = Float32(ec["placement"]["min_ion_drift_cm"])
    maximum_cm == Float32(det["geometry"]["bulk_reference_length_m"] * 100) ||
        error("Electron placement length disagrees with detector card")
    minimum_cm < maximum_cm || error("Invalid drift range")
    mkpath(dirname(output))
    started = time_ns()
    h5open(input, "r") do source
        ids = read(source["events/event_id"])
        ids = ids[1:min(limit, length(ids))]
        h5open(output, "w") do target
            attributes(target)["schema"] = "nausicaa.three_cloud_prototype.v1"
            target["metadata/sample"] = sample
            target["metadata/carrier"] = carrier
            target["metadata/source_path"] = input
            target["metadata/source_sha256"] = open(io -> bytes2hex(sha256(io)), input)
            target["metadata/pressure_bar"] = pressure
            target["metadata/helium_fraction"] =
                carrier == "pure_xe" ? 0f0 : Float32(ec["gas"]["helium_fraction"])
            target["metadata/l_max_cm"] = maximum_cm
            target["metadata/min_ion_drift_cm"] = minimum_cm
            target["metadata/nbins"] = Int64(tc["nbins"])
            target["metadata/nsigma"] = Float32(tc["nsigma"])
            target["metadata/w_eV"] = w
            target["metadata/base_seed"] = Int64(card["run"]["seed"])
            group = create_group(target, "events")
            for (index, raw_id) in enumerate(ids)
                id = Int64(raw_id)
                hits, _ = read_event(source, id)
                counts = ion_counts(hits, w).sampled
                location_seed = event_seed(input, id, card["run"]["seed"], "position")
                ion_cm = sample_drift_cm(Xoshiro(location_seed), minimum_cm, maximum_cm)
                electron_cm = maximum_cm - ion_cm
                drift_cm = carrier == "ion" ? ion_cm : electron_cm
                if carrier == "ion"
                    sigma_t = sigma_l = sigma_t_ion_mm(Float32(det["gas"]["temperature_K"]),
                        drift_cm, Float32(det["bias"]["bulk_field_V_m"]) / 100f0)
                else
                    coeff = ec["diffusion"][carrier]
                    scale = sqrt(drift_cm / pressure)
                    sigma_t = Float32(coeff["d_T_mm_sqrtbar_per_sqrtcm"]) * scale
                    sigma_l = Float32(coeff["d_L_mm_sqrtbar_per_sqrtcm"]) * scale
                end
                seed = event_seed(input, id, card["run"]["seed"], carrier)
                cloud = diffuse_xyz_image_mc(Xoshiro(seed), hits, counts;
                    sigma_t_mm=sigma_t, sigma_l_mm=sigma_l,
                    nbins=Int(tc["nbins"]), nsigma=Float32(tc["nsigma"]),
                    ionisation_energy_eV=w)
                g = create_group(group, string(id))
                for field in ("x", "y", "z", "ions")
                    values = cloud.bins[!, Symbol(field)]
                    ds = create_dataset(g, field, eltype(values), (length(values),);
                        chunk=(max(1, min(length(values), 4096)),), shuffle=true, deflate=3)
                    write(ds, values)
                end
                for field in (:x_edges_mm, :y_edges_mm, :z_edges_mm,
                              :sampled_ions, :retained_ions, :outside_ions)
                    g[string(field)] = getproperty(cloud, field)
                end
                g["ion_drift_cm"] = ion_cm
                g["electron_drift_cm"] = electron_cm
                g["sigma_t_mm"] = sigma_t
                g["sigma_l_mm"] = sigma_l
                g["location_seed"] = string(location_seed)
                g["diffusion_seed"] = string(seed)
                index % 20 == 0 && (println("$sample/$carrier $index/$(length(ids))"); flush(stdout))
            end
            target["metadata/event_count"] = Int64(length(ids))
            target["metadata/complete"] = true
        end
    end
    println("completed $sample/$carrier in ", round((time_ns()-started)/1e9; digits=2), " s")
    println("peak_rss_bytes=", Sys.maxrss())
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 4 || error("Expected SAMPLE CARRIER LIMIT OUTPUT")
    propagate_stream(ARGS[1], ARGS[2], parse(Int, ARGS[3]), abspath(ARGS[4]))
end
