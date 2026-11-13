#!/usr/bin/env julia
# Reduce saved detector counts; no transport is rerun.
include("run_framed_detector_response.jl")
include("../../src/analysis/gaussian_peak.jl")
using .GaussianPeak
using Distributions: Normal, cdf, pdf
function reduce_energy(path)
    target=joinpath(dirname(path),"energy_controls.h5");tmp=target*".tmp"
    budget=8*1024^2
    GuardedFiles.space_check(dirname(path),budget,1024^3)
    try
      h5open(path,"r") do f
        c=TOML.parse(read(f["card"],String));ids=read(f["events/event_number"])
        s=event_settings(c,EventPlacement(first(ids),75f0,75f0))
        W=Float64(s.W_eV)/1000
        photon_scale=W/(Float64(s.el_photons_per_electron)*Float64(s.optical_acceptance))
        L=read(f["events/electron_drift_cm"])
        correction=exp.(Float64.(L).*10 ./ c["electron_drift_velocity_mm_us"] ./1000 ./ c["electron_lifetime_ms"])
        h5open(tmp,"w") do out
            out["schema"]="itaca.framed_detector_energy.v1"
            out["input_sha256"]=bytes2hex(open(sha256,path))
            out["code_sha256"]=bytes2hex(open(sha256,@__FILE__))
            out["fit_code_sha256"]=bytes2hex(open(sha256,joinpath(ROOT,"src/analysis/gaussian_peak.jl")))
            out["event_number"]=ids;out["electron_drift_cm"]=L;out["attachment_correction"]=correction
            out["card"]=read(f["card"],String)
            out["ion_keV_per_count"]=W;out["photon_keV_per_count_before_attachment_correction"]=photon_scale
            for arm in ("ions","photons")
                inside=read(f["$arm/events/collected"]);outside=read(f["$arm/events/loss"]);total=read(f["$arm/events/total"])
                @assert inside+outside==total
                countrows=[(inside=inside[i],outside=outside[i],total=total[i]) for i in eachindex(ids)]
                histograms(create_group(out,"$arm/counts"),countrows,c["histogram_bins"])
                corrected_counts=[(inside=inside[i]*(arm=="ions" ? 1. : correction[i]),outside=outside[i]*(arm=="ions" ? 1. : correction[i]),total=total[i]*(arm=="ions" ? 1. : correction[i])) for i in eachindex(ids)]
                histograms(create_group(out,"$arm/corrected_counts"),corrected_counts,c["histogram_bins"])
                means=read(f["$arm/events/mean"]);maxima=read(f["$arm/events/maximum"])
                corrected_signal=[(collected=corrected_counts[i].inside,mean=means[i]*(arm=="ions" ? 1. : correction[i]),maximum=maxima[i]*(arm=="ions" ? 1. : correction[i])) for i in eachindex(ids)]
                histograms(create_group(out,"$arm/corrected_signal"),corrected_signal,c["histogram_bins"])

                scale=arm=="ions" ? W : photon_scale
                raw=[(inside_keV=inside[i]*scale,outside_keV=outside[i]*scale,total_keV=total[i]*scale) for i in eachindex(ids)]
                corrected=[(inside_keV=inside[i]*scale*(arm=="ions" ? 1. : correction[i]),outside_keV=outside[i]*scale*(arm=="ions" ? 1. : correction[i]),total_keV=total[i]*scale*(arm=="ions" ? 1. : correction[i])) for i in eachindex(ids)]
                histograms(create_group(out,"$arm/raw_energy"),raw,c["histogram_bins"])
                histograms(create_group(out,"$arm/corrected_energy"),corrected,c["histogram_bins"])
                columns(create_group(out,"$arm/events_raw"),raw);columns(create_group(out,"$arm/events_corrected"),corrected)
                for key in (:inside_keV,:total_keV)
                    x=[r[key] for r in corrected]
                    fit=fit_peak(x)
                    group=create_group(out,"$arm/peak_fits/$key")
                    for k in keys(fit);group[string(k)]=fit[k];end
                    for (label,width) in (("narrow",2.),("wide",3.))
                        r=fit_peak(x;window_sigma=width)
                        for k in keys(r);group["$label/"*string(k)]=r[k];end
                    end
                    edges=collect(range(fit.centre-5fit.scale,fit.centre+5fit.scale;length=101))
                    counts=zeros(Int64,100)
                    for v in x
                        first(edges)<=v<=last(edges)&&(counts[min(searchsortedlast(edges,v),100)]+=1)
                    end
                    group["plot_edges"]=edges;group["plot_counts"]=counts
                    curve=collect(range(fit.lower,fit.upper;length=300));dist=Normal(fit.mu,fit.sigma)
                    norm=fit.entries/(cdf(dist,fit.upper)-cdf(dist,fit.lower))
                    group["curve_x"]=curve;group["curve_y"]=norm.*pdf.(dist,curve).*(edges[2]-edges[1])
                    println("$arm $key peak: sigma=$(fit.sigma) ± $(fit.sigma_error) keV; mu=$(fit.mu); window [$(fit.lower),$(fit.upper)]")
                end
                x=[r.outside_keV for r in corrected]
                edges=collect(range(0.,10.;length=51));counts=zeros(Int64,50)
                for v in x;0<=v<=10&&(counts[min(searchsortedlast(edges,v),50)]+=1);end
                group=create_group(out,"$arm/outside_display")
                group["edges"]=edges;group["counts"]=counts;group["above_range"]=count(>(10),x)
                group["entries"]=length(x);group["mean"]=mean(x)
                for (label,rows) in (("raw",raw),("corrected",corrected))
                    x=[r.total_keV for r in rows]
                    out["$arm/$(label)_relative_std"]=std(x;corrected=false)/mean(x)
                    println("$arm $label total: mean $(mean(x)) keV; std $(std(x;corrected=false)) keV; relative std $(100std(x;corrected=false)/mean(x)) %")
                end
            end
            idrow=Dict(id=>i for (i,id) in enumerate(ids))
            for examples in filter(k->haskey(f,k),("examples","largest_losses","frame_checks")),id in keys(f[examples])
                factor=correction[idrow[parse(Int64,id)]]
                for name in ("xy_count","xz_count","outside_xy_count","outside_xz_count")
                    out["photons/example_counts/$examples/$id/$name"]=Float64.(read(f["$examples/$id/photons/$name"])).*factor
                end
            end
            out["expected_survival_min"]=1/maximum(correction);out["expected_survival_max"]=1/minimum(correction)
            out["photon_total_vs_drift_correlation"]=cor(Float64.(read(f["photons/events/total"])),Float64.(L))
            println("Photon total vs drift correlation: ",read(out["photon_total_vs_drift_correlation"]))
            println("Survival range: ",extrema(1 ./ correction))
            flush(out);GuardedFiles.check_size(tmp,budget)
        end
      end
      mv(tmp,target;force=true)
    catch
      isfile(tmp)&&rm(tmp);rethrow()
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    path=isempty(ARGS) ? joinpath(ROOT,"results/readout/bb0nu_full_ion3mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42/detector_response/controls.h5") : abspath(ARGS[1])
    reduce_energy(path)
end
