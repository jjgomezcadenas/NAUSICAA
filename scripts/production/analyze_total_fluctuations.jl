#!/usr/bin/env julia
# Reduce saved 100-event detector histories to event-total fluctuation observables.
using HDF5, TOML, Statistics, SHA

const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis

"""Read a scalar HDF5 dataset, preserving its numerical or string value."""
scalar(group,name)=read(group[name])

"""Return the predicted pair count and integer-compatible Fano variance from source hits."""
function pair_moments(hits::EventHits,W_eV::Float64,F::Float64)
    means=Float64.(hits.energy_MeV).*(1e6/W_eV)
    fractions=means.-floor.(means)
    variances=max.(F.*means,fractions.*(1 .- fractions))
    return sum(means),sum(variances)
end

"""Return variance of a photon total converted to the pair-count scale.

The terms are shared-pair Fano variance, electron attachment, conditional EL
variance, and binomial optical collection. Event-specific drift survival is
applied before comparing photon totals across the sampled TPC positions.
"""
function photon_equivalent_variance(μ::Float64,v_pairs::Float64,
        survival::Float64,Y::Float64,J::Float64,A::Float64)
    return v_pairs + μ*(1-survival)/survival +
        J*μ/(survival*Y) + μ*(1-A)/(survival*Y*A)
end

"""Read one completed pilot summary and write aligned per-event numerical results."""
function analyze(summary_path::String,output_path::String)
    h5open(summary_path,"r") do summary
        String(read(attributes(summary)["schema"]))=="nausicaa.detector_pilot_summary.v1" ||
            error("Unexpected pilot summary schema")
        Bool(scalar(summary["metadata"],"complete")) || error("Incomplete pilot summary")
        group=summary["bb0nu"]
        n=Int(scalar(group,"events"))
        r=TOML.parse(String(scalar(summary["metadata"],"realisation_card")))
        e=TOML.parse(String(scalar(summary["metadata"],"electron_card")))
        W=Float64(r["w_eV"]);F=Float64(r["fano_factor"])
        J=Float64(e["j_variance_ratio"]);A=Float64(scalar(summary["derived"],"optical_acceptance"))
        Yxe=Float64(scalar(summary["derived/pure_xe"],"yield_photons_per_e"))
        Yhe=Float64(scalar(summary["derived/xe_he15"],"yield_photons_per_e"))
        speed=Float64(e["drift_speed_mm_us"]);lifetime=Float64(e["lifetime_ms"])
        source_path=String(scalar(group,"realisation_path"))
        ion_path=String(scalar(group,"ions_path"))
        xe_path=String(scalar(group,"pure_xe_path"))
        he_path=String(scalar(group,"xe_he15_path"))
        event_id=Vector{Int64}(undef,n)
        expected_pairs=zeros(Float64,n);pair_variance=zeros(Float64,n)
        ions=zeros(Int64,n);survivors=zeros(Int64,n)
        drift_cm=zeros(Float64,n);survival=zeros(Float64,n)
        xe_photons=zeros(Int64,n);he_photons=zeros(Int64,n)
        xe_equiv=zeros(Float64,n);he_equiv=zeros(Float64,n)
        xe_variance=zeros(Float64,n);he_variance=zeros(Float64,n)
        hit_path=Ref("")
        h5open(source_path,"r") do real
            h5open(ion_path,"r") do ion
                h5open(xe_path,"r") do xe
                    h5open(he_path,"r") do he
                        input_path=String(scalar(real["metadata"],"source_path"))
                        hit_path[]=input_path
                        h5open(input_path,"r") do input
                            ids=read(input["events/event_id"])[1:n]
                            for (i,id) in enumerate(ids)
                                key=string(id);event_id[i]=id
                                hits,_=read_event(input,Int64(id))
                                expected_pairs[i],pair_variance[i]=pair_moments(hits,W,F)
                                rr=real["events/$key"];ii=ion["events/$key"]
                                xx=xe["events/$key"];hh=he["events/$key"]
                                ions[i]=Int64(scalar(rr,"initial_pairs"))
                                ions[i]==sum(read(ii["count"])) || error("Ion ledger failed at $key")
                                survivors[i]=Int64(scalar(xx,"surviving_electrons"))
                                survivors[i]==Int64(scalar(hh,"surviving_electrons")) ||
                                    error("Gas branches disagree on attachment at $key")
                                drift_cm[i]=Float64(scalar(rr,"electron_drift_cm"))
                                survival[i]=exp(-10*drift_cm[i]/(1000*speed*lifetime))
                                xe_photons[i]=sum(read(xx["count"]))
                                he_photons[i]=sum(read(hh["count"]))
                                xe_equiv[i]=xe_photons[i]/(survival[i]*Yxe*A)
                                he_equiv[i]=he_photons[i]/(survival[i]*Yhe*A)
                                xe_variance[i]=photon_equivalent_variance(
                                    expected_pairs[i],pair_variance[i],survival[i],Yxe,J,A)
                                he_variance[i]=photon_equivalent_variance(
                                    expected_pairs[i],pair_variance[i],survival[i],Yhe,J,A)
                            end
                        end
                    end
                end
            end
        end
        mkpath(dirname(output_path))
        h5open(output_path,"w") do out
            attributes(out)["schema"]="nausicaa.total_fluctuations.v1"
            m=create_group(out,"metadata")
            m["source_summary"]=abspath(summary_path)
            m["source_summary_sha256"]=bytes2hex(sha256(read(summary_path)))
            m["source_hits"]=hit_path[]
            m["sample"]="bb0nu";m["events"]=n
            m["W_eV"]=W;m["fano_factor"]=F;m["j_variance_ratio"]=J
            m["optical_acceptance"]=A;m["xe_yield_photons_per_e"]=Yxe
            m["he_yield_photons_per_e"]=Yhe
            m["electron_lifetime_ms"]=lifetime;m["electron_speed_mm_us"]=speed
            m["normalization"]="photons/(event survival * mean EL yield * optical acceptance)"
            m["pair_variance_model"]="sum per hit max(F*mu, fractional(mu)*(1-fractional(mu)))"
            m["photon_variance_model"]="pair + attachment + EL + optical; on pair-equivalent scale"
            d=create_group(out,"events")
            for (name,data) in (("event_id",event_id),("expected_pairs",expected_pairs),
                    ("pair_variance",pair_variance),("ions",ions),("survivors",survivors),
                    ("electron_drift_cm",drift_cm),("survival",survival),
                    ("xe_photons",xe_photons),("he_photons",he_photons),
                    ("xe_pair_equivalent",xe_equiv),("he_pair_equivalent",he_equiv),
                    ("expected_energy_MeV",expected_pairs.*W./1e6),
                    ("ion_energy_MeV",ions.*W./1e6),
                    ("xe_energy_MeV",xe_equiv.*W./1e6),
                    ("he_energy_MeV",he_equiv.*W./1e6),
                    ("ion_energy_residual_keV",(ions.-expected_pairs).*W./1e3),
                    ("xe_energy_residual_keV",(xe_equiv.-expected_pairs).*W./1e3),
                    ("he_energy_residual_keV",(he_equiv.-expected_pairs).*W./1e3),
                    ("xe_pair_variance",xe_variance),("he_pair_variance",he_variance))
                d[name]=data
            end
            s=create_group(out,"summary")
            s["source_mean_pairs"]=mean(expected_pairs)
            s["source_sd_pairs"]=std(expected_pairs)
            s["source_mean_energy_MeV"]=mean(expected_pairs)*W/1e6
            s["mean_ions"]=mean(ions)
            s["sd_ions"]=std(ions)
            s["sd_ion_residual"]=std(ions.-expected_pairs)
            s["predicted_ion_sd"]=sqrt(mean(pair_variance))
            s["ion_relative_sd_percent"]=100*std(ions.-expected_pairs)/mean(expected_pairs)
            s["ion_predicted_relative_sd_percent"]=100*sqrt(mean(pair_variance))/mean(expected_pairs)
            s["ion_energy_residual_sd_keV"]=std((ions.-expected_pairs).*W./1e3)
            s["ion_predicted_energy_sd_keV"]=sqrt(mean(pair_variance))*W/1e3
            s["mean_pair_variance"]=mean(pair_variance)
            s["mean_attachment_variance"]=mean(expected_pairs.*(1 .- survival)./survival)
            s["mean_survival"]=mean(survival)
            s["min_survival"]=minimum(survival)
            s["max_survival"]=maximum(survival)
            for (gas,photons,equivalent,variances) in
                    (("xe",xe_photons,xe_equiv,xe_variance),
                     ("he",he_photons,he_equiv,he_variance))
                s["$(gas)_mean_photons"]=mean(photons)
                s["$(gas)_sd_photons"]=std(photons)
                s["$(gas)_mean_pair_equivalent"]=mean(equivalent)
                s["$(gas)_sd_pair_equivalent"]=std(equivalent)
                s["$(gas)_sd_residual"]=std(equivalent.-expected_pairs)
                s["$(gas)_predicted_sd"]=sqrt(mean(variances))
                s["$(gas)_relative_sd_percent"]=100*std(equivalent.-expected_pairs)/mean(expected_pairs)
                s["$(gas)_predicted_relative_sd_percent"]=100*sqrt(mean(variances))/mean(expected_pairs)
                s["$(gas)_energy_residual_sd_keV"]=std((equivalent.-expected_pairs).*W./1e3)
                s["$(gas)_predicted_energy_sd_keV"]=sqrt(mean(variances))*W/1e3
                s["$(gas)_residual_drift_correlation"]=cor(drift_cm,equivalent.-expected_pairs)
                Y=gas=="xe" ? Yxe : Yhe
                s["$(gas)_mean_el_variance"]=mean(J.*expected_pairs./(survival.*Y))
                s["$(gas)_mean_optical_variance"]=mean(
                    expected_pairs.*(1-A)./(survival.*Y.*A))
            end
        end
    end
    return output_path
end

length(ARGS)==2 || error("Usage: analyze_total_fluctuations.jl SUMMARY.h5 OUTPUT.h5")
println(analyze(ARGS[1],ARGS[2]))
