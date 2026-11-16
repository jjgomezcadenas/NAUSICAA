#!/usr/bin/env julia
# Classical separation of two samples by the smaller endpoint energy E2.
#   julia --project=. scripts/production/prepare_species_comparison.jl SIGNAL_events.h5 BACKGROUND_events.h5 OUTPUT_DIR
# The inputs are the statistics files of prepare_batch_statistics.jl. For each
# track and for a cut E2 > threshold, saves the signal efficiency and the
# background acceptance against the threshold, the area under that curve and
# the background acceptance at fixed signal efficiencies, for three sets of
# events: all accepted, those with one component in that track, and "selected":
# inside the ROI and passing the track selection of that track (the lattice
# selection `itrack_selected` for the itrack, one component for the etrack).
# The ROI is one FWHM of the signal electron-energy peak around Qbb
# (config/detector/energy.toml).
using HDF5,SHA,Statistics,TOML

"""ROI bounds (keV): one FWHM of a Gaussian of the given sigma, centred on Qbb."""
roi_bounds(qbb_keV::Real,sigma_keV::Real)=(qbb_keV-2.3548sigma_keV/2,qbb_keV+2.3548sigma_keV/2)

"""Events passing the track selection of a track: the lattice selection of the itrack, one component in the etrack."""
track_selected(d::AbstractDict,arm::AbstractString)=arm=="itrack"&&haskey(d,"itrack_selected") ? d["itrack_selected"].==1 : d["$(arm)_components"].==1

"""Efficiencies of the cut `value > threshold` for signal and background values.

Events with a non-finite value fail the cut. Thresholds are the sorted distinct
finite values of both samples, preceded by -Inf. The area is the integral of
the signal efficiency over the background acceptance (trapezoids); 0.5 means no
separation.
"""
function cut_curve(signal::AbstractVector,background::AbstractVector)
    s=sort(filter(isfinite,Float64.(signal)));b=sort(filter(isfinite,Float64.(background)))
    thresholds=vcat(-Inf,sort(unique(vcat(s,b))))
    es=[(length(s)-searchsortedlast(s,t))/length(signal) for t in thresholds]
    eb=[(length(b)-searchsortedlast(b,t))/length(background) for t in thresholds]
    area=sum((eb[i]-eb[i+1])*(es[i]+es[i+1])/2 for i in 1:length(thresholds)-1)
    (;thresholds,signal_efficiency=es,background_acceptance=eb,area)
end

"""Threshold, signal efficiency and background acceptance at the loosest cut whose signal efficiency does not exceed `target`."""
function working_point(c,target::Real)
    i=findfirst(<=(target),c.signal_efficiency);i===nothing&&return (threshold=NaN,signal_efficiency=NaN,background_acceptance=NaN,background_error=NaN)
    (threshold=c.thresholds[i],signal_efficiency=c.signal_efficiency[i],background_acceptance=c.background_acceptance[i])
end

function prepare_species_comparison(signal_path,background_path,out_dir)
    mkpath(out_dir);target=joinpath(out_dir,"comparison.h5");isfile(target)&&error("Preserve existing comparison: $target")
    read_events(p)=h5open(p,"r") do f;Dict(k=>read(f["events"][k]) for k in keys(f["events"])),read(f["species"]),read(f["run_label"]),read(f["track_analysis"]);end
    s,sname,slabel,scard=read_events(signal_path);b,bname,blabel,bcard=read_events(background_path)
    scard==bcard||error("The two samples use different track settings")
    qbb=TOML.parsefile(joinpath(@__DIR__,"../../config/detector/energy.toml"))["qbb_keV"];sigma=h5open(f->read(f["energy_fits/electron_energy_keV/sigma"]),signal_path,"r")
    lower,upper=roi_bounds(qbb,sigma);in_roi(d)=(d["electron_energy_keV"].>lower).&(d["electron_energy_keV"].<upper)
    tmp=target*".tmp"
    try
        h5open(tmp,"w") do f
            f["schema"]="itaca.species_comparison.v1";f["signal/species"]=sname;f["signal/run_label"]=slabel;f["background/species"]=bname;f["background/run_label"]=blabel
            f["signal/events"]=length(s["event_number"]);f["background/events"]=length(b["event_number"])
            f["signal/file"]=relpath(signal_path,out_dir);f["background/file"]=relpath(background_path,out_dir)
            f["signal/sha256"]=bytes2hex(open(sha256,signal_path));f["background/sha256"]=bytes2hex(open(sha256,background_path));f["track_analysis"]=scard
            f["roi/lower_keV"]=lower;f["roi/upper_keV"]=upper;f["roi/qbb_keV"]=qbb;f["roi/signal_sigma_keV"]=sigma
            for (side,d) in (("signal",s),("background",b))
                f["flow/$side/accepted"]=length(d["event_number"]);f["flow/$side/roi"]=count(in_roi(d))
                for arm in ("itrack","etrack");f["flow/$side/$(arm)_selected"]=count(in_roi(d).&track_selected(d,arm));end
            end
            f["itrack_selection"]=haskey(s,"itrack_selected")&&haskey(b,"itrack_selected") ? "lattice" : "one_component"
            for arm in ("itrack","etrack"), (selection,keep) in (("all",(d)->trues(length(d["event_number"]))),("one_component",(d)->d["$(arm)_components"].==1),("selected",(d)->in_roi(d).&track_selected(d,arm)))
                ks=keep(s);kb=keep(b);c=cut_curve(s["$(arm)_E2_keV"][ks],b["$(arm)_E2_keV"][kb]);g="E2_cut/$arm/$selection"
                # The curve has one point per event; save it thinned to at most 2000 points
                step=max(1,length(c.thresholds)÷2000);idx=unique(vcat(1:step:length(c.thresholds),length(c.thresholds)))
                f["$g/threshold_keV"]=c.thresholds[idx];f["$g/signal_efficiency"]=c.signal_efficiency[idx];f["$g/background_acceptance"]=c.background_acceptance[idx]
                f["$g/area"]=c.area;f["$g/signal_events"]=count(ks);f["$g/background_events"]=count(kb)
                f["$g/signal_selected_fraction"]=count(ks)/length(ks);f["$g/background_selected_fraction"]=count(kb)/length(kb)
                targets=[0.9,0.8,0.7,0.6,0.5];w=[working_point(c,t) for t in targets]
                f["$g/working_points/target_signal_efficiency"]=targets;f["$g/working_points/threshold_keV"]=[x.threshold for x in w]
                f["$g/working_points/signal_efficiency"]=[x.signal_efficiency for x in w];f["$g/working_points/background_acceptance"]=[x.background_acceptance for x in w]
                # Binomial standard error of the background acceptance
                f["$g/working_points/background_acceptance_error"]=[sqrt(max(x.background_acceptance*(1-x.background_acceptance),0)/count(kb)) for x in w]
            end
        end
        Base.rename(tmp,target)
    catch
        isfile(tmp)&&rm(tmp);rethrow()
    end
    println("COMPARISON $target");target
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==3||error("Usage: prepare_species_comparison.jl SIGNAL_events.h5 BACKGROUND_events.h5 OUTPUT_DIR")
    prepare_species_comparison(abspath(ARGS[1]),abspath(ARGS[2]),abspath(ARGS[3]))
end
