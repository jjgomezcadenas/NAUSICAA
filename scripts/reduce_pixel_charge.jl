#!/usr/bin/env julia
using HDF5, SHA, Statistics

"""Replace a saved numerical reduction while preserving raw transport and pixel histories."""
function pixel_save(f,path,value)
    haskey(f,path) && HDF5.delete_object(f,path)
    f[path]=value isa BitArray ? Array(value) : value
end

"""Rebuild uncut pulse distributions from saved pixel/frame counts.

A pulse is an occupied cell with positive ions. Inclusive contains all pulses;
event mean/maximum use all pulses. Compute each event's true minimum over all
pulses, then include that event in the minimum histogram only if minimum > n_min.
Save exact event pulse sums/counts and independent masks for pulses and minimum
histogram acceptance. Empty events have zero summaries and no histogram entries.
No transport or collection is repeated. Use common edges across both species.
"""
function reduce_pixel_charge(output)
    maxima=Dict(name=>1.0 for name in ("inclusive","mean","maximum","minimum"));cut=0
    for sample in ("bb0nu","electrons")
        h5open(joinpath(output,sample*".h5"),"r+") do f
            ids=read(f["events/event_id"]);rows=Dict(id=>i for (i,id) in enumerate(ids))
            sums=zeros(Int64,length(ids));counts=zeros(Int64,length(ids));hi=zeros(Int64,length(ids))
            lo=fill(typemax(Int64),length(ids));cut=read(f["settings/n_min"])
            for (id,n) in zip(read(f["cells/event_id"]),read(f["cells/n_ions"]))
                if n>0
                    i=rows[id];sums[i]+=n;counts[i]+=1;hi[i]=max(hi[i],n);lo[i]=min(lo[i],n)
                end
            end
            valid=counts.>0;lo[.!valid].=0
            means=Float32[counts[i]>0 ? Float32(sums[i])/Float32(counts[i]) : 0f0 for i in eachindex(ids)]
            minimum_mask=valid .& (lo.>cut)
            @assert sums==read(f["events/collected_ions"])
            for (name,value) in (("mean_ions",means),("maximum_ions",hi),("minimum_ions",lo),
                ("pulse_ion_sum",sums),("pulse_count",counts),("has_pulse",valid),
                ("minimum_histogram_accepted",minimum_mask))
                pixel_save(f,"events/"*name,value)
            end
            for key in ("has_cell_above_threshold","cells_above_threshold")
                haskey(f,"events/"*key) && HDF5.delete_object(f,"events/"*key)
            end
            for key in ("events_with_charge_above_threshold","cells_above_threshold")
                haskey(f,"summary/"*key) && HDF5.delete_object(f,"summary/"*key)
            end
            pixel_save(f,"summary/events_with_charge",count(valid))
            pixel_save(f,"summary/pulse_count",sum(counts))
            pixel_save(f,"summary/events_in_minimum_histogram",count(minimum_mask))
            pixel_save(f,"settings/threshold_scope","Minimum histogram event selection only: true event minimum over positive cells > n_min")
            pixel_save(f,"provenance/mean_reducer_sha256",bytes2hex(sha256(read(@__FILE__))))
            pixel_save(f,"provenance/mean_definition","Uncut event pulse-ion sum/pulse count; pulse is a positive pixel-frame cell. Minimum histogram alone requires true event minimum > n_min.")
            maxima["inclusive"]=max(maxima["inclusive"],maximum(hi;init=1))
            maxima["maximum"]=max(maxima["maximum"],maximum(hi;init=1))
            maxima["mean"]=max(maxima["mean"],maximum(means;init=1f0))
            maxima["minimum"]=max(maxima["minimum"],maximum(lo[minimum_mask];init=cut+1))
        end
    end
    for sample in ("bb0nu","electrons")
        h5open(joinpath(output,sample*".h5"),"r+") do f
            valid=Bool.(read(f["events/has_pulse"]));minimum_mask=Bool.(read(f["events/minimum_histogram_accepted"]))
            ions=read(f["cells/n_ions"])
            for (name,values,step) in (("inclusive",ions[ions.>0],10),("mean",read(f["events/mean_ions"])[valid],2),
                ("maximum",read(f["events/maximum_ions"])[valid],50),("minimum",read(f["events/minimum_ions"])[minimum_mask],1))
                lower=name=="minimum" ? cut+.5 : .5
                edges=collect(lower:step:lower+step*max(1,ceil(Int,(maxima[name]-lower)/step)))
                counts=zeros(Int64,length(edges)-1)
                for value in values;counts[min(searchsortedlast(edges,value),length(counts))]+=1;end
                pixel_save(f,"histograms/"*name*"/edges_ions",edges)
                pixel_save(f,"histograms/"*name*"/counts",counts)
                pixel_save(f,"histograms/"*name*"/underflow",0);pixel_save(f,"histograms/"*name*"/overflow",0)
            end
        end
    end
end

"""Save unbinned distribution statistics and uncut collection retention.

Event maximum/mean use all events with pulses; minimum statistics use only events
whose true minimum is greater than n_min, matching that histogram. Empty selected
distributions store zero entries and NaN statistics, presented as unavailable.
Retention is all collected ions/produced ions. The minimum histogram cut removes
no charge and creates no threshold loss. Float32 statistics are checked against
Float64 references; exact Int64 charge ledgers are preserved.
"""
function summarize_pixel_charge(output)
    for sample in ("bb0nu","electrons")
        h5open(joinpath(output,sample*".h5"),"r+") do f
            valid=Bool.(read(f["events/has_pulse"]))
            produced=read(f["events/produced_ions"]);collected=read(f["events/collected_ions"])
            outside=read(f["events/outside_ions"])
            @assert all(produced.>0) && produced-collected==outside
            denominator=sum(produced);ratios=Float32.(collected)./Float32.(produced)
            haskey(f,"ion_retention") && HDF5.delete_object(f,"ion_retention")
            haskey(f,"ion_losses") && HDF5.delete_object(f,"ion_losses")
            for (key,value) in (("retention_per_event",ratios),("total_pixel_ions",sum(collected)),
                ("total_event_ions",denominator),("aggregate_retention",Float32(sum(collected))/Float32(denominator)),
                ("aggregate_loss",Float32(sum(outside))/Float32(denominator)),
                ("mean_event_retention",mean(ratios)),("median_event_retention",median(ratios)),
                ("max_Float64_error",maximum(abs.(Float64.(ratios).-Float64.(collected)./Float64.(produced)))))
                f["ion_retention/all_pixels/"*key]=value
            end
            f["ion_losses/outside_histogram_ions"]=sum(outside)
            f["ion_losses/total_event_ions"]=denominator
            f["ion_losses/outside_histogram_fraction"]=Float32(sum(outside))/Float32(denominator)
            for name in ("maximum","minimum","mean")
                mask=name=="minimum" ? Bool.(read(f["events/minimum_histogram_accepted"])) : valid
                values=Float32.(read(f["events/"*name*"_ions"])[mask])
                average=isempty(values) ? NaN32 : mean(values);middle=isempty(values) ? NaN32 : median(values)
                error64=isempty(values) ? NaN : abs(Float64(average)-mean(Float64.(values)))
                for (key,value) in (("entries",length(values)),("mean_ions",average),("median_ions",middle),("mean_Float64_error_ions",error64))
                    pixel_save(f,"distribution_statistics/"*name*"/"*key,value)
                end
            end
            pixel_save(f,"provenance/distribution_statistics_sha256",bytes2hex(sha256(read(@__FILE__))))
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    summary_only="--summary-only" in ARGS
    overrides=filter(x->startswith(x,"--minimum-cut="),ARGS)
    paths=filter(x->x!="--summary-only" && !startswith(x,"--minimum-cut="),ARGS)
    output=isempty(paths) ? joinpath(dirname(@__DIR__),"results/pixel_charge") : abspath(only(paths))
    if !isempty(overrides)
        cut=parse(Int64,split(only(overrides),'=')[2]);cut>=0 || error("Minimum cut must be nonnegative")
        for sample in ("bb0nu","electrons")
            h5open(joinpath(output,sample*".h5"),"r+") do f
                pixel_save(f,"settings/n_min",cut)
                pixel_save(f,"provenance/minimum_cut_override",cut)
            end
        end
    end
    summary_only || reduce_pixel_charge(output)
    summarize_pixel_charge(output)
    println("Saved uncut pulse histograms; minimum histogram alone selects event minimum > n_min.")
end
