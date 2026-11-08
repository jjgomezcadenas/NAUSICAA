#!/usr/bin/env julia
using HDF5, TOML, Statistics, SHA
const ROOT=normpath(joinpath(@__DIR__,"../.."))
resolve(p)=isabspath(p) ? p : joinpath(ROOT,p)

"""Summarise values on the positive-input support, preserving signed response.

Returns the number of selected channel–bucket cells, total, mean and maximum.
An empty support gives zero count and total and NaN mean/maximum.
"""
function response_metrics(input, values)
    size(input)==size(values) || throw(DimensionMismatch("Pulse arrays must align"))
    selected=Float64.(values[input .> 0])
    all(isfinite,selected) || error("Nonfinite response")
    isempty(selected) && return (0.0,0.0,NaN,NaN)
    (Float64(length(selected)),sum(selected),mean(selected),maximum(selected))
end

"""Write per-event statistics and their ensemble distributions to HDF5.

Each histogram has 80 linear bins including the complete finite range. Save
ensemble means and medians for presentation without further numerical reduction.
"""
function save_series(out,name,ids,rows,unit)
    g=create_group(out,name);g["event_id"]=ids;g["unit"]=unit
    for (i,key) in enumerate(("occupied_cells","total","mean","maximum"))
        values=[r[i] for r in rows];g[key]=values
        v=filter(isfinite,values)
        isempty(v) && continue
        g[key*"_ensemble_mean"]=mean(v);g[key*"_ensemble_median"]=median(v)
        lo,hi=extrema(v);lo==hi && (lo-=0.5;hi+=0.5)
        edges=collect(range(lo,hi;length=81));counts=zeros(Int64,80)
        for x in v;counts[clamp(searchsortedlast(edges,x),1,80)]+=1;end
        g[key*"_edges"]=edges;g[key*"_histogram"]=counts
    end
end

"""Give overlaid distributions identical bin edges, retaining all finite values."""
function shared_histograms!(out)
    families=[["ions_$(p)mm_$s" for p in (2,3,4) for s in ("input","corrected")]]
    append!(families,[["$(g)_$s" for g in ("pure_xe","xe_he15")] for s in ("photons","detected","sensor","readout")])
    for family in families,metric in ("mean","maximum")
        values=[filter(isfinite,read(out[name*"/"*metric])) for name in family]
        lo,hi=extrema(vcat(values...));lo==hi && (lo-=0.5;hi+=0.5)
        edges=collect(range(lo,hi;length=81))
        for (name,v) in zip(family,values)
            counts=zeros(Int64,80)
            for x in v;counts[clamp(searchsortedlast(edges,x),1,80)]+=1;end
            write(out[name*"/"*metric*"_edges"],edges)
            write(out[name*"/"*metric*"_histogram"],counts)
        end
    end
end

"""Stream completed CSA and SiPM files into paired event-level pulse summaries.

Ion support is counts > 0; electron support is incident photons > 0. Photon
keys are joined exactly to SiPM response keys. Pair checksums and drift lengths
must agree across all branches. The source card and code digest are retained.
"""
function main(cardpath;limit=typemax(Int))
    card=TOML.parsefile(cardpath)
    # Check inputs before opening the saved reduction for replacement.
    sources=[resolve(source["path"]) for source in card["ions"]]
    push!(sources,resolve(card["sipm_manifest"]))
    missing=filter(p->!isfile(p),sources)
    if !isempty(missing)
        advice=get(card,"historical_sources_retired",false) ?
            "Historical response files were retired. Use make response-overview-present for saved plots or make response-viewer for the current matched production." :
            "Supply completed response inputs before reducing."
        throw(ArgumentError(advice*" Missing input: "*first(missing)))
    end
    manifest=TOML.parsefile(resolve(card["sipm_manifest"]))
    for shard in manifest["shards"]
        shard["sample"]==card["sample"] || continue
        for key in ("sensor","asic")
            isfile(resolve(shard[key])) || throw(ArgumentError("Missing $key response: $(shard[key]); saved reduction preserved"))
        end
    end
    outdir=resolve(card["output"]);mkpath(outdir)
    ledger=Dict{Int,Tuple{String,Float64,Float64}}()
    function check(g,id)
        signature=(String(read(g["pair_checksum"])),Float64(read(g["ion_drift_cm"])),Float64(read(g["electron_drift_cm"])))
        haskey(ledger,id) && ledger[id]!=signature && error("Event pairing mismatch $id")
        ledger[id]=signature
    end
    h5open(joinpath(outdir,"response_overview.h5"),"w") do out
        attributes(out)["schema"]="nausicaa.response_overview.v1"
        out["card_toml"]=read(cardpath,String);out["reducer_sha256"]=bytes2hex(sha256(read(@__FILE__)))
        out["sample"]=card["sample"]
        for source in card["ions"]
            p=source["pitch_mm"];println("Reducing ions $p mm");flush(stdout)
            h5open(resolve(source["path"]),"r") do f
                read(f["metadata/complete"]) || error("Incomplete ion source")
                ids=sort(parse.(Int,keys(f["events"])))[1:min(limit,length(f["events"]))]
                rowsin=[];rowsout=[]
                for id in ids
                    g=f["events/$id"];check(g,id);input=read(g["counts"])
                    push!(rowsin,response_metrics(input,input))
                    push!(rowsout,response_metrics(input,read(g["recovered_ions"])))
                end
                save_series(out,"ions_$(p)mm_input",ids,rowsin,"ions")
                save_series(out,"ions_$(p)mm_corrected",ids,rowsout,"ions")
            end
        end
        manifest=TOML.parsefile(resolve(card["sipm_manifest"]))
        out["sipm_manifest_toml"]=read(resolve(card["sipm_manifest"]),String)
        manifest["complete"] || error("Incomplete SiPM manifest")
        for gas in ("pure_xe","xe_he15")
            println("Reducing electrons $gas");flush(stdout)
            ids=Int[];rows=[[] for _ in 1:4]
            for shard in manifest["shards"]
                shard["gas"]==gas && shard["sample"]==card["sample"] || continue
                length(ids)>=limit && break
                h5open(resolve(shard["sensor"]),"r") do sensor
                    h5open(resolve(shard["asic"]),"r") do asic
                        h5open(String(read(sensor["metadata/source_path"])),"r") do photons
                            for id in shard["event_ids"]
                                length(ids)>=limit && break
                                s=sensor["events/$id"];a=asic["events/$id"];g=photons["events/$id"]
                                check(s,id);check(a,id);check(g,id)
                                sk=collect(zip(read(s["pixel_ix"]),read(s["pixel_iy"]),read(s["bucket"])))
                                ak=collect(zip(read(a["pixel_ix"]),read(a["pixel_iy"]),read(a["bucket"])))
                                sk==ak || error("Sensor/ASIC key mismatch")
                                lookup=Dict(k=>i for (i,k) in enumerate(sk))
                                pk=zip(read(g["pixel_ix"]),read(g["pixel_iy"]),read(g["bucket"]))
                                indices=[lookup[k] for k in pk];counts=read(g["count"])
                                push!(ids,id);push!(rows[1],response_metrics(counts,counts))
                                for (j,group,key) in ((2,s,"detected_pe"),(3,s,"sensor_pe_equiv"),(4,a,"readout_pe_equiv"))
                                    push!(rows[j],response_metrics(counts,read(group[key])[indices]))
                                end
                            end
                        end
                    end
                end
            end
            for (i,stage) in enumerate(("photons","detected","sensor","readout"))
                save_series(out,"$(gas)_$stage",ids,rows[i],i==1 ? "photons" : "PE")
            end
        end
        shared_histograms!(out)
        out["validated_paired_events"]=length(ledger);out["complete"]=true
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    main(isempty(ARGS) ? resolve("config/specs/runs/response_overview.toml") : ARGS[1];limit=length(ARGS)>1 ? parse(Int,ARGS[2]) : typemax(Int))
end
