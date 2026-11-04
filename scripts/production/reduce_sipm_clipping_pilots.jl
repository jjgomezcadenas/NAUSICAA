#!/usr/bin/env julia
using HDF5,TOML,Statistics,SHA
module ExistingMetrics
include("reduce_electron_operating_points.jl")
end
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const OUT=joinpath(ROOT,"results/sipm_statistical_production_pilot")

"""Save full-distribution RMS/mean and a histogram on supplied common edges."""
function save_distribution(parent,key,values,edges)
    all(isfinite,values) && !isempty(values) || error("Finite, populated distribution required")
    g=create_group(parent,key);g["values"]=values;g["edges"]=edges
    counts=zeros(Int64,length(edges)-1)
    for value in values
        first(edges)<=value<=last(edges) || error("Histogram range excludes an event")
        counts[clamp(searchsortedlast(edges,value),1,length(counts))]+=1
    end
    g["counts"]=counts;g["mean"]=mean(values);g["median"]=median(values)
    g["relative_rms_percent"]=100std(values;corrected=false)/mean(values)
end

"""Compare charge loss and event-total spread for the same source events at E/P=2,2.5,3.

Use the existing expected-dark and known-attachment correction. Divide by each
source event's deposited MeV, retaining the forward SiPM nonlinearity. Save raw
and corrected totals, clipping fractions with explicit denominators, event-level
retention and histograms. This is a response-spread study of the 100-event pilot.
"""
function compare_clipping()
    timings=Dict("ep2"=>joinpath(OUT,"ep2/timing.h5"),
                 "ep2p5"=>joinpath(OUT,"timing.h5"),"ep3"=>joinpath(OUT,"timing.h5"))
    energies=h5open(joinpath(ROOT,"events/0nubb_hitdf.h5"),"r") do f
        Dict(zip(Int64.(read(f["events/event_id"])),Float64.(read(f["events/energy_MeV"]))))
    end
    reference=Dict{Tuple{String,Int64},Tuple{String,Float64}}()
    results=Dict{Tuple{String,String},Any}()
    for tag in ("ep2","ep2p5","ep3")
        timing=timings[tag]
        manifestpath=h5open(timing,"r") do f
            read(f["metadata/complete"]) || error("Incomplete pilot")
            String(read(f["$tag/manifest"]))
        end
        manifest=TOML.parsefile(manifestpath);manifest["complete"] || error("Incomplete response")
        for gas in ("pure_xe","xe_he15")
            ids=Int64[];qs=Float64[];qa=Float64[];bs=Float64[];ba=Float64[];ret=Float64[]
            clips=Int64[];acquired=Int64[];positive=Int64[];positive_clips=Int64[];loss=Float64[]
            ecard="";scard="";acard="";limit_pe=NaN
            for shard in manifest["shards"]
                shard["gas"]==gas || continue
                h5open(shard["sensor"],"r") do sf
                    h5open(shard["asic"],"r") do af
                        scard=String(read(sf["metadata/sensor_card_toml"]))
                        acard=String(read(af["metadata/electronics_card_toml"]))
                        sensor=TOML.parse(scard);dt=read(sf["metadata/bucket_width_s"])
                        limit_pe=read(af["metadata/current_limit_A"])*dt/read(af["metadata/charge_per_pe_C"])
                        h5open(String(read(sf["metadata/source_path"])),"r") do pf
                            ecard=String(read(pf["metadata/branch_card_toml"]));electron=TOML.parse(ecard)
                            for id in shard["event_ids"]
                                sg=sf["events/$id"];ag=af["events/$id"];pg=pf["events/$id"]
                                identity=(String(read(sg["pair_checksum"])),Float64(read(sg["electron_drift_cm"])))
                                if tag=="ep2";reference[(gas,id)]=identity
                                else;reference[(gas,id)]==identity || error("Primary/placement comparison mismatch")
                                end
                                read(sg["pair_checksum"])==read(ag["pair_checksum"])==read(pg["pair_checksum"]) || error("Stage pair mismatch")
                                q=read(sg["sensor_pe_equiv"]);a=read(ag["readout_pe_equiv"])
                                sat=Bool.(read(ag["saturated"]))
                                lookup=Dict(k=>i for (i,k) in enumerate(zip(read(sg["pixel_ix"]),read(sg["pixel_iy"]),read(sg["bucket"]))))
                                photon_keys=zip(read(pg["pixel_ix"]),read(pg["pixel_iy"]),read(pg["bucket"]))
                                positions=[lookup[k] for (k,n) in zip(photon_keys,read(pg["count"])) if n>0]
                                totalq=sum(q);totala=sum(a)
                                cq,_=ExistingMetrics.corrected_total(totalq,length(q),dt,sensor,electron,identity[2])
                                ca,_=ExistingMetrics.corrected_total(totala,length(q),dt,sensor,electron,identity[2])
                                push!(ids,id);push!(qs,totalq);push!(qa,totala);push!(bs,cq/energies[id]);push!(ba,ca/energies[id])
                                push!(ret,totala/totalq);push!(clips,count(sat));push!(acquired,length(q))
                                push!(positive,length(positions));push!(positive_clips,count(sat[positions]))
                                push!(loss,read(ag["clip_loss_pe_equiv"]))
                            end
                        end
                    end
                end
            end
            allunique(ids) && length(ids)==100 || error("Pilot requires the same 100 distinct events")
            results[(tag,gas)]=(;ids,qs,qa,bs,ba,ret,clips,acquired,positive,positive_clips,loss,
                ecard,scard,acard,manifestpath,timing,limit_pe)
        end
    end
    h5open(joinpath(OUT,"clipping_comparison.h5.tmp"),"w") do out
        meta=create_group(out,"metadata");meta["complete"]=false
        meta["reducer_source"]=read(@__FILE__,String)
        meta["correction_source"]=read(joinpath(@__DIR__,"reduce_electron_operating_points.jl"),String)
        for gas in ("pure_xe","xe_he15")
            normalized=vcat([results[(t,gas)].ba./mean(results[(t,gas)].ba) for t in ("ep2","ep2p5","ep3")]...)
            lo,hi=extrema(normalized);edges=collect(range(min(0.,lo),max(2.,hi);length=61))
            for tag in ("ep2","ep2p5","ep3")
                r=results[(tag,gas)];g=create_group(out,"$tag/$gas")
                electron=TOML.parse(r.ecard)
                for (k,v) in pairs((;event_id=r.ids,sensor_total_pe=r.qs,asic_total_pe=r.qa,
                    corrected_sensor_pe_per_MeV=r.bs,corrected_asic_pe_per_MeV=r.ba,
                    clipped_cells=r.clips,acquired_cells=r.acquired,photon_positive_cells=r.positive,
                    clipped_photon_positive_cells=r.positive_clips,clip_loss_pe=r.loss,
                    event_count=length(r.ids),events_with_clipping=count(>(0),r.clips),
                    events_with_clipping_percent=100count(>(0),r.clips)/length(r.ids),
                    charge_retained_percent=100sum(r.qa)/sum(r.qs),
                    clipped_charge_percent=100sum(r.loss)/sum(r.qs),
                    clipped_acquired_cells_percent=100sum(r.clips)/sum(r.acquired),
                    clipped_photon_positive_cells_percent=100sum(r.positive_clips)/sum(r.positive),
                    sensor_response_rms_percent=100std(r.bs;corrected=false)/mean(r.bs),
                    asic_response_rms_percent=100std(r.ba;corrected=false)/mean(r.ba),
                    reduced_field_kV_cm_bar=electron["reduced_field_kV_cm_bar"],J=electron["j_variance_ratio"],bucket_limit_pe=r.limit_pe,
                    electron_card=r.ecard,sensor_card=r.scard,electronics_card=r.acard,manifest=r.manifestpath))
                    g[string(k)]=v
                end
                save_distribution(g,"event_retention",r.ret,collect(range(0.,1.;length=41)))
                save_distribution(g,"normalized_asic_response",r.ba./mean(r.ba),edges)
            end
        end
        write(meta["complete"],true)
    end
    mv(joinpath(OUT,"clipping_comparison.h5.tmp"),joinpath(OUT,"clipping_comparison.h5");force=true)
    for gas in ("pure_xe","xe_he15"),tag in ("ep2","ep2p5","ep3")
        r=results[(tag,gas)]
        @show tag gas 100sum(r.qa)/sum(r.qs) 100std(r.bs;corrected=false)/mean(r.bs) 100std(r.ba;corrected=false)/mean(r.ba)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;compare_clipping();end
