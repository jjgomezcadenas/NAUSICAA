#!/usr/bin/env julia
using HDF5,TOML,Statistics,SHA
module Metrics
include("reduce_electron_operating_points.jl")
end
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const OUT=joinpath(ROOT,"results/sipm_tpb_pilot")

"""Reduce matched TPB pilots into pulse and attachment-corrected energy observables.

Photon-positive cells define pulse means/maxima at every stage. Event sums
include all acquired cells and signed electronics noise. Correct known electron
attachment and normalize by the deposited MC energy; subtract stationary dark
baseline from sensor/ASIC totals. Preserve the uncorrected recovery nonlinearity
so reported widths measure the forward response of the detector configuration.
"""
function reduce_tpb()
    energies=h5open(joinpath(ROOT,"events/0nubb_hitdf.h5"),"r") do f
        Dict(zip(read(f["events/event_id"]),Float64.(read(f["events/energy_MeV"]))))
    end
    refs=Dict{Tuple{String,Int64},Tuple{String,Float64}}()
    temp=joinpath(OUT,"analysis.h5.tmp")
    h5open(temp,"w") do out
        out["metadata/complete"]=false
        out["metadata/source"]=read(@__FILE__,String)
        out["metadata/model_note"]=read(joinpath(ROOT,"docs/sipm_tpb_pilot.md"),String)
        for setting in ("tau12p5_n512","tau12p5_n1024","tau25p0_n512")
            timing=joinpath(OUT,setting,"timing.h5")
            isfile(timing) || continue
            h5open(timing,"r") do tf
                read(tf["metadata/complete"]) || error("Incomplete timing")
                for tag in sort(filter(!=("metadata"),collect(keys(tf))))
                    manifest=String(read(tf["$tag/manifest"]))
                    run=TOML.parsefile(manifest)
                    for gas in ("pure_xe","xe_he15")
                        rows=Dict(k=>Float64[] for k in ("event_id","incident_total","detected_total","sensor_total","asic_total",
                            "incident_corrected","detected_corrected","sensor_corrected","asic_corrected",
                            "event_mean_sensor","event_max_sensor","event_mean_asic","event_max_asic","clip_loss",
                            "clipped_cells","positive_cells","clipped_positive_cells","clipped_events"))
                        pulse=Float64[];ids=Int64[];detectedpulse=Float64[];checks=String[]
                        scard="";ecard="";acard="";dt=0.;q0=0.;limit=0.;acceptance=0.;parent=""
                        for shard in run["shards"]
                            shard["gas"]==gas || continue
                            h5open(shard["sensor"],"r") do sf
                                h5open(shard["asic"],"r") do af
                                    read(sf["metadata/complete"]) && read(af["metadata/complete"]) || error("Incomplete response")
                                    scard=String(read(sf["metadata/sensor_card_toml"]));sensor=TOML.parse(scard)
                                    acard=String(read(af["metadata/electronics_card_toml"]))
                                    dt=Float64(read(sf["metadata/bucket_width_s"]));q0=Float64(read(sf["metadata/charge_per_pe_C"]))
                                    limit=Float64(read(af["metadata/current_limit_A"]))*dt/q0
                                    parent=String(read(sf["metadata/source_path"]))
                                    h5open(parent,"r") do pf
                                        ecard=String(read(pf["metadata/branch_card_toml"]));electron=TOML.parse(ecard)
                                        acceptance=Float64(read(pf["metadata/square_acceptance"]))
                                        for id in shard["event_ids"]
                                            s=sf["events/$id"];a=af["events/$id"];p=pf["events/$id"]
                                            identity=(String(read(s["pair_checksum"])),Float64(read(s["electron_drift_cm"])))
                                            get!(refs,(gas,id),identity)==identity || error("Primary/placement mismatch")
                                            read(p["pair_checksum"])==read(s["pair_checksum"])==read(a["pair_checksum"]) || error("Stage mismatch")
                                            push!(ids,id);push!(checks,identity[1])
                                            q=Float64.(read(s["sensor_pe_equiv"]));v=Float64.(read(a["readout_pe_equiv"]))
                                            d=read(s["detected_pe"]);sat=Bool.(read(a["saturated"]))
                                            lookup=Dict(k=>i for (i,k) in enumerate(zip(read(s["pixel_ix"]),read(s["pixel_iy"]),read(s["bucket"]))))
                                            positions=[lookup[k] for k in zip(read(p["pixel_ix"]),read(p["pixel_iy"]),read(p["bucket"]))]
                                            isempty(positions) && error("No positive input pulses")
                                            append!(pulse,q[positions]);append!(detectedpulse,d[positions])
                                            factor=exp(identity[2]*10/electron["drift_speed_mm_us"]/1000/electron["lifetime_ms"])/energies[id]
                                            cq,_=Metrics.corrected_total(sum(q),length(q),dt,sensor,electron,identity[2])
                                            ca,_=Metrics.corrected_total(sum(v),length(q),dt,sensor,electron,identity[2])
                                            vals=(;event_id=id,incident_total=sum(read(p["count"])),detected_total=sum(d),sensor_total=sum(q),asic_total=sum(v),
                                                incident_corrected=sum(read(p["count"]))*factor,detected_corrected=sum(d)*factor,
                                                sensor_corrected=cq/energies[id],asic_corrected=ca/energies[id],
                                                event_mean_sensor=mean(q[positions]),event_max_sensor=maximum(q[positions]),
                                                event_mean_asic=mean(v[positions]),event_max_asic=maximum(v[positions]),
                                                clip_loss=read(a["clip_loss_pe_equiv"]),clipped_cells=count(sat),positive_cells=length(positions),
                                                clipped_positive_cells=count(sat[positions]),clipped_events=any(sat))
                                            for (k,x) in pairs(vals);push!(rows[string(k)],x);end
                                        end
                                    end
                                end
                            end
                        end
                        ids==collect(Int64,0:99) || error("Expected matched 100-event pilot")
                        ep=TOML.parse(ecard)["reduced_field_kV_cm_bar"]
                        g=create_group(out,"$setting/$(ep==2 ? "ep2" : "ep2p5")/$gas")
                        for (k,v) in rows;g[k]=v;end
                        g["pair_checksum"]=checks
                        for (k,v) in pairs((;sensor_card=scard,electron_card=ecard,electronics_card=acard,manifest,
                            photon_path=parent,bucket_limit_pe=limit,square_acceptance=acceptance,
                            elapsed_s=read(tf["$tag/total_s"]),
                            charge_retained_percent=100sum(rows["asic_total"])/sum(rows["sensor_total"]),
                            clip_loss_percent=100sum(rows["clip_loss"])/sum(rows["sensor_total"]),
                            sensor_charge_retained_percent=100sum(rows["sensor_total"])/sum(rows["detected_total"]),
                            events_clipped_percent=100mean(rows["clipped_events"]),
                            positive_cells_clipped_percent=100sum(rows["clipped_positive_cells"])/sum(rows["positive_cells"])))
                            g[string(k)]=v
                        end
                        for k in ("incident_corrected","detected_corrected","sensor_corrected","asic_corrected",
                            "detected_total","sensor_total","asic_total","event_mean_sensor","event_max_sensor","event_mean_asic","event_max_asic")
                            Metrics.distribution(g,"$(k)_distribution",rows[k])
                        end
                        Metrics.distribution(g,"sensor_pulses",pulse);Metrics.distribution(g,"detected_pulses",detectedpulse)
                        positive=filter(>(0),pulse)
                        edges=exp.(collect(range(log(minimum(positive)),log(maximum(positive));length=101)))
                        edges[1]=minimum(positive);edges[end]=maximum(positive)
                        counts=zeros(Int64,100)
                        for x in positive;counts[clamp(searchsortedlast(edges,x),1,100)]+=1;end
                        g["sensor_pulses/log_edges"]=edges;g["sensor_pulses/log_counts"]=counts
                        g["sensor_pulses/zero_count"]=count(iszero,pulse)
                        sum(counts)+count(iszero,pulse)==length(pulse) || error("Pulse histogram ledger")
                        println("$setting E/P=$ep $gas: clip loss=$(round(read(g["clip_loss_percent"]);digits=4))%, ASIC RMS=$(round(read(g["asic_corrected_distribution/relative_rms_percent"]);digits=3))%")
                    end
                end
            end
        end
        if haskey(out,"tau12p5_n1024")
            for point in ("ep2","ep2p5"),gas in ("pure_xe","xe_he15")
                a=out["tau12p5_n512/$point/$gas"];b=out["tau12p5_n1024/$point/$gas"]
                g=create_group(out,"convergence/$point/$gas")
                delta=100*(read(b["sensor_total_distribution/mean"])/read(a["sensor_total_distribution/mean"])-1)
                loss=read(b["clip_loss_percent"])-read(a["clip_loss_percent"])
                g["mean_sensor_change_percent"]=delta;g["clip_loss_change_percentage_points"]=loss
                g["passed"]=abs(delta)<.1 && abs(loss)<.02
                read(g["passed"]) || error("Time convergence exceeds 0.1% charge or 0.02 percentage points clipping")
            end
        end
        write(out["metadata/complete"],true)
    end
    mv(temp,joinpath(OUT,"analysis.h5");force=true)
end
if abspath(PROGRAM_FILE)==@__FILE__;reduce_tpb();end
