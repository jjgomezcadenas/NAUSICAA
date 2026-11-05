#!/usr/bin/env julia
# Reduce saved ELM photon/linear-SiPM/ASIC histories; physics stays in Julia.
module ELMReduction
using HDF5,TOML,Statistics,SHA
include("sipm_production.jl")
const PROD=SiPMProduction
const ROOT=PROD.ROOT
const OUT=joinpath(ROOT,"results/elm_pilot")
const EA=PROD.EA
put(g,d)=PROD.put(g,d)
rms_percent(x)=100std(x;corrected=false)/mean(x)

"""Summarize 100 matched events with distinct pulse/event/charge clipping denominators.

Pulse distributions use the photon-positive SiPM/bucket cells, preserving the
same selection before and after electronics. Event totals sum all acquired
signed cells. Energy-resolution diagnostics subtract the expected dark count,
correct known attachment survival, and divide by the deposited MC energy.
The unclipped ASIC sum is reconstructed exactly using its saved clipping loss.
"""
function main()
    files=TOML.parsefile(joinpath(OUT,"files.toml"))
    manifest=TOML.parsefile(files["response_manifest"])
    run=PROD.settings("config/specs/runs/elm_bb0nu_ep2_pilot100.toml")
    branches=PROD.sources(run)
    settings=run.pilot
    q0=PROD.P.charge_per_pe_C(settings.sensor)
    limit=settings.electronics.current_limit_A*settings.electronics.bucket_width_s/q0
    output=joinpath(OUT,"elm_pilot.h5")
    h5open(output*".tmp","w") do out
        attributes(out)["schema"]="nausicaa.elm_clipping_analysis.v1"
        meta=create_group(out,"metadata")
        put(meta,(;complete=false,source_manifest=files["response_manifest"],
            source_manifest_sha256=PROD.sha(files["response_manifest"]),
            reduction_sha256=PROD.sha(@__FILE__),photon_index=files["photon_index"],
            photon_wall_s=files["photon_wall_s"],response_wall_s=files["response_wall_s"],
            total_wall_s=files["total_wall_s"],sensor_card=read(settings.cp,String),
            electronics_card=read(settings.ep,String),sample="bb0nu"))
        p=create_group(out,"parameters")
        put(p,(;pde=settings.sensor.pde,gain=settings.sensor.gain,
            dark_rate_Hz=settings.sensor.dark_rate_Hz,bucket_us=settings.electronics.bucket_width_s*1e6,
            noise_rms_pe=settings.electronics.noise_rms_pe,current_limit_mA=settings.electronics.current_limit_A*1000,
            clip_threshold_pe=limit))
        for branch in branches
            gas=branch.src.gas
            g=create_group(out,gas)
            columns=Dict{String,Vector{Float64}}(k=>Float64[] for k in
                ("event_id","energy_MeV","electron_drift_cm","ion_drift_cm","survival",
                 "photons_total","detected_total","sensor_total","unclipped_total","readout_total",
                 "clip_loss","clip_loss_percent","input_cells","clipped_input_cells","clipped_acquired_cells",
                 "acquired_cells","dark_expected","event_mean_pe","event_max_pe",
                 "event_mean_after_pe","event_max_after_pe","detected_corrected","unclipped_corrected",
                 "readout_corrected","photons_corrected"))
            inclusive=Float64[];inclusive_after=Float64[]
            quiet=Float64[];pair_checks=0;ledger_checks=0
            pde_trials=Int64(0);pde_successes=Int64(0)
            source=branch.src.path
            source_hits=h5open(source,"r") do ph
                String(read(ph["metadata/source_path"]))
            end
            h5open(source_hits,"r") do hitsfile
              h5open(source,"r") do ph
                phmeta=ph["metadata"]
                String(read(phmeta["photon_band"]))=="blue_tpb" || error("Expected blue photons")
                card=TOML.parse(String(read(phmeta["branch_card_toml"])))
                put(g,(;photon_path=source,photon_sha256=PROD.sha(source),
                    detector_card=String(read(phmeta["branch_card_toml"])),
                    optical_acceptance=read(phmeta["optical_acceptance"])))
                for entry in manifest["shards"]
                    entry["gas"]==gas || continue
                    put(g,(;sensor_path=entry["sensor"],electronics_path=entry["asic"]))
                    h5open(entry["sensor"],"r") do sensor
                      h5open(entry["asic"],"r") do asic
                        for id in entry["event_ids"]
                            pg=ph["events/$id"];sg=sensor["events/$id"];ag=asic["events/$id"]
                            String(read(sg["response_model"]))=="linear" || error("Expected linear sensor")
                            all(String(read(x["pair_checksum"]))==String(read(pg["pair_checksum"])) for x in (sg,ag)) || error("Pair mismatch")
                            pair_checks+=1
                            h=EA.read_photon_history(pg)
                            isapprox(h.grid.pitch_mm,5f0)&&isapprox(h.grid.bucket_width_mm,.25f0) || error("ELM grid")
                            pe=Float64.(read(sg["detected_pe"]));dark=Float64.(read(sg["dark_pulses"]))
                            sensor_pe=read(sg["sensor_pe_equiv"]);measured=read(ag["readout_pe_equiv"])
                            sat=read(ag["saturated"]);loss=read(ag["clip_loss_pe_equiv"])
                            @assert all(isapprox.(sensor_pe,pe+dark;atol=1e-10))
                            @assert all(measured.<=limit+1e-8)
                            @assert sum(h.cells.count)==read(sg["incident_photons_total"])
                            coordinates=Dict((x,y,b)=>i for (i,(x,y,b)) in enumerate(zip(read(sg["pixel_ix"]),read(sg["pixel_iy"]),read(sg["bucket"]))))
                            selected=[coordinates[(x,y,b)] for (x,y,b) in zip(h.cells.pixel_ix,h.cells.pixel_iy,h.cells.bucket)]
                            @assert all(pe[selected].<=h.cells.count)
                            @assert sum(pe[selected])==sum(pe)
                            ledger_checks+=1
                            append!(inclusive,pe[selected]);append!(inclusive_after,measured[selected])
                            append!(quiet,measured[findall(iszero,sensor_pe)])
                            hits,_=EA.read_event(hitsfile,Int64(id));energy=sum(Float64.(hits.energy_MeV))
                            drift=Float64(h.placement.electron_drift_cm)
                            survival=exp(-(drift*10/card["drift_speed_mm_us"]*.001)/card["lifetime_ms"])
                            dark_expected=length(pe)*settings.sensor.dark_rate_Hz*branch.src.dt
                            detector=sum(pe);raw=sum(measured)+loss;after=sum(measured)
                            photon_total=sum(h.cells.count)
                            pde_trials+=photon_total;pde_successes+=Int64(detector)
                            values=(;event_id=id,energy_MeV=energy,electron_drift_cm=drift,
                                ion_drift_cm=h.placement.ion_drift_cm,survival,
                                photons_total=photon_total,detected_total=detector,sensor_total=sum(sensor_pe),
                                unclipped_total=raw,readout_total=after,clip_loss=loss,
                                clip_loss_percent=100loss/sum(sensor_pe),input_cells=length(selected),
                                clipped_input_cells=count(identity,sat[selected]),clipped_acquired_cells=count(identity,sat),
                                acquired_cells=length(pe),dark_expected,
                                event_mean_pe=mean(pe[selected]),event_max_pe=maximum(pe[selected]),
                                event_mean_after_pe=mean(measured[selected]),event_max_after_pe=maximum(measured[selected]),
                                detected_corrected=detector/survival/energy,
                                photons_corrected=photon_total/survival/energy,
                                unclipped_corrected=(raw-dark_expected)/survival/energy,
                                readout_corrected=(after-dark_expected)/survival/energy)
                            for (key,value) in pairs(values);push!(columns[string(key)],Float64(value));end
                        end
                      end
                    end
                end
              end
            end
            for (key,value) in columns;g[key]=value;end
            g["inclusive_pe",compress=4]=inclusive
            g["inclusive_after_pe",compress=4]=inclusive_after
            npe=length(inclusive);nevents=length(columns["event_id"])
            put(g,(;event_count=nevents,pair_checks,ledger_checks,
                clipped_pulse_percent=100sum(columns["clipped_input_cells"])/sum(columns["input_cells"]),
                clipped_event_percent=100count(>(0),columns["clipped_acquired_cells"])/nevents,
                charge_lost_percent=100sum(columns["clip_loss"])/sum(columns["sensor_total"]),
                inclusive_mean_pe=mean(inclusive),inclusive_median_pe=median(inclusive),
                inclusive_max_pe=maximum(inclusive),event_mean_average_pe=mean(columns["event_mean_pe"]),
                event_mean_median_pe=median(columns["event_mean_pe"]),
                event_max_average_pe=mean(columns["event_max_pe"]),
                event_max_median_pe=median(columns["event_max_pe"]),
                detected_mean_total_pe=mean(columns["detected_total"]),
                photons_mean_total=mean(columns["photons_total"]),
                corrected_detected_rms_percent=rms_percent(columns["detected_corrected"]),
                corrected_unclipped_rms_percent=rms_percent(columns["unclipped_corrected"]),
                corrected_readout_rms_percent=rms_percent(columns["readout_corrected"]),
                quiet_noise_mean_pe=mean(quiet),quiet_noise_rms_pe=std(quiet;corrected=false),
                pde_fraction=pde_successes/pde_trials,
                pde_mean_z=(pde_successes-pde_trials*settings.sensor.pde)/sqrt(pde_trials*settings.sensor.pde*(1-settings.sensor.pde))))
            @assert abs(read(g["pde_mean_z"]))<6
            @assert abs(mean(quiet))<6settings.electronics.noise_rms_pe/sqrt(length(quiet))
            @assert abs(std(quiet;corrected=false)-settings.electronics.noise_rms_pe)<6settings.electronics.noise_rms_pe/sqrt(2length(quiet))
            println("$gas: $(round(read(g["clipped_event_percent"]);digits=2))% events clipped; $(round(read(g["charge_lost_percent"]);digits=5))% charge lost; max $(maximum(inclusive)) PE")
        end
        write(meta["complete"],true)
    end
    mv(output*".tmp",output;force=true)
    println(output)
end
end
if abspath(PROGRAM_FILE)==@__FILE__;ELMReduction.main();end
