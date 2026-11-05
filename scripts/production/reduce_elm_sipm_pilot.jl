#!/usr/bin/env julia
"""Reduce ELM finite-cell SiPM response into pulse, occupancy and energy diagnostics."""
module ELMSiPMReduction
using HDF5,TOML,Statistics
include("sipm_production.jl")
const PROD=SiPMProduction
const ROOT=PROD.ROOT
const OUT=joinpath(ROOT,"results/elm_sipm_pilot")
put(g,d)=PROD.put(g,d)
rms_percent(x)=100std(x;corrected=false)/mean(x)

"""Electron reference bounds from shared placement or an explicit legacy override."""
function electron_reference_bounds(e,r)
    chamber=Float64(r["chamber_length_cm"])
    if haskey(e,"placement")
        p=e["placement"]
        return (Float64(p["minimum_electron_drift_cm"]),Float64(p["maximum_electron_drift_cm"]),chamber)
    end
    (chamber-Float64(get(r,"max_ion_drift_cm",chamber)),chamber-Float64(r["min_ion_drift_cm"]),chamber)
end

"""Compare matched photon histories for two recovery constants and both gases.

Input-positive SiPM/bucket cells define all pulse comparisons. Optical triggers
count PDE successes; optical charge is measured in fully recovered avalanche
units. Distinct occupancy is the conditional expected number of cells hit by
optical triggers, divided by the illuminated cell count. Event totals preserve
signed noise over the complete acquisition window. Known attachment survival
and deposited MC energy normalize every resolution diagnostic.
"""
function main(;output_dir=OUT,report_name="elm_sipm_pilot",reference_results=nothing,store_inclusive=true)
    files=TOML.parsefile(joinpath(output_dir,"files.toml"))
    output=joinpath(output_dir,report_name*".h5")
    h5open(output*".tmp","w") do out
        attributes(out)["schema"]="nausicaa.elm_sipm_analysis.v1"
        meta=create_group(out,"metadata")
        put(meta,(;complete=false,files_index=read(joinpath(output_dir,"files.toml"),String),
            reduction_sha256=PROD.sha(@__FILE__),photon_wall_s=files["photon_wall_s"],
            response_wall_s=files["response_wall_s"],total_wall_s=files["total_wall_s"]))
        for tag in sort!(collect(keys(files["responses"])))
          entry=files["responses"][tag];run=PROD.settings(entry["card"])
          settings=run.pilot;manifest=TOML.parsefile(entry["manifest"])
          q0=PROD.P.charge_per_pe_C(settings.sensor)
          limit=settings.electronics.current_limit_A*settings.electronics.bucket_width_s/q0
          M=settings.c["illuminated_microcell_count"]
          tg=create_group(out,"tau_"*tag);p=create_group(tg,"parameters")
          put(p,(;pde=settings.sensor.pde,gain=settings.sensor.gain,
              dark_rate_Hz=settings.sensor.dark_rate_Hz,recovery_ns=settings.sensor.recovery_s*1e9,
              physical_microcells=settings.sensor.microcells,illuminated_microcells=M,
              bucket_us=settings.electronics.bucket_width_s*1e6,noise_rms_pe=settings.electronics.noise_rms_pe,
              current_limit_mA=settings.electronics.current_limit_A*1000,clip_threshold_pe=limit,
              response_wall_s=entry["elapsed_s"],sensor_card=read(settings.cp,String),
              electronics_card=read(settings.ep,String),source_manifest=entry["manifest"]))
          for branch in PROD.sources(run)
            gas=branch.src.gas;g=create_group(tg,gas)
            columns=Dict{String,Vector{Float64}}()
            inclusive=Dict(k=>Float64[] for k in ("detected_pe","optical_charge_pe","after_pe","occupancy_percent"))
            quiet_n=0;quiet_sum=0.;quiet_sumsq=0.;checks=0;repairs=0;limits=0;maxrepair=0.;batches=0;pde_trials=Int64(0);pde_successes=Int64(0)
            source=branch.src.path
            source_hits=h5open(source,"r") do ph;String(read(ph["metadata/source_path"]));end
            h5open(source_hits,"r") do hitsfile
              h5open(source,"r") do ph
                card=TOML.parse(String(read(ph["metadata/branch_card_toml"])))
                realization=TOML.parse(String(read(ph["metadata/realisation_card_toml"])))
                min_drift,max_drift,chamber=electron_reference_bounds(card,realization)
                @assert String(read(ph["metadata/photon_band"]))=="blue_tpb"
                put(g,(;photon_path=source,photon_sha256=PROD.sha(source),
                    detector_card=String(read(ph["metadata/branch_card_toml"])),
                    optical_acceptance=read(ph["metadata/optical_acceptance"]),
                    minimum_electron_drift_cm=min_drift,maximum_electron_drift_cm=max_drift,chamber_length_cm=chamber))
                matching_shards=filter(x->x["gas"]==gas,manifest["shards"])
                put(g,(;sensor_paths=[x["sensor"] for x in matching_shards],
                    electronics_paths=[x["asic"] for x in matching_shards]))
                for shard in matching_shards
                    shard["gas"]==gas || continue
                    if !haskey(g,"sensor_path")
                        put(g,(;sensor_path=shard["sensor"],electronics_path=shard["asic"]))
                    end
                    h5open(shard["sensor"],"r") do sensor
                      h5open(shard["asic"],"r") do asic
                        for id in shard["event_ids"]
                            pg=ph["events/$id"];sg=sensor["events/$id"];ag=asic["events/$id"]
                            @assert String(read(sg["response_model"]))=="statistical"
                            @assert all(String(read(x["pair_checksum"]))==String(read(pg["pair_checksum"])) for x in (sg,ag))
                            h=PROD.EA.read_photon_history(pg)
                            pe=Float64.(read(sg["detected_pe"]));dark=Float64.(read(sg["dark_pulses"]))
                            optical=read(sg["optical_charge_pe_equiv"]);dq=read(sg["dark_charge_pe_equiv"])
                            sensor_pe=read(sg["sensor_pe_equiv"]);measured=read(ag["readout_pe_equiv"])
                            distinct=read(sg["expected_distinct_optical_cells"])
                            sat=read(ag["saturated"]);loss=read(ag["clip_loss_pe_equiv"])
                            @assert all(isapprox.(sensor_pe,optical+dq;atol=1e-8))
                            @assert all(optical.<=pe.+1e-8) && all(dq.<=dark.+1e-8)
                            @assert all(distinct.<=min.(pe,M)) && all(measured.<=limit+1e-8)
                            @assert sum(h.cells.count)==read(sg["incident_photons_total"])
                            coordinates=Dict((x,y,b)=>i for (i,(x,y,b)) in enumerate(zip(read(sg["pixel_ix"]),read(sg["pixel_iy"]),read(sg["bucket"]))))
                            selected=[coordinates[(x,y,b)] for (x,y,b) in zip(h.cells.pixel_ix,h.cells.pixel_iy,h.cells.bucket)]
                            @assert all(pe[selected].<=h.cells.count) && sum(pe[selected])==sum(pe)
                            checks+=1;repairs+=read(sg["statistical_state_repairs"]);limits+=read(sg["statistical_variance_limits"])
                            maxrepair=max(maxrepair,read(sg["maximum_statistical_state_repair"]));batches+=read(sg["statistical_batches"])
                            if store_inclusive
                                append!(inclusive["detected_pe"],pe[selected]);append!(inclusive["optical_charge_pe"],optical[selected])
                                append!(inclusive["after_pe"],measured[selected]);append!(inclusive["occupancy_percent"],100distinct[selected]/M)
                            end
                            quiet_values=measured[findall(iszero,sensor_pe)]
                            quiet_n+=length(quiet_values);quiet_sum+=sum(quiet_values);quiet_sumsq+=sum(abs2,quiet_values)
                            hits,_=PROD.EA.read_event(hitsfile,Int64(id));energy=sum(Float64.(hits.energy_MeV))
                            drift=Float64(h.placement.electron_drift_cm)
                            @assert min_drift<=drift<=max_drift && isapprox(drift+h.placement.ion_drift_cm,chamber;atol=1e-4)
                            survival=exp(-(drift*10/card["drift_speed_mm_us"]*.001)/card["lifetime_ms"])
                            dark_expected=length(pe)*settings.sensor.dark_rate_Hz*branch.src.dt
                            detector=sum(pe);raw=sum(measured)+loss;after=sum(measured);oq=sum(optical)
                            photon_total=sum(h.cells.count);sipm_loss=detector-oq
                            pde_trials+=photon_total;pde_successes+=Int64(detector)
                            values=(;event_id=id,energy_MeV=energy,electron_drift_cm=drift,
                                ion_drift_cm=h.placement.ion_drift_cm,survival,photons_total=photon_total,
                                detected_total=detector,optical_charge_total=oq,sensor_total=sum(sensor_pe),
                                unclipped_total=raw,readout_total=after,clip_loss=loss,sipm_loss,
                                sipm_loss_percent=100sipm_loss/detector,clip_loss_percent=100loss/sum(sensor_pe),
                                input_cells=length(selected),clipped_input_cells=count(identity,sat[selected]),
                                clipped_acquired_cells=count(identity,sat),acquired_cells=length(pe),dark_expected,
                                event_mean_photons=mean(h.cells.count),event_max_photons=maximum(h.cells.count),
                                event_mean_pe=mean(pe[selected]),event_max_pe=maximum(pe[selected]),
                                event_mean_sensor_pe=mean(optical[selected]),event_max_sensor_pe=maximum(optical[selected]),
                                event_mean_after_pe=mean(measured[selected]),event_max_after_pe=maximum(measured[selected]),
                                event_mean_occupancy_percent=mean(100distinct[selected]/M),
                                event_max_occupancy_percent=maximum(100distinct[selected]/M),
                                event_max_expected_distinct_cells=maximum(distinct[selected]),
                                detected_corrected=detector/survival/energy,
                                sensor_corrected=oq/survival/energy,
                                unclipped_corrected=(raw-dark_expected)/survival/energy,
                                readout_corrected=(after-dark_expected)/survival/energy)
                            for (key,value) in pairs(values)
                                push!(get!(columns,string(key),Float64[]),Float64(value))
                            end
                        end
                      end
                    end
                end
              end
            end
            for (key,value) in columns;g[key]=value;end
            if store_inclusive
                for (key,value) in inclusive;g["inclusive_"*key,compress=4]=value;end
            end
            quiet_mean=quiet_sum/quiet_n;quiet_rms=sqrt(max(0.,quiet_sumsq/quiet_n-quiet_mean^2))
            nevents=length(columns["event_id"])
            put(g,(;event_count=nevents,ledger_checks=checks,state_repairs=repairs,variance_limits=limits,
                maximum_state_repair=maxrepair,statistical_batches=batches,
                maximum_state_repair_fraction=maxrepair/min(M,settings.sensor.microcells-M),
                clipped_pulse_percent=100sum(columns["clipped_input_cells"])/sum(columns["input_cells"]),
                clipped_event_percent=100count(>(0),columns["clipped_acquired_cells"])/nevents,
                asic_charge_lost_percent=100sum(columns["clip_loss"])/sum(columns["sensor_total"]),
                sipm_charge_lost_percent=100sum(columns["sipm_loss"])/sum(columns["detected_total"]),
                event_mean_average_pe=mean(columns["event_mean_pe"]),event_max_average_pe=mean(columns["event_max_pe"]),
                detected_mean_total_pe=mean(columns["detected_total"]),
                largest_detected_pulse=maximum(columns["event_max_pe"]),
                largest_sensor_pulse=maximum(columns["event_max_sensor_pe"]),
                largest_expected_occupancy_percent=maximum(columns["event_max_occupancy_percent"]),
                corrected_detected_rms_percent=rms_percent(columns["detected_corrected"]),
                corrected_sensor_rms_percent=rms_percent(columns["sensor_corrected"]),
                corrected_unclipped_rms_percent=rms_percent(columns["unclipped_corrected"]),
                corrected_readout_rms_percent=rms_percent(columns["readout_corrected"]),
                drift_sipm_loss_correlation=cor(columns["electron_drift_cm"],columns["sipm_loss_percent"]),
                quiet_noise_mean_pe=quiet_mean,quiet_noise_rms_pe=quiet_rms,quiet_cells=quiet_n,
                pde_fraction=pde_successes/pde_trials,
                pde_mean_z=(pde_successes-pde_trials*settings.sensor.pde)/sqrt(pde_trials*settings.sensor.pde*(1-settings.sensor.pde))))
            for key in ("event_mean_photons","event_max_photons","event_mean_pe","event_max_pe","event_mean_sensor_pe","event_max_sensor_pe",
                        "event_mean_after_pe","event_max_after_pe","acquired_cells","input_cells")
                put(g,Dict("average_"*key=>mean(columns[key]),"median_"*key=>median(columns[key])))
            end
            retained=columns["clipped_acquired_cells"].==0
            put(g,(;retained_events=count(retained),retained_fraction=mean(retained)))
            for stage in ("detected","sensor","unclipped","readout")
                g["retained_$(stage)_rms_percent"]=rms_percent(columns["$(stage)_corrected"][retained])
            end
            g["retained_event_mask"]=Vector{Bool}(retained)
            @assert abs(read(g["pde_mean_z"]))<6 && limits==0
            # Repairs affect only S=sum(r_i^2), not sampled charge or R. For
            # one subsequent trigger Var(q/q0)=S/M-(R/M)^2; hence |delta S|/M
            # bounds its variance change, using the smaller population here.
            repair_rms=sqrt(read(g["maximum_state_repair_fraction"]))
            put(g,(;maximum_repair_single_trigger_rms_pe=repair_rms,
                maximum_repair_noise_fraction=repair_rms/settings.electronics.noise_rms_pe))
            @assert repair_rms<.01settings.electronics.noise_rms_pe && repairs/max(1,batches)<1e-6
            @assert abs(quiet_mean)<6settings.electronics.noise_rms_pe/sqrt(quiet_n)
            @assert abs(quiet_rms-settings.electronics.noise_rms_pe)<6settings.electronics.noise_rms_pe/sqrt(2quiet_n)
            println("$tag ns / $gas: SiPM loss $(read(g["sipm_charge_lost_percent"]))%; ASIC loss $(read(g["asic_charge_lost_percent"]))%; clipped events $(read(g["clipped_event_percent"]))%; corrected RMS before/SiPM/ASIC $(read(g["corrected_detected_rms_percent"])), $(read(g["corrected_sensor_rms_percent"])), $(read(g["corrected_readout_rms_percent"]))%")
          end
        end
        if reference_results!==nothing
            reference=create_group(out,"reference")
            put(reference,(;path=reference_results,sha256=PROD.sha(reference_results)))
            h5open(reference_results,"r") do old
                for tag in sort!(collect(keys(files["responses"])))
                    branch=create_group(reference,"tau_"*tag)
                    for gas in ("pure_xe","xe_he15")
                        haskey(out["tau_"*tag],gas) || continue
                        og=old["tau_$tag/$gas"];ng=out["tau_$tag/$gas"]
                        @assert read(og["event_id"])==read(ng["event_id"])
                        @assert read(og["electron_drift_cm"])==read(ng["electron_drift_cm"])
                        rg=create_group(branch,gas)
                        for key in ("event_id","electron_drift_cm","event_count","clipped_event_percent",
                                    "asic_charge_lost_percent","sipm_charge_lost_percent","event_mean_average_pe",
                                    "event_max_average_pe","largest_detected_pulse","largest_sensor_pulse",
                                    "corrected_detected_rms_percent","corrected_sensor_rms_percent",
                                    "corrected_unclipped_rms_percent","corrected_readout_rms_percent",
                                    "clip_loss_percent","sipm_loss_percent","event_mean_pe","event_max_pe")
                            rg[key]=read(og[key])
                        end
                        keep=read(og["clipped_acquired_cells"]).==0
                        put(rg,(;retained_events=count(keep),retained_fraction=mean(keep)))
                        for stage in ("detected","sensor","unclipped","readout")
                            rg["retained_$(stage)_rms_percent"]=rms_percent(read(og["$(stage)_corrected"])[keep])
                        end
                    end
                end
            end
        end
        write(meta["complete"],true)
    end
    mv(output*".tmp",output;force=true);println(output)
end
end
if abspath(PROGRAM_FILE)==@__FILE__;ELMSiPMReduction.main();end
