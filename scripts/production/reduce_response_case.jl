# Join completed ion CSA and SiPM response branches into numerical results.
using HDF5, TOML, SHA, Statistics
const ROOT=normpath(joinpath(@__DIR__,"../.."))
resolve(path)=normpath(isabspath(path) ? path : joinpath(ROOT,path))
sha(path)=open(io->bytes2hex(sha256(io)),path)
str(g,key)=String(read(g[key]))
require(ok,message)=ok || throw(ArgumentError(message))

function corrected_summary(ionmanifest,parent)
    root=joinpath(dirname(ionmanifest),"execution/pitch3/corrected")
    require(isdir(root),"Corrected CSA study directory is absent")
    matches=String[]
    for dir in readdir(root;join=true)
        path=joinpath(dir,"summary.h5");isfile(path) || continue
        h5open(path,"r") do h
            haskey(h,"bb0nu_p3") && String(read(attributes(h["bb0nu_p3"])["response_path"]))==parent &&
                push!(matches,path)
        end
    end
    require(length(matches)==1,"Expected one matching corrected CSA summary")
    only(matches)
end

"""Count event values in fixed bins, including the upper endpoint."""
function histogram(values,edges)
    counts=zeros(Int64,length(edges)-1)
    for x in values
        i=clamp(searchsortedlast(edges,x),1,length(counts))
        counts[i]+=1
    end
    require(sum(counts)==length(values),"Histogram ledger")
    counts
end

function join_case(ionmanifest,sipmmanifest,output)
    ion=TOML.parsefile(ionmanifest);sipm=TOML.parsefile(sipmmanifest)
    require(ion["schema"]=="nausicaa.ion_response_run.v1" && ion["complete"],"Ion manifest")
    require(sipm["schema"]=="nausicaa.sipm_production_manifest.v1" && sipm["complete"],"SiPM manifest")
    responses=ion["responses"]
    require(length(responses)==1 && responses[1]["sample"]=="bb0nu" &&
            responses[1]["pitch_mm"]==3,"Expected bb0nu 3 mm ion branch")
    response=only(responses);ionpath=response["path"]
    require(isfile(ionpath) && sha(ionpath)==response["sha256"],"Ion response content")
    csummary=corrected_summary(ionmanifest,response["parent_corrected_response"])
    ssumm=sipm["summary"]
    require(isfile(ssumm) && sha(ssumm)==sipm["summary_sha256"],"SiPM summary content")
    ionrun=TOML.parsefile(resolve(ion["run_card"]))
    sipmrun=TOML.parsefile(resolve(sipm["run_card"]))
    pilotrun=TOML.parsefile(resolve(sipmrun["pilot_card"]))
    mkpath(dirname(output));tmp=output*".tmp";isfile(tmp)&&rm(tmp)
    h5open(tmp,"w") do out
        attributes(out)["schema"]="nausicaa.response_case_analysis.v1"
        m=create_group(out,"metadata")
        m["complete"]=false;m["ion_manifest"]=abspath(ionmanifest)
        m["ion_manifest_sha256"]=sha(ionmanifest)
        m["sipm_manifest"]=abspath(sipmmanifest)
        m["sipm_manifest_sha256"]=sha(sipmmanifest)
        m["ion_response"]=ionpath;m["ion_response_sha256"]=response["sha256"]
        m["sipm_summary"]=ssumm;m["sipm_summary_sha256"]=sha(ssumm)
        m["sample"]="bb0nu";m["ion_pitch_mm"]=3.0;m["ion_frame_ms"]=30.0
        m["electron_bucket_us"]=0.25;m["reducer_sha256"]=sha(@__FILE__)
        m["ion_card_toml"]=read(resolve(ionrun["response_card"]),String)
        m["ion_reconstruction_card_toml"]=read(resolve(ionrun["reconstruction_card"]),String)
        m["sipm_card_toml"]=read(resolve(pilotrun["sipm_card"]),String)
        m["asic_card_toml"]=read(resolve(pilotrun["electronics_card"]),String)
        iongroup=create_group(out,"ions")
        h5open(csummary,"r") do h
            g=h["bb0nu_p3"];a=attributes(g)
            ids=Int64.(read(g["event_id"]))
            require(length(ids)==8000 && ids==sort!(unique(ids)),"Ion event set")
            iongroup["event_id"]=ids
            for key in ("initial_ions","event_mean_snr","noiseless_charge_bias_fraction",
                        "noisy_charge_bias_fraction","compression_deficit_fraction",
                        "negative_estimates_per_event","pads_per_event","positive_cells")
                iongroup[key]=read(g[key])
            end
            iongroup["snr_edges"]=read(g["snr/edges"])
            iongroup["snr_counts"]=read(g["snr/count"])
            iongroup["noisy_bias_edges"]=read(g["noisy_bias/edges"])
            iongroup["noisy_bias_counts"]=read(g["noisy_bias/count"])
            iongroup["enc_ions_rms"]=Float64(read(a["ENC_e_rms"]))
            iongroup["mean_event_snr"]=mean(read(g["event_mean_snr"]))
            iongroup["mean_noiseless_bias_fraction"]=mean(read(g["noiseless_charge_bias_fraction"]))
            iongroup["mean_noisy_bias_fraction"]=mean(read(g["noisy_charge_bias_fraction"]))
            iongroup["compression50_event_fraction"]=Float64(read(g["fraction_events_with_compression"])[3])
            iongroup["input_ion_total"]=sum(Int64.(read(g["initial_ions"])))
            require(read(iongroup["input_ion_total"])==response["input_ion_total"],"Ion total ledger")
        end
        pair_ids=Dict{Int64,Tuple{String,Int64}}()
        h5open(ionpath,"r") do h
            for id in read(iongroup["event_id"])
                g=h["events/$id"]
                pair_ids[id]=(str(g,"pair_checksum"),Int64(read(g["initial_pairs"])))
            end
        end
        h5open(ssumm,"r") do sh
            for gas in ("pure_xe","xe_he15")
                src=sh["bb0nu/$gas"];ids=Int64.(read(src["event_id"]))
                require(ids==read(iongroup["event_id"]),"Ion/SiPM event ID mismatch: $gas")
                group=create_group(out,"photons/$gas")
                for key in ("event_id","incident_photons","detected_pe","dark_pulses",
                            "sensor_pe_equiv","readout_pe_equiv","clip_loss_pe_equiv",
                            "clipped_cells","acquired_cells")
                    group[key]=read(src[key])
                end
                clipped=Int64.(read(src["clipped_cells"]))
                acquired=Int64.(read(src["acquired_cells"]))
                sensor=Float64.(read(src["sensor_pe_equiv"]))
                detected=Int64.(read(src["detected_pe"]))
                dark=Int64.(read(src["dark_pulses"]))
                loss=Float64.(read(src["clip_loss_pe_equiv"]))
                require(all(acquired.>0) && all(clipped.<=acquired) &&
                        all(detected.<=read(src["incident_photons"])),"SiPM event ledgers")
                recovery=1 .-sensor./(detected.+dark)
                clip_fraction=clipped./acquired
                group["recovery_loss_fraction"]=recovery
                group["clipped_cell_fraction"]=clip_fraction
                group["incident_total"]=sum(Int64.(read(src["incident_photons"])))
                group["detected_total"]=sum(detected)
                group["dark_total"]=sum(dark)
                group["sensor_pe_total"]=sum(sensor)
                group["readout_pe_total"]=sum(Float64.(read(src["readout_pe_equiv"])))
                group["clipped_cell_total"]=sum(clipped)
                group["acquired_cell_total"]=sum(acquired)
                group["clipped_cell_fraction_total"]=sum(clipped)/sum(acquired)
                group["events_with_clipping_fraction"]=count(>(0),clipped)/length(ids)
                group["clip_loss_pe_total"]=sum(loss)
                group["clip_loss_per_sensor_pe_fraction"]=sum(loss)/sum(sensor)
                group["mean_recovery_loss_fraction"]=mean(recovery)
                for (name,values) in (("recovery",recovery),("clipped",clip_fraction))
                    edges=collect(range(0.0,stop=max(maximum(values)*1.01,1e-6),length=81))
                    group["$(name)_edges"]=edges
                    group["$(name)_counts"]=histogram(values,edges)
                end
            end
        end
        # Compare all gas-specific pair identities with the selected ion branch.
        for gas in ("pure_xe","xe_he15")
            entries=filter(e->e["sample"]=="bb0nu" && e["gas"]==gas,sipm["shards"])
            ids=Int64[]
            photon_index=h5open(ssumm,"r") do h
                str(h["metadata"],"photon_index")
            end
            sourcepath=h5open(photon_index,"r") do h
                str(h["bb0nu"],"$(gas)_path")
            end
            h5open(sourcepath,"r") do photons
                for entry in entries
                    path=entry["asic"]
                    require(isfile(path) && sha(path)==entry["asic_sha256"],"ASIC shard content")
                    h5open(path,"r") do h
                        require(Bool(read(h["metadata/complete"])),"ASIC shard incomplete")
                        for id in entry["event_ids"]
                            g=h["events/$id"]
                            require(str(g,"pair_checksum")==pair_ids[id][1] &&
                                    Int64(read(photons["events/$id/initial_pairs"]))==pair_ids[id][2],
                                    "Ion/SiPM primary-pair mismatch")
                            push!(ids,id)
                        end
                    end
                end
            end
            require(ids==read(iongroup["event_id"]),"ASIC shard event coverage: $gas")
        end
        write(m["complete"],true)
    end
    mv(tmp,output;force=true)
    output
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==3 || error("Usage: reduce_response_case.jl ION_MANIFEST SIPM_MANIFEST OUTPUT.h5")
    println(join_case(ARGS[1],ARGS[2],ARGS[3]))
end
