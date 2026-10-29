#!/usr/bin/env julia
# Reduce six completed grid histories into per-event and ensemble HDF5 summaries.
using HDF5, TOML
const ROOT=normpath(joinpath(@__DIR__,".."))
const STREAMS=[(sample,carrier) for sample in ("bb0nu","electrons")
               for carrier in ("ion","pure_xe","xe_he")]

function summarize(directory::String)
    summary=joinpath(directory,"grid_readout_summary.h5")
    isfile(summary) || error("Missing parallel-run summary: $summary")
    h5open(summary,"r+") do target
        read(attributes(target)["schema"])=="nausicaa.direct_grid_readout_summary.v2" ||
            error("Expected direct-grid run summary v2")
        n=Int(read(target["events_per_stream"]))
        haskey(target,"stream_event_totals") && HDF5.delete_object(target,"stream_event_totals")
        haskey(target,"summary_complete") && HDF5.delete_object(target,"summary_complete")
        for (sample,carrier) in STREAMS
            path=joinpath(directory,"$(sample)_$(carrier).h5")
            h5open(path,"r") do file
                read(attributes(file)["schema"])=="nausicaa.direct_grid_readout.v2" ||
                    error("Expected direct-grid stream v2: $path")
                read(file["metadata/complete"])==true || error("Incomplete grid stream: $path")
                Int(read(file["metadata/event_count"]))==n || error("Event count mismatch: $path")
                if sample=="bb0nu" && carrier=="ion"
                    for key in ("electron_card_toml","ion_card_toml",
                                "detector_card_toml","code_paths","code_sha256")
                        haskey(target,key) && HDF5.delete_object(target,key)
                        target[key]=read(file["metadata/"*key])
                    end
                else
                    for key in ("electron_card_toml","ion_card_toml",
                                "detector_card_toml","code_paths","code_sha256")
                        read(file["metadata/"*key])==read(target[key]) ||
                            error("Input provenance differs between grid streams: $key")
                    end
                end
                carriers=Vector{Int64}(undef,n)
                sampled=Vector{Int64}(undef,n)
                outside=Vector{Int64}(undef,n)
                signals=Vector{Int64}(undef,n)
                pulses=Vector{Int64}(undef,n)
                for index in 1:n
                    id=index-1
                    g=file["events/"*string(id)]
                    carriers[index]=read(g["input_carriers"])
                    sampled[index]=read(g["sampled_carriers"])
                    outside[index]=read(g["outside_carriers"])
                    sampled[index]==carriers[index]+outside[index] ||
                        error("Diffusion ledger mismatch at $sample/$carrier event $id")
                    signals[index]=read(g["total_signal"])
                    pulses[index]=length(g["signal_count"])
                    if carrier=="ion"
                        signals[index]==carriers[index] || error("Ion loss at $sample event $id")
                    end
                end
                base="stream_event_totals/$(sample)_$(carrier)/"
                target[base*"event_id"]=Int64.(0:n-1)
                target[base*"input_carriers"]=carriers
                target[base*"sampled_carriers"]=sampled
                target[base*"outside_carriers"]=outside
                target[base*"signal_count"]=signals
                target[base*"occupied_bucket_count"]=pulses
                target[base*"sum_input_carriers"]=sum(carriers)
                target[base*"sum_sampled_carriers"]=sum(sampled)
                target[base*"sum_outside_carriers"]=sum(outside)
                target[base*"sum_signal_count"]=sum(signals)
                target[base*"sum_occupied_bucket_count"]=sum(pulses)
                target[base*"mean_signal_per_input_carrier"]=
                    Float64(sum(signals))/Float64(sum(carriers))
                target[base*"pitch_mm"]=read(file["metadata/pitch_mm"])
                target[base*"bucket_mm"]=read(file["metadata/bucket_mm"])
                target[base*"signal_unit"]=read(file["metadata/signal_unit"])
                target[base*"time_unit"]=read(file["metadata/time_unit"])
                target[base*"source_sha256"]=read(file["metadata/source_sha256"])
                target[base*"el_photons_per_electron"]=read(file["metadata/el_photons_per_electron"])
                target[base*"pe_probability_per_photon"]=read(file["metadata/pe_probability_per_photon"])
                target[base*"sipm_pde"]=Float32(TOML.parse(read(file["metadata/electron_card_toml"]))["readout"]["sipm_pde"])
                target[base*"axial_acceptance"]=carrier=="ion" ? 0f0 :
                    Float32(read(file["metadata/pe_probability_per_photon"]))/
                    read(target[base*"sipm_pde"])
                println("$sample/$carrier: ",sum(signals)," ",read(file["metadata/signal_unit"]),
                        "; ",sum(pulses)," occupied buckets")
            end
        end
        target["summary_complete"]=true
    end
    return summary
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("Expected grid readout directory")
    summarize(abspath(ARGS[1]))
end
