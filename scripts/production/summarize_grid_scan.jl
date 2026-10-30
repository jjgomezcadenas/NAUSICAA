#!/usr/bin/env julia
# Collect the completed 2/3/4 mm ion readouts and shared electron histories.
include(joinpath(@__DIR__,"run_ion_readout.jl"))

"""Read one HDF5 scalar as an owning Julia string."""
saved_string(group,name)=String(read(group[name]))

"""Verify one ion-only readout against the full-sample reference."""
function check_ion_file(path::String,expected_id::String,source_sha::String,
        stream_identity::String,n::Int,expected_pairs::Int64)
    return h5open(path,"r") do file
        String(read(attributes(file)["schema"]))=="nausicaa.ion_pad_file.pilot.v1" ||
            error("Ion file schema mismatch: $path")
        m=file["metadata"]
        Bool(read(m["complete"])) && Int(read(m["event_count"]))==n ||
            error("Ion file incomplete: $path")
        saved_string(m,"card_hash")==expected_id &&
            saved_string(m,"source_sha256")==source_sha &&
            saved_string(m,"stream_identity")==stream_identity ||
            error("Ion file provenance mismatch: $path")
        pairs=Int64(read(m["total_pairs"]));ions=Int64(read(m["total_ions"]))
        pairs==expected_pairs && ions==pairs || error("Ion charge mismatch: $path")
        return (;pairs,ions,elapsed=Float64(read(m["elapsed_s"])),bytes=filesize(path))
    end
end

"""Write one numerical HDF5 index for two species, six ion grids and four
gas-specific photon files. The full 3 mm summary provides electron totals;
each ion-only file has already verified every event's shared pair checksum.
"""
function summarize(card2::String,card3::String,card4::String,full_summary::String)
    specs=Dict(p=>load_specs(ROOT,card) for (p,card) in
        ((2,card2),(3,card3),(4,card4)))
    for p in (2,4)
        canonical(specs[p].realisation.data)==canonical(specs[3].realisation.data) ||
            error("The pair realisation differs at pitch $p")
        canonical(specs[p].electrons.data)==canonical(specs[3].electrons.data) ||
            error("Electron settings differ at pitch $p")
    end
    ion_hash=ion_replay_code_hash()
    h5open(full_summary,"r") do full
        String(read(attributes(full)["schema"]))=="nausicaa.detector_pilot_summary.v1" ||
            error("Full summary schema mismatch")
        Bool(read(full["metadata/complete"])) || error("Full summary incomplete")
        canonical(TOML.parse(saved_string(full["metadata"],"realisation_card")))==
            canonical(specs[3].realisation.data) || error("Full realisation card differs")
        canonical(TOML.parse(saved_string(full["metadata"],"ion_card")))==
            canonical(specs[3].ions.data) || error("Full 3 mm ion card differs")
        canonical(TOML.parse(saved_string(full["metadata"],"electron_card")))==
            canonical(specs[3].electrons.data) || error("Full electron card differs")
        run_id=saved_string(full["metadata"],"run_id")
        output=joinpath(specs[3].results_dir,"grid_scan_summary_$run_id.h5")
        temporary=output*".tmp"
        isfile(temporary) && rm(temporary)
        h5open(temporary,"w") do out
            attributes(out)["schema"]="nausicaa.detector_grid_scan.v1"
            m=create_group(out,"metadata")
            m["complete"]=false;m["run_id"]=run_id
            m["full_summary"]=abspath(full_summary)
            m["full_summary_sha256"]=content_hash(full_summary;file=true)
            m["realisation_card"]=read(specs[3].realisation.path,String)
            m["electron_card"]=read(specs[3].electrons.path,String)
            m["ion_card_2mm"]=read(specs[2].ions.path,String)
            m["ion_card_3mm"]=read(specs[3].ions.path,String)
            m["ion_card_4mm"]=read(specs[4].ions.path,String)
            for sample in ("bb0nu","single_electrons")
                reference=full[sample]
                n=Int(read(reference["events"]))
                n==8000 || error("Expected 8000 events for $sample")
                Bool(read(reference["ion_charge_conserved"])) &&
                    Bool(read(reference["photon_ledgers_valid"])) ||
                    error("Full branch ledger failed for $sample")
                pairs=Int64(read(reference["pairs"]))
                group=create_group(out,sample)
                group["events"]=n;group["pairs"]=pairs
                group["survivors"]=Int64(read(reference["survivors"]))
                group["source_sha256"]=saved_string(reference,"source_sha256")
                group["realisation_id"]=saved_string(reference,"realisation_id")
                group["realisation_path"]=saved_string(reference,"realisation_path")
                group["joint_run_elapsed_s"]=Float64(read(reference["elapsed_s"]))
                ions_group=create_group(group,"ions")
                for pitch in (2,3,4)
                    card=specs[pitch].ions.data
                    pitch==Int(card["pad_pitch_mm"]) || error("Pitch card mismatch")
                    path=pitch==3 ? saved_string(reference,"ions_path") :
                        ion_paths(specs[pitch],sample,ion_hash).ions
                    isfile(path) || error("Missing ion output: $path")
                    result=if pitch==3
                        (;pairs,ions=Int64(read(reference["ions"])),
                          elapsed=NaN,bytes=filesize(path))
                    else
                        resolved=ion_paths(specs[pitch],sample,ion_hash)
                        resolved.realisation_id==saved_string(reference,"realisation_id") ||
                            error("Ion grid uses a different realisation")
                        check_ion_file(path,resolved.ion_id,
                            saved_string(reference,"source_sha256"),
                            saved_string(reference,"stream_identity"),n,pairs)
                    end
                    result.ions==pairs || error("Ion charge mismatch at pitch $pitch")
                    g=create_group(ions_group,string(pitch))
                    g["path"]=path;g["pad_pitch_mm"]=card["pad_pitch_mm"]
                    g["bucket_width_ms"]=card["bucket_width_ms"]
                    g["ions"]=result.ions;g["file_bytes"]=result.bytes
                    g["elapsed_s"]=result.elapsed
                end
                photons=create_group(group,"photons")
                for (gas,key,count_key) in (("pure_xe","pure_xe_path","pure_xe_incident"),
                                             ("xe_he15","xe_he15_path","xe_he15_incident"))
                    path=saved_string(reference,key)
                    isfile(path) || error("Missing photon output: $path")
                    g=create_group(photons,gas)
                    g["path"]=path;g["incident_photons"]=Int64(read(reference[count_key]))
                    g["file_bytes"]=filesize(path)
                end
            end
            write(m["complete"],true)
        end
        mv(temporary,output;force=true)
        println(output)
        return output
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==4 || error("Usage: summarize_grid_scan.jl PITCH2.toml PITCH3.toml PITCH4.toml FULL_SUMMARY.h5")
    summarize(ARGS...)
end
