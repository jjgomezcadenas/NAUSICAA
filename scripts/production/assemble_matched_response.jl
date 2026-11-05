#!/usr/bin/env julia
module MatchedResponse
using HDF5,TOML,SHA
const ROOT=normpath(joinpath(@__DIR__,"../.."))
sha(path)=bytes2hex(sha256(read(path)))
function write_toml(path,data)
    mkpath(dirname(path));open(io->TOML.print(io,data),path*".tmp","w");mv(path*".tmp",path;force=true)
end

"""Validate shared identities/placement across six histories and publish one graph.

The common realization, ion collection, photon collection, corrected CSA,
SiPM charge and ASIC output must agree exactly on event ID, pair checksum and
both drift distances. The viewer card uses existing response file schemas.
"""
function main(kind,sipm_manifest)
    kind in ("smoke","full") || error("Choose smoke or full")
    out=joinpath(ROOT,"results",kind=="smoke" ? "elm_matched_smoke" : "elm_matched_production")
    detector=TOML.parsefile(joinpath(out,"detector_files.toml"))
    ion_manifest=TOML.parsefile(joinpath(out,"ion_csa/latest.toml"))["manifest"]
    ionrun=TOML.parsefile(ion_manifest);sipm=TOML.parsefile(sipm_manifest)
    @assert ionrun["complete"] && sipm["complete"] && length(ionrun["responses"])==1
    ions=ionrun["responses"][1]["path"]
    ph=h5open(detector["detector_index"],"r") do f;String(read(f["bb0nu/xe_he15_path"]));end
    paths=[detector["common_realisation"],detector["ions"],ph,ions]
    handles=HDF5.File[];ids=Int64[];li=Float32[];le=Float32[];pairs=Int64[];checksums=String[]
    try
        append!(handles,[h5open(p,"r") for p in paths])
        canonical=sort!(parse.(Int64,collect(keys(handles[1]["events"]))))
        @assert all(sort!(parse.(Int64,collect(keys(h["events"]))))==canonical for h in handles)
        for shard in sipm["shards"]
            @assert shard["gas"]=="xe_he15" && shard["sample"]=="bb0nu"
            h5open(shard["sensor"],"r") do sensor
              h5open(shard["asic"],"r") do asic
                for id in shard["event_ids"]
                    groups=[h["events/$id"] for h in vcat(handles,[sensor,asic])]
                    for key in ("event_id","pair_checksum","ion_drift_cm","electron_drift_cm")
                        values=[read(g[key]) for g in groups]
                        all(==(first(values)),values) || error("Matched response mismatch: event $id field $key")
                    end
                    a=Float32(read(groups[1]["ion_drift_cm"]));b=Float32(read(groups[1]["electron_drift_cm"]))
                    @assert 5<=a<=145 && 5<=b<=145 && isapprox(a+b,150;atol=1f-4)
                    push!(ids,id);push!(li,a);push!(le,b);push!(pairs,read(groups[1]["initial_pairs"]))
                    push!(checksums,String(read(groups[1]["pair_checksum"])))
                end
              end
            end
        end
        @assert sort(ids)==canonical && allunique(ids)
    finally
        foreach(close,handles)
    end
    validation=joinpath(out,"matching_validation.h5")
    h5open(validation*".tmp","w") do f
        attributes(f)["schema"]="nausicaa.matched_response_validation.v1"
        f["complete"]=true;f["event_id"]=ids;f["ion_drift_cm"]=li;f["electron_drift_cm"]=le
        f["initial_pairs"]=pairs;f["pair_checksum"]=checksums;f["validated_branches"]=6
    end
    mv(validation*".tmp",validation;force=true)
    manifest=joinpath(out,"manifest.toml")
    viewer=joinpath(ROOT,"config/specs/runs",kind=="smoke" ? "elm_matched_smoke_viewer.toml" : "elm_matched_response_viewer.toml")
    write_toml(manifest,Dict("schema"=>"nausicaa.matched_response_run.v1","complete"=>true,
        "sample"=>"bb0nu","event_count"=>length(ids),"gases"=>["xe_he15"],
        "common_realisation"=>paths[1],"common_sha256"=>sha(paths[1]),
        "detector_index"=>detector["detector_index"],"ion_collection"=>paths[2],"photon_collection"=>ph,
        "ion_manifest"=>ion_manifest,"ion_manifest_sha256"=>sha(ion_manifest),
        "ion_response"=>ions,"sipm_manifest"=>sipm_manifest,"sipm_manifest_sha256"=>sha(sipm_manifest),
        "matching_validation"=>validation,"matching_validation_sha256"=>sha(validation),
        "viewer_card"=>viewer,"chamber_length_cm"=>150.,"minimum_drift_cm"=>5.,"maximum_drift_cm"=>145.))
    write_toml(viewer,Dict("sample"=>"bb0nu","sipm_manifest"=>relpath(sipm_manifest,ROOT),
        "production_manifest"=>relpath(manifest,ROOT),"label"=>"ELM · E/P 2 · 4 mm · 5–145 cm",
        "default_gas"=>"xe_he","default_pitch_mm"=>4,"chamber_length_cm"=>150.,
        "id"=>kind=="smoke" ? "elm_matched_smoke" : "elm_matched",
        "ions"=>[Dict("pitch_mm"=>4,"path"=>relpath(ions,ROOT))]))
    println("Validated $(length(ids)) events across six matched histories")
    println("Combined manifest: $manifest");println("Viewer card: $viewer");manifest
end
end
if abspath(PROGRAM_FILE)==@__FILE__;MatchedResponse.main(ARGS[1],abspath(ARGS[2]));end
