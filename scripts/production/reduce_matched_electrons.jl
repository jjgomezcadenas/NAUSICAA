#!/usr/bin/env julia
module MatchedElectronAnalysis
using HDF5,TOML,SHA
include("reduce_elm_sipm_pilot.jl")
const ROOT=ELMSiPMReduction.ROOT
sha(path)=bytes2hex(sha256(read(path)))

"""Reduce or validate cached electron observables independently of the ion arm.

The cache binds the completed SiPM manifest, detector index and both reduction
sources by content hash. Its timing describes detector plus electron processing;
the paired report reads the separate full-workflow timing record.
"""
function main(kind="full")
    out=joinpath(ROOT,"results",kind=="smoke" ? "elm_matched_smoke" : "elm_matched_production")
    detector=TOML.parsefile(joinpath(out,"detector_files.toml"))
    response=TOML.parsefile(joinpath(out,"sipm_files.toml"))
    output=joinpath(out,"electron_analysis.h5")
    manifest_sha=sha(response["manifest"]);index_sha=sha(detector["detector_index"])
    reducer_sha=sha(joinpath(@__DIR__,"reduce_elm_sipm_pilot.jl"));wrapper_sha=sha(@__FILE__)
    if isfile(output)
        valid=h5open(output,"r") do f
            m=f["metadata"]
            all(haskey(m,k) for k in ("source_manifest_sha256","input_index_sha256","wrapper_sha256")) &&
                read(m["complete"]) && read(m["source_manifest_sha256"])==manifest_sha &&
                read(m["input_index_sha256"])==index_sha && read(m["wrapper_sha256"])==wrapper_sha &&
                read(m["reduction_sha256"])==reducer_sha
        end
        valid && (println("Reused validated electron analysis: $output");return output)
    end
    files=Dict("complete"=>true,"photon_index"=>detector["detector_index"],
        "photon_wall_s"=>detector["detector_wall_s"],"response_wall_s"=>response["response_wall_s"],
        "total_wall_s"=>detector["detector_wall_s"]+response["response_wall_s"],
        "responses"=>Dict("12p5"=>Dict("manifest"=>response["manifest"],
            "card"=>"config/specs/runs/elm_matched_sipm_$kind.toml","elapsed_s"=>response["response_wall_s"])))
    open(io->TOML.print(io,files),joinpath(out,"files.toml"),"w")
    ELMSiPMReduction.main(output_dir=out,report_name="electron_analysis",store_inclusive=false)
    h5open(output,"r+") do f
        m=f["metadata"];m["source_manifest_sha256"]=manifest_sha;m["input_index_sha256"]=index_sha
        m["wrapper_sha256"]=wrapper_sha
    end
    output
end
end
if abspath(PROGRAM_FILE)==@__FILE__;MatchedElectronAnalysis.main(isempty(ARGS) ? "full" : only(ARGS));end
