#!/usr/bin/env julia

using TOML,HDF5
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/ion_signal/ion_csa_pipeline.jl"))
const P=IonCSAProduction
"""Finish a byte-identical interrupted publication copy and register the saved study.

Usage: julia --project=. scripts/production/recover_ion_csa_publication.jl ORIGINAL_CARD REGISTRATION_CARD

The original corrected HDF5 and interrupted temporary copy must have identical
SHA256 checksums. Numerical arrays remain unchanged. The production metadata and
validation follow IonCSAProduction.publish; registration uses the existing public
pipeline with a saved corrected_study. Record recovery time separately from the
fresh numerical run, whose branch timing survives the publication failure.
"""
function main(original_card,registration_card)
    began=time_ns();s=P.settings(ROOT,original_card)
    @assert length(s.branches)==1
    pre=TOML.parsefile(joinpath(s.dir,"preflight.toml"))
    @assert s.codehash==pre["production_code_hash"]
    for key in ("actual_deck","actual_model","kernel_waveform","nonlinear_calibration")
        @assert P.sha(pre[key*"_path"])==pre[key*"_sha256"]
    end
    pitch=only(s.branches)["pitch_mm"]
    execution=joinpath(s.dir,"execution","pitch$pitch")
    choices=filter(p->isfile(joinpath(p,"summary.h5")),readdir(joinpath(execution,"corrected");join=true))
    study=only(choices);parent=only(P.study_responses(s,study))
    ledger=P.validate_response(parent);P.check_selection(s,parent,ledger,pitch)
    parenthash=P.sha(parent)
    nodeid=P.S.content_hash(parenthash*s.codehash*P.CHANNEL_POLICY*P.ABSENT_POLICY*P.PRODUCTION_SCHEMA)
    target=joinpath(dirname(dirname(parent)),nodeid,"ion_response.h5");tmp=target*".tmp"
    @assert !isfile(target) && isfile(tmp)
    P.GuardedFiles.check_size(tmp,filesize(parent))
    @assert P.sha(tmp)==parenthash "Interrupted copy differs from corrected source"
    # The original temporary is verified evidence; publication uses the bounded copier.
    P.publish(s,parent,ledger)
    rm(tmp)
    P.write_toml(joinpath(dirname(target),"production_manifest.toml"),Dict("schema"=>"nausicaa.ion_response_node.v1",
        "response_id"=>nodeid,"response_file"=>target,"sha256"=>P.sha(target),"parent_response"=>parent,
        "parent_sha256"=>parenthash,"channel_selection"=>P.CHANNEL_POLICY,"absent_channel_policy"=>P.ABSENT_POLICY,"complete"=>true))
    reg=TOML.parsefile(P.resolve(ROOT,registration_card))
    @assert only(reg["branches"])["corrected_study"]==relpath(study,ROOT)
    manifest=P.run(ROOT,registration_card)
    elapsed=(time_ns()-began)/1e9
    branch=TOML.parsefile(joinpath(execution,"timing.toml"))
    log=joinpath(ROOT,"work","matched_pitch$(pitch)_csa.log")
    calibration_s=parse(Float64,only(match(r"covariance_setup_s=([0-9.]+)",read(log,String)).captures))
    summary=copy(branch)
    merge!(summary,Dict("complete"=>true,"manifest"=>manifest,"events"=>ledger.event_count,"pitch_mm"=>pitch,
        "fresh_filtered_outputs"=>pre["fresh_filtered_outputs"],"publication_recovered"=>true,
        "recovery_wall_s"=>elapsed,"event_work_s"=>branch["numerical_s"],
        "event_work_seconds_per_event"=>branch["numerical_s"]/ledger.event_count,
        "execution_directory"=>execution,"original_numerical_run"=>s.dir,
        "recovery_script_sha256"=>P.sha(@__FILE__),"calibration_setup_s"=>calibration_s,
        "calibration_timing_source"=>relpath(log,ROOT)*" FRESH PRODUCTION line"))
    P.write_toml(joinpath(dirname(manifest),"production_timing.toml"),summary)
    println("Recovered publication: $manifest\nRecovery: $elapsed s\nFresh numerical work: $(summary["event_work_s"]) s; $(summary["event_work_seconds_per_event"]) s/event")
end
main(ARGS...)
