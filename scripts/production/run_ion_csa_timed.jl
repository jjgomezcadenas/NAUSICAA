#!/usr/bin/env julia
# Run one fresh CSA production with calibration and presentation timed separately.
entry_started=time_ns()
using TOML,HDF5,SHA
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/ion_signal/ion_csa_pipeline.jl"))
const P=IonCSAProduction
sha(path)=open(io->bytes2hex(sha256(io)),path)

"""
    main(card)

Run a single-pitch, production-mode CSA card against fresh filtered outputs.
Save source/circuit preflight hashes and covariance setup time separately from
the complete response calculation. The calculation retains all selected events
and signed buckets. The final timing summary reports seconds per event excluding
calibration setup and Python presentation, plus the complete measured wall time.
Existing filtered outputs cause an error so cache reuse cannot masquerade as a
fresh-production benchmark. The run manifest and library validate output ledgers.
"""
function main(card)
    loaded_s=(time_ns()-entry_started)/1e9;s=P.settings(ROOT,card)
    @assert length(s.branches)==1 && s.c["mode"]=="production"
    pitch=only(s.branches)["pitch_mm"];base=joinpath(s.dir,"execution","pitch$pitch")
    cp=P.write_toml(joinpath(base,"resolved_response.toml"),s.c)
    raw=P.write_toml(joinpath(base,"filtered_run.toml"),Dict("schema"=>"nausicaa.ion_csa_run.v1",
        "response_card"=>relpath(cp,ROOT),"ion_index"=>s.r["ion_index"],"samples"=>s.r["samples"],
        "pitches_mm"=>[pitch],"events_per_sample"=>s.r["events_per_sample"],"selection"=>s.r["selection"],
        "results_directory"=>relpath(joinpath(base,"filtered"),ROOT),"report_directory"=>relpath(joinpath(base,"filtered_report"),ROOT)))
    probe=P.R.settings(ROOT,raw)
    @assert all(x->!isfile(x.response),probe.sources) "Filtered response already exists; this is not a fresh production"
    caldir=P.resolve(ROOT,s.c["calibration_directory"]);variant=s.c["variant"]
    inputs=Dict("actual_deck"=>joinpath(caldir,"decks",variant,"30ms_single_N100_r2.cir"),
                "actual_model"=>joinpath(caldir,"decks",variant,"corrected.inc"),
                "kernel_waveform"=>joinpath(caldir,"variants",variant*".h5"),
                "nonlinear_calibration"=>P.resolve(ROOT,s.c["nonlinear_calibration"]))
    provenance=Dict{String,Any}("card"=>relpath(s.path,ROOT),"script_sha256"=>sha(@__FILE__),
        "production_code_hash"=>s.codehash,"fresh_filtered_outputs"=>true,
        "source_paths"=>[x.path for x in probe.sources],"source_sha256"=>[x.sourcehash for x in probe.sources])
    for (key,path) in inputs;provenance[key*"_path"]=path;provenance[key*"_sha256"]=sha(path);end
    P.write_toml(joinpath(s.dir,"preflight.toml"),provenance)
    maxframes=64
    for source in probe.sources
        h5open(source.path,"r") do h
            for id in source.selected
                buckets=read(h["events/$id/bucket"])
                maxframes=max(maxframes,maximum(buckets)-minimum(buckets)+probe.c["quiet_before_frames"]+ceil(Int,probe.c["recovery_tail_s"]/source.frame)+1)
            end
        end
    end
    before=time_ns();P.R.calibration(probe,maxframes);calibration_s=(time_ns()-before)/1e9
    println("FRESH PRODUCTION pitch=$pitch mm events=$(s.r["events_per_sample"]) covariance_setup_s=$calibration_s");flush(stdout)
    before=time_ns();manifest=P.run(ROOT,card);run_s=(time_ns()-before)/1e9
    @assert P.settings(ROOT,card).codehash==s.codehash "CSA code changed during production"
    for (key,path) in inputs;@assert sha(path)==provenance[key*"_sha256"];end
    timing=TOML.parsefile(joinpath(base,"timing.toml"));work_s=run_s-timing["presentation_s"]
    n=s.r["events_per_sample"]*length(s.r["samples"])
    summary=Dict{String,Any}("complete"=>true,"manifest"=>manifest,"events"=>n,"pitch_mm"=>pitch,
        "fresh_filtered_outputs"=>true,"module_loading_s"=>loaded_s,"calibration_setup_s"=>calibration_s,
        "maximum_frames"=>maxframes,"response_run_s"=>run_s,"event_work_s"=>work_s,
        "event_work_seconds_per_event"=>work_s/n,"complete_wall_s"=>(time_ns()-entry_started)/1e9)
    merge!(summary,timing)
    P.write_toml(joinpath(s.dir,"production_timing.toml"),summary)
    println("COMPLETE pitch=$pitch mm: $work_s s event work; $(work_s/n) s/event; total $(summary["complete_wall_s"]) s");flush(stdout)
    manifest
end
main(only(ARGS))
