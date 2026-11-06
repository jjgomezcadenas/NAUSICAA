#!/usr/bin/env julia
# Collect a new ion grid while preserving the completed matched realization.
module MatchedIonPitch
using HDF5,TOML,SHA
module Ions
include("run_ion_readout.jl")
end
const ROOT=Ions.ROOT
sha(path)=open(io->bytes2hex(sha256(io)),path)

"""
    main(pitch)

Replay the existing 8,000-event bb0nu common realization into a 2 or 3 mm
ion grid, with the corresponding 20 or 30 ms frame. The generation vertex,
primary-pair checksums, shared diffusion stream and complementary 5–145 cm
drift placement are retained. Save a CSA-compatible index after checking
every event against the common realization and completed 4 mm ion readout.
The electron files and common ledger are read only.
"""
function main(pitch::Int)
    pitch in (2,3) || error("Choose pitch 2 or 3 mm")
    card="config/specs/runs/elm_matched_ions_pitch$pitch.toml"
    spec=Ions.load_specs(ROOT,card);code=Ions.ion_replay_code_hash()
    paths=Ions.ion_paths(spec,"bb0nu",code)
    out=joinpath(ROOT,"results/ion_csa_nonlinear/matched_pitch$pitch");mkpath(out)
    write(joinpath(out,".gitignore"),"*.h5\n*.tmp\n")
    common_sha=sha(paths.realisation)
    reference=h5open(h->String(read(h["bb0nu/ions/4/path"])),joinpath(ROOT,"results/elm_matched_production/detector_index.h5"),"r")
    started=time_ns();ions=Ions.run_ion_sample(spec,"bb0nu",code);elapsed=(time_ns()-started)/1e9
    @assert sha(paths.realisation)==common_sha
    count=0;ledger=Int64(0)
    h5open(paths.realisation,"r") do common
    h5open(reference,"r") do four
    h5open(ions,"r") do current
        @assert Set(keys(common["events"]))==Set(keys(four["events"]))==Set(keys(current["events"]))
        for id in keys(common["events"])
            for key in ("event_id","ion_drift_cm","electron_drift_cm","initial_pairs","pair_checksum")
                @assert read(common["events/$id/$key"])==read(four["events/$id/$key"])==read(current["events/$id/$key"])
            end
            g=current["events/$id"];li=read(g["ion_drift_cm"]);le=read(g["electron_drift_cm"])
            @assert 5<=li<=145 && 5<=le<=145 && isapprox(li+le,150;atol=1e-4)
            @assert sum(read(g["count"]))==read(g["initial_pairs"])
            count+=1;ledger+=read(g["initial_pairs"])
        end
    end
    end
    end
    @assert count==8000
    target=joinpath(out,"detector_index.h5")
    h5open(target*".tmp","w") do h
        h["bb0nu/ions/$pitch/path"]=ions
        h["metadata/common_realisation_path"]=paths.realisation
        h["metadata/common_realisation_sha256"]=common_sha
        h["metadata/reference_4mm_path"]=reference
        h["metadata/ion_path"]=ions;h["metadata/ion_sha256"]=sha(ions)
        h["metadata/event_count"]=count;h["metadata/total_ions"]=ledger
        h["metadata/placement_validation_passed"]=true;h["metadata/complete"]=true
    end
    mv(target*".tmp",target;force=true)
    open(joinpath(out,"collection.toml"),"w") do io
        TOML.print(io,Dict("complete"=>true,"pitch_mm"=>pitch,"bucket_ms"=>10pitch,"events"=>count,
            "total_ions"=>ledger,"collection_wall_s"=>elapsed,"ion_path"=>ions,"ion_sha256"=>sha(ions),
            "common_realisation"=>paths.realisation,"common_sha256"=>common_sha,"detector_index"=>target,
            "card"=>card,"script_sha256"=>sha(@__FILE__),"ion_replay_code_sha256"=>code))
    end
    println("Matched ion index: $target");target
end
end
if abspath(PROGRAM_FILE)==@__FILE__;MatchedIonPitch.main(parse(Int,only(ARGS)));end
