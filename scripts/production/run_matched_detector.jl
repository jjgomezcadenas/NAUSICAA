#!/usr/bin/env julia
module MatchedDetector
using HDF5,TOML,SHA
module Common
include("run_common_realisation.jl")
end
module Ions
include("run_ion_readout.jl")
end
module Electrons
include("run_electron_250ns.jl")
end
include(joinpath(@__DIR__,"../../src/ion_signal/guarded_files.jl"))
const ROOT=Common.ROOT

"""Run shared placement, both carrier collections, and validate every paired event.

Writes a joint detector index compatible with the existing ion CSA and SiPM
runners. Its `bb0nu/ions/4/path` and `bb0nu/xe_he15_path` branches point to
immutable physical histories. Pair checksums, source IDs and both drift values
must match the common realization exactly before the index is marked complete.
"""
function main(kind="smoke")
    kind in ("smoke","full") || error("Choose smoke or full")
    card="config/specs/runs/elm_matched_detector_$kind.toml"
    spec=Common.load_specs(ROOT,card)
    !haskey(spec.electrons.data,"placement") || error("Matched production uses common placement")
    out=joinpath(ROOT,"results",kind=="smoke" ? "elm_matched_smoke" : "elm_matched_production")
    mkpath(out);started=time();println("Matched detector PID $(getpid())");flush(stdout)
    common=Common.run_common(spec,"bb0nu")
    ions=Ions.run_ion_sample(spec,"bb0nu",Ions.ion_replay_code_hash())
    index=Electrons.main_electrons(card)
    target=joinpath(out,"detector_index.h5");GuardedFiles.checked_copy(index,target*".tmp";force=true)
    h5open(target*".tmp","r+") do idx
        idx["bb0nu/ions/4/path"]=ions
        ph=String(read(idx["bb0nu/xe_he15_path"]))
        h5open(common,"r") do c
          h5open(ions,"r") do i
            h5open(ph,"r") do p
              @assert Set(keys(c["events"]))==Set(keys(i["events"]))==Set(keys(p["events"]))
              for id in keys(c["events"])
                for key in ("event_id","ion_drift_cm","electron_drift_cm","initial_pairs","pair_checksum")
                  @assert read(c["events/$id/$key"])==read(i["events/$id/$key"])==read(p["events/$id/$key"])
                end
                li=read(c["events/$id/ion_drift_cm"]);le=read(c["events/$id/electron_drift_cm"])
                @assert 5<=li<=145 && 5<=le<=145 && isapprox(li+le,150;atol=1e-4)
              end
              idx["metadata/matched_events"]=length(keys(c["events"]))
              idx["metadata/placement_validation_passed"]=true
            end
          end
        end
        idx["metadata/common_realisation_path"]=common
        idx["metadata/ion_path"]=ions
        idx["metadata/detector_wall_s"]=time()-started
    end
    mv(target*".tmp",target;force=true)
    open(joinpath(out,"detector_files.toml"),"w") do io
        TOML.print(io,Dict("complete"=>true,"detector_index"=>target,"common_realisation"=>common,
            "ions"=>ions,"detector_wall_s"=>time()-started,"detector_card"=>card))
    end
    println("Matched detector index: $target");target
end
end
if abspath(PROGRAM_FILE)==@__FILE__;MatchedDetector.main(isempty(ARGS) ? "smoke" : only(ARGS));end
