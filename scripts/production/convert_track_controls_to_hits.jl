#!/usr/bin/env julia
include("run_hit_production.jl")
"""Convert the existing all-event track controls to compact hits, without simulation.
Preserve stored hit values and pre-cut energy estimates. Read acquisition frames
from the checksum-verified selection controls recorded by the original run.
"""
function convert_track_controls(source,target)
 h5open(source,"r") do f
  rc=TOML.parse(read(f["card"],String));dc=TOML.parse(read(f["cards/detector_card"],String))
  selection=joinpath(ROOT,rc["selection_controls"])
  bytes2hex(open(sha256,selection))==read(f["detector_controls_sha256"],String)||error("Selection changed")
  cards=Dict("run"=>read(f["card"],String))
  for key in keys(f["cards"]);cards[key]=read(f["cards/$key"],String);end
  hashes=Dict("source_controls"=>bytes2hex(open(sha256,source)),"selection"=>bytes2hex(open(sha256,selection)),"converter"=>bytes2hex(open(sha256,@__FILE__)))
  meta=HitFileMetadata(dc["run_label"],"bb0nu",Xe_85_He_15,cards,hashes)
  tc=rc["track_analysis"];settings=HitOutputSettings(tc["ion_hit_cut_keV"],tc["electron_hit_cut_keV"],4,65536,64*1024^2)
  w=open_hit_output(target,meta,settings)
  try
   h5open(selection,"r") do selected
    g=selected["events"];lookup=Dict(id=>i for (i,id) in enumerate(read(g["event_number"])))
    ids=read(f["events/event_number"])
    for (j,id) in enumerate(ids)
     row=lookup[id];value(key)=read(g[key])[row]
     frame=EventFrame(mc_vertex_frame,ntuple(k->value("frame_$("xyz"[k])_min"),3),ntuple(k->value("frame_$("xyz"[k])_max"),3))
     sigma(arm)=read(f["hit_histograms/$id/cut_statistics/$arm/topology/sigma_t_mm"])
     info=HitEventInfo(id,value("ion_drift_cm"),value("electron_drift_cm"),frame,sigma("ions"),sigma("pe"))
     energy=EventEnergy((read(f["events/$key"])[j] for key in ("ions_total","pe_total","ion_energy_keV","electron_energy_keV","invalid_ion_buckets","invalid_pe_buckets"))...)
     hits(arm)=SavedHits((read(f["track_examples/$id/$arm/hits/$key"]) for key in ("x_mm","y_mm","z_mm","energy_keV"))...)
     write_hit_event(w,HitEvent(info,energy,hits("itrack"),hits("etrack")))
    end
   end
   finish_hit_output(w)
  catch;abort_hit_output(w);rethrow();end
 end
 println("CONVERTED $target bytes=$(filesize(target))")
end
if abspath(PROGRAM_FILE)==@__FILE__
 convert_track_controls(abspath(ARGS[1]),abspath(ARGS[2]))
end
