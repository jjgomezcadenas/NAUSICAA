include("../production/run_framed_electronics.jl")
"""Benchmark both initial states on identical photon vectors and SiPM settings.
Only the initial recovery mode changes. Full SiPM and ASIC calls are timed after
warm-up, including all zero-input channels and dark activity."""
function benchmark_initialisation(;comparison=:initialisation)
 comparison in (:initialisation,:substeps)||error("Unknown benchmark comparison")
 BLAS.set_num_threads(1)
 ec=TOML.parsefile(joinpath(ROOT,"config/readout/bb0nu_electronics_pilot_100.toml"))
 dc=TOML.parsefile(joinpath(ROOT,ec["detector_card"]))
 base=joinpath(ROOT,"results/readout",dc["run_label"]);source=joinpath(base,"detector_response/controls.h5")
 outdir=joinpath(base,comparison==:initialisation ? "electronics_sipm_initialisation" : "electronics_sipm_substeps64");mkpath(outdir);target=joinpath(outdir,"controls.h5")
 isfile(target)&&error("Preserve existing benchmark")
 s=prepare_electronics_run(ec).sipm
 before=SiPMReadoutSettings(s.sipm,s.illuminated_microcells,SiPMStatisticalOptions(128,s.statistical.sparse_trigger_limit;initial_state=comparison==:initialisation ? :stationary_dark : :fully_recovered),s.asic)
 after=SiPMReadoutSettings(s.sipm,s.illuminated_microcells,SiPMStatisticalOptions(comparison==:initialisation ? 128 : 64,s.statistical.sparse_trigger_limit;initial_state=:fully_recovered),s.asic)
 ids,ps,fs=h5open(source,"r") do f
  e=f["events"];ids=read(e["event_number"])[1:10]
  ps=[EventPlacement(id,read(e["ion_drift_cm"])[j],read(e["electron_drift_cm"])[j]) for (j,id) in enumerate(ids)]
  fs=[EventFrame(mc_vertex_frame,ntuple(k->read(e["frame_$("xyz"[k])_min"])[j],3),ntuple(k->read(e["frame_$("xyz"[k])_max"])[j],3)) for j in eachindex(ids)]
  ids,ps,fs
 end
 rows=NamedTuple[]
 h5open(joinpath(ROOT,dc["hits"]),"r") do f
  rowmap=Dict(id=>i for (i,id) in enumerate(read(f["events/event_id"])))
  for (number,id) in enumerate(ids)
   hits,info=read_mc_event(f,id;event_row=rowmap[id]);hits=hits_in_frame(center_hits(hits,info),fs[number])
   rng(tag)=Xoshiro(seedfor(dc["seed"],id,tag));erng(tag)=Xoshiro(seedfor(ec["seed"],id,tag))
   pairs=sample_fano_pairs(rng("pairs"),hits,Float32(dc["W_eV"]),Float32(dc["fano_factor"]))
   photons=sim_electron_detector_response(rng("electrons"),hits,pairs,ps[number],fs[number],event_settings(dc,ps[number]);max_vector_bytes=dc["max_vector_bytes"])
   bt=Float64[];at=Float64[];bb=Int[];ab=Int[]
   for mode in (before,after);sim_sipm_readout(erng("sipm"),erng("asic"),photons,mode);end
   for repeat in 1:3
    # Alternate execution order to avoid assigning a systematic warm-cache advantage.
    for (label,mode) in (isodd(repeat) ? ((:old,before),(:new,after)) : ((:new,after),(:old,before)))
     result=@timed sim_sipm_readout(erng("sipm"),erng("asic"),photons,mode)
     if label==:old;push!(bt,result.time);push!(bb,result.bytes)
     else;push!(at,result.time);push!(ab,result.bytes);end
    end
   end
   push!(rows,(event_number=id,channels=size(photons.photons.counts,2),buckets=size(photons.photons.counts,1),before_s=median(bt),after_s=median(at),before_allocated_bytes=median(bb),after_allocated_bytes=median(ab)))
   println("Event $id: $(round(median(bt);digits=4)) → $(round(median(at);digits=4)) s");flush(stdout)
  end
 end
 budget=4*1024^2;GuardedFiles.space_check(outdir,budget,1024^3);tmp=target*".tmp"
 try
  h5open(tmp,"w") do f
   f["schema"]="itaca.sipm_timing_benchmark.v1";f["comparison"]=String(comparison);
   f["before_substeps"]=before.statistical.substeps;f["after_substeps"]=after.statistical.substeps
   f["before_initial_state"]=String(before.statistical.initial_state);f["after_initial_state"]=String(after.statistical.initial_state);columns(create_group(f,"events"),rows)
   f["summary/events"]=length(rows);f["summary/before_s"]=mean(r.before_s for r in rows);f["summary/after_s"]=mean(r.after_s for r in rows)
   f["summary/speedup"]=mean(r.before_s for r in rows)/mean(r.after_s for r in rows)
   f["summary/before_MiB"]=mean(r.before_allocated_bytes for r in rows)/1024^2;f["summary/after_MiB"]=mean(r.after_allocated_bytes for r in rows)/1024^2
   f["summary/mean_initial_draws_removed"]=comparison==:initialisation ? mean(r.channels for r in rows)*s.sipm.microcells : 0.
   f["sipm_card"]=read(joinpath(ROOT,ec["sipm_card"]),String);f["asic_card"]=read(joinpath(ROOT,ec["asic_card"]),String)
   f["source_sha256"]=bytes2hex(open(sha256,source));f["code_sha256"]=bytes2hex(open(sha256,joinpath(ROOT,"src/ion_signal/sipm_statistical.jl")))
   flush(f);GuardedFiles.check_size(tmp,budget)
  end
  mv(tmp,target)
 catch;isfile(tmp)&&rm(tmp);rethrow();end
 println(target)
end
if abspath(PROGRAM_FILE)==@__FILE__
 benchmark_initialisation()
end
