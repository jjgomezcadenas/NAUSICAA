#!/usr/bin/env julia
const PROCESS_STARTED=time()
include("../production/run_framed_reconstruction.jl")

"""Measure the complete in-memory readout of 100 accepted events after one warm-up.
Save stage times and response checks; do not save event vectors or perform ZS.
Existing detector controls supply the selected events, positions and frames.
"""
function benchmark_chain()
 BLAS.set_num_threads(1)
 rcpath=joinpath(ROOT,"config/readout/bb0nu_reconstruction_100.toml")
 rc=TOML.parsefile(rcpath);ec=TOML.parsefile(joinpath(ROOT,rc["electronics_card"]))
 dc=TOML.parsefile(joinpath(ROOT,ec["detector_card"]))
 base=joinpath(ROOT,"results/readout",dc["run_label"])
 source=joinpath(base,"detector_response/controls.h5")
 outdir=joinpath(base,"readout_timing_100");mkpath(outdir)
 target=joinpath(outdir,"controls.h5");isfile(target)&&error("Preserve existing benchmark")
 input_started=time();hp=joinpath(ROOT,dc["hits"]);sha=bytes2hex(open(sha256,hp))
 ids,ps,fs=h5open(source,"r") do f
  read(f["source_sha256"])==sha||error("Source changed")
  e=f["events"];ids=read(e["event_number"])[1:100]
  ps=[EventPlacement(id,read(e["ion_drift_cm"])[j],read(e["electron_drift_cm"])[j]) for (j,id) in enumerate(ids)]
  fs=[EventFrame(mc_vertex_frame,ntuple(k->read(e["frame_$("xyz"[k])_min"])[j],3),ntuple(k->read(e["frame_$("xyz"[k])_max"])[j],3)) for j in eachindex(ids)]
  ids,ps,fs
 end
 input_setup_s=time()-input_started
 setup=@timed prepare_electronics_run(ec);prepared=setup.value;workspace=ER.CSAWorkspace()
 c=rc["ions"];settings=IonReconstructionSettings(c["baseline_ms"],c["average_ms"],c["nonlinear_below_V"],c["max_iterations"],c["convergence_ions"])
 es=ElectronReconstructionSettings(prepared.sipm.sipm,prepared.sipm.illuminated_microcells,rc["electrons"]["subtract_mean_dark"])
 rows=NamedTuple[];warmup_s=0.;loop_started=time()
 h5open(hp,"r") do input
  rowmap=Dict(id=>j for (j,id) in enumerate(read(input["events/event_id"])))
  for number in vcat(1,collect(1:100))
   warmup=isempty(rows)&&warmup_s==0.;id=ids[number];t=time()
   loading=@timed begin
    h,info=read_mc_event(input,id;event_row=rowmap[id]);hits_in_frame(center_hits(h,info),fs[number])
   end
   h=loading.value;p=ps[number];frame=fs[number]
   rng(tag)=Xoshiro(seedfor(dc["seed"],id,tag));erng(tag)=Xoshiro(seedfor(ec["seed"],id,tag))
   pairs_setup=@timed (event_settings(dc,p),sample_fano_pairs(rng("pairs"),h,Float32(dc["W_eV"]),Float32(dc["fano_factor"])))
   ds,pairs=pairs_setup.value
   iondet=@timed sim_ion_detector_response(rng("ions"),h,pairs,p,frame,ds;max_vector_bytes=dc["max_vector_bytes"])
   eledet=@timed sim_electron_detector_response(rng("electrons"),h,pairs,p,frame,ds;max_vector_bytes=dc["max_vector_bytes"])
   ions=iondet.value;photons=eledet.value
   ep=@timed electronics_settings(prepared,ions.ions.time);s=ep.value
   csa=@timed ER.sim_ion_electronics_prepared(erng("csa_noise"),erng("adc_noise"),ions,s.ion,workspace)
   adc,ic=csa.value
   sipm=@timed sim_sipm_readout(erng("sipm"),erng("asic"),photons,s.sipm);asic,ac=sipm.value
   prep=@timed begin
    cal=IonReconstructionCalibration(s.ion.adc,s.ion.csa.current_pulse,s.ion.csa.linear_response,s.ion.csa.nonlinear_response,s.ion.csa.noise)
    SR.prepare_ion_reconstruction(adc.time,ions.ions.time,cal,settings)
   end
   ireco=@timed SR.reconstruct_ions(adc,prep.value);ri=ireco.value
   ereco=@timed reco_electron_signal(asic,es)
   rebin=@timed rebin_electron_signal(ereco.value,Int64(rc["electrons"]["rebin_factor"]))
   event_s=time()-t
   row=(event_number=id,input_s=loading.time,pairs_settings_s=pairs_setup.time,ion_detector_s=iondet.time,electron_detector_s=eledet.time,
    electronics_preparation_s=ep.time,csa_adc_s=csa.time,sipm_asic_s=sipm.time,ion_reconstruction_preparation_s=prep.time,
    ion_reconstruction_s=ireco.time,electron_reconstruction_s=ereco.time,electron_rebin_s=rebin.time,event_s=event_s,
    gc_s=sum(x.gctime for x in (loading,pairs_setup,iondet,eledet,ep,csa,sipm,prep,ireco,ereco,rebin)),
    ion_channels=size(adc.counts,2),sipm_channels=size(asic.response_pe,2),nonlinear_pads=count(ri.checks.nonlinear),failed_pads=count(.!ri.checks.converged),invalid_pe=count(.!rebin.value.valid))
   if warmup;warmup_s=event_s;println("Warm-up: $warmup_s s")
   else;push!(rows,row);println("Event $number/100 ($id): $(round(event_s;digits=3)) s; CSA $(round(csa.time;digits=3)); SiPM $(round(sipm.time;digits=3)); ion reco $(round(ireco.time;digits=3))")
   end
   flush(stdout)
  end
 end
 loop_s=time()-loop_started
 budget=4*1024^2;GuardedFiles.space_check(outdir,budget,1024^3);tmp=target*".tmp"
 try
  h5open(tmp,"w") do f
   f["schema"]="itaca.readout_timing.v1";columns(create_group(f,"events"),rows)
   f["warmup_s"]=warmup_s;f["input_setup_s"]=input_setup_s;f["electronics_setup_s"]=setup.time;f["loop_s"]=loop_s
   f["elapsed_before_write_s"]=time()-PROCESS_STARTED;f["source_sha256"]=sha
   f["cards/reconstruction"]=read(rcpath,String);f["cards/electronics"]=read(joinpath(ROOT,rc["electronics_card"]),String)
   for key in ["detector_card","csa_card","sipm_card","asic_card"];f["cards/$key"]=read(joinpath(ROOT,ec[key]),String);end
   f["calibration_sha256"]=bytes2hex(open(sha256,joinpath(ROOT,ec["nonlinear_calibration"])))
   files=vcat([@__FILE__,joinpath(ROOT,"scripts/production/framed_electronics_settings.jl")],[joinpath(dir,n) for sub in ["src/reconstruction","src/electronics","src/detector_response","src/ion_signal"] for (dir,_,names) in walkdir(joinpath(ROOT,sub)) for n in names if endswith(n,".jl")])
   f["code_files"]=relpath.(files,ROOT);f["code_sha256"]=[bytes2hex(open(sha256,p)) for p in files]
   flush(f);GuardedFiles.check_size(tmp,budget)
  end
  mv(tmp,target)
 catch;isfile(tmp)&&rm(tmp);rethrow();end
 println("FINISHED $target; elapsed $(time()-PROCESS_STARTED) s")
end
benchmark_chain()
