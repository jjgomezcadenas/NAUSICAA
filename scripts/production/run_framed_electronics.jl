#!/usr/bin/env julia
include("run_framed_detector_response.jl")
include("../../src/electronics/ElectronicsResponse.jl")
using .ElectronicsResponse
include("framed_electronics_settings.jl")

"""Run electronics on selected accepted MC events; retain numerical controls and example vectors."""
function run_electronics(cardpath)
 card=TOML.parsefile(cardpath);c=TOML.parsefile(joinpath(ROOT,card["detector_card"]))
 base=joinpath(ROOT,"results/readout",c["run_label"]);source=joinpath(base,"detector_response/controls.h5")
 outdir=joinpath(base,get(card,"output_subdir","electronics"));mkpath(outdir);target=joinpath(outdir,"controls.h5")
 isfile(target)&&error("Preserve existing controls: $target")
 hp=joinpath(ROOT,c["hits"])
 source_hash=bytes2hex(open(sha256,hp))
 ids,placements,frames=h5open(source,"r") do f
  read(f["source_sha256"])==source_hash||error("Detector source changed")
  ev=f["events"];ids=read(ev["event_number"])[1:card["accepted_events"]]
  ps=[EventPlacement(id,read(ev["ion_drift_cm"])[j],read(ev["electron_drift_cm"])[j]) for (j,id) in enumerate(ids)]
  fs=[EventFrame(mc_vertex_frame,ntuple(k->read(ev["frame_$("xyz"[k])_min"])[j],3),ntuple(k->read(ev["frame_$("xyz"[k])_max"])[j],3)) for j in eachindex(ids)]
  ids,ps,fs
 end
 rows=NamedTuple[];current_exceedances=NamedTuple[];examples=Dict{Int64,Any}();started=time()
 prepared=prepare_electronics_run(card);workspace=ER.CSAWorkspace()
 h5open(hp,"r") do input
  inputids=read(input["events/event_id"]);rowmap=Dict(id=>j for (j,id) in enumerate(inputids))
  for (number,id) in enumerate(ids)
   h,info=read_mc_event(input,id;event_row=rowmap[id]);h=center_hits(h,info);p=placements[number];frame=frames[number]
   selected=hits_in_frame(h,frame);settings=event_settings(c,p)
   rng(tag)=Xoshiro(seedfor(c["seed"],id,tag));erng(tag)=Xoshiro(seedfor(card["seed"],id,tag))
   pairs=sample_fano_pairs(rng("pairs"),selected,Float32(c["W_eV"]),Float32(c["fano_factor"]))
   t=time()
   ions=sim_ion_detector_response(rng("ions"),selected,pairs,p,frame,settings;max_vector_bytes=c["max_vector_bytes"])
   photons=sim_electron_detector_response(rng("electrons"),selected,pairs,p,frame,settings;max_vector_bytes=c["max_vector_bytes"])
   detector_s=time()-t;ion_started=time();t=time()
   preparation=@timed electronics_settings(prepared,ions.ions.time);s=preparation.value;setup_s=time()-t
   calculation=@timed ER.sim_ion_electronics_prepared(erng("csa_noise"),erng("adc_noise"),ions,s.ion,workspace)
   adc,ic=calculation.value;ion_s=calculation.time;spice_s=0.
   at=ion_adc_time(ions.ions.time,s.ion.csa.current_pulse,s.ion.acquisition,s.ion.adc)
   t=time();exceedances=0;maximum_current_ratio=0.
   for j in axes(ions.ions.counts,2)
    all(iszero,view(ions.ions.counts,:,j))&&continue
    peak=workspace.peak_current_A[j];ratio=peak/s.nonlinear.maximum_current_A
    maximum_current_ratio=max(maximum_current_ratio,ratio)
    if ratio>1
     exceedances+=1
     push!(current_exceedances,(event_number=id,channel=j,pixel_ix=ions.ions.channels.ix[j],pixel_iy=ions.ions.channels.iy[j],peak_current_A=peak,calibration_limit_A=s.nonlinear.maximum_current_A,ratio))
    end
   end
   current_check_s=time()-t
   ion_total_s=time()-ion_started
   ion_allocated_bytes=preparation.bytes+calculation.bytes
   ion_prepared_bytes=Base.summarysize(prepared.cache)+Base.summarysize(workspace)
   t=time();asic,ec=sim_sipm_readout(erng("sipm"),erng("asic"),photons,s.sipm);sipm_s=time()-t
   push!(rows,(event_number=id,detector_s,setup_s,current_check_s,ion_s,ion_total_s,ion_allocated_bytes,ion_prepared_bytes,spice_s,sipm_s,peak_current_exceedances=exceedances,maximum_current_ratio,ion_start_ms=adc.time.start_ms,electron_start_us=asic.time.start_us,ion_channels=size(adc.counts,2),ion_samples=size(adc.counts,1),sipm_channels=size(asic.response_pe,2),sipm_buckets=size(asic.response_pe,1),
    ion_vector_bytes=sizeof(adc.counts),sipm_vector_bytes=sizeof(asic.response_pe)+sizeof(asic.clipped),linear_channels=ic.linear_channels,nonlinear_channels=ic.nonlinear_channels,spice_channels=ic.spice_channels,
    adc_below=ic.adc_below_range,adc_above=ic.adc_above_range,detected_pe=ec.detected_pe,dark_pulses=ec.dark_pulses,sipm_charge_pe=ec.sipm_charge_pe,
    asic_clipped_buckets=ec.asic_clipped_buckets,asic_charge_lost_pe=ec.asic_charge_lost_pe,asic_sum_pe=sum(Float64,asic.response_pe)))
   if number<=get(card,"example_events",5)
   examples[id]=Dict{String,Any}()
   for (name,counts) in (("ions",ions.ions.counts),("sipm",photons.photons.counts))
    totals=vec(sum(counts;dims=1));chosen=unique(vcat(argmax(totals),something(findfirst(iszero,totals),argmin(totals))))
    for j in chosen
     ex=Dict{String,Any}("channel"=>j,"input"=>counts[:,j])
     if name=="ions"
      wave=ER.pad_voltage(view(counts,:,j),ions.ions.time,s.ion,at)
      ex["time_ms"]=1000wave.time_s;ex["voltage_V"]=wave.voltage_V;ex["adc_counts"]=adc.counts[:,j]
      ex["adc_midpoint_V"]=s.ion.adc.minimum_V.+(Float64.(adc.counts[:,j]).+.5).*((s.ion.adc.maximum_V-s.ion.adc.minimum_V)/2^s.ion.adc.bits)
      ex["input_time_ms"]=[(ions.ions.time.start_ms-at.start_ms)+(k-.5)*ions.ions.time.bucket_ms for k in axes(counts,1)]
      ex["current_time_ms"]=1000wave.current.time_s;ex["current_A"]=wave.current.current_A
     else
      ex["time_us"]=(collect(1:size(counts,1)).-.5).*photons.photons.time.bucket_us
      ex["response_pe"]=asic.response_pe[:,j];ex["clipped"]=UInt8.(asic.clipped[:,j])
     end
     examples[id]["$(name)_$j"]=ex
    end
   end
   end
   println("Event $id: ion $(round(ion_total_s;digits=3)) s; SiPM $(round(sipm_s;digits=3)) s; nonlinear $(ic.nonlinear_channels)/$(size(adc.counts,2)); current exceedances $exceedances; ADC range $(ic.adc_below_range+ic.adc_above_range); ASIC clipped $(ec.asic_clipped_buckets)");flush(stdout)
  end
 end
 # Explicit bound: retained example arrays, scalar controls and HDF5 metadata.
 data_bytes=sum(sizeof(v) for d in values(examples) for e in values(d) for v in values(e))
 budget=8*1024^2+2data_bytes;GuardedFiles.space_check(outdir,budget,1024^3);tmp=target*".tmp"
 try
  h5open(tmp,"w") do f
   f["schema"]="itaca.framed_electronics_controls.v1";f["noise_sampler"]="prepared_circulant_gaussian_v2";f["card"]=read(cardpath,String);f["source_sha256"]=source_hash
   f["detector_controls_sha256"]=bytes2hex(open(sha256,source));f["total_s"]=time()-started
   f["ion_setup_s"]=prepared.ion_setup_s;f["sipm_setup_s"]=prepared.sipm_setup_s
   f["linear_method"]="exact_piecewise_linear_current_recurrence"
   for key in ["detector_card","csa_card","sipm_card","asic_card"]
    f["cards/$key"]=read(joinpath(ROOT,card[key]),String)
   end
   f["nonlinear_calibration_sha256"]=bytes2hex(open(sha256,joinpath(ROOT,card["nonlinear_calibration"])))
   files=vcat([@__FILE__,joinpath(@__DIR__,"framed_electronics_settings.jl")],
       [joinpath(dir,n) for (dir,_,names) in walkdir(joinpath(ROOT,"src/electronics")) for n in names if endswith(n,".jl")])
   append!(files,[joinpath(ROOT,"src/ion_signal",name) for name in ["ion_csa_nonlinear.jl","ion_csa_nonlinear_io.jl","sipm_response.jl","sipm_statistical.jl","sipm_electronics.jl","bucket_noise_covariance.jl","averaging_scan.jl","csa_resources.jl"]])
   csa_card=TOML.parsefile(joinpath(ROOT,card["csa_card"]))
   for (key,value) in CSAResources.verify(ROOT,csa_card);f["resources/$key"]=value;end
   f["code_files"]=relpath.(files,ROOT);f["code_sha256"]=[bytes2hex(open(sha256,p)) for p in files]
   f["current_limit_policy"]="record_only";f["nonlinear_voltage_limit_policy"]="stop";f["automatic_spice"]=false
   columns(create_group(f,"current_exceedances"),current_exceedances)
   columns(create_group(f,"events"),rows);histograms(create_group(f,"histograms"),rows,20)
   for (id,d) in examples,(name,e) in d,(key,value) in e
    f["examples/$id/$name/$key"]=value
   end
   validation=TOML.parsefile(joinpath(ROOT,"work/framed_electronics/validation.toml"))
   for (key,value) in validation;f["validation/$key"]=value;end
   flush(f);GuardedFiles.check_size(tmp,budget)
  end
  mv(tmp,target)
 catch
  isfile(tmp)&&rm(tmp);rethrow()
 end
 println("FINISHED: $target")
end
if abspath(PROGRAM_FILE)==@__FILE__
 run_electronics(isempty(ARGS) ? joinpath(ROOT,"config/readout/bb0nu_electronics_pilot.toml") : abspath(only(ARGS)))
end
