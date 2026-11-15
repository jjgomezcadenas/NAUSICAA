#!/usr/bin/env julia
include("run_framed_reconstruction.jl")
include("../../src/vh1/VH1Readout.jl")
using .VH1Readout
using Statistics, LinearAlgebra

"""Compare VH1 and the present CSA on the same 100 saved detector selections."""
function run_vh1_pilot(cardpath)
 BLAS.set_num_threads(1)
 card=TOML.parsefile(cardpath)
 ec=TOML.parsefile(joinpath(ROOT,card["electronics_card"]))
 dc=TOML.parsefile(joinpath(ROOT,ec["detector_card"]))
 source=joinpath(ROOT,card["selection_controls"])
 outdir=joinpath(ROOT,"results/readout",dc["run_label"],card["output_subdir"])
 mkpath(outdir);target=joinpath(outdir,"controls.h5")
 isfile(target)&&error("Preserve existing pilot: $target")
 hp=joinpath(ROOT,dc["hits"]);source_sha=bytes2hex(open(sha256,hp))
 ids,placements,frames=h5open(source,"r") do f
  read(f["source_sha256"])==source_sha||error("Source hits changed")
  ev=f["events"];ids=read(ev["event_number"])[1:card["accepted_events"]]
  ps=[EventPlacement(id,read(ev["ion_drift_cm"])[j],read(ev["electron_drift_cm"])[j]) for (j,id) in enumerate(ids)]
  fs=[EventFrame(mc_vertex_frame,ntuple(k->read(ev["frame_$("xyz"[k])_min"])[j],3),ntuple(k->read(ev["frame_$("xyz"[k])_max"])[j],3)) for j in eachindex(ids)]
  ids,ps,fs
 end
 setup=@timed prepare_electronics_run(ec);prepared=setup.value;workspace=ER.CSAWorkspace()
 vhsetup=@timed load_vh1_settings(ROOT,joinpath(ROOT,card["vh1_card"]);minimum_samples=card["minimum_samples"])
 vh=vhsetup.value
 c=card["ions"];is=IonReconstructionSettings(c["baseline_ms"],c["average_ms"],c["nonlinear_below_V"],c["max_iterations"],c["convergence_ions"])
 rows=NamedTuple[];examples=Dict{String,Any}();plans=Dict{Int,Any}();started=time()
 h5open(hp,"r") do input
  rowmap=Dict(id=>j for (j,id) in enumerate(read(input["events/event_id"])))
  for (number,id) in enumerate(ids)
   event_started=time();h,info=read_mc_event(input,id;event_row=rowmap[id]);h=center_hits(h,info)
   p=placements[number];selected=hits_in_frame(h,frames[number]);ds=event_settings(dc,p)
   rng(tag)=Xoshiro(seedfor(dc["seed"],id,tag));erng(tag)=Xoshiro(seedfor(ec["seed"],id,tag))
   pairs=sample_fano_pairs(rng("pairs"),selected,Float32(dc["W_eV"]),Float32(dc["fano_factor"]))
   det=@timed sim_ion_detector_response(rng("ions"),selected,pairs,p,frames[number],ds;max_vector_bytes=dc["max_vector_bytes"])
   ions=det.value
   old=@timed begin
    s=electronics_settings(prepared,ions.ions.time)
    adc,_=ER.sim_ion_electronics_prepared(erng("csa_noise"),erng("adc_noise"),ions,s.ion,workspace)
    cal=IonReconstructionCalibration(s.ion.adc,s.ion.csa.current_pulse,s.ion.csa.linear_response,s.ion.csa.nonlinear_response,s.ion.csa.noise)
    plan=SR.prepare_ion_reconstruction(adc.time,ions.ions.time,cal,is)
    SR.reconstruct_ions(adc,plan),plan
   end
   ri_old,oldplan=old.value
   simulation=@timed sim_vh1_electronics(erng("vh1_noise"),erng("vh1_adc"),ions,vh)
   response=simulation.value
   bucket_count=Int(ions.ions.time.bucket_count)
   prep=@timed get!(plans,bucket_count) do
    prepare_vh1_reconstruction(response.time,ions.ions.time,vh)
   end
   plan=prep.value
   fit=@timed reco_vh1_ions(response,plan);ri=fit.value
   empty=vec(all(iszero,ions.ions.counts;dims=1))
   quiet=vec(all(ri.valid;dims=1)).&empty
   quiet_old=vec(all(ri_old.valid;dims=1)).&empty
   quiet_sums=vec(sum(Float64.(ri.ions[:,quiet]);dims=1))
   quiet_old_sums=vec(sum(Float64.(ri_old.ions[:,quiet_old]);dims=1))
   predicted=sum(plan.ion_covariance)
   signal_cells=(ions.ions.counts .> 0).&ri.valid
   signal_count=count(signal_cells)
   signal_sum=sum(Float64,ri.ions[signal_cells];init=0.)
   input_signal_sum=sum(ions.ions.counts[signal_cells])
   bucket_noise=sqrt.(diag(plan.ion_covariance))
   snr_sum=sum(Float64(ri.ions[i,j])/bucket_noise[i] for j in axes(ri.ions,2), i in axes(ri.ions,1) if signal_cells[i,j];init=0.)
   push!(rows,(event_number=id,ion_buckets=bucket_count,ion_pads=size(ions.ions.counts,2),ions_input=sum(ions.ions.counts),
    signal_pad_buckets=signal_count,signal_reconstructed_ions=signal_sum,signal_input_ions=input_signal_sum,signal_snr_sum=snr_sum,
    vh1_ions_reconstructed=sum(Float64,ri.ions[ri.valid];init=0.),old_ions_reconstructed=sum(Float64,ri_old.ions[ri_old.valid];init=0.),
    vh1_ions_on_invalid_pads=sum(ions.ions.counts[.!ri.valid]),old_ions_on_invalid_pads=sum(ions.ions.counts[.!ri_old.valid]),
    vh1_saturated_pads=count(response.saturated),vh1_adc_limit_pads=count(ri.checks.adc_at_limit),vh1_invalid_pads=count(.!ri.checks.converged),old_invalid_pads=count(.!ri_old.checks.converged),
    vh1_bucket_residual_rms=sqrt(mean(abs2,Float64.(ri.ions[ri.valid]).-Float64.(ions.ions.counts[ri.valid]))),
    old_bucket_residual_rms=sqrt(mean(abs2,Float64.(ri_old.ions[ri_old.valid]).-Float64.(ions.ions.counts[ri_old.valid]))),
    vh1_quiet_pads=length(quiet_sums),vh1_quiet_sum=sum(quiet_sums),vh1_quiet_sumsq=sum(abs2,quiet_sums),
    old_quiet_pads=length(quiet_old_sums),old_quiet_sum=sum(quiet_old_sums),old_quiet_sumsq=sum(abs2,quiet_old_sums),
    vh1_predicted_quiet_variance=predicted,vh1_predicted_bucket_noise=sqrt(mean(diag(plan.ion_covariance))),
    detector_s=det.time,old_chain_s=old.time,vh1_simulation_s=simulation.time,vh1_plan_s=prep.time,vh1_reconstruction_s=fit.time,event_s=time()-event_started))
   if number<=card["example_events"]
    chosen=unique([argmax(vec(sum(ions.ions.counts;dims=1))),something(findfirst(empty),1)])
    ex=Dict{String,Any}()
    for j in chosen
     clean=ER.prepare_bucket_current(vh.current_pulse,ions.ions.time.bucket_ms)
     grid=VH1Readout.internal_grid(response.time,vh)
     cur=ER.prepared_ion_current(vec(ions.ions.counts[:,j]),ions.ions.time,clean,grid)
     readings=zeros(size(response.codes,1));VH1Readout.frame_readings!(readings,cur,response.time,vh)
     ex["pad_$j"]=Dict("input_ions"=>ions.ions.counts[:,j],"vh1_codes"=>response.codes[:,j],"vh1_noiseless_V"=>readings,
      "vh1_reconstructed_ions"=>ri.ions[:,j],"old_reconstructed_ions"=>ri_old.ions[:,j])
    end
    examples[string(id)]=ex
   end
   println("VH1 event $number/$(length(ids)) MC $id: saturated $(count(response.saturated)), invalid $(count(.!ri.checks.converged)), time $(round(rows[end].event_s;digits=3)) s");flush(stdout)
  end
 end
 budget=16*1024^2+3Base.summarysize(examples)+3Base.summarysize(rows)
 GuardedFiles.space_check(outdir,budget,1024^3);tmp=target*".tmp"
 try
  h5open(tmp,"w") do f
   f["schema"]="itaca.vh1_pilot.v1";f["card"]=read(cardpath,String);f["vh1_card"]=read(joinpath(ROOT,card["vh1_card"]),String)
   f["electronics_card"]=read(joinpath(ROOT,card["electronics_card"]),String)
   f["source_sha256"]=source_sha;f["selection_sha256"]=bytes2hex(open(sha256,source))
   f["setup_old_s"]=setup.time;f["setup_vh1_s"]=vhsetup.time;f["total_s"]=time()-started
   columns(create_group(f,"events"),rows);write_control_tree(create_group(f,"examples"),examples)
   flush(f);GuardedFiles.check_size(tmp,budget)
  end
  Base.rename(tmp,target)
 catch;isfile(tmp)&&rm(tmp);rethrow();end
 println("FINISHED $target")
end
if abspath(PROGRAM_FILE)==@__FILE__;run_vh1_pilot(abspath(ARGS[1]));end
