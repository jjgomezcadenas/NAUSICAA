#!/usr/bin/env julia
include("run_framed_electronics.jl")
include("../../src/reconstruction/SignalReconstruction.jl")
using .SignalReconstruction
include("reconstruction_controls.jl")
const SR=SignalReconstruction
"""Run the same accepted events through electronics and reconstruction in memory.
Save per-event numerical ZS controls and selected vectors, without event products
or suppression. The source detector controls fix event positions and frames."""
function run_reconstruction(cardpath;limit=0)
 BLAS.set_num_threads(1)
 rc=TOML.parsefile(cardpath);ec=TOML.parsefile(joinpath(ROOT,rc["electronics_card"]))
 dc=TOML.parsefile(joinpath(ROOT,ec["detector_card"]))
 base=joinpath(ROOT,"results/readout",dc["run_label"])
 source=joinpath(base,"detector_response/controls.h5")
 outdir=limit>0 ? joinpath(ROOT,"work/reconstruction","pilot$limit") : joinpath(base,rc["output_subdir"])
 mkpath(outdir);target=joinpath(outdir,"controls.h5");isfile(target)&&error("Preserve existing controls: $target")
 hp=joinpath(ROOT,dc["hits"]);sha=bytes2hex(open(sha256,hp))
 ids,placements,frames=h5open(source,"r") do f
  read(f["source_sha256"])==sha||error("Detector source changed")
  ev=f["events"];n=limit>0 ? limit : rc["accepted_events"];ids=read(ev["event_number"])[1:n]
  ps=[EventPlacement(id,read(ev["ion_drift_cm"])[j],read(ev["electron_drift_cm"])[j]) for (j,id) in enumerate(ids)]
  fs=[EventFrame(mc_vertex_frame,ntuple(k->read(ev["frame_$("xyz"[k])_min"])[j],3),ntuple(k->read(ev["frame_$("xyz"[k])_max"])[j],3)) for j in eachindex(ids)]
  ids,ps,fs
 end
 setup=@timed prepare_electronics_run(ec);prepared=setup.value;workspace=ER.CSAWorkspace()
 c=rc["ions"];is=IonReconstructionSettings(c["baseline_ms"],c["average_ms"],c["nonlinear_below_V"],c["max_iterations"],c["convergence_ions"])
 es=ElectronReconstructionSettings(prepared.sipm.sipm,prepared.sipm.illuminated_microcells,rc["electrons"]["subtract_mean_dark"])
 factor=Int64(rc["electrons"]["rebin_factor"]);zs=rc["zs_controls"]
 cuts_i=collect(0.:zs["ion_threshold_step"]:zs["ion_threshold_maximum"])
 cuts_e=collect(0.:zs["pe_threshold_step"]:zs["pe_threshold_maximum"])
 rows=NamedTuple[];event_controls=Dict{String,Any}();examples=Dict{String,Any}();started=time()
 h5open(hp,"r") do input
  rowmap=Dict(id=>j for (j,id) in enumerate(read(input["events/event_id"])))
  for (number,id) in enumerate(ids)
   event_started=time()
   h,info=read_mc_event(input,id;event_row=rowmap[id]);h=center_hits(h,info);p=placements[number];frame=frames[number]
   selected=hits_in_frame(h,frame);ds=event_settings(dc,p)
   rng(tag)=Xoshiro(seedfor(dc["seed"],id,tag));erng(tag)=Xoshiro(seedfor(ec["seed"],id,tag))
   pairs=sample_fano_pairs(rng("pairs"),selected,Float32(dc["W_eV"]),Float32(dc["fano_factor"]))
   det=@timed begin
    ions=sim_ion_detector_response(rng("ions"),selected,pairs,p,frame,ds;max_vector_bytes=dc["max_vector_bytes"])
    photons=sim_electron_detector_response(rng("electrons"),selected,pairs,p,frame,ds;max_vector_bytes=dc["max_vector_bytes"])
    ions,photons
   end
   ions,photons=det.value
   electronics=@timed begin
    s=electronics_settings(prepared,ions.ions.time)
    adc,ic=ER.sim_ion_electronics_prepared(erng("csa_noise"),erng("adc_noise"),ions,s.ion,workspace)
    asic,ac=sim_sipm_readout(erng("sipm"),erng("asic"),photons,s.sipm)
    s,adc,asic
   end
   s,adc,asic=electronics.value
   cal=IonReconstructionCalibration(s.ion.adc,s.ion.csa.current_pulse,s.ion.csa.linear_response,s.ion.csa.nonlinear_response,s.ion.csa.noise)
   preparation=@timed SR.prepare_ion_reconstruction(adc.time,ions.ions.time,cal,is);plan=preparation.value
   ionfit=@timed SR.reconstruct_ions(adc,plan);ri=ionfit.value
   electronfit=@timed rebin_electron_signal(reco_electron_signal(asic,es),factor);pe=electronfit.value
   ei=grouped_input(photons.photons.counts,factor).*s.sipm.sipm.pde
   control_time=@timed begin
    ci=zs_controls(ri.ions,ri.valid,ions.ions.counts,cuts_i,zs["histogram_bins"])
    ce=zs_controls(pe.pe,pe.valid,ei,cuts_e,zs["histogram_bins"])
    event_controls[string(id)]=Dict("ions"=>ci,"pe"=>ce)
   end
   # Noise in fully empty pads: test the complete acquisition and inverse.
   empty_pads=vec(all(iszero,ions.ions.counts;dims=1));quiet=ri.ions[:,empty_pads]
   quiet_sums=vec(sum(Float64.(quiet);dims=1))
   push!(rows,(event_number=id,detector_s=det.time,electronics_s=electronics.time,
    ion_preparation_s=preparation.time,ion_reconstruction_s=ionfit.time,electron_reconstruction_s=electronfit.time,controls_s=control_time.time,
    reconstruction_allocated_bytes=preparation.bytes+ionfit.bytes+electronfit.bytes,prepared_bytes=Base.summarysize(plan),
    ion_channels=size(ri.ions,2),ion_buckets=size(ri.ions,1),sipm_channels=size(pe.pe,2),electron_buckets=size(pe.pe,1),
    nonlinear_pads=count(ri.checks.nonlinear),failed_pads=count(.!ri.checks.converged),adc_limit_pads=count(ri.checks.adc_at_limit),
    invalid_pe_buckets=count(.!pe.valid),clipped_pe_buckets=count(pe.asic_clipped),
    ions_input=sum(ions.ions.counts),ions_reconstructed=sum(Float64,ri.ions[ri.valid];init=0.),ions_input_invalid=sum(ions.ions.counts[.!ri.valid]),
    maximum_photons_per_bucket=maximum(photons.photons.counts),expected_pe=sum(ei),reconstructed_pe=sum(Float64,pe.pe[pe.valid];init=0.),expected_pe_invalid=sum(ei[.!pe.valid];init=0.),
    quiet_pad_count=length(quiet_sums),quiet_sum_squares=sum(abs2,quiet_sums),quiet_sum=sum(quiet_sums),
    predicted_quiet_sum_variance=sum(plan.ion_covariance),predicted_bucket_variance=mean(diag(plan.ion_covariance)),
    event_s=time()-event_started))
   if number<=rc["example_events"]
    chosen=unique(vcat(argmax(vec(sum(ions.ions.counts;dims=1))),something(findfirst(empty_pads),1)))
    ex=Dict{String,Any}()
    for j in chosen
     volts=adc_voltages(IonADCResponse(id,adc.channels,adc.time,adc.counts[:,j:j]),cal.adc)
     measured=[mean(volts[plan.groups[1],:])-mean(volts[g,:]) for g in plan.groups[2:end]]
     q=Float64.(ri.ions[:,j]);fit=ri.checks.converged[j] ? (ri.checks.nonlinear[j] ? SR.nonlinear_prediction(plan,q).voltage : plan.response*q) : fill(NaN,length(measured))
     ex["pad_$j"]=Dict("input"=>ions.ions.counts[:,j],"reconstructed"=>ri.ions[:,j],"time_ms"=>(collect(1:plan.time.bucket_count).-.5).*plan.time.bucket_ms,
       "measured_drop_V"=>measured,"fitted_drop_V"=>fit,"average_time_ms"=>[1000mean(plan.times[g]) for g in plan.groups[2:end]])
    end
    examples[string(id)]=ex
   end
   println("Event $id: ion prep $(round(preparation.time;digits=3)) s; ion fit $(round(ionfit.time;digits=3)) s; PE $(round(electronfit.time;digits=3)) s; nonlinear $(count(ri.checks.nonlinear)); failed $(count(.!ri.checks.converged)); invalid PE $(count(.!pe.valid))");flush(stdout)
  end
 end
 total_s=time()-started
 # Histograms, scans and five example-vector sets only: bound output independently of source size.
 budget=16*1024^2+3Base.summarysize(event_controls)+3Base.summarysize(examples)
 GuardedFiles.space_check(outdir,budget,1024^3);tmp=target*".tmp"
 try
  h5open(tmp,"w") do f
   f["schema"]="itaca.reconstruction_controls.v1";f["ion_reconstruction_method"]="linear_voltage_inverse_and_direct_nonlinear_integral";f["card"]=read(cardpath,String);f["electronics_card"]=read(joinpath(ROOT,rc["electronics_card"]),String)
   f["source_sha256"]=sha;f["detector_controls_sha256"]=bytes2hex(open(sha256,source))
   write_control_tree(create_group(f,"electron_inverse_validation"),TOML.parsefile(joinpath(ROOT,"work/reconstruction/electron_validation.toml")))
   f["setup_s"]=setup.time;f["total_s"]=total_s
   columns(create_group(f,"events"),rows)
   write_control_tree(create_group(f,"zs"),event_controls);write_control_tree(create_group(f,"examples"),examples)
   for key in ["detector_card","csa_card","sipm_card","asic_card"];f["cards/$key"]=read(joinpath(ROOT,ec[key]),String);end
   f["calibration_sha256"]=bytes2hex(open(sha256,joinpath(ROOT,ec["nonlinear_calibration"])))
   files=vcat([@__FILE__,joinpath(@__DIR__,"reconstruction_controls.jl"),joinpath(@__DIR__,"framed_electronics_settings.jl")],
    [joinpath(dir,n) for sub in ["src/reconstruction","src/electronics","src/detector_response"] for (dir,_,names) in walkdir(joinpath(ROOT,sub)) for n in names if endswith(n,".jl")])
   f["code_files"]=relpath.(files,ROOT);f["code_sha256"]=[bytes2hex(open(sha256,p)) for p in files]
   flush(f);GuardedFiles.check_size(tmp,budget)
  end
  mv(tmp,target)
 catch;isfile(tmp)&&rm(tmp);rethrow();end
 println("FINISHED: $target")
end
if abspath(PROGRAM_FILE)==@__FILE__
 run_reconstruction(isempty(ARGS) ? joinpath(ROOT,"config/readout/bb0nu_reconstruction_100.toml") : abspath(ARGS[1]);limit=length(ARGS)>1 ? parse(Int,ARGS[2]) : 0)
end
