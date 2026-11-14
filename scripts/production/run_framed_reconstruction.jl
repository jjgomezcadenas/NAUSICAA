#!/usr/bin/env julia
include("run_framed_electronics.jl")
include("../../src/reconstruction/SignalReconstruction.jl")
using .SignalReconstruction
include("reconstruction_controls.jl")
include("track_controls.jl")
const SR=SignalReconstruction
"""Run the same accepted events through electronics and reconstruction in memory.
Save per-event numerical ZS controls and selected vectors, without event products
or suppression. The source detector controls fix event positions and frames."""
function run_reconstruction(cardpath;limit=0)
 BLAS.set_num_threads(1)
 rc=TOML.parsefile(cardpath);ec=TOML.parsefile(joinpath(ROOT,rc["electronics_card"]))
 dc=TOML.parsefile(joinpath(ROOT,ec["detector_card"]))
 base=joinpath(ROOT,"results/readout",dc["run_label"])
 source=haskey(rc,"selection_controls") ? joinpath(ROOT,rc["selection_controls"]) : joinpath(base,"detector_response/controls.h5")
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
 track_mode=haskey(rc,"track_analysis");track_results=Dict{String,Any}();track_examples=Dict{String,Any}()
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
    if track_mode
     charge_estimates=event_charge_estimates(ri,pe,p,ds,s.sipm.sipm.pde)
     ihits=reconstruct_hits(ri,p,ds;keV_per_count=charge_estimates.ion_keV_per_count)
     ehits=reconstruct_hits(pe,p,ds;charge_correction=charge_estimates.attachment_correction,keV_per_count=charge_estimates.pe_keV_per_count)
     tc=rc["track_analysis"];edges=collect(0.:tc["hit_energy_bin_keV"]:tc["hit_energy_max_keV"])
     ion_cut=get(tc,"ion_hit_cut_keV",0.);electron_cut=get(tc,"electron_hit_cut_keV",0.)
     cut_stats=Dict("ions"=>Dict{String,Any}(string(k)=>v for (k,v) in Base.pairs(hit_cut_statistics(ri.ions,ri.valid,ions.ions.counts,charge_estimates.ion_keV_per_count,ion_cut))),
      "pe"=>Dict{String,Any}(string(k)=>v for (k,v) in Base.pairs(hit_cut_statistics(pe.pe,pe.valid,ei,charge_estimates.attachment_correction*charge_estimates.pe_keV_per_count,electron_cut))))
     ihits_before=ihits;ehits_before=ehits
     if get(tc,"check_connectivity",false)
      cut_stats["ions"]["topology"]=hit_cut_topology(ri.ions,ri.valid,ions.ions.counts,ri.channels,ihits_before,charge_estimates.ion_keV_per_count,ion_cut,ds.ion_transport.sigma_t_mm,get(tc,"ion_graph_k",tc["graph_k"]),get(tc,"ion_edge_sigma_scale",tc["edge_sigma_scale"]))
      cut_stats["pe"]["topology"]=hit_cut_topology(pe.pe,pe.valid,ei,pe.channels,ehits_before,charge_estimates.attachment_correction*charge_estimates.pe_keV_per_count,electron_cut,ds.electron_transport.sigma_t_mm,get(tc,"electron_graph_k",tc["graph_k"]),tc["edge_sigma_scale"])
     end
     if haskey(tc,"graph_k_scan")
      for arm in ("ions","pe");cut_stats[arm]["topology_scan"]=Dict{String,Any}();end
      for graph_k in tc["graph_k_scan"]
       if graph_k==get(tc,"ion_graph_k",tc["graph_k"])
        cut_stats["ions"]["topology_scan"]["k_$graph_k"]=cut_stats["ions"]["topology"]
       else
        cut_stats["ions"]["topology_scan"]["k_$graph_k"]=hit_cut_topology(ri.ions,ri.valid,ions.ions.counts,ri.channels,ihits_before,charge_estimates.ion_keV_per_count,ion_cut,ds.ion_transport.sigma_t_mm,graph_k,get(tc,"ion_edge_sigma_scale",tc["edge_sigma_scale"]);compute_endpoints=false)
       end
       if graph_k==get(tc,"electron_graph_k",tc["graph_k"])
        cut_stats["pe"]["topology_scan"]["k_$graph_k"]=cut_stats["pe"]["topology"]
       else
        cut_stats["pe"]["topology_scan"]["k_$graph_k"]=hit_cut_topology(pe.pe,pe.valid,ei,pe.channels,ehits_before,charge_estimates.attachment_correction*charge_estimates.pe_keV_per_count,electron_cut,ds.electron_transport.sigma_t_mm,graph_k,tc["edge_sigma_scale"];compute_endpoints=false)
       end
      end
     end
     ihits=reconstruct_hits(ri,p,ds;keV_per_count=charge_estimates.ion_keV_per_count,minimum_energy_keV=ion_cut)
     ehits=reconstruct_hits(pe,p,ds;charge_correction=charge_estimates.attachment_correction,keV_per_count=charge_estimates.pe_keV_per_count,minimum_energy_keV=electron_cut)
     track_results[string(id)]=Dict("ions"=>control_histogram(ihits.energy_keV,edges),"pe"=>control_histogram(ehits.energy_keV,edges),
      "ions_before"=>control_histogram(ihits_before.energy_keV,edges),"pe_before"=>control_histogram(ehits_before.energy_keV,edges),"cut_statistics"=>cut_stats)
     if get(tc,"save_all_track_hits",false) || number%tc["example_every"]==0
      track_examples[string(id)]=Dict("ordinal"=>number,
       "itrack"=>hit_projections(ihits,dc["ion_pitch_mm"],dc["ion_bucket_ms"]*dc["ion_drift_velocity_mm_ms"];z_centres=SR.hit_z_centres(ri,p,ds)),
       "etrack"=>hit_projections(ehits,dc["EL_pitch_mm"],dc["electron_bucket_us"]*factor*dc["electron_drift_velocity_mm_us"];z_centres=SR.hit_z_centres(pe,p,ds)))
      if get(tc,"save_all_track_hits",false)
       for (arm,hits) in (("itrack",ihits),("etrack",ehits))
        track_examples[string(id)][arm]["hits"]=Dict("x_mm"=>hits.x_mm,"y_mm"=>hits.y_mm,"z_mm"=>hits.z_mm,"energy_keV"=>hits.energy_keV)
       end
      end
     end
    else
     ci=channel_histograms(ri.ions,ri.valid,ions.ions.counts,zs["histogram_bins"])
     ce=channel_histograms(pe.pe,pe.valid,ei,zs["histogram_bins"])
     event_controls[string(id)]=Dict("ions"=>ci,"pe"=>ce)
    end
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
   if track_mode
    rows[end]=merge(rows[end],charge_estimates,(itrack_hits=length(ihits.charge),etrack_hits=length(ehits.charge)))
   end
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
 track_summary=Dict{String,Any}()
 if track_mode
  for key in ("ions_total","pe_total","ion_energy_keV","electron_energy_keV")
   values=filter(isfinite,[getproperty(r,Symbol(key)) for r in rows])
   track_summary[key]=control_histogram(values,histogram_edges(values,80))
  end
  for arm in ("ions","pe");track_summary[arm*"_hit_energy_keV"]=combine_histograms([track_results[string(id)][arm] for id in ids]);end
 end
 budget=16*1024^2+3Base.summarysize(event_controls)+3Base.summarysize(examples)+3Base.summarysize(track_results)+3Base.summarysize(track_examples)
 GuardedFiles.space_check(outdir,budget,1024^3);tmp=target*".tmp"
 try
  h5open(tmp,"w") do f
   f["schema"]=track_mode ? "itaca.track_controls.v1" : "itaca.reconstruction_controls.v3";f["ion_reconstruction_method"]="linear_voltage_inverse_and_direct_nonlinear_integral";f["card"]=read(cardpath,String);f["electronics_card"]=read(joinpath(ROOT,rc["electronics_card"]),String)
   f["source_sha256"]=sha;f["detector_controls_sha256"]=bytes2hex(open(sha256,source))
   write_control_tree(create_group(f,"electron_inverse_validation"),TOML.parsefile(joinpath(ROOT,"work/reconstruction/electron_validation.toml")))
   f["setup_s"]=setup.time;f["total_s"]=total_s
   columns(create_group(f,"events"),rows)
   if track_mode
    write_control_tree(create_group(f,"track_summary"),track_summary)
    write_control_tree(create_group(f,"hit_histograms"),track_results)
    write_control_tree(create_group(f,"track_examples"),track_examples)
   end
   write_control_tree(create_group(f,"zs"),event_controls);write_control_tree(create_group(f,"examples"),examples)
   for key in ["detector_card","csa_card","sipm_card","asic_card"];f["cards/$key"]=read(joinpath(ROOT,ec[key]),String);end
   f["calibration_sha256"]=bytes2hex(open(sha256,joinpath(ROOT,ec["nonlinear_calibration"])))
   files=vcat([@__FILE__,joinpath(@__DIR__,"reconstruction_controls.jl"),joinpath(@__DIR__,"track_controls.jl"),joinpath(@__DIR__,"framed_electronics_settings.jl")],
    [joinpath(dir,n) for sub in ["src/reconstruction","src/electronics","src/detector_response"] for (dir,_,names) in walkdir(joinpath(ROOT,sub)) for n in names if endswith(n,".jl")])
   f["code_files"]=relpath.(files,ROOT);f["code_sha256"]=[bytes2hex(open(sha256,p)) for p in files]
   flush(f);GuardedFiles.check_size(tmp,budget)
  end
  mv(tmp,target)
 catch;isfile(tmp)&&rm(tmp);rethrow();end
 println("FINISHED: $target")
end
if abspath(PROGRAM_FILE)==@__FILE__
 run_reconstruction(isempty(ARGS) ? joinpath(ROOT,"config/readout/bb0nu_reconstruction_bucket_histograms_100.toml") : abspath(ARGS[1]);limit=length(ARGS)>1 ? parse(Int,ARGS[2]) : 0)
end
