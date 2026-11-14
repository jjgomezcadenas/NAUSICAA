include("run_framed_detector_response.jl")
include("reconstruction_controls.jl")
"""Reduce reconstruction controls for summary figures without repeating simulation."""
function reduce_reconstruction(path)
 target=joinpath(dirname(path),"plot_controls.h5");tmp=target*".tmp";budget=8*1024^2
 GuardedFiles.space_check(dirname(path),budget,1024^3)
 try
  h5open(path,"r") do f
   data=Dict{String,Any}();ev=f["events"];ids=read(ev["event_number"]);n=length(ids)
   summary=Dict{String,Any}("events"=>n,"failed_pads"=>sum(read(ev["failed_pads"])),"nonlinear_pads"=>sum(read(ev["nonlinear_pads"])),
    "adc_limit_pads"=>sum(read(ev["adc_limit_pads"])),"invalid_pe_buckets"=>sum(read(ev["invalid_pe_buckets"])),"clipped_pe_buckets"=>sum(read(ev["clipped_pe_buckets"])),
    "events_with_invalid_ions"=>count(>(0),read(ev["failed_pads"])),"events_with_invalid_pe"=>count(>(0),read(ev["invalid_pe_buckets"])),
    "maximum_photons_per_bucket"=>maximum(read(ev["maximum_photons_per_bucket"])),
    "maximum_prepared_MiB"=>maximum(read(ev["prepared_bytes"]))/1024^2)
   for name in ["ion_preparation_s","ion_reconstruction_s","electron_reconstruction_s","electronics_s","detector_s","controls_s","event_s","reconstruction_allocated_bytes"]
    v=read(ev[name]);summary[name*"_mean"]=mean(v[2:end]);summary[name*"_median"]=median(v[2:end])
   end
   summary["quiet_pad_variance_ratio"]=(sum(read(ev["quiet_sum_squares"]))-sum(read(ev["quiet_sum"]))^2/sum(read(ev["quiet_pad_count"])))/sum(read(ev["quiet_pad_count"]).*read(ev["predicted_quiet_sum_variance"]))
   summary["quiet_pad_mean_total"]=sum(read(ev["quiet_sum"]))/sum(read(ev["quiet_pad_count"]))
   data["summary"]=summary
   for (arm,input,reco,invalid) in [("ions","ions_input","ions_reconstructed","ions_input_invalid"),("pe","expected_pe","reconstructed_pe","expected_pe_invalid")]
    ref=Float64.(read(ev[input]));bad=Float64.(read(ev[invalid]));reconstructed=read(ev[reco]);good=ref.-bad
    residual=(reconstructed.-good)./good.*100
    d=Dict{String,Any}("event_relative_residual_percent"=>control_histogram(residual,histogram_edges(residual,25)),
       "relative_residual_percent"=>residual,"invalid_input_fraction"=>sum(bad)/sum(ref))
    firstgroup=f["zs/$(first(ids))/$arm"];cuts=read(firstgroup["thresholds"])
    noise=zeros(length(cuts));quiet=zeros(length(cuts));loss=zeros(length(cuts));noise_n=0;quiet_n=0;input_sum=0.
    event_loss=zeros(length(cuts),n);event_noise=zeros(length(cuts),n)
    for (j,id) in enumerate(ids)
     g=f["zs/$id/$arm"];nn=read(g["noise_cells"]);qn=read(g["quiet_cells"]);total=read(g["input_total"])
     nv=read(g["noise_above"]);qv=read(g["quiet_above"]);lost=read(g["input_below"])
     noise.+=nv;quiet.+=qv;loss.+=lost;noise_n+=nn;quiet_n+=qn;input_sum+=total
     event_loss[:,j]=lost./total;event_noise[:,j]=nn>0 ? nv./nn : fill(NaN,length(cuts))
    end
    d["thresholds"]=cuts;d["noise_fraction"]=noise./noise_n;d["quiet_fraction"]=quiet./quiet_n
    d["input_lost_fraction"]=loss./input_sum;d["noise_cells_per_event"]=noise./n
    d["loss_p16"]=[quantile(row,.16) for row in eachrow(event_loss)];d["loss_p84"]=[quantile(row,.84) for row in eachrow(event_loss)]
    d["noise_cells"]=noise_n;d["quiet_cells"]=quiet_n
    data[arm]=d
   end
   exstats=Dict{String,Any}()
   for event in keys(f["examples"])
    exstats[event]=Dict{String,Any}()
    for name in keys(f["examples/$event"])
     g=f["examples/$event/$name"]
     exstats[event][name]=Dict("input_sum"=>sum(read(g["input"])),"reconstructed_sum"=>sum(read(g["reconstructed"])),
      "voltage_residual_rms_V"=>sqrt(mean(abs2,read(g["measured_drop_V"])-read(g["fitted_drop_V"]))))
    end
   end
   data["examples"]=exstats
   h5open(tmp,"w") do out
    write_control_tree(out,data);out["source_sha256"]=bytes2hex(open(sha256,path));out["code_sha256"]=bytes2hex(open(sha256,@__FILE__))
    flush(out);GuardedFiles.check_size(tmp,budget)
   end
  end
  mv(tmp,target;force=true)
 catch;isfile(tmp)&&rm(tmp);rethrow();end
end
if abspath(PROGRAM_FILE)==@__FILE__
 p=isempty(ARGS) ? joinpath(ROOT,"results/readout/bb0nu_full_ion3mm_20ms_el4mm_mainframe20mm_XeHe15_E300_tau30_seed42/reconstruction/controls.h5") : abspath(only(ARGS))
 reduce_reconstruction(p)
end
