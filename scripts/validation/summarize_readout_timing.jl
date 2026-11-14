using HDF5, Statistics
path=ARGS[1]
h5open(path,"r+") do f
 haskey(f,"summary")&&error("Summary exists")
 e=f["events"];g=create_group(f,"summary")
 names=["input_s","pairs_settings_s","ion_detector_s","electron_detector_s","electronics_preparation_s","csa_adc_s","sipm_asic_s","ion_reconstruction_preparation_s","ion_reconstruction_s","electron_reconstruction_s","electron_rebin_s","event_s","gc_s"]
 values=Dict(n=>read(e[n]) for n in names)
 values["detector_total_s"]=values["pairs_settings_s"]+values["ion_detector_s"]+values["electron_detector_s"]
 values["ion_electronics_total_s"]=values["electronics_preparation_s"]+values["csa_adc_s"]
 values["ion_reconstruction_total_s"]=values["ion_reconstruction_preparation_s"]+values["ion_reconstruction_s"]
 values["electron_reconstruction_total_s"]=values["electron_reconstruction_s"]+values["electron_rebin_s"]
 values["other_s"]=values["event_s"]-sum(values[n] for n in names[1:11])
 for (n,v) in sort(collect(values))
  s=create_group(g,n)
  for (key,x) in [("mean",mean(v)),("median",median(v)),("p90",quantile(v,.9)),("maximum",maximum(v)),("total",sum(v)),("first_half_mean",mean(v[1:50])),("second_half_mean",mean(v[51:100]))];s[key]=x;end
  println("$n: mean $(round(1000mean(v);digits=3)) ms; median $(round(1000median(v);digits=3)) ms; p90 $(round(1000quantile(v,.9);digits=3)) ms")
 end
 for n in ["failed_pads","nonlinear_pads","invalid_pe"];g[n*"_total"]=sum(read(e[n]));end
end
