using HDF5,Statistics
include("../../src/ion_signal/guarded_files.jl")
path=ARGS[1];out=joinpath(dirname(path),"noise_summary.h5")
isfile(out)&&error("Preserve existing noise summary")
budget=1024^2
GuardedFiles.space_check(dirname(out),budget,1024^3)
tmp=out*".tmp"
try
h5open(path,"r") do f
 e=f["events"];energy=read(e["ion_energy_keV"]);good=isfinite.(energy)
 W=read(e["ion_keV_per_count"]);input=read(e["ions_input"]).*W
 variance=read(e["predicted_quiet_sum_variance"]);channels=read(e["ion_channels"])
 h5open(tmp,"w") do o
  o["input_ion_energy_std_keV"]=std(input;corrected=false)
  o["ion_residual_std_keV"]=std(energy[good]-input[good];corrected=false)
  o["predicted_linear_noise_std_keV"]=sqrt(mean((variance.*channels.*W.^2)[good]))
  o["mean_ion_channels"]=mean(channels);o["mean_ion_buckets"]=mean(read(e["ion_buckets"]))
  o["empty_pad_sum_rms_ions"]=sqrt(sum(read(e["quiet_sum_squares"]))/sum(read(e["quiet_pad_count"])))
  o["predicted_empty_pad_sum_rms_ions"]=sqrt(mean(variance))
 end
end

GuardedFiles.check_size(tmp,budget)
mv(tmp,out)
catch
 isfile(tmp)&&rm(tmp)
 rethrow()
end
