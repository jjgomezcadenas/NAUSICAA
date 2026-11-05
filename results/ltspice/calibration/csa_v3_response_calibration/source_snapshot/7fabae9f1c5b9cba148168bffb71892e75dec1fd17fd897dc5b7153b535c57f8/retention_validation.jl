module RetentionValidation
using HDF5,TOML,SHA,Test
sha(p)=open(io->bytes2hex(sha256(io)),p)

"""Measure zero-input numerical drift in volts and calibrated ion units before setting a history tolerance."""
function run(R,root,card)
 c=R.settings(root,card);dir=joinpath(root,c["output"]);file=joinpath(dir,"baseline_validation.h5")
 if isfile(file)
  h5open(file,"r") do h
   @assert read(attributes(h)["analysis_sha256"])==sha(joinpath(dir,"analysis.h5"))
  end
 else
 @testset "Zero-input drift in recovered ion units" begin
 h5open(file,"w") do out
  attributes(out)["analysis_sha256"]=sha(joinpath(dir,"analysis.h5"))
  attributes(out)["validator_sha256"]=sha(@__FILE__)
  for variant in c["variants"]
   name=variant["name"];source=joinpath(dir,"variants",name*".h5")
   h5open(source,"r") do raw
   h5open(joinpath(dir,"analysis.h5"),"r") do analysis
    baseline=only(read(raw["op/V(out)"]))
    g=create_group(out,name);attributes(g)["source_sha256"]=sha(source)
    for frame in c["frames_s"]
     key="$(round(Int,1000frame))ms";cal=analysis[name*"/"*key]
     gain=only(read(cal["gain_V_per_ion"]));noise=only(read(cal["noise_rms_V"]))
     alpha=read(attributes(cal)["alpha"]);errors=Float64[];charge_errors=Float64[]
     for r in (1,2)
      z=raw[key*"/single_N0_r$r"];t=read(z["time"]).-c["prehistory_s"]
      voltage=read(z["V(out)"]).-baseline
      push!(errors,maximum(abs.(voltage)))
      push!(charge_errors,abs(R.A.waveform_gain(t,voltage,c["averaging_s"],frame,alpha))/gain)
     end
     bound=1e-7
     @test maximum(errors)<bound
     @test bound/gain<=.005
     @test maximum(charge_errors)<.001
     @test maximum(errors)/noise<1e-4
     target=create_group(g,key);a=attributes(target)
     a["maximum_baseline_error_V"]=maximum(errors)
     a["maximum_frame_error_e"]=maximum(charge_errors)
     a["equivalent_baseline_error_e"]=maximum(errors)/gain
     a["fraction_of_noise_rms"]=maximum(errors)/noise
     a["validated_limit_V"]=bound
     println(name," ",key," baseline_error_V=",maximum(errors)," frame_error_e=",maximum(charge_errors));flush(stdout)
    end
   end
   end
  end
 end
 end
 end
 # Keep existing strict cards unchanged; record the bounded tolerance in a new card.
 h5open(file,"r") do h
  for variant in c["variants"]
   name=variant["name"]
   error=maximum(read(attributes(h[name*"/$(round(Int,1000f))ms"])["maximum_baseline_error_V"]) for f in c["frames_s"])
   kernel_error=h5open(joinpath(dir,"variants",name*".h5"),"r") do raw
    base=only(read(raw["op/V(out)"]))
    maximum(maximum(abs.(read(raw["$(round(Int,1000f))ms/single_N0_r2/V(out)"]).-base)) for f in c["frames_s"])
   end
   if kernel_error>=1e-9
    original=joinpath(dir,"cards",name*"_viewer.toml");v=TOML.parsefile(original)
    v["baseline_validation_max_error_V"]=1e-7
    v["baseline_validation_source"]=relpath(file,root)
    v["baseline_validation_sha256"]=sha(file)
    v["worker_script"]="src/pad_history_viewer/retention_worker.jl"
    v["worker_sha256"]=sha(joinpath(root,v["worker_script"]))
    v["kernel_override_sha256"]=sha(joinpath(root,"src/ion_signal/retention_compression.jl"))
    target=joinpath(dir,"cards",name*"_validated_worker_viewer.toml")
    open(io->TOML.print(io,v),target,"w")
   end
  end
 end
end
end
