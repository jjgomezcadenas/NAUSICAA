module PitchExercise
using HDF5,TOML,SHA,Test
include("discharge_calibration.jl")
const D=DischargeCalibration
const B=D.B
const A=D.A
sha(p)=bytes2hex(sha256(read(p)))
settings(root)=TOML.parsefile(joinpath(root,"config/pitch_4mm_40ms.toml"))
function inputs(root)
 c=settings(root);out=joinpath(root,c["output"]);mkpath(out);file=joinpath(out,"input.h5")
 isfile(file) && error("Preserving existing input: $file")
 d=TOML.parsefile(joinpath(root,c["detector"]));a=d["geometry"]["pad_width_m"]/2;g=d["geometry"]["gap_m"]
 v=d["gas"]["bulk_speed_m_s"]/d["bias"]["bulk_field_V_m"]*(d["bias"]["ic_voltage_V"]-d["bias"]["pad_voltage_V"])/g
 h5open(file,"w") do h
  at=attributes(h)
  for (k,x) in ("config_sha256"=>sha(joinpath(root,"config/pitch_4mm_40ms.toml")),"detector_sha256"=>sha(joinpath(root,c["detector"])),"generator_sha256"=>sha(@__FILE__),"waveform_source_sha256"=>sha(joinpath(@__DIR__,"bucket_study.jl")),"single_ion_charge_C"=>B.E*(1-B.phi(g,a)),"pad_width_m"=>2a,"gap_m"=>g,"gap_speed_m_s"=>v,"pitch_m"=>c["pitch_m"],"frame_s"=>c["frame_s"],"mean_ions"=>c["mean_ions"])
   at[k]=x
  end
  at["assumptions"]="0.5 mm central square pad in infinite grounded plane; homogeneous gas half-space; weighting potential one on pad, zero on surrounding plane and at infinity. Central straight approach through 5 mm at constant 3 m/s. Uniform mean arrivals over 40 ms. Finite-path induced charge retained. Same pad, external capacitance and CSA as the 20 ms paper. Pitch 4 mm and mean 400 ions are user design inputs; focusing and event charge redistribution are not recalculated."
  for (r,step) in enumerate(c["max_steps_s"])
   t,i=B.pwl(c["frame_s"],a,g,v,c["prehistory_s"],c["stop_s"],step)
   group=create_group(h,"r$r");group["time_s"]=t;group["single_A_per_ion"]=i
  end
 end
end
function analyze(root)
 c=settings(root);out=joinpath(root,c["output"]);frame=c["frame_s"];width=c["averaging_s"];pre=c["prehistory_s"]
 h5open(joinpath(out,"raw.h5"),"r") do raw
  @assert read(attributes(raw)["config_sha256"])==sha(joinpath(root,"config/pitch_4mm_40ms.toml"))
  h5open(joinpath(out,"analysis.h5"),"w") do h
   at=attributes(h)
   for k in keys(attributes(raw));at[k]=read(attributes(raw)[k]);end
   at["raw_sha256"]=sha(joinpath(out,"raw.h5"));at["analysis_sha256"]=sha(@__FILE__)
   at["filter_sha256"]=sha(joinpath(@__DIR__,"averaging_scan.jl"))
   at["reference_sha256"]=sha(joinpath(root,c["reference"]))
   for k in ("pitch_m","frame_s","averaging_s","mean_ions","reference_mean_ions");at[k]=c[k];end
   fit(n,r;lo=c["fit_start_s"],hi=c["fit_stop_s"])=begin
    t,v=D.waveform(raw,n,r);D.fit_tail(t.-pre,v,lo,hi,c["fit_sample_step_s"])
   end
   tau,err,ts,y,yfit=fit(100,2);tauc,_,_,_,_=fit(100,1);tau10,_,_,_,_=fit(10,2)
   localtaus=[first(fit(100,2;lo=x[1],hi=x[2])) for x in c["fit_checks_s"]]
   alpha=exp(-frame/tau)
   signal(t,v;shift=0.)=A.waveform_gain(t.-pre.-shift,v,width,frame,alpha)
   t,v=D.waveform(raw,1,2);tc,vc=D.waveform(raw,1,1);gain=signal(t,v)
   f=read(raw["noise_400/frequency"]);asd=read(raw["noise_400/V(onoise)"])
   fc=read(raw["noise_200/frequency"]);ac=read(raw["noise_200/V(onoise)"])
   sigma=A.noise_rms(f,asd,width,frame,alpha)
   terr=abs(signal(tc,vc)/gain-1);serr=abs(A.noise_rms(fc,ac,width,frame,alpha)/sigma-1)
   qerr=abs(A.noise_rms(f,asd,width,frame,alpha;order=16)/sigma-1)
   for (k,x) in ("tau_s"=>tau,"alpha"=>alpha,"tail_fit_relative_error"=>err,"tau_timestep_error"=>abs(tauc/tau-1),"tau_10_vs_100_error"=>abs(tau10/tau-1),"transient_error"=>terr,"spectral_error"=>serr,"quadrature_error"=>qerr,"gain_V_per_ion"=>gain,"noise_rms_V"=>sigma,"ENC_e_rms"=>sigma/gain,"previous_frame_fraction"=>signal(t,v;shift=frame)/gain,"following_frame_fraction"=>-D.average(t.-pre,v,0.,width)/gain)
    at[k]=x
   end
   ns=filter(>(1),c["counts"]);snrs=Float64[];gains=Float64[];charges=Float64[];errors=Float64[]
   for n in ns
    tn,vn=D.waveform(raw,n,2);tnc,vnc=D.waveform(raw,n,1);s=signal(tn,vn)
    push!(snrs,s/sigma);push!(gains,s/(n*gain));push!(charges,B.integral(tn,read(raw["single_N$(n)_r2/I(Is1)"])))
    push!(errors,abs(signal(tnc,vnc)/s-1))
   end
   for (k,x) in ("counts"=>ns,"snr"=>snrs,"gain_relative_to_small_signal"=>gains,"input_charge_C"=>charges,"finite_count_refinement_error"=>errors,"local_tau_s"=>localtaus,"tail_time_s"=>ts,"tail_magnitude_V"=>y,"tail_fit_V"=>yfit)
    h[k]=x
   end
   n=c["mean_ions"];j=findfirst(==(n),ns);tn,vn=D.waveform(raw,n,2)
   at["mean_signal_V"]=snrs[j]*sigma;at["mean_snr"]=snrs[j];at["mean_input_charge_C"]=charges[j]
   current=read(raw["single_N$(n)_r2/I(Is1)"]);at["mean_peak_current_A"]=maximum(current)
   h["time_s"]=tn.-pre;h["voltage_V"]=vn;h["current_A"]=current
   h["snr_targets"]=c["snr_targets"];h["linear_ions_at_target"]=c["snr_targets"].*(sigma/gain)
   h["frequency_Hz"]=f;h["output_asd_V_sqrtHz"]=asd
   names=sort(filter(x->occursin(r"^V\(m\w+\)$",x),keys(raw["noise_400"])))
   vars=[A.noise_rms(f,max.(read(raw["noise_400"][name]),1e-150),width,frame,alpha)^2 for name in names]
   h["device_names"]=names;h["device_variance_fraction"]=vars./sigma^2
   at["noise_variance_closure_error"]=abs(sum(vars)/sigma^2-1)
   ref=create_group(h,"reference_20ms")
   h5open(joinpath(root,c["reference"]),"r") do old
    for k in ("frame_s","gain_V_per_ion","noise_rms_V","ENC_e_rms","tau_s");attributes(ref)[k]=read(attributes(old)[k]);end
    gg=read(attributes(old)["gain_V_per_ion"]);ss=read(attributes(old)["noise_rms_V"])
    attributes(ref)["mean_ions"]=c["reference_mean_ions"];attributes(ref)["mean_signal_V"]=c["reference_mean_ions"]*gg
    attributes(ref)["mean_snr"]=c["reference_mean_ions"]*gg/ss
   end
   @assert gain>0 && sigma>0 && max(terr,serr,maximum(errors))<1e-3 && qerr<1e-5
   @assert maximum(abs.(charges./(ns.*read(attributes(raw)["single_ion_charge_C"])).-1))<2e-5
   @assert abs(tauc/tau-1)<1e-3 && abs(tau10/tau-1)<.01
   @assert abs(sum(vars)/sigma^2-1)<1e-3
   println("40 ms: gain=",gain," noise=",sigma," ENC=",sigma/gain," mean S/N=",snrs[j]," tau=",tau)
  end
 end
end
function checks(root)
 c=settings(root);out=joinpath(root,c["output"])
 @testset "4 mm pitch / 40 ms exercise" begin
  h5open(joinpath(out,"analysis.h5"),"r") do h
   a=attributes(h);get(k)=read(a[k]);q=get("single_ion_charge_C")
   @test get("frame_s")==.04
   @test get("pitch_m")==.004
   @test get("mean_ions")==400
   @test get("tail_fit_relative_error")<1e-4
   @test get("tau_timestep_error")<1e-3
   @test get("tau_10_vs_100_error")<.01
   for k in ("transient_error","spectral_error","noise_variance_closure_error");@test get(k)<1e-3;end
   @test get("quadrature_error")<1e-5
   @test maximum(read(h["finite_count_refinement_error"]))<1e-3
   @test maximum(abs.(read(h["input_charge_C"])./(read(h["counts"]).*q).-1))<2e-5
   @test maximum(abs.(read(h["gain_relative_to_small_signal"]).-1))<.01
   @test maximum(abs.(read(h["local_tau_s"])./get("tau_s").-1))<1e-3
   @test get("mean_snr")≈get("mean_signal_V")/get("noise_rms_V")
   @test get("ENC_e_rms")≈get("noise_rms_V")/get("gain_V_per_ion")
   h5open(joinpath(out,"input.h5"),"r") do inp
    for r in (1,2)
     t=read(inp["r$r/time_s"]);i=read(inp["r$r/single_A_per_ion"])
     @test B.integral(t,i)≈q rtol=1e-9
     # Independent uniform-time quadrature of the analytical current.
     dt=2e-7;tt=collect(-.002:dt:.042)
     vals=B.current.(tt,.04,.00025,.005,3.)
     @test B.integral(tt,vals)≈q rtol=1e-5
    end
   end
  end
 end
end
end
