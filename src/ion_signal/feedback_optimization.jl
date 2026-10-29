module FeedbackOptimization
using HDF5,TOML,SHA,Test
include("discharge_calibration.jl")
const D=DischargeCalibration
const B=D.B
const A=D.A
sha(p)=bytes2hex(sha256(read(p)))
"""Extend the noiseless pre-pulse baseline for centred acquisition windows."""
function wave(group,n,r,pre,widthmax)
 t,v=D.waveform(group,n,r);t.-=pre
 if first(t)>-widthmax/2
  @assert first(v)==0.
  t=vcat(-widthmax/2-.0001,t);v=vcat(0.,v)
 end
 t,v
end
function analyze(root;names=nothing)
 cp=joinpath(root,"config/feedback_optimization.toml");c=TOML.parsefile(cp);dir=joinpath(root,c["output"])
 h5open(joinpath(dir,"analysis.h5"),"w") do out
  attributes(out)["config_sha256"]=sha(cp);attributes(out)["analysis_sha256"]=sha(@__FILE__)
  attributes(out)["assumptions"]="Same central 0.5 mm square-pad benchmark and 1.25 pF external load as paper; finite 5 mm path, 3 m/s, uniform mean arrivals. Unit pad weighting potential; surrounding infinite plane and infinity at zero; homogeneous gas. Stationary 0.1 Hz to 1 MHz noise. Only feedback dimensions and gate bias varied. No charge renormalization. Signal-free prehistory extended as zero for centred windows longer than the archived prehistory."
  names_all=vcat([v["name"] for v in c["variants"]],sort([splitext(f)[1] for f in readdir(dir) if endswith(f,"_retuned.h5") || endswith(f,"_rounded.h5")]))
  for name in names_all
   names===nothing || name in names || continue
   path=joinpath(dir,name*".h5");isfile(path) || continue
   h5open(path,"r") do raw
    vg=create_group(out,name);at=attributes(vg)
    for k in keys(attributes(raw));at[k]=read(attributes(raw)[k]);end
    at["raw_sha256"]=sha(path)
    if haskey(at,"rejected_operating_point");println(name," rejected DC");return;end
    haskey(raw,"40ms/single_N2000_r2") || (at["incomplete"]=true;return)
    f=read(raw["noise_400/frequency"]);asd=read(raw["noise_400/V(onoise)"])
    fc=read(raw["noise_200/frequency"]);ac=read(raw["noise_200/V(onoise)"])
    at["power_W"]=-1.8*only(read(raw["op/I(Vdd)"]))
    for (frame,meanions,inputfile) in zip(c["frames_s"],c["mean_ions"],c["input_sources"])
     key="$(round(Int,frame*1000))ms";g=create_group(vg,key);a=attributes(g)
     a["frame_s"]=frame;a["mean_ions"]=meanions
     src=raw[key];pre=c["prehistory_s"];maxwidth=maximum(c["averaging_s"])
     t100,v100=wave(src,100,2,pre,maxwidth);tco,vco=wave(src,100,1,pre,maxwidth)
     t10,v10=wave(src,10,2,pre,maxwidth)
     try
      tau,fiterr,_,_,_=D.fit_tail(t100,v100,frame+.001,.12,.0002)
      tauc,_,_,_,_=D.fit_tail(tco,vco,frame+.001,.12,.0002)
      tau10,_,_,_,_=D.fit_tail(t10,v10,frame+.001,.12,.0002)
      for (k,v) in ("tau_s"=>tau,"fit_error"=>fiterr,"tau_refinement_error"=>abs(tauc/tau-1),"tau_10_vs_100_error"=>abs(tau10/tau-1));a[k]=v;end
      alpha=exp(-frame/tau);a["alpha"]=alpha
      tm,vm=wave(src,meanions,2,pre,maxwidth);tmc,vmc=wave(src,meanions,1,pre,maxwidth)
      tmax,vmax=wave(src,2000,2,pre,maxwidth)
      a["minimum_output_2000_V"]=minimum(read(src["single_N2000_r2/V(out)"]))
      a["maximum_output_2000_V"]=maximum(read(src["single_N2000_r2/V(out)"]))
      gain=[];noise=[];enc=[];snr=[];transient=[];spectral=[];quadrature=[];linearity=[];prev=[];next=[];sig=[]
      for width in c["averaging_s"]
       measure(t,v;shift=0.)=A.waveform_gain(t.-shift,v,width,frame,alpha)
       gg=measure(t100,v100)/100;ss=A.noise_rms(f,asd,width,frame,alpha)
       sm=measure(tm,vm)
       push!(gain,gg);push!(noise,ss);push!(enc,ss/gg);push!(snr,sm/ss);push!(sig,sm)
       push!(transient,max(abs(measure(tco,vco)/(100gg)-1),abs(measure(tmc,vmc)/sm-1)))
       push!(spectral,abs(A.noise_rms(fc,ac,width,frame,alpha)/ss-1))
       push!(quadrature,abs(A.noise_rms(f,asd,width,frame,alpha;order=16)/ss-1))
       push!(linearity,max(abs(measure(tmax,vmax)/(2000gg)-1),abs(sm/(meanions*gg)-1)))
       push!(prev,measure(t100,v100;shift=frame)/(100gg))
       push!(next,-D.average(t100,v100,0.,width)/(100gg))
      end
      for (k,v) in ("width_s"=>c["averaging_s"],"gain_V_per_ion"=>gain,"noise_rms_V"=>noise,"ENC_e_rms"=>enc,"mean_snr"=>snr,"mean_signal_V"=>sig,"transient_error"=>transient,"spectral_error"=>spectral,"quadrature_error"=>quadrature,"max_gain_change"=>linearity,"previous_frame_fraction"=>prev,"following_frame_fraction"=>next)
       g[k]=Float64.(v)
      end
      h5open(joinpath(root,inputfile),"r") do inp
       q=read(attributes(inp)["single_ion_charge_C"]);errors=Float64[]
       for r in (1,2),n in filter(>(0),c["counts"])
        d=src["single_N$(n)_r$r"]
        push!(errors,abs(B.integral(read(d["time"]),read(d["I(Is1)"]))/(n*q)-1))
       end
       a["charge_error"]=maximum(errors)
      end
      valid=[gain[j]>0 && max(transient[j],spectral[j])<1e-3 && quadrature[j]<1e-5 && linearity[j]<c["max_gain_change"] && max(abs(prev[j]),abs(next[j]))<c["max_frame_leakage_fraction"] && read(a["tau_refinement_error"])<1e-3 && read(a["tau_10_vs_100_error"])<.01 && read(a["charge_error"])<2e-5 && read(a["minimum_output_2000_V"])>c["minimum_output_V"] && read(a["maximum_output_2000_V"])<c["maximum_output_V"] for j in eachindex(gain)]
      g["valid"]=Int.(valid)
      passing=findall(valid)
      if !isempty(passing)
       j=passing[argmin(Float64.(enc[passing]))];a["best_index"]=j
       println(name," ",key," tau=",tau," best width=",c["averaging_s"][j]," ENC=",enc[j]," mean SNR=",snr[j]," leakage=",next[j])
       device_names=sort(filter(n->occursin(r"^V\(m\w+\)$",n),keys(raw["noise_400"])))
       vars=[A.noise_rms(f,max.(read(raw["noise_400"][n]),1e-150),c["averaging_s"][j],frame,alpha)^2 for n in device_names]
       g["device_names"]=device_names;g["device_variance_fraction"]=vars./noise[j]^2;g["device_enc_e"]=sqrt.(vars)./gain[j]
       a["noise_variance_closure_error"]=abs(sum(vars)/noise[j]^2-1)
      else
       a["no_valid_setting"]=true;println(name," ",key," no valid setting; ENC=",enc," gainChange=",linearity," fitErr=",fiterr," tauChange=",tau10/tau-1)
      end
      g["time_s"]=tm;g["mean_voltage_V"]=vm;g["mean_current_A"]=read(src["single_N$(meanions)_r2/I(Is1)"])
      # Extended time contains one extra prehistory point; keep the acquisition time for current.
      g["current_time_s"]=read(src["single_N$(meanions)_r2/time"]).-pre
     catch e
      a["analysis_failure"]=sprint(showerror,e);println(name," ",key," failed: ",sprint(showerror,e))
     end
    end
   end
  end
 end
end
function checks(root)
 c=TOML.parsefile(joinpath(root,"config/feedback_optimization.toml"));dir=joinpath(root,c["output"])
 @testset "Feedback and averaging optimization" begin
  h5open(joinpath(dir,"analysis.h5"),"r") do h
   for name in keys(h)
    vg=h[name];a=attributes(vg)
    if haskey(a,"acquisition_failure")
     @test haskey(a,"incomplete")
     continue
    end
    if haskey(a,"rejected_operating_point");@test !(c["minimum_output_V"]<read(a["baseline_V"])<c["maximum_output_V"]);continue;end
    @test !haskey(a,"incomplete")
    for key in ("20ms","40ms")
     g=vg[key];at=attributes(g)
     @test !haskey(at,"analysis_failure")
     @test read(at["charge_error"])<2e-5
     @test read(at["tau_refinement_error"])<1e-3
     @test read(at["tau_10_vs_100_error"])<.01
     @test maximum(read(g["transient_error"]))<1e-3
     @test maximum(read(g["spectral_error"]))<1e-3
     @test maximum(read(g["quadrature_error"]))<1e-5
     if haskey(at,"best_index")
      j=read(at["best_index"])
      @test read(g["valid"])[j]==1
      @test read(g["ENC_e_rms"])[j]≈read(g["noise_rms_V"])[j]/read(g["gain_V_per_ion"])[j]
      @test read(at["noise_variance_closure_error"])<1e-3
     end
    end
   end
   for (key,path) in (("20ms","results/ltspice/discharge_calibration/20ms/paper_readout.h5"),("40ms","results/ltspice/pitch_4mm_40ms/analysis.h5"))
    h5open(joinpath(root,path),"r") do old
     g=h["reference/$key"]
     @test first(read(g["ENC_e_rms"]))≈read(attributes(old)["ENC_e_rms"]) rtol=.002
    end
   end
  end
 end
end
end
