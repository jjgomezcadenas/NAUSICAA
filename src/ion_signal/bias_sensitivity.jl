module BiasSensitivity
using HDF5,TOML,SHA
include("discharge_calibration.jl")
const D=DischargeCalibration
const A=D.A
const B=D.B
function analyze(root)
 cp=joinpath(root,"config/bias_sensitivity.toml");c=TOML.parsefile(cp);dir=joinpath(root,c["output"])
 h5open(joinpath(dir,"analysis.h5"),"w") do out
  attributes(out)["config_sha256"]=D.sha(cp);attributes(out)["analysis_sha256"]=D.sha(@__FILE__)
  frame=c["frame_s"];width=c["averaging_s"];pre=c["prehistory_s"]
  attributes(out)["frame_s"]=frame;attributes(out)["averaging_s"]=width
  h5open(joinpath(dir,"points.h5"),"r") do points
   for name in keys(points)
    path=joinpath(dir,name*".h5");isfile(path) || continue
    h5open(path,"r") do h
     g=create_group(out,name);a=attributes(g)
     for k in keys(attributes(h));a[k]=read(attributes(h)[k]);end
     a["raw_sha256"]=D.sha(path)
     a["baseline_V"]=only(read(h["op/V(out)"]))
     if haskey(h,"ac")
      ratio=read(h["ac/V(out)"])./read(h["ac/V(in)"])
      a["gain_ratio_1Hz"]=abs(first(ratio))
     end
     haskey(h,"single_N1000_r1") || return
     r=haskey(h,"single_N1000_r2") ? 2 : 1
     t,v=D.waveform(h,100,r);t.-=pre
     try
      tau,err,ts,y,yfit=D.fit_tail(t,v,c["fit_start_s"],c["fit_stop_s"],c["fit_sample_step_s"])
      a["tau_s"]=tau;a["fit_max_relative_error"]=err
      alpha=exp(-frame/tau)
      gain=A.waveform_gain(t,v,width,frame,alpha)/100
      f=read(h["noise_400/frequency"]);asd=read(h["noise_400/V(onoise)"])
      sigma=A.noise_rms(f,asd,width,frame,alpha)
      t1000,v1000=D.waveform(h,1000,r);t1000.-=pre
      signal=A.waveform_gain(t1000,v1000,width,frame,alpha)
      a["gain_V_per_ion"]=gain;a["noise_rms_V"]=sigma;a["ENC_e_rms"]=sigma/gain
      a["snr_1000"]=signal/sigma;a["gain1000_relative_change"]=signal/(1000gain)-1
      actual=read(h["single_N1000_r$r/V(out)"])
      a["minimum_output_1000_V"]=minimum(actual);a["maximum_output_1000_V"]=maximum(actual)
      a["max_input_shift_1000_V"]=maximum(abs.(read(h["single_N1000_r$r/V(in)"]).-only(read(h["op/V(in)"]))))
      q=B.integral(t1000,read(h["single_N1000_r$r/I(Is1)"]))
      a["charge_relative_error"]=q/(1000read(attributes(h)["single_ion_charge_C"]))-1
      @assert abs(read(a["charge_relative_error"]))<1e-8
      g["time_s"]=t1000;g["voltage_1000_V"]=v1000
      if r==2
       tc,vc=D.waveform(h,100,1);tc.-=pre
       a["signal_refinement_error"]=abs(A.waveform_gain(tc,vc,width,frame,alpha)/(100gain)-1)
       fc=read(h["noise_200/frequency"]);sc=read(h["noise_200/V(onoise)"])
       a["noise_refinement_error"]=abs(A.noise_rms(fc,sc,width,frame,alpha)/sigma-1)
       a["quadrature_error"]=abs(A.noise_rms(f,asd,width,frame,alpha;order=16)/sigma-1)
       @assert read(a["signal_refinement_error"])<1e-3
       @assert read(a["noise_refinement_error"])<1e-3
       @assert read(a["quadrature_error"])<1e-5
      end
      println(name," Vout=",read(a["baseline_V"])," gain=",gain*1e6," ENC=",sigma/gain," SNR1000=",signal/sigma," nonlinear=",read(a["gain1000_relative_change"]))
     catch e
      a["analysis_failure"]=sprint(showerror,e)
      println(name," calculation failed: ",sprint(showerror,e))
     end
    end
   end
  end
  for family in ("paper","v2")
   nominal=family*"_m4_6"
   haskey(out,nominal) && haskey(attributes(out[nominal]),"ENC_e_rms") || continue
   base=read(attributes(out[nominal])["ENC_e_rms"])
   gainbase=read(attributes(out[nominal])["gain_V_per_ion"])
   for name in keys(out)
    g=out[name];a=attributes(g)
    read(a["family"])==family && haskey(a,"ENC_e_rms") || continue
    a["ENC_relative_to_nominal"]=read(a["ENC_e_rms"])/base
    a["gain_relative_to_nominal"]=read(a["gain_V_per_ion"])/gainbase
    a["within_five_percent"]=read(a["ENC_relative_to_nominal"])<=1.05
    a["within_ten_percent"]=read(a["ENC_relative_to_nominal"])<=1.10
    a["useful_plateau"]=read(a["ENC_relative_to_nominal"])<=1.05 && abs(read(a["gain_relative_to_nominal"])-1)<=.10 && abs(read(a["gain1000_relative_change"]))<=.01
   end
  end
  ranges=create_group(out,"ranges")
  attributes(ranges)["criterion"]="ENC <=105% of nominal, gain within10% of nominal, gain change100 to1000 ions <=1%"
  for (family,kind) in (("paper","m3"),("paper","m4"),("v2","m4"))
   names=sort(filter(n->startswith(n,"$(family)_$(kind)_"),collect(keys(out)));by=n->read(attributes(out[n])["offset_V"]))
   center=only(findall(n->read(attributes(out[n])["offset_V"])==0,names))
   lo=center;hi=center
   while lo>1 && read(attributes(out[names[lo-1]])["useful_plateau"]);lo-=1;end
   while hi<length(names) && read(attributes(out[names[hi+1]])["useful_plateau"]);hi+=1;end
   g=create_group(ranges,"$(family)_$(kind)")
   attributes(g)["low_offset_V"]=read(attributes(out[names[lo]])["offset_V"])
   attributes(g)["high_offset_V"]=read(attributes(out[names[hi]])["offset_V"])
  end
 end
end
end
