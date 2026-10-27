module PaperReadout
using HDF5,TOML,SHA
include("discharge_calibration.jl")
const D=DischargeCalibration
function analyze(root)
 cfg,configpath=D.configuration(root);dir=joinpath(root,cfg["output"])
 rawpath=joinpath(dir,"raw.h5");analysispath=joinpath(dir,"analysis.h5")
 h5open(rawpath,"r") do raw
  h5open(analysispath,"r") do saved
   @assert read(attributes(saved)["raw_sha256"])==D.sha(rawpath)
   @assert read(attributes(saved)["config_sha256"])==D.sha(configpath)
   frame=cfg["frame_s"];width=cfg["averaging_s"];pre=cfg["prehistory_s"]
   calibrated=saved["calibrated"];tau=read(attributes(calibrated)["tau_s"]);alpha=exp(-frame/tau)
   t,v=D.waveform(raw,1,2);tc,vc=D.waveform(raw,1,1);t.-=pre;tc.-=pre
   f=read(raw["noise_400/frequency"]);asd=read(raw["noise_400/V(onoise)"])
   fc=read(raw["noise_200/frequency"]);ac=read(raw["noise_200/V(onoise)"])
   widths=[180.,200.,250.,300.,350.].*1e-6
   gains=[D.A.waveform_gain(t,v,w,frame,alpha) for w in widths]
   sigmas=[D.A.noise_rms(f,asd,w,frame,alpha) for w in widths]
   qerr=[abs(D.A.noise_rms(f,asd,w,frame,alpha;order=16)/s-1) for (w,s) in zip(widths,sigmas)]
   serr=[abs(D.A.noise_rms(fc,ac,w,frame,alpha)/s-1) for (w,s) in zip(widths,sigmas)]
   terr=[abs(D.A.waveform_gain(tc,vc,w,frame,alpha)/g-1) for (w,g) in zip(widths,gains)]
   @assert maximum(qerr)<1e-6 && maximum(serr)<1e-3 && maximum(terr)<1e-3
   index=only(findall(w->isapprox(w,width;rtol=1e-12),widths))
   @assert isapprox(gains[index],read(attributes(calibrated)["gain_V_per_ion"]);rtol=1e-10)
   @assert isapprox(sigmas[index],read(attributes(calibrated)["noise_rms_V"]);rtol=1e-10)
   h5open(joinpath(dir,"paper_readout.h5"),"w") do out
    a=attributes(out)
    for k in keys(attributes(calibrated));a[k]=read(attributes(calibrated)[k]);end
    a["frame_s"]=frame;a["averaging_s"]=width;a["frequency_low_Hz"]=first(f);a["frequency_high_Hz"]=last(f)
    a["raw_sha256"]=D.sha(rawpath);a["calibration_sha256"]=D.sha(analysispath)
    a["analysis_sha256"]=D.sha(@__FILE__);a["config_sha256"]=D.sha(configpath)
    a["filter_sha256"]=D.sha(joinpath(@__DIR__,"averaging_scan.jl"))
    a["assumptions"]=read(attributes(raw)["assumptions"])
    a["calibration_ions"]=cfg["calibration_ions"];a["fit_start_s"]=cfg["fit_start_s"];a["fit_stop_s"]=cfg["fit_stop_s"]
    out["time_s"]=t;out["voltage_V_per_ion"]=v;out["current_A_per_ion"]=read(raw["single_N1_r2/I(Is1)"])
    out["sample_time_s"]=[0.,frame];out["sample_voltage_V_per_ion"]=[D.average(t,v,x,width) for x in (0.,frame)]
    out["frequency_Hz"]=f;out["output_asd_V_sqrtHz"]=asd
    out["counts"]=read(calibrated["counts"]);out["snr"]=read(calibrated["snr"])
    scan=create_group(out,"averaging")
    scan["width_s"]=widths;scan["gain_V_per_ion"]=gains;scan["noise_rms_V"]=sigmas;scan["ENC_e_rms"]=sigmas./gains
    scan["quadrature_error"]=qerr;scan["spectral_error"]=serr;scan["transient_error"]=terr
    attributes(scan)["selected_index"]=index
    names=sort(filter(n->occursin(r"^V\(m\w+\)$",n),keys(raw["noise_400"])))
    vars=[D.A.noise_rms(f,max.(read(raw["noise_400"][n]),1e-150),width,frame,alpha)^2 for n in names]
    @assert abs(sum(vars)/sigmas[index]^2-1)<1e-3
    out["device_names"]=names;out["device_variance_fraction"]=vars./sigmas[index]^2
    flickers=filter(n->occursin("1overf",n),keys(raw["noise_400"]))
    @assert !isempty(flickers)
    flickervar=sum(D.A.noise_rms(f,max.(read(raw["noise_400"][n]),1e-150),width,frame,alpha)^2 for n in flickers)
    a["flicker_variance_fraction"]=flickervar/sigmas[index]^2
    @assert read(a["flicker_variance_fraction"])>0
   end
   for (w,g,s) in zip(widths,gains,sigmas)
    println(w*1e6," us: gain=",g*1e6," uV/ion, ENC=",s/g)
   end
  end
 end
end
end
