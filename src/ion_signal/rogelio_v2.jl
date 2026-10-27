module RogelioV2
using HDF5,TOML,SHA
include("discharge_calibration.jl")
const D=DischargeCalibration
const A=D.A
const B=D.B
sha(p)=bytes2hex(sha256(read(p)))
function diagnostics(root)
  dir=joinpath(root,"results/ltspice/rogelio_v2")
  h5open(joinpath(dir,"diagnostics.h5"),"r") do raw
    h5open(joinpath(dir,"diagnostic_analysis.h5"),"w") do out
      attributes(out)["raw_sha256"]=sha(joinpath(dir,"diagnostics.h5"))
      for name in ("paper","v2_native","v2_matched","v2_retuned")
        haskey(raw,name*"_ac") || continue
        g=create_group(out,name);a=attributes(g)
        ac=raw[name*"_ac"];f=real.(read(ac["frequency"]))
        gain=read(ac["V(out)"])./read(ac["V(in)"])
        g["frequency_Hz"]=f;g["voltage_gain"]=gain
        a["gain_magnitude_1Hz"]=abs(first(gain))
        log=read(joinpath(dir,"decks/diagnostics/$(name)_device.run.txt"),String)
        lines=split(log,'\n');names=String[]
        for line in lines
          startswith(line,"Name:") && (names=split(line)[2:end])
          if startswith(line,"Gds:") && "m4" in names
            gd=parse(Float64,split(line)[findfirst(==("m4"),names)+1]);a["M4_gds_S"]=gd;a["M4_ro_Ohm"]=1/gd
          end
          if startswith(line,"Gm:") && "m4" in names
            a["M4_gm_S"]=parse(Float64,split(line)[findfirst(==("m4"),names)+1])
          end
        end
        a["baseline_V"]=only(read(raw[name*"_device/V(out)"]))
        println(name," |Vout/Vin|(1Hz)=",abs(first(gain))," M4 ro=",read(a["M4_ro_Ohm"])," Vout=",read(a["baseline_V"]))
      end
    end
  end
end
function analyze(root;config="config/rogelio_v2.toml")
    cp=joinpath(root,config);c=TOML.parsefile(cp)
    dir=joinpath(root,c["output"]);path=joinpath(dir,"raw.h5")
    h5open(path,"r") do raw
      @assert read(attributes(raw)["config_sha256"])==sha(cp)
      h5open(joinpath(dir,"analysis.h5"),"w") do out
        at=attributes(out);at["raw_sha256"]=sha(path);at["config_sha256"]=sha(cp)
        at["analysis_sha256"]=sha(@__FILE__);at["filter_sha256"]=sha(joinpath(@__DIR__,"averaging_scan.jl"))
        at["assumptions"]=read(attributes(raw)["assumptions"])
        frame=c["frame_s"];width=c["averaging_s"];pre=c["prehistory_s"]
        at["frame_s"]=frame;at["averaging_s"]=width
        if haskey(raw,"reproduction")
          rg=create_group(out,"reproduction")
          r=raw["reproduction"]
          for name in ("original","corrected","supplied_noise")
            f=read(r[name]["frequency"])
            asd=read(r[name][name=="supplied_noise" ? "v(onoise)" : "V(onoise)"])
            g=create_group(rg,name);g["frequency_Hz"]=f;g["output_asd_V_sqrtHz"]=asd
          end
          supplied=rg["supplied_noise"];f=read(supplied["frequency_Hz"]);s=read(supplied["output_asd_V_sqrtHz"])
          for name in ("original","corrected")
            ff=read(rg[name]["frequency_Hz"]);ss=read(rg[name]["output_asd_V_sqrtHz"])
            @assert length(ff)==length(f) && maximum(abs.(ff./f.-1))<1e-10
            attributes(rg[name])["max_relative_difference_supplied"]=maximum(abs.(ss./s.-1))
            attributes(rg[name])["asd_1Hz_V_sqrtHz"]=first(ss)
          end
        end
        for variant in c["variants"]
          name=variant["name"];h=raw[name];g=create_group(out,name);a=attributes(g)
          for (k,v) in variant;a[k]=v;end
          t,v=D.waveform(h,1,2);tc,vc=D.waveform(h,1,1)
          t.-=pre;tc.-=pre
          fit(n,r)=begin
            tt,vv=D.waveform(h,n,r)
            D.fit_tail(tt.-pre,vv,c["fit_start_s"],c["fit_stop_s"],c["fit_sample_step_s"])
          end
          tau,fiterr,ts,y,yfit=fit(100,2);tauco,_,_,_,_=fit(100,1)
          tau10,_,_,_,_=fit(10,2)
          alpha=exp(-frame/tau)
          a["tau_s"]=tau;a["tail_fit_max_relative_error"]=fiterr
          a["tau_timestep_relative_error"]=abs(tauco/tau-1)
          a["tau_amplitude_relative_change"]=tau10/tau-1
          f=read(h["noise_400/frequency"]);asd=read(h["noise_400/V(onoise)"])
          fc=read(h["noise_200/frequency"]);ac=read(h["noise_200/V(onoise)"])
          gain=A.waveform_gain(t,v,width,frame,alpha);sigma=A.noise_rms(f,asd,width,frame,alpha)
          a["gain_V_per_ion"]=gain;a["noise_rms_V"]=sigma;a["ENC_e_rms"]=sigma/gain
          a["linear_N_snr5"]=floor(Int,5sigma/gain)+1
          a["transient_relative_error"]=abs(A.waveform_gain(tc,vc,width,frame,alpha)/gain-1)
          a["spectral_relative_error"]=abs(A.noise_rms(fc,ac,width,frame,alpha)/sigma-1)
          a["quadrature_relative_error"]=abs(A.noise_rms(f,asd,width,frame,alpha;order=16)/sigma-1)
          op=h["op"];a["baseline_V"]=only(read(op["V(out)"]))
          a["input_baseline_V"]=only(read(op["V(in)"]))
          a["M5_Vds_V"]=only(read(op["V(vd6)"]))-only(read(op["V(out)"]))
          a["power_W"]=-1.8*only(read(op["I(Vdd)"]))
          counts=sort([parse(Int,match(r"^single_N(\d+)_r2$",k)[1]) for k in keys(h) if occursin(r"^single_N\d+_r2$",k) && parse(Int,match(r"^single_N(\d+)_r2$",k)[1])>0])
          snrs=Float64[];gains=Float64[];taus=Float64[];errors=Float64[];charges=Float64[]
          for n in counts
            tn,vn=D.waveform(h,n,2);tnc,vnc=D.waveform(h,n,1);tn.-=pre;tnc.-=pre
            signal=A.waveform_gain(tn,vn,width,frame,alpha)
            push!(snrs,signal/sigma);push!(gains,signal/n)
            push!(errors,abs(A.waveform_gain(tnc,vnc,width,frame,alpha)/signal-1))
            q=B.integral(tn,read(h["single_N$(n)_r2/I(Is1)"]))
            push!(charges,q/(n*read(attributes(raw)["single_ion_charge_C"]))-1)
            τ,_,_,_,_=D.fit_tail(tn,vn,c["fit_start_s"],c["fit_stop_s"],c["fit_sample_step_s"])
            push!(taus,τ)
          end
          g["counts"]=counts;g["snr"]=snrs;g["gain_by_count_V_per_ion"]=gains
          a["max_gain_relative_change_to_2000"]=maximum(abs.(gains./gain.-1))
          g["tau_by_count_s"]=taus;g["signal_refinement_error"]=errors;g["charge_relative_error"]=charges
          j=findfirst(>(5),snrs)
          if j!==nothing && j>1
            cross=counts[j-1]+(5-snrs[j-1])*(counts[j]-counts[j-1])/(snrs[j]-snrs[j-1])
            lo=floor(Int,cross);g["suggested_counts"]=collect(max(1,lo-1):lo+2)
            if counts[j]-counts[j-1]==1
              a["verified_N_snr5"]=counts[j]
              a["threshold_below_snr"]=snrs[j-1];a["threshold_above_snr"]=snrs[j]
              @assert maximum(errors[j-1:j])+read(a["spectral_relative_error"])+read(a["quadrature_relative_error"]) < minimum(abs.(snrs[j-1:j]./5 .-1))/5
            end
          end
          g["time_s"]=t;g["voltage_V_per_ion"]=v;g["frequency_Hz"]=f;g["output_asd_V_sqrtHz"]=asd
          g["tail_time_s"]=ts;g["tail_magnitude_V"]=y;g["tail_fit_V"]=yfit
          names=sort(filter(n->occursin(r"^V\(m\w+\)$",n),keys(h["noise_400"])))
          vars=[A.noise_rms(f,max.(read(h["noise_400"][n]),1e-150),width,frame,alpha)^2 for n in names]
          g["device_names"]=names;g["device_variance_fraction"]=vars./sigma^2
          a["device_variance_sum_relative_error"]=sum(vars)/sigma^2-1
          @assert maximum(abs.(charges))<2e-5
          @assert maximum(errors)<1e-3
          @assert read(a["spectral_relative_error"])<1e-3 && read(a["quadrature_relative_error"])<1e-5
          @assert abs(read(a["device_variance_sum_relative_error"]))<1e-3
          println(name," Vout=",read(a["baseline_V"])," tau=",tau," fiterr=",fiterr," gain uV=",gain*1e6," ENC=",sigma/gain," counts=",counts," SNR=",snrs," gains=",gains)
          println("  refinement=",errors," charge=",charges," tau10/100=",tau10/tau," noise refinement=",read(a["spectral_relative_error"]))
        end
        h5open(joinpath(root,c["reference"]),"r") do ref
          g=create_group(out,"paper_reference")
          for k in keys(attributes(ref["calibrated"]));attributes(g)[k]=read(attributes(ref["calibrated"])[k]);end
          g["time_s"]=read(ref["single_time_s"]);g["voltage_V_per_ion"]=read(ref["single_voltage_V"])
          g["frequency_Hz"]=read(ref["frequency_Hz"]);g["output_asd_V_sqrtHz"]=read(ref["output_asd_V_sqrtHz"])
          g["counts"]=read(ref["calibrated/counts"]);g["snr"]=read(ref["calibrated/snr"])
          attributes(g)["source_sha256"]=sha(joinpath(root,c["reference"]))
        end
      end
    end
end
end
