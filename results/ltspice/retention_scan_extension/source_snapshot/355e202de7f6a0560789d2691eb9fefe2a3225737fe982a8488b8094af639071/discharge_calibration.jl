module DischargeCalibration
using HDF5,TOML,SHA,LinearAlgebra
include("bucket_study.jl")
const B=BucketStudy
const A=B.AveragingScan
sha(path)=bytes2hex(sha256(read(path)))
function configuration(root)
    file=joinpath(root,"config/discharge_calibration.toml")
    TOML.parsefile(file),file
end
function inputs(root)
    cfg,file=configuration(root);dir=joinpath(root,cfg["output"]);mkpath(dir)
    path=joinpath(dir,"input.h5");isfile(path) && error("Preserving existing $path")
    d=TOML.parsefile(joinpath(root,cfg["detector"]))
    a=d["geometry"]["pad_width_m"]/2;g=d["geometry"]["gap_m"]
    v=d["gas"]["bulk_speed_m_s"]/d["bias"]["bulk_field_V_m"]*(d["bias"]["ic_voltage_V"]-d["bias"]["pad_voltage_V"])/g
    h5open(path,"w") do h
        at=attributes(h);at["config_sha256"]=sha(file);at["generator_sha256"]=sha(@__FILE__)
        at["waveform_source_sha256"]=sha(joinpath(@__DIR__,"bucket_study.jl"))
        at["detector_sha256"]=sha(joinpath(root,cfg["detector"]))
        at["single_ion_charge_C"]=B.E*(1-B.phi(g,a))
        at["assumptions"]="Central square pad of side 0.5 mm in an infinite grounded plane, homogeneous gas half-space; pad weighting potential one, others zero, zero at infinity. Constant 3 m/s central approach over 5 mm; finite-path charge preserved. Mean uniform collection times over each 20 ms bucket."
        for (r,step) in enumerate(cfg["max_steps_s"])
            t,i=B.pwl(cfg["frame_s"],a,g,v,cfg["prehistory_s"],cfg["stop_s"],step)
            group=create_group(h,"r$r");group["time_s"]=t;group["single_A_per_ion"]=i
            # Translate the charge-preserving single-packet PWL by exact grid intervals.
            train=zeros(length(t));centres=t[2:end-1];base=i[2:end-1]
            width=centres[2]-centres[1]
            for (k,n) in enumerate(cfg["train_counts"])
                shift=round(Int,(k-1)*cfg["frame_s"]/width)
                @assert isapprox(shift*width,(k-1)*cfg["frame_s"];atol=1e-12)
                if n>0;train[2+shift:end-1].+=n.*base[1:end-shift];end
            end
            @assert abs(B.integral(t,train)/(sum(cfg["train_counts"])*B.E*(1-B.phi(g,a)))-1)<1e-9
            group["train_A"]=train
        end
    end
end
"""Fit an exponential to a baseline-subtracted, current-free tail on a uniform time grid."""
function fit_tail(t,v,start,stop,step)
    ts=collect(start:step:stop);y=-[B.interp(t,v,x) for x in ts]
    @assert all(y.>0) && first(t)<=start<stop<=last(t)
    X=hcat(ones(length(ts)),ts.-start);c=X\log.(y)
    @assert c[2]<0 "Tail is not decaying"
    fit=exp.(X*c);tau=-1/c[2]
    tau,maximum(abs.(fit./y.-1)),ts,y,fit
end
function average(t,v,center,width)
    @assert first(t)<=center-width/2 && center+width/2<=last(t)
    lo=center-width/2;hi=center+width/2
    inside=findall(x->lo<x<hi,t)
    times=vcat(lo,t[inside],hi);values=vcat(B.interp(t,v,lo),v[inside],B.interp(t,v,hi))
    B.integral(times,values)/width
end
function waveform(h,n,r)
    d=h["single_N$(n)_r$r"];b=h["single_N0_r$r"]
    t=read(d["time"]);tb=read(b["time"]);vb=read(b["V(out)"])
    t,read(d["V(out)"]).-[B.interp(tb,vb,x) for x in t]
end
function analyze(root)
    cfg,file=configuration(root);dir=joinpath(root,cfg["output"])
    frame=cfg["frame_s"];width=cfg["averaging_s"];pre=cfg["prehistory_s"]
    h5open(joinpath(dir,"raw.h5"),"r") do h
        @assert read(attributes(h)["config_sha256"])==sha(file)
        h5open(joinpath(dir,"analysis.h5"),"w") do out
            at=attributes(out)
            for k in ("frame_s","averaging_s","calibration_ions","old_tau_s");at[k]=cfg[k];end
            at["raw_sha256"]=sha(joinpath(dir,"raw.h5"));at["analysis_sha256"]=sha(@__FILE__)
            at["config_sha256"]=sha(file)
            at["filter_sha256"]=sha(joinpath(@__DIR__,"averaging_scan.jl"))
            at["assumptions"]=read(attributes(h)["assumptions"])
            at["fit_start_s"]=cfg["fit_start_s"];at["fit_stop_s"]=cfg["fit_stop_s"]
            ncal=cfg["calibration_ions"]
            fit(n,r)=begin
                t,v=waveform(h,n,r)
                fit_tail(t.-pre,v,cfg["fit_start_s"],cfg["fit_stop_s"],cfg["fit_sample_step_s"])
            end
            tau,err,ts,y,yfit=fit(ncal,2);tauc,_,_,_,_=fit(ncal,1)
            tau10,_,_,_,_=fit(10,2)
            at["tau_s"]=tau;at["alpha"]=exp(-frame/tau)
            at["tail_fit_max_relative_error"]=err
            at["tau_timestep_relative_error"]=abs(tauc/tau-1)
            at["tau_10_vs_100_relative_change"]=tau10/tau-1
            @assert abs(tauc/tau-1)<1e-3
            @assert abs(tau10/tau-1)<.01 "Calibration depends on amplitude"
            out["tail_time_s"]=ts;out["tail_magnitude_V"]=y;out["tail_fit_V"]=yfit
            windows=cfg["fit_checks_s"];localtaus=Float64[]
            tcal,vcal=waveform(h,ncal,2)
            for bounds in windows
                localtau,_,_,_,_=fit_tail(tcal.-pre,vcal,bounds[1],bounds[2],cfg["fit_sample_step_s"])
                push!(localtaus,localtau)
            end
            out["fit_window_start_s"]=[x[1] for x in windows];out["fit_window_stop_s"]=[x[2] for x in windows];out["local_tau_s"]=localtaus
            f=read(h["noise_400/frequency"]);asd=read(h["noise_400/V(onoise)"])
            fc=read(h["noise_200/frequency"]);ac=read(h["noise_200/V(onoise)"])
            t1,v1=waveform(h,1,2);tc,vc=waveform(h,1,1)
            counts=sort([parse(Int,match(r"^single_N(\d+)_r2$",k)[1]) for k in keys(h) if occursin(r"^single_N\d+_r2$",k) && parse(Int,match(r"^single_N(\d+)_r2$",k)[1])>1])
            for (name,tau_used) in (("old",cfg["old_tau_s"]),("calibrated",tau))
                group=create_group(out,name);a=attributes(group);alpha=exp(-frame/tau_used)
                signal(t,v;shift=0.)=A.waveform_gain(t.-pre.-shift,v,width,frame,alpha)
                gain=signal(t1,v1);sigma=A.noise_rms(f,asd,width,frame,alpha)
                terr=abs(signal(tc,vc)/gain-1)
                serr=abs(A.noise_rms(fc,ac,width,frame,alpha)/sigma-1)
                qerr=abs(A.noise_rms(f,asd,width,frame,alpha;order=16)/sigma-1)
                @assert max(terr,serr)<1e-3 && qerr<1e-5
                for (key,value) in ("tau_s"=>tau_used,"alpha"=>alpha,"gain_V_per_ion"=>gain,"noise_rms_V"=>sigma,"ENC_e_rms"=>sigma/gain,"linear_N_snr5"=>floor(Int,5sigma/gain)+1,"transient_error"=>terr,"spectral_error"=>serr,"quadrature_error"=>qerr)
                    a[key]=value
                end
                a["previous_frame_residual_fraction"]=signal(t1,v1;shift=frame)/gain
                # Future packet begins inducing before collection at the right boundary.
                a["following_frame_fraction"]=-average(t1.-pre,v1,0.,width)/gain
                snrs=Float64[];errors=Float64[];taus=Float64[];residuals=Float64[]
                for n in counts
                    tn,vn=waveform(h,n,2);tnc,vnc=waveform(h,n,1)
                    d=signal(tn,vn);er=abs(signal(tnc,vnc)/d-1)
                    @assert er<1e-3
                    expected=read(attributes(h)["single_ion_charge_C"])*n
                    @assert abs(B.integral(tn,read(h["single_N$(n)_r2/I(Is1)"]))/expected-1)<2e-5
                    tnfit,_,_,_,_=fit_tail(tn.-pre,vn,cfg["fit_start_s"],cfg["fit_stop_s"],cfg["fit_sample_step_s"])
                    push!(taus,tnfit);push!(residuals,signal(tn,vn;shift=frame)/d)
                    push!(snrs,d/sigma);push!(errors,er)
                end
                group["counts"]=counts;group["snr"]=snrs;group["refinement_error"]=errors
                group["tail_tau_by_count_s"]=taus;group["previous_residual_by_count"]=residuals
                pass=counts[snrs.>5];fail=counts[snrs.<=5]
                if !isempty(pass) && !isempty(fail) && minimum(pass)-maximum(fail)==1
                    n=minimum(pass)
                    for nn in (n-1,n)
                        j=findfirst(==(nn),counts)
                        @assert errors[j]+serr+qerr<abs(snrs[j]/5-1)/10 "Threshold needs refinement"
                    end
                    a["verified_N_snr5"]=n
                end
                d=h["train_r2"];b=h["single_N0_r2"]
                tt=read(d["time"]);tb=read(b["time"]);vb=read(b["V(out)"])
                vt=read(d["V(out)"]).-[B.interp(tb,vb,x) for x in tt]
                dco=h["train_r1"];bco=h["single_N0_r1"]
                tco=read(dco["time"]);tbc=read(bco["time"]);vbc=read(bco["V(out)"])
                vco=read(dco["V(out)"]).-[B.interp(tbc,vbc,x) for x in tco]
                estimates=[signal(tt,vt;shift=(k-1)*frame)/gain for k in eachindex(cfg["train_counts"])]
                estimatesco=[signal(tco,vco;shift=(k-1)*frame)/gain for k in eachindex(cfg["train_counts"])]
                group["train_true_counts"]=cfg["train_counts"];group["train_estimated_counts"]=estimates
                group["train_refinement_ions"]=abs.(estimatesco.-estimates)
                @assert maximum(abs.(estimatesco.-estimates))<.05
                group["train_error_ions"]=estimates.-cfg["train_counts"]
                if name=="old"
                    h5open(joinpath(root,cfg["historical"],"analysis.h5"),"r") do old
                        a["historical_enc_relative_change"]=sigma/gain/read(attributes(old)["ENC_e_rms"])-1
                        @assert abs(read(a["historical_enc_relative_change"]))<1e-3
                    end
                end
                # Sensitivity to fitting-window choice propagated through signal and noise.
                group["local_fit_ENC_e"]=[A.noise_rms(f,asd,width,frame,exp(-frame/ta))/A.waveform_gain(t1.-pre,v1,width,frame,exp(-frame/ta)) for ta in localtaus]
                println(name," tau=",tau_used," ENC=",sigma/gain," gain=",gain," threshold linear=",read(a["linear_N_snr5"])," residual=",read(a["previous_frame_residual_fraction"])," counts=",counts," snr=",snrs," train=",estimates)
            end
            out["single_time_s"]=t1.-pre;out["single_voltage_V"]=v1
            out["frequency_Hz"]=f;out["output_asd_V_sqrtHz"]=asd
            println("calibration tau=",tau," fit error=",err," timestep=",abs(tauc/tau-1)," amplitude=",tau10/tau-1," local taus=",localtaus)
        end
    end
end
end
