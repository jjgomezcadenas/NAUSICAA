module BucketStudy
using HDF5, TOML, SHA, LinearAlgebra, Printf
include("averaging_scan.jl")
using .AveragingScan
const E = 1.602176634e-19
integral(t,v)=sum(diff(t).*(v[1:end-1].+v[2:end])./2)
phi(z,a)=z<=0 ? 1. : 2/pi*atan(a^2/(z*sqrt(z^2+2a^2)))
"""Finite-gap current averaged over uniform collection times 0..frame, per ion."""
function current(t,frame,a,g,v)
    cumulative(u)=E*(phi(clamp(-v*u,0,g),a)-phi(g,a))
    (cumulative(t)-cumulative(t-frame))/frame
end
function gauss(order=8)
    b=[k/sqrt(4k^2-1) for k in 1:order-1]
    eig=eigen(SymTridiagonal(zeros(order),b))
    eig.values,2 .* eig.vectors[1,:].^2
end
"""Bin-integrated current converted to charge-preserving centre-based PWL."""
function pwl(frame,a,g,v,pre,stop,step; train=false)
    edges=collect(range(-pre,stop,length=ceil(Int,(stop+pre)/step)+1))
    nodes,weights=gauss(); offsets=train ? [0.,frame] : [0.]
    boundaries=sort!(unique([x+y for x in (-g/v,0.,frame-g/v,frame) for y in offsets]))
    currents=zeros(length(edges)-1)
    for j in eachindex(currents)
        splits=vcat(edges[j],filter(x->edges[j]<x<edges[j+1],boundaries),edges[j+1])
        q=0.
        for k in 1:length(splits)-1
            mid=(splits[k]+splits[k+1])/2; half=(splits[k+1]-splits[k])/2
            q+=half*sum(w*sum(current(mid+half*n-o,frame,a,g,v) for o in offsets) for (n,w) in zip(nodes,weights))
        end
        currents[j]=q/(edges[j+1]-edges[j])
    end
    t=vcat(first(edges),(edges[1:end-1].+edges[2:end])./2,last(edges)).+pre
    i=vcat(first(currents),currents,last(currents))
    @assert abs(integral(t,i)/(length(offsets)*E*(1-phi(g,a)))-1)<1e-9
    t,i
end
function settings(root)
    cfg=TOML.parsefile(joinpath(root,"config/bucket_study.toml"))
    d=TOML.parsefile(joinpath(root,cfg["detector"]))
    a=d["geometry"]["pad_width_m"]/2;g=d["geometry"]["gap_m"]
    v=d["gas"]["bulk_speed_m_s"]/d["bias"]["bulk_field_V_m"]*(d["bias"]["ic_voltage_V"]-d["bias"]["pad_voltage_V"])/g
    cfg,a,g,v
end
function inputs(root)
    cfg,a,g,v=settings(root)
    for frame in cfg["bucket_sizes_s"]
        dir=joinpath(root,cfg["output"],"$(round(Int,frame*1000))ms");mkpath(dir)
        target=joinpath(dir,"input.h5")
        isfile(target) && error("Preserving existing $target; choose a new output directory for a new configuration")
        h5open(target,"w") do h
            at=attributes(h)
            for (k,x) in ("frame_s"=>frame,"averaging_s"=>cfg["averaging_s"],"correction_tau_s"=>cfg["correction_tau_s"],"prehistory_s"=>cfg["prehistory_s"],"pad_halfwidth_m"=>a,"gap_m"=>g,"gap_speed_m_s"=>v,"single_ion_charge_C"=>E*(1-phi(g,a)))
                at[k]=x
            end
            at["assumptions"]="Central square pad in infinite grounded plane, uniform dielectric half-space; central straight trajectory at constant speed in finite 5 mm gap; uniform mean arrival density. Unit pad weighting potential, zero elsewhere and at infinity. No charge renormalization or collection impulse."
            at["config_sha256"]=bytes2hex(sha256(read(joinpath(root,"config/bucket_study.toml"))))
            at["detector_sha256"]=bytes2hex(sha256(read(joinpath(root,cfg["detector"]))))
            at["generator_sha256"]=bytes2hex(sha256(read(@__FILE__)))
            for (r,step) in enumerate(cfg["max_steps_s"]), train in (false,true)
                t,i=pwl(frame,a,g,v,cfg["prehistory_s"],cfg["post_frames"]*frame,step;train)
                group=create_group(h,"$(train ? "train" : "ions")_r$r")
                group["time_s"]=t;group["current_A_per_ion"]=i
                attributes(group)["max_step_s"]=step
            end
        end
        println(target)
    end
end
function interp(t,v,x)
    j=clamp(searchsortedlast(t,x),1,length(t)-1)
    v[j]+(v[j+1]-v[j])*(x-t[j])/(t[j+1]-t[j])
end
function wave(h,kind,n,r)
    d=h["$(kind)_N$(n)_r$r"];b=h["ions_N0_r$r"]
    t=read(d["time"]);v=read(d["V(out)"]);tb=read(b["time"]);vb=read(b["V(out)"])
    t,v.-[interp(tb,vb,x) for x in t]
end
function analyze(root)
 cfg,a,g,v=settings(root)
 for frame in cfg["bucket_sizes_s"]
  dir=joinpath(root,cfg["output"],"$(round(Int,frame*1000))ms")
  h5open(joinpath(dir,"raw.h5"),"r") do h
   h5open(joinpath(dir,"analysis.h5"),"w") do out
    at=attributes(out)
    for k in keys(attributes(h));at[k]=read(attributes(h)[k]);end
    at["raw_sha256"]=bytes2hex(sha256(read(joinpath(dir,"raw.h5"))))
    at["analysis_sha256"]=bytes2hex(sha256(read(@__FILE__)))
    at["filter_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,"averaging_scan.jl"))))
    T=cfg["averaging_s"];pre=cfg["prehistory_s"];alpha=exp(-frame/cfg["correction_tau_s"])
    gainat(t,v;shift=0.)=waveform_gain(t.-pre.-shift,v,T,frame,alpha)
    f=read(h["noise_400/frequency"]);asd=read(h["noise_400/V(onoise)"])
    fc=read(h["noise_200/frequency"]);ac=read(h["noise_200/V(onoise)"])
    sigma=noise_rms(f,asd,T,frame,alpha)
    serr=abs(noise_rms(fc,ac,T,frame,alpha)/sigma-1)
    qerr=abs(noise_rms(f,asd,T,frame,alpha;order=16)/sigma-1)
    t,v=wave(h,"ions",1,2);tc,vc=wave(h,"ions",1,1)
    gain=gainat(t,v);terr=abs(gainat(tc,vc)/gain-1)
    expected=E*(1-phi(g,a));qerrsignal=abs(integral(t,read(h["ions_N1_r2/I(Is1)"]))/expected-1)
    @assert gain>0 && max(terr,serr)<1e-3 && qerr<1e-5 && qerrsignal<2e-5
    # Isolated packet's response in the following empty frame, relative to its own frame.
    residual=gainat(t,v;shift=frame)/gain
    at["alpha"]=alpha;at["gain_V_per_ion"]=gain;at["noise_rms_V"]=sigma;at["ENC_e_rms"]=sigma/gain
    at["linear_N_snr5"]=floor(Int,5sigma/gain)+1
    at["previous_frame_residual_fraction"]=residual
    at["ideal_estimator_charge_fraction"]=3e-15*(frame/cfg["correction_tau_s"])/(1-alpha)*gain/expected
    at["transient_error"]=terr;at["spectral_error"]=serr;at["quadrature_error"]=qerr;at["charge_error"]=qerrsignal
    counts=sort([parse(Int,match(r"^ions_N(\d+)_r2$",k)[1]) for k in keys(h) if occursin(r"^ions_N\d+_r2$",k) && parse(Int,match(r"^ions_N(\d+)_r2$",k)[1])>1])
    snrs=Float64[];ratios=Float64[];errors=Float64[]
    for n in counts
        tn,vn=wave(h,"ions",n,2);tnc,vnc=wave(h,"ions",n,1)
        d=gainat(tn,vn);err=abs(gainat(tnc,vnc)/d-1)
        @assert err<1e-3
        @assert abs(integral(tn,read(h["ions_N$(n)_r2/I(Is1)"]))/(n*expected)-1)<2e-5
        push!(snrs,d/sigma);push!(ratios,d/(n*gain));push!(errors,err)
    end
    out["counts"]=counts;out["snr"]=snrs;out["gain_relative_to_small_signal"]=ratios;out["finite_count_refinement_error"]=errors
    passing=counts[snrs.>5];failing=counts[snrs.<=5]
    if !isempty(passing) && !isempty(failing) && minimum(passing)-maximum(failing)==1
        threshold=minimum(passing)
        for n in (threshold-1,threshold)
            j=findfirst(==(n),counts)
            @assert errors[j]+serr+qerr<abs(snrs[j]/5-1)/10 "Threshold requires finer convergence"
        end
        at["verified_N_snr5"]=threshold
    end
    nt=cfg["train_count"]
    tt,vt=wave(h,"train",nt,2);ttc,vtc=wave(h,"train",nt,1)
    ti,vi=wave(h,"ions",nt,2)
    trainsecond=gainat(tt,vt;shift=frame);isolated=gainat(ti,vi)
    at["train_second_frame_snr"]=trainsecond/sigma
    at["train_second_vs_isolated_fraction"]=trainsecond/isolated-1
    at["train_refinement_error"]=abs(gainat(ttc,vtc;shift=frame)/trainsecond-1)
    @assert read(at["train_refinement_error"])<1e-3
    names=sort(filter(n->occursin(r"^V\(m\w+\)$",n),keys(h["noise_400"])))
    vars=[noise_rms(f,max.(read(h["noise_400"][n]),1e-150),T,frame,alpha)^2 for n in names]
    @assert abs(sum(vars)/sigma^2-1)<1e-3
    out["device_names"]=names;out["device_variance_fraction"]=vars./sigma^2
    out["time_s"]=ti.-pre;out["voltage_V"]=vi;out["current_A"]=read(h["ions_N$(nt)_r2/I(Is1)"])
    out["train_time_s"]=tt.-pre;out["train_voltage_V"]=vt
    out["frequency_Hz"]=f;out["output_asd_V_sqrtHz"]=asd
    if frame==.01
        h5open(joinpath(root,"results/ltspice/m3_scan/analysis.h5"),"r") do old
            ref=attributes(old[cfg["noise_group"]]);ratio=gain/read(ref["frame_gain_V_per_ion"])-1
            at["historical_10ms_gain_relative_change"]=ratio
            @assert abs(ratio)<1e-3
        end
    end
    println("frame=",frame," ENC=",sigma/gain," gain=",gain," threshold=",read(at["linear_N_snr5"])," previous=",residual," counts=",counts," SNR=",snrs)
   end
  end
 end
end
end
