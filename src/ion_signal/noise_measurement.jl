module NoiseMeasurement
using HDF5, SHA
include("averaging_scan.jl")
const A=AveragingScan

"""Compare output-noise RMS for three explicitly defined measurements."""
function analyze(root)
    reference=joinpath(root,"results/ltspice/bias_sensitivity/analysis.h5")
    p=h5open(reference,"r") do h
        a=attributes(h["paper_rounded_4"])
        (gain=read(a["gain_V_per_ion"]),tau=read(a["tau_s"]))
    end
    frame=.02; width=.0003; alpha=exp(-frame/p.tau)
    out=joinpath(root,"results/ltspice/noise_measurement/analysis.h5");mkpath(dirname(out))
    h5open(out,"w") do h
        a=attributes(h);a["frame_s"]=frame;a["average_s"]=width;a["alpha"]=alpha
        a["gain_V_per_ion"]=p.gain;a["tau_s"]=p.tau
        a["analysis_sha256"]=bytes2hex(sha256(read(@__FILE__)))
        a["reference_sha256"]=bytes2hex(sha256(read(reference)))
        a["filter_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,"averaging_scan.jl"))))
        a["rogelio_reported_peak_V"]=27.8e-6;a["rogelio_reported_noise_V"]=.0108
        a["reported_injected_charge_e"]=.096e-15*.0017/1.602176634e-19
        for (name,path) in [("paper_rounded","results/ltspice/bias_sensitivity/paper_rounded_4.h5"),("v3","results/ltspice/noise_measurement/raw.h5")]
            isfile(joinpath(root,path)) || continue
            group=create_group(h,name);attributes(group)["source"]=path
            attributes(group)["sha256"]=bytes2hex(sha256(read(joinpath(root,path))))
            h5open(joinpath(root,path),"r") do raw
                for density in (200,400)
                    f=read(raw["noise_$density/frequency"]);asd=read(raw["noise_$density/V(onoise)"])
                    g=create_group(group,string(density))
                    for (label,w,al) in [("raw",0.,0.),("average",width,0.),("difference",width,alpha)]
                        sigma=A.noise_rms(f,asd,w,frame,al)
                        attributes(g)[label*"_rms_V"]=sigma
                        attributes(g)[label*"_quadrature_relative_error"]=abs(A.noise_rms(f,asd,w,frame,al;order=12)/sigma-1)
                    end
                    rawrms=read(attributes(g)["raw_rms_V"]);av=read(attributes(g)["average_rms_V"]);d=read(attributes(g)["difference_rms_V"])
                    attributes(g)["averaged_covariance_V2"]=(av^2*(1+alpha^2)-d^2)/(2alpha)
                    attributes(g)["averaged_noise_correlation"]=(av^2*(1+alpha^2)-d^2)/(2alpha*av^2)
                    attributes(g)["noise_over_reference_gain_e"]=d/p.gain
                    attributes(g)["raw_over_reported_peak_e"]=rawrms/27.8e-6
                    if name=="v3"
                        text=read(attributes(raw["noise_$density"])["log"])
                        m=match(r"int_sigma: INTEG\(v\(onoise\)\)=([\d.eE+-]+)",text)
                        @assert !isnothing(m)
                        attributes(g)["ltspice_measured_rms_V"]=parse(Float64,m[1])
                    end
                    println(name," ",density,": raw=",rawrms," average=",av," difference=",d," ENC(using rounded gain)=",d/p.gain)
                end
            end
        end
    end
    out
end
end
