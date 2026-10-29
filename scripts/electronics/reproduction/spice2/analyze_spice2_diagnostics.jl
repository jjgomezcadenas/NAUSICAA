using HDF5, SHA
include(joinpath(@__DIR__,"../../../../src/ion_signal/averaging_scan.jl"))
using .AveragingScan
root=normpath(joinpath(@__DIR__,"../../../.."));dir=joinpath(root,"results/ltspice/spice2_diagnostics")
integral(t,v)=sum(diff(t).*(v[1:end-1].+v[2:end])./2)
function interp(t,v,x)
    j=clamp(searchsortedlast(t,x),1,length(t)-1)
    v[j]+(v[j+1]-v[j])*(x-t[j])/(t[j+1]-t[j])
end
function wave(g,kind,r)
    d=g["$(kind)_N1_r$r"];b=g["$(kind)_N0_r$r"]
    t=read(d["time"]);v=read(d["V(out)"]);tb=read(b["time"]);vb=read(b["V(out)"])
    t,v.-[interp(tb,vb,x) for x in t]
end
function primitive(t,v)
    c=vcat(0.,cumsum(diff(t).*(v[1:end-1].+v[2:end])./2))
    function at(x)
        @assert first(t)<=x<=last(t)
        j=min(searchsortedlast(t,x),length(t)-1);dx=x-t[j]
        c[j]+v[j]*dx+(v[j+1]-v[j])/(t[j+1]-t[j])*dx^2/2
    end
end
function injection_gain(t,v,T;step=.25e-6)
    p=primitive(t,v)
    average(x)=(p(x)-p(x-T))/T
    times=collect(.005:step:.012)
    response=[-(average(x)-average(x-2T)) for x in times]
    j=argmax(response)
    response[j],times[j]
end
h5open(joinpath(dir,"raw.h5"),"r") do h
 h5open(joinpath(dir,"analysis.h5"),"w") do o
    attributes(o)["raw_sha256"]=bytes2hex(sha256(read(joinpath(dir,"raw.h5"))))
    attributes(o)["script_sha256"]=bytes2hex(sha256(read(@__FILE__)))
    for k in keys(attributes(h));attributes(o)[k]=read(attributes(h)[k]);end
    for tag in ("cap_0.75pF","cap_1.25pF","cap_3pF")
        g=h[tag];out=create_group(o,tag);a=attributes(out)
        a["external_capacitance_F"]=read(attributes(g)["external_capacitance_F"])
        f=read(g["noise_400/frequency"]);asd=read(g["noise_400/V(onoise)"])
        fc=read(g["noise_200/frequency"]);ac=read(g["noise_200/V(onoise)"])
        t,v=wave(g,"ions",2);tc,vc=wave(g,"ions",1)
        T=300e-6;frame=.01;alpha=exp(-frame/.3)
        gain=waveform_gain(t.-.002,v,T,frame,alpha)
        sigma=noise_rms(f,asd,T,frame,alpha)
        terr=abs(waveform_gain(tc.-.002,vc,T,frame,alpha)/gain-1)
        serr=abs(noise_rms(fc,ac,T,frame,alpha)/sigma-1)
        qerr=abs(noise_rms(f,asd,T,frame,alpha;order=16)/sigma-1)
        e=1.602176634e-19;half=.00025;gap=.005
        expected=e*(1-2/pi*atan(half^2/(gap*sqrt(gap^2+2half^2))))
        chargeerr=abs(integral(t,read(g["ions_N1_r2/I(Is1)"]))/expected-1)
        @assert max(terr,serr)<1e-3 && qerr<1e-5 && chargeerr<2e-5
        for (k,x) in ("frame_gain_V_per_ion"=>gain,"frame_noise_rms_V"=>sigma,"frame_enc_e"=>sigma/gain,"transient_error"=>terr,"spectral_error"=>serr,"quadrature_error"=>qerr,"charge_error"=>chargeerr)
            a[k]=x
        end
        names=sort(filter(n->occursin(r"^V\(m\w+\)$",n),keys(g["noise_400"])))
        filtered(n)=noise_rms(f,max.(read(g["noise_400"][n]),1e-150),T,frame,alpha)
        variances=[filtered(n)^2 for n in names]
        @assert abs(sum(variances)/sigma^2-1)<1e-3
        out["device_names"]=names;out["device_variance_fraction"]=variances./sigma^2
        flickers=filter(n->occursin("1overf",n),keys(g["noise_400"]))
        a["flicker_variance_fraction"]=sum(filtered(n)^2 for n in flickers)/sigma^2
        println(tag," frame ENC = ",sigma/gain,"; fractions ",collect(zip(names,variances./sigma^2)))
        if tag=="cap_1.25pF"
            ti,vi=wave(g,"injection",2);tic,vic=wave(g,"injection",1)
            @assert abs(integral(ti,read(g["injection_N1_r2/I(Is1)"]))/e-1)<2e-5
            gframe=waveform_gain(ti.-.002,vi,T,frame,alpha)
            a["injection_frame_gain_V_per_e"]=gframe;a["injection_frame_enc_e"]=sigma/gframe
            for width in (180.,200.,250.,300.,350.).*1e-6
                gaini,peak=injection_gain(ti,vi,width)
                gainc,_=injection_gain(tic,vic,width)
                gainfine,_=injection_gain(ti,vi,width;step=.125e-6)
                sigi=noise_rms(f,asd,width,2width,1.)
                se=abs(noise_rms(fc,ac,width,2width,1.)/sigi-1)
                qe=abs(noise_rms(f,asd,width,2width,1.;order=16)/sigi-1)
                te=abs(gainc/gaini-1);pe=abs(gainfine/gaini-1)
                @assert max(se,te,pe)<1e-3 && qe<1e-5
                row=create_group(out,"injection_$(round(Int,width*1e6))us")
                for (k,x) in ("rise_s"=>width,"flat_s"=>width,"gain_V_per_e"=>gaini,"noise_rms_V"=>sigi,"enc_e"=>sigi/gaini,"peak_time_s"=>peak,"transient_error"=>te,"peak_grid_error"=>pe,"spectral_error"=>se,"quadrature_error"=>qe)
                    attributes(row)[k]=x
                end
                println("Injection ",width*1e6," us: ENC = ",sigi/gaini)
            end
            println("Injection with 10 ms boundary estimator: ENC = ",sigma/gframe)
            old=h5open(joinpath(root,"results/ltspice/spice2_enc/analysis.h5"),"r") do old
                read(attributes(old)["enc_electrons_rms"])
            end
            @assert abs(sigma/gain/old-1)<1e-9
        end
    end
 end
end
