module M4Scan
using HDF5,SHA
include("averaging_scan.jl")
using .AveragingScan
integral(t,v)=sum(diff(t).*(v[1:end-1].+v[2:end])./2)
function interp(t,v,x)
    j=clamp(searchsortedlast(t,x),1,length(t)-1)
    v[j]+(v[j+1]-v[j])*(x-t[j])/(t[j+1]-t[j])
end
function wave(g,kind,n,r)
    d=g["$(kind)_N$(n)_r$r"];b=g["$(kind)_N0_r$r"]
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
    p=primitive(t,v);average(x)=(p(x)-p(x-T))/T
    times=collect(.005:step:.012)
    response=[-(average(x)-average(x-2T)) for x in times]
    j=argmax(response);response[j],times[j]
end
function rise_time(t,v)
    y=-v;jp=argmax(y);peak=y[jp]
    crossing(f)=begin
        j=findfirst(i->t[i]>=.005 && y[i]>=f*peak,1:jp)
        @assert j!==nothing && j>1
        t[j-1]+(f*peak-y[j-1])/(y[j]-y[j-1])*(t[j]-t[j-1])
    end
    crossing(.9)-crossing(.1)
end

function analyze(root; study="m4_scan", reference="spice2_enc", reference_group=nothing, gate_node="V(vncasc)")
 dir=joinpath(root,"results/ltspice",study)
 h5open(joinpath(dir,"raw.h5"),"r") do h
  h5open(joinpath(dir,"analysis.h5"),"w") do o
    for k in keys(attributes(h));attributes(o)[k]=read(attributes(h)[k]);end
    attributes(o)["raw_sha256"]=bytes2hex(sha256(read(joinpath(dir,"raw.h5"))))
    attributes(o)["analysis_sha256"]=bytes2hex(sha256(read(@__FILE__)))
    attributes(o)["filter_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,"averaging_scan.jl"))))
    tags=sort(filter(k->startswith(k,"scale_") && haskey(h[k],"noise_400"),keys(h));by=k->parse(Int,split(k,'_')[2]))
    for tag in tags
        g=h[tag];out=create_group(o,tag);a=attributes(out)
        for k in keys(attributes(g));a[k]=read(attributes(g)[k]);end
        f=read(g["noise_400/frequency"]);asd=read(g["noise_400/V(onoise)"])
        fc=read(g["noise_200/frequency"]);ac=read(g["noise_200/V(onoise)"])
        t,v=wave(g,"ions",1,2);tc,vc=wave(g,"ions",1,1)
        T=300e-6;frame=.01;alpha=exp(-frame/.3)
        gain=waveform_gain(t.-.002,v,T,frame,alpha)
        sigma=noise_rms(f,asd,T,frame,alpha)
        terr=abs(waveform_gain(tc.-.002,vc,T,frame,alpha)/gain-1)
        serr=abs(noise_rms(fc,ac,T,frame,alpha)/sigma-1)
        qerr=abs(noise_rms(f,asd,T,frame,alpha;order=16)/sigma-1)
        e=1.602176634e-19;half=.00025;gap=.005
        expected=e*(1-2/pi*atan(half^2/(gap*sqrt(gap^2+2half^2))))
        chargeerr=abs(integral(t,read(g["ions_N1_r2/I(Is1)"]))/expected-1)
        @assert gain>0 && max(terr,serr)<1e-3 && qerr<1e-5 && chargeerr<2e-5
        ikind="injection"
        inum=haskey(g,"injection_N100_r2") ? 100 : 1
        ir=sort([parse(Int,split(k,"_r")[end]) for k in keys(g) if startswith(k,ikind*"_N$(inum)_r")])
        @assert length(ir)>=2
        irfine=ir[end];ircoarse=ir[end-1]
        ti,vi=wave(g,ikind,inum,irfine);tic,vic=wave(g,ikind,inum,ircoarse)
        vi=vi./inum;vic=vic./inum
        a["injection_calibration_electrons"]=inum
        a["injection_max_timestep_s"]=1e-6/irfine
        a["injection_acquisition_kind"]=ikind
        @assert abs(integral(ti,read(g["$(ikind)_N$(inum)_r$(irfine)/I(Is1)"]))/(inum*e)-1)<2e-5
        gi,peak=injection_gain(ti,vi,T);gic,_=injection_gain(tic,vic,T)
        gif,_=injection_gain(ti,vi,T;step=.125e-6)
        sigi=noise_rms(f,asd,T,2T,1.)
        if inum==100
            t10,v10=wave(g,ikind,10,2);g10,_=injection_gain(t10,v10./10,T)
            a["injection_gain_10_vs_100_relative"]=abs(g10/gi-1)
            @assert abs(g10/gi-1)<1e-3 "Injection amplitude linearity: tag=$tag ratio=$(g10/gi)"
        end
        iterr=abs(gic/gi-1);perr=abs(gif/gi-1)
        iserr=abs(noise_rms(fc,ac,T,2T,1.)/sigi-1)
        iqerr=abs(noise_rms(f,asd,T,2T,1.;order=16)/sigi-1)
        @assert max(iterr,perr,iserr)<1e-3 && iqerr<1e-5 "Injection convergence: tag=$tag transient=$iterr peak=$perr spectral=$iserr quadrature=$iqerr"
        for (k,x) in ("frame_gain_V_per_ion"=>gain,"frame_noise_rms_V"=>sigma,"frame_enc_e"=>sigma/gain,"linear_N_snr5"=>floor(Int,5sigma/gain)+1,"transient_error"=>terr,"spectral_error"=>serr,"quadrature_error"=>qerr,"charge_error"=>chargeerr,"injection_enc_e"=>sigi/gi,"injection_gain_V_per_e"=>gi,"injection_peak_s"=>peak,"injection_rise_10_90_s"=>rise_time(ti,vi),"injection_transient_error"=>iterr,"peak_grid_error"=>perr,"injection_spectral_error"=>iserr,"injection_quadrature_error"=>iqerr)
            a[k]=x
        end
        names=sort(filter(n->occursin(r"^V\(m\w+\)$",n),keys(g["noise_400"])))
        filtered(n)=noise_rms(f,max.(read(g["noise_400"][n]),1e-150),T,frame,alpha)
        variances=[filtered(n)^2 for n in names]
        @assert abs(sum(variances)/sigma^2-1)<1e-3
        out["device_names"]=names;out["device_variance_fraction"]=variances./sigma^2;out["device_enc_e"]=sqrt.(variances)./gain
        flickers=filter(n->occursin("1overf",n),keys(g["noise_400"]))
        a["flicker_variance_fraction"]=sum(filtered(n)^2 for n in flickers)/sigma^2
        counts=sort([parse(Int,match(r"^ions_N(\d+)_r2$",k)[1]) for k in keys(g) if occursin(r"^ions_N\d+_r2$",k)])
        snrs=Float64[];ratios=Float64[];refinements=Float64[]
        for n in filter(>(1),counts)
            tn,vn=wave(g,"ions",n,2);tnc,vnc=wave(g,"ions",n,1)
            gn=waveform_gain(tn.-.002,vn,T,frame,alpha)
            refinement=abs(waveform_gain(tnc.-.002,vnc,T,frame,alpha)/gn-1)
            @assert refinement<1e-3
            push!(refinements,refinement)
            @assert abs(integral(tn,read(g["ions_N$(n)_r2/I(Is1)"]))/(n*expected)-1)<2e-5
            push!(snrs,gn/sigma);push!(ratios,gn/(n*gain))
        end
        out["checked_counts"]=filter(>(1),counts);out["checked_snr"]=snrs;out["gain_relative_to_small_signal"]=ratios
        out["finite_count_refinement_error"]=refinements
        checked=filter(>(1),counts);passing=checked[snrs.>5];failing=checked[snrs.<=5]
        if !isempty(passing) && !isempty(failing) && minimum(passing)-maximum(failing)==1
            threshold=minimum(passing);a["verified_N_snr5"]=threshold
            for n in (threshold-1,threshold)
                j=findfirst(==(n),checked)
                @assert refinements[j]+serr+qerr<abs(snrs[j]/5-1)/10
            end
            println("Verified threshold: ",threshold," counts=",checked," snrs=",snrs)
        end
        op=g["op"];lookup(name)=only(read(op[only(filter(k->lowercase(k)==lowercase(name),keys(op)))]))
        for (k,name) in ("dc_output_V"=>"V(out)","dc_input_V"=>"V(in)","dc_med_V"=>"V(med)","M3_current_A"=>"Id(M3)","M4_current_A"=>"Id(M4)","supply_current_A"=>"I(Vdd)")
            a[k]=lookup(name)
        end
        a["power_W"]=-1.8*lookup("I(Vdd)")
        haskey(a,"gate_V") || (a["gate_V"]=lookup(gate_node))
        baseattrs=attributes(o["scale_1"])
        a["quiescent_output_error_V"]=lookup("V(out)")-read(baseattrs["dc_output_V"])
        a["supply_current_relative_change"]=lookup("I(Vdd)")/read(baseattrs["supply_current_A"])-1
        a["M3_current_relative_change"]=lookup("Id(M3)")/read(baseattrs["M3_current_A"])-1
        out["ion_time_s"]=t.-.002;out["ion_voltage_V"]=v
        out["injection_time_s"]=ti.-.005;out["injection_voltage_V"]=vi
        out["frequency_Hz"]=f;out["output_asd_V_sqrtHz"]=asd
        println(tag," ENC=",sigma/gain," injection=",sigi/gi," gain=",gain," rise=",rise_time(ti,vi)," ratios=",ratios," fractions=",collect(zip(names,variances./sigma^2)))
        if tag=="scale_1"
            h5open(joinpath(root,"results/ltspice",reference,"analysis.h5"),"r") do old
                ref=reference_group===nothing ? old : old[reference_group]
                field=reference_group===nothing ? "enc_electrons_rms" : "frame_enc_e"
                @assert abs(sigma/gain/read(attributes(ref)[field])-1)<1e-9
            end
        end
    end
  end
 end
end
end
