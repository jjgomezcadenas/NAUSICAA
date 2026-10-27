# Recalculate calibrated frame charge and stationary noise with the revised card.
using HDF5,SHA
include(joinpath(@__DIR__,"../src/ion_signal/averaging_scan.jl"))
using .AveragingScan
root=normpath(joinpath(@__DIR__,".."));dir=joinpath(root,"results/ltspice/spice2_enc")
raw=joinpath(dir,"raw.h5")
function interpolate(t,v,x)
    j=clamp(searchsortedlast(t,x),1,length(t)-1)
    v[j]+(v[j+1]-v[j])*(x-t[j])/(t[j+1]-t[j])
end
integral(t,v)=sum(diff(t).*(v[1:end-1].+v[2:end])./2)
h5open(raw,"r") do h
    a=attributes(h);frame=read(a["frame_s"]);T=read(a["averaging_s"]);tau=read(a["correction_tau_s"]);alpha=exp(-frame/tau)
    function wave(n,r)
        local d,b,t,v,tb,vb,delta
        d=h["ions_N$(n)_r$(r)"];b=h["ions_N0_r$(r)"]
        t=read(d["time"]);v=read(d["V(out)"]);tb=read(b["time"]);vb=read(b["V(out)"])
        delta=v.-[interpolate(tb,vb,x) for x in t]
        return t.-.002,delta
    end
    f=read(h["noise_400/frequency"]);asd=read(h["noise_400/V(onoise)"])
    fc=read(h["noise_200/frequency"]);ac=read(h["noise_200/V(onoise)"])
    t,v=wave(1,2);tc,vc=wave(1,1)
    gain=waveform_gain(t,v,T,frame,alpha);sigma=noise_rms(f,asd,T,frame,alpha);enc=sigma/gain
    @assert gain>0
    qerr=abs(noise_rms(f,asd,T,frame,alpha;order=16)/sigma-1)
    serr=abs(noise_rms(fc,ac,T,frame,alpha)/sigma-1)
    terr=abs(waveform_gain(tc,vc,T,frame,alpha)/gain-1)
    @assert qerr<1e-5 && serr<1e-3 && terr<1e-3
    names=sort(filter(n->occursin(r"^V\(m\w+\)$",n),keys(h["noise_400"])))
    flicker=sort(filter(n->occursin("1overf",n),keys(h["noise_400"])))
    @assert length(names)==length(flicker)==7
    spectrum_closure=maximum(abs.(sum(read(h["noise_400"][n]).^2 for n in names)./asd.^2 .-1))
    @assert spectrum_closure<1e-10
    # Integrate each source with precisely the same measurement filter.
    filtered(name)=noise_rms(f,max.(read(h["noise_400"][name]),1e-150),T,frame,alpha)
    device_sigma=[filtered(n) for n in names]; flicker_sigma=[filtered(n) for n in flicker]
    thermal_names=filter(n->occursin(r"^V\(m\w+\.(id|rd|rs)\)$",n),keys(h["noise_400"]))
    thermal_sigma=sqrt(sum(filtered(n)^2 for n in thermal_names))
    total_flicker=sqrt(sum(flicker_sigma.^2))
    integration_closure=abs((total_flicker^2+thermal_sigma^2)/sigma^2-1)
    @assert integration_closure<1e-3
    # Validate finite-path charge against independent analytical weighting potential.
    halfwidth=.00025;gap=.005;e=1.602176634e-19
    expected=e*(1-2/pi*atan(halfwidth^2/(gap*sqrt(gap^2+2halfwidth^2))))
    qinput=integral(read(h["input_time_s"]),read(h["input_current_per_ion_A"]))
    @assert abs(qinput/expected-1)<1e-9
    matches=[match(r"^ions_N(\d+)_r2$",k) for k in keys(h)]
    counts=sort([parse(Int,m[1]) for m in matches if m!==nothing && parse(Int,m[1])>1])
    snrs=Float64[]; refinements=Float64[]; charge_errors=Float64[]
    for n in vcat([1],counts)
        for r in (1,2)
            d=h["ions_N$(n)_r$(r)"]
            push!(charge_errors,abs(integral(read(d["time"]),read(d["I(Is1)"]))/(n*expected)-1))
        end
        if n>1
            d=waveform_gain(wave(n,2)...,T,frame,alpha)
            push!(snrs,d/sigma);push!(refinements,abs(waveform_gain(wave(n,1)...,T,frame,alpha)/d-1))
        end
    end
    @assert maximum(charge_errors)<2e-5
    @assert isempty(refinements) || maximum(refinements)<1e-3
    @assert all(diff(snrs).>0)
    candidate=floor(Int,5enc)+1;pending=Int[];threshold=0
    if isempty(counts)
        pending=[candidate-1,candidate]
    else
        fail=counts[snrs.<=5];pass=counts[snrs.>5]
        if isempty(pass);pending=[maximum(counts)+max(2,ceil(Int,.01candidate))]
        elseif isempty(fail);pending=[max(2,minimum(counts)-max(2,ceil(Int,.01candidate)))]
        elseif minimum(pass)-maximum(fail)==1;threshold=minimum(pass)
        else;pending=[(minimum(pass)+maximum(fail))÷2]
        end
    end
    if threshold>0
        for n in (threshold-1,threshold)
            j=findfirst(==(n),counts)
            @assert refinements[j]+serr+qerr < abs(snrs[j]/5-1)/10
            for r in (1,2)
                value=waveform_gain(wave(n,r)...,T,frame,alpha)/sigma
                @assert (value>5)==(n==threshold)
            end
        end
    end
    previous=joinpath(root,"results/ltspice/averaging_scan.h5")
    previous_gain,previous_enc=h5open(previous,"r") do old
        j=argmin(abs.(read(old["averaging_s"]).-T))
        read(old["gain_V_per_ion"])[j],read(old["enc_electrons_rms"])[j]
    end
    @assert abs(gain/previous_gain-1)<1e-9
    widths=[180.,200.,250.,300.,350.]*1e-6
    scan_gain=[waveform_gain(t,v,w,frame,alpha) for w in widths]
    scan_sigma=[noise_rms(f,asd,w,frame,alpha) for w in widths]
    h5open(joinpath(dir,"analysis.h5"),"w") do o
        attrs=attributes(o)
        for key in keys(a);attrs[key]=read(a[key]);end
        attrs["previous_analysis_sha256"]=bytes2hex(sha256(read(previous)))
        attrs["raw_sha256"]=bytes2hex(sha256(read(raw)))
        attrs["analysis_script_sha256"]=bytes2hex(sha256(read(@__FILE__)))
        attrs["filter_script_sha256"]=bytes2hex(sha256(read(joinpath(root,"src/ion_signal/averaging_scan.jl"))))
        for (key,value) in ("gain_V_per_ion"=>gain,"noise_rms_V"=>sigma,"enc_electrons_rms"=>enc,
            "flicker_rms_V"=>total_flicker,"thermal_rms_V"=>thermal_sigma,
            "flicker_variance_fraction"=>total_flicker^2/sigma^2,"flicker_enc_electrons_rms"=>total_flicker/gain,"thermal_enc_electrons_rms"=>thermal_sigma/gain,"previous_enc_electrons_rms"=>previous_enc,"enc_ratio_to_previous"=>enc/previous_enc,"linear_threshold_ions"=>candidate,
            "verified_threshold_ions"=>threshold,"quadrature_error"=>qerr,"spectral_error"=>serr,
            "transient_error"=>terr,"charge_error"=>maximum(charge_errors),"noise_variance_closure"=>spectrum_closure,
            "filtered_variance_closure"=>integration_closure)
            attrs[key]=value
        end
        o["pending_counts"]=pending;o["checked_counts"]=counts;o["checked_snr"]=snrs;o["count_refinement_error"]=refinements
        o["device_variance_fraction"]=device_sigma.^2 ./sigma^2;o["device_enc_electrons_rms"]=device_sigma./gain;
        o["device_names"]=names;o["device_noise_rms_V"]=device_sigma
        o["flicker_names"]=flicker;o["device_flicker_rms_V"]=flicker_sigma
        o["averaging_s"]=widths;o["scan_gain_V_per_ion"]=scan_gain;o["scan_noise_rms_V"]=scan_sigma;o["scan_enc_electrons_rms"]=scan_sigma./scan_gain
        o["physical_time_s"]=t;o["one_ion_voltage_V"]=v
    end
    println("gain=$gain V/ion; sigma=$sigma V; ENC=$enc electrons RMS; flicker fraction=$(total_flicker^2/sigma^2)")
    println("counts=$counts; S/N=$snrs; verified threshold=$threshold; next=$pending")
end
