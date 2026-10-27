# Validate spectra and reduce benchmark numbers in Julia; preserve SI units.
using HDF5, SHA
root=normpath(joinpath(@__DIR__,"..")); dir=joinpath(root,"results/ltspice/spice2_reproduction")
raw=joinpath(dir,"raw.h5"); out=joinpath(dir,"analysis.h5")
h5open(raw,"r") do input
    h5open(out,"w") do result
        attributes(result)["raw_sha256"]=bytes2hex(sha256(read(raw)))
        attributes(result)["script_sha256"]=bytes2hex(sha256(read(@__FILE__)))
        attributes(result)["conditions"]=read(attributes(input)["conditions"])
        for case in ("updated_1Hz","updated_0p1Hz","previous_0p1Hz")
            src=input[case]; f=read(src["frequency"]); total=read(src["V(onoise)"])
            names=filter(n->occursin(r"^V\(m\w+\)$",n),keys(src))
            flicker=filter(n->occursin("1overf",n),keys(src))
            @assert length(names)==length(flicker)==7
            @assert all(diff(f).>0) && all(total.>0)
            summed=sum(read(src[n]).^2 for n in names)
            closure=maximum(abs.(summed./total.^2 .-1))
            @assert closure<1e-10
            # Each device total must also equal the sum of its mechanisms.
            component_error=0.
            for name in names
                prefix=name[1:end-1]*"."
                components=filter(n->startswith(n,prefix),keys(src))
                dev=read(src[name]); mech=sum(read(src[n]).^2 for n in components)
                component_error=max(component_error,maximum(abs.(mech./dev.^2 .-1)))
            end
            @assert component_error<1e-10
            flicker_power=sum(read(src[n]).^2 for n in flicker)
            if startswith(case,"updated")
                @assert all(maximum(read(src[n]))>0 for n in flicker)
            else
                @assert all(iszero,flicker_power)
            end
            targets=first(f)<.5 ? [.1,1.,10.,100.,1e3,1e6,1e9] : [1.,10.,100.,1e3,1e6,1e9]
            indices=[argmin(abs.(f.-x)) for x in targets]
            dest=create_group(result,case)
            dest["frequency_Hz"]=f[indices]
            dest["output_ASD_V_sqrtHz"]=total[indices]
            dest["input_ASD_A_sqrtHz"]=read(src["inoise"])[indices]
            dest["gain_ohm"]=read(src["gain"])[indices]
            dest["flicker_ASD_V_sqrtHz"]=sqrt.(flicker_power[indices])
            dest["flicker_variance_fraction"]=flicker_power[indices]./total[indices].^2
            dest["M1_flicker_ASD_V_sqrtHz"]=read(src["V(m1.1overf)"])[indices]
            dest["M1_channel_ASD_V_sqrtHz"]=read(src["V(m1.id)"])[indices]
            attributes(dest)["device_variance_closure"]=closure
            attributes(dest)["mechanism_variance_closure"]=component_error
            println(case)
            for j in indices
                println("$(f[j]) Hz: total $(total[j]) V/sqrt(Hz), flicker fraction $(flicker_power[j]/total[j]^2)")
            end
        end
        # Noise parameter changes should leave DC and signal transfer unchanged.
        op_error=maximum(abs.(read(input["updated_op"][n]).-read(input["previous_op"][n]))[1] for n in keys(input["updated_op"]))
        gain_error=maximum(abs.(read(input["updated_0p1Hz/gain"])./read(input["previous_0p1Hz/gain"]).-1))
        @assert op_error<1e-12 && gain_error<1e-10
        attributes(result)["maximum_op_absolute_change"]=op_error
        attributes(result)["maximum_gain_relative_change"]=gain_error
        # Independently started sweeps agree at common frequencies.
        fwide=read(input["updated_0p1Hz/frequency"]); fnarrow=read(input["updated_1Hz/frequency"])
        ix=[argmin(abs.(fwide.-f)) for f in fnarrow]
        sweep_error=maximum(abs.(read(input["updated_0p1Hz/V(onoise)"])[ix]./read(input["updated_1Hz/V(onoise)"]).-1))
        @assert sweep_error<1e-8
        attributes(result)["overlapping_sweep_relative_difference"]=sweep_error
        dest=create_group(result,"operating_point")
        for n in keys(input["updated_op"]);dest[n]=read(input["updated_op"][n]);end
    end
end
println(out)
