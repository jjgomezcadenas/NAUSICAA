using HDF5,Statistics,SHA
"""Compare saved independent ensembles with their measured sampling errors.

Mean and variance z-scores use the quadrature sum of Monte Carlo errors. The
numerical review budget is 0.05% mean bias and 5% variance discrepancy, or four
Monte Carlo standard errors when larger. These budgets describe this pilot.
"""
function compare(path)
    h5open(path,"r+") do f
        read(f["metadata/complete"]) || error("Incomplete benchmark")
        haskey(f,"comparison") && delete_object(f,"comparison")
        out=create_group(f,"comparison")
        for category in ("synthetic","events"),name in keys(f[category])
            input=f["$category/$name"];ref=input["microscopic"]
            for key in filter(k->startswith(k,"steps_"),keys(input))
                s=input[key];g=create_group(out,"$(category)_$(name)_$key")
                rmean=read(ref["mean"]);sm=read(s["mean"])
                rv=read(ref["variance"]);sv=read(s["variance"])
                meanse=sqrt.(read(ref["mean_se"]).^2+read(s["mean_se"]).^2)
                varse=sqrt.(read(ref["variance_se"]).^2+read(s["variance_se"]).^2)
                covref=read(ref["adjacent_covariance"]);covs=read(s["adjacent_covariance"])
                covse=sqrt.(read(ref["adjacent_covariance_se"]).^2+read(s["adjacent_covariance_se"]).^2)
                covbudget=4covse.+.05sqrt.(rv[1:end-1].*rv[2:end])
                totalrel=read(s["total_mean"])/read(ref["total_mean"])-1
                vr=read(s["total_variance"])/read(ref["total_variance"])
                reftime=haskey(ref,"median_seconds_per_replica") ? read(ref["median_seconds_per_replica"]) : read(ref["elapsed_s"])/read(ref["repetitions"])
                stime=haskey(s,"median_seconds_per_replica") ? read(s["median_seconds_per_replica"]) : read(s["elapsed_s"])/read(s["repetitions"])
                speed=reftime/stime
                passmean=all(abs.(sm-rmean).<=max.(.0005abs.(rmean),4meanse).+1e-10)
                passvar=all(abs.(sv-rv).<=max.(.05rv,4varse).+1e-10)
                passcov=all(abs.(covs-covref).<=covbudget.+1e-10)
                g["total_mean_bias_percent"]=100totalrel;g["total_variance_ratio"]=vr;g["speedup"]=speed
                g["total_mean_bias_se_percent"]=100sqrt(read(ref["total_mean_se"])^2+read(s["total_mean_se"])^2)/read(ref["total_mean"])
                g["total_variance_ratio_se"]=sqrt(read(s["total_variance_se"])^2+vr^2*read(ref["total_variance_se"])^2)/read(ref["total_variance"])
                g["reference_seconds"]=reftime
                g["statistical_seconds"]=stime
                g["mean_pass"]=passmean;g["variance_pass"]=passvar;g["covariance_pass"]=passcov
                g["state_repairs"]=read(s["state_repairs"]);g["variance_limits"]=read(s["variance_limits"])
                g["maximum_state_repair"]=read(s["maximum_state_repair"]);g["statistical_batches"]=read(s["batches"])
                if haskey(s,"detected_totals")
                    for (label,model) in (("reference",ref),("statistical",s))
                        n=read(model["detected_totals"]);v=var(n;corrected=false)
                        se=sqrt(max(0.,mean((n.-mean(n)).^4)-v^2)/length(n))
                        predicted=read(model["binomial_variance"])
                        g["$(label)_detected_variance_ratio"]=v/predicted
                        g["$(label)_detected_variance_ratio_se"]=se/predicted
                        # A fixed, common control-variate coefficient changes
                        # estimator precision, not its expected charge variance.
                        q=read(model["charge_totals"]);nr=length(q)
                        score=(q .- mean(q)).^2 .* (nr/(nr-1)) .-
                            0.5 .* ((n .- read(model["binomial_mean"])).^2 .- predicted)
                        g["$(label)_controlled_variance"]=mean(score)
                        g["$(label)_controlled_variance_se"]=std(score)/sqrt(nr)
                    end
                    g["variance_control_coefficient"]=.5
                    cr=read(g["reference_controlled_variance"]);cs=read(g["statistical_controlled_variance"])
                    ratio=cs/cr
                    g["controlled_variance_ratio"]=ratio
                    g["controlled_variance_ratio_se"]=sqrt(read(g["statistical_controlled_variance_se"])^2+
                        ratio^2*read(g["reference_controlled_variance_se"])^2)/cr
                end
                println(category,"/",name," ",key," bias=",round(100totalrel;digits=5),"% varratio=",round(vr;digits=3)," speedup=",round(speed;digits=2)," pass=",(passmean,passvar,passcov)," repairs=",read(s["state_repairs"])," limits=",read(s["variance_limits"]))
            end
        end
        haskey(f,"default_summary") && delete_object(f,"default_summary")
        summary=create_group(f,"default_summary")
        selected=[out[k] for k in keys(out) if endswith(k,"steps_128")]
        eventrows=[out[k] for k in keys(out) if startswith(k,"events_") && endswith(k,"steps_128")]
        summary["substeps"]=128
        summary["state_repairs"]=sum(read(g["state_repairs"]) for g in selected)
        summary["statistical_batches"]=sum(read(g["statistical_batches"]) for g in selected)
        summary["maximum_state_repair"]=maximum(read(g["maximum_state_repair"]) for g in selected)
        summary["variance_limits"]=sum(read(g["variance_limits"]) for g in selected)
        summary["max_event_mean_bias_percent"]=maximum(abs(read(g["total_mean_bias_percent"])) for g in eventrows)
        summary["min_event_speedup"]=minimum(read(g["speedup"]) for g in eventrows)
        summary["max_event_speedup"]=maximum(read(g["speedup"]) for g in eventrows)
        meta=f["metadata"]
        haskey(meta,"reduction") && delete_object(meta,"reduction")
        provenance=create_group(meta,"reduction")
        provenance["source"]=read(@__FILE__,String)
        provenance["sha256"]=bytes2hex(sha256(read(@__FILE__)))
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    compare(isempty(ARGS) ? "results/sipm_statistical_validation/validation.h5" : only(ARGS))
end
