#!/usr/bin/env julia
using HDF5,TOML,Statistics,Random,SHA,Printf
include("sipm_pipeline.jl")
const ROOT=normpath(joinpath(@__DIR__,"../.."));const P=SiPMPipeline.P
const OUT=joinpath(ROOT,"results/sipm_statistical_validation")

"""Return first/second moments, sampling errors and successive-bucket covariance."""
function moments(values)
    n=size(values,2);totals=vec(sum(values;dims=1));mu=vec(mean(values;dims=2))
    vars=vec(var(values;dims=2,corrected=false));centered=values.-mu
    variance_se=sqrt.(max.(0.,vec(mean(centered.^4;dims=2)).-vars.^2)./n)
    covs=[mean(centered[i,:].*centered[i+1,:]) for i in 1:length(mu)-1]
    covse=[std(centered[i,:].*centered[i+1,:];corrected=false)/sqrt(n) for i in 1:length(mu)-1]
    sigma=std(totals;corrected=false);v=sigma^2
    (;mean=mu,variance=vars,mean_se=sqrt.(vars./n),variance_se,
      adjacent_covariance=covs,adjacent_covariance_se=covse,
      total_mean=mean(totals),total_variance=v,total_mean_se=sigma/sqrt(n),
      total_variance_se=sqrt(max(0.,mean((totals.-mean(totals)).^4)-v^2)/n),
      total_relative_rms=sigma/mean(totals))
end

"""Benchmark independent channel replicas, preserving fixed incident counts."""
function ensemble(sensor,counts,repetitions,model;steps=64)
    result=zeros(length(counts),repetitions);repair=0;limits=0;maxrepair=0.;batches=0
    elapsed=@elapsed for i in 1:repetitions
        rng=Xoshiro(1000+i)
        trace=model=="microscopic" ? P.simulate_channel(rng,counts,250e-9,sensor) :
            P.simulate_channel_statistical(rng,counts,250e-9,sensor;options=P.StatisticalOptions(steps))
        result[:,i]=trace.charge_C./P.charge_per_pe_C(sensor)
        if model!="microscopic"
            repair+=trace.state_repairs;limits+=trace.variance_limits
            maxrepair=max(maxrepair,trace.maximum_state_repair);batches+=trace.batches
        end
    end
    (;moments(result)...,elapsed_s=elapsed,repetitions,state_repairs=repair,
      variance_limits=limits,maximum_state_repair=maxrepair,batches)
end

"""Store a numerical benchmark row with independently regenerable provenance."""
function save(g,results)
    for (k,v) in pairs(results);g[string(k)]=v;end
end

"""Validate fixed-count transient histories at three resolutions and real photon histories.

Synthetic ensembles test means, variances, adjacent covariance and total charge.
Real events use both gases and both new E/P photon files; SiPM production jobs
are not resumed. The first event is repeated to test the total-charge spread.
"""
function main(;repetitions=2000,event_repetitions=256)
    mkpath(OUT);write(joinpath(OUT,".gitignore"),"*.h5\n")
    physical=SiPMPipeline.settings(ROOT,"config/specs/runs/sipm_response_ep3.toml")
    s=P.Sensor(1.,physical.sensor.microcells,physical.sensor.gain,physical.sensor.recovery_s,0.)
    cases=[("few",Int64[0,16,0,16,0]),("moderate",Int64[0,1000,1000,0,1000]),
        ("bright",Int64[0,20000,20000,0,20000]),("very_bright",Int64[0,100000,100000,0,100000]),
        ("transient",Int64[1000,50000,1000,0,0,10000])]
    # Compile both paths before timing.
    ensemble(s,[100],1,"microscopic");ensemble(s,[100],1,"statistical")
    h5open(joinpath(OUT,"validation.h5"),"w") do out
        attributes(out)["schema"]="nausicaa.sipm_statistical_validation.v1"
        meta=create_group(out,"metadata");meta["complete"]=false
        meta["sensor_card_toml"]=read(physical.cp,String)
        for path in ("src/ion_signal/sipm_statistical.jl","src/ion_signal/sipm_response.jl",relpath(@__FILE__,ROOT))
            g=create_group(meta,replace(path,"/"=>"_"));g["source"]=read(joinpath(ROOT,path),String);g["sha256"]=bytes2hex(sha256(read(joinpath(ROOT,path))))
        end
        for (name,counts) in cases
            println("Synthetic $name");flush(stdout)
            g=create_group(out,"synthetic/$name");g["incident_counts"]=counts
            nrep=name=="few" ? 20000 : repetitions
            save(create_group(g,"microscopic"),ensemble(s,counts,nrep,"microscopic"))
            for steps in (32,64,128)
                save(create_group(g,"steps_$steps"),ensemble(s,counts,nrep,"statistical";steps))
            end
        end
        dark_sensor=P.Sensor(0.,s.microcells,s.gain,s.recovery_s,physical.sensor.dark_rate_Hz)
        g=create_group(out,"synthetic/dark");counts=zeros(Int64,20);g["incident_counts"]=counts
        save(create_group(g,"microscopic"),ensemble(dark_sensor,counts,repetitions,"microscopic"))
        for steps in (32,64,128);save(create_group(g,"steps_$steps"),ensemble(dark_sensor,counts,repetitions,"statistical";steps));end
        for tag in ("ep2p5","ep3")
            run=TOML.parsefile(joinpath(ROOT,"config/specs/runs/sipm_response_$tag.toml"))
            for gas in ("pure_xe","xe_he15")
                println("Real event $tag/$gas");flush(stdout)
                path=h5open(joinpath(ROOT,run["photon_index"]),"r") do f;String(read(f["bb0nu/$(gas)_path"]));end
                g=create_group(out,"events/$(tag)_$gas");g["photon_path"]=path
                g["photon_sha256"]=SiPMPipeline.sha(path)
                history=h5open(path,"r") do f;SiPMPipeline.EA.read_photon_history(f["events/0"]);end
                grid=SiPMPipeline.acquisition(history,physical.c["acquisition"])
                g["incident_photons"]=sum(grid.counts);g["event_id"]=0
                for (model,steps) in (("microscopic",0),("statistical",64),("statistical",128))
                    totals=zeros(1,event_repetitions);repairs=0;limits=0;maxrepair=0.;batches=0
                    times=zeros(event_repetitions)
                    elapsed=@elapsed for rep in 1:event_repetitions
                        times[rep]=@elapsed begin
                        for j in axes(grid.counts,2)
                            rng=Xoshiro(100000rep+j)
                            trace=model=="microscopic" ? P.simulate_channel(rng,view(grid.counts,:,j),250e-9,physical.sensor) :
                                P.simulate_channel_statistical(rng,view(grid.counts,:,j),250e-9,physical.sensor;options=P.StatisticalOptions(steps))
                            totals[1,rep]+=sum(trace.charge_C)/P.charge_per_pe_C(physical.sensor)
                            if model!="microscopic"
                                repairs+=trace.state_repairs;limits+=trace.variance_limits
                                maxrepair=max(maxrepair,trace.maximum_state_repair);batches+=trace.batches
                            end
                        end
                        end
                    end
                    name=model=="microscopic" ? model : "steps_$steps"
                    save(create_group(g,name),(;moments(totals)...,elapsed_s=elapsed,
                        median_seconds_per_replica=median(times),replica_seconds=times,repetitions=event_repetitions,
                        state_repairs=repairs,variance_limits=limits,maximum_state_repair=maxrepair,batches))
                end
            end
        end
        write(meta["complete"],true)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
