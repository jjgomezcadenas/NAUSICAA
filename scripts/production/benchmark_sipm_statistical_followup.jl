#!/usr/bin/env julia

include("benchmark_sipm_statistical.jl")

"""Repeat a complete fixed-photon event with independently seeded channels.

Returns charge moments, detected-count moments, the analytical binomial count
variance, raw event totals and timing/feasibility diagnostics.
"""
function event_replicas(grid,sensor,model;repetitions=512,steps=128)
    q=zeros(1,repetitions);detected=zeros(repetitions);times=zeros(repetitions)
    repairs=0;limits=0;maxrepair=0.;batches=0
    elapsed=@elapsed for rep in 1:repetitions
        times[rep]=@elapsed for j in axes(grid.counts,2)
            rng=Xoshiro(1_000_000_000_000+100000rep+j)
            t=model=="microscopic" ? P.simulate_channel(rng,view(grid.counts,:,j),250e-9,sensor) :
                P.simulate_channel_statistical(rng,view(grid.counts,:,j),250e-9,sensor;options=P.StatisticalOptions(steps))
            q[1,rep]+=sum(t.charge_C)/P.charge_per_pe_C(sensor)
            detected[rep]+=sum(t.detected_pe)
            if model!="microscopic"
                repairs+=t.state_repairs;limits+=t.variance_limits
                maxrepair=max(maxrepair,t.maximum_state_repair);batches+=t.batches
            end
        end
    end
    (;moments(q)...,elapsed_s=elapsed,median_seconds_per_replica=median(times),replica_seconds=times,
      repetitions,state_repairs=repairs,variance_limits=limits,maximum_state_repair=maxrepair,batches,
      charge_totals=vec(q),detected_totals=detected,detected_mean=mean(detected),
      detected_variance=var(detected;corrected=false),
      binomial_mean=sum(grid.counts)*sensor.pde,binomial_variance=sum(grid.counts)*sensor.pde*(1-sensor.pde))
end

"""Independent-seed follow-up of event variance and explicit binomial-PDE tests.

Append these checks to the completed primary benchmark. The saved first event
at E/P=2.5 in pure Xe had the largest variance discrepancy in the initial sample;
512 new replicas separate that fluctuation from a reproducible model bias.
"""
function followup()
    settings=SiPMPipeline.settings(ROOT,"config/specs/runs/sipm_response_ep2p5.toml")
    sensor=P.Sensor(.3,settings.sensor.microcells,settings.sensor.gain,settings.sensor.recovery_s,0.)
    h5open(joinpath(OUT,"validation.h5"),"r+") do out
        read(out["metadata/complete"]) || error("Primary benchmark must complete first")
        write(out["metadata/complete"],false)
        meta=out["metadata"]
        haskey(meta,"followup") && delete_object(meta,"followup")
        provenance=create_group(meta,"followup")
        provenance["source"]=read(@__FILE__,String);provenance["sha256"]=SiPMPipeline.sha(@__FILE__)
        for (name,counts) in (("pde_moderate",Int64[0,1000,1000,0,1000]),("pde_bright",Int64[0,100000,100000,0,100000]))
            println("PDE check $name");flush(stdout)
            key="synthetic/$name";haskey(out,key) && delete_object(out,key)
            g=create_group(out,key);g["incident_counts"]=counts
            save(create_group(g,"microscopic"),ensemble(sensor,counts,10000,"microscopic"))
            save(create_group(g,"steps_128"),ensemble(sensor,counts,10000,"statistical";steps=128))
        end
        path=h5open(joinpath(ROOT,settings.run["photon_index"]),"r") do f
            String(read(f["bb0nu/pure_xe_path"]))
        end
        history=h5open(path,"r") do f;SiPMPipeline.EA.read_photon_history(f["events/0"]);end
        grid=SiPMPipeline.acquisition(history,settings.c["acquisition"])
        key="events/ep2p5_pure_xe_independent";haskey(out,key) && delete_object(out,key)
        g=create_group(out,key);g["photon_path"]=path;g["photon_sha256"]=SiPMPipeline.sha(path)
        g["incident_photons"]=sum(grid.counts);g["event_id"]=0
        for model in ("microscopic","statistical")
            println("Independent real-event check $model");flush(stdout)
            key=model=="microscopic" ? model : "steps_128"
            save(create_group(g,key),event_replicas(grid,settings.sensor,model))
        end
        write(meta["complete"],true)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;followup();end
