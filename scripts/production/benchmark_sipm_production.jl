#!/usr/bin/env julia
using HDF5,TOML,Statistics,Printf,Dates
include("sipm_production.jl")
const Prod=SiPMProduction
const ROOT=Prod.ROOT
const OUT=joinpath(ROOT,"results/sipm_statistical_production_pilot")

"""Verify the saved sensor/ASIC cell correspondence, charge ledger and quiet-cell noise.

Returns per-event PE summaries and aggregate diagnostics. Quiet cells are chosen
by zero input sensor charge, preserving the independent ASIC noise distribution.
All selected event IDs and signed pixel/time indices must agree between stages.
"""
function validate_response(sensor_path,asic_path,ids,s)
    repairs=0;limits=0;maxrepair=0.;batches=Int64(0);cells=Int64(0);clipped=Int64(0)
    noise_n=Int64(0);noise_sum=0.;noise_sumsq=0.
    sensor_totals=Float64[];asic_totals=Float64[];max_pulses=Float64[]
    h5open(sensor_path,"r") do sf
        h5open(asic_path,"r") do af
            read(sf["metadata/complete"]) && read(af["metadata/complete"]) || error("Incomplete shard")
            Int64.(read(sf["metadata/event_ids"]))==ids==Int64.(read(af["metadata/event_ids"])) || error("Event IDs")
            String(read(af["metadata/source_sha256"]))==Prod.sha(sensor_path) || error("ASIC parent hash")
            limit=s.pilot.electronics.current_limit_A*s.pilot.electronics.bucket_width_s/Prod.P.charge_per_pe_C(s.pilot.sensor)
            for id in ids
                sg=sf["events/$id"];ag=af["events/$id"]
                String(read(sg["response_model"]))=="statistical" || error("Wrong sensor model")
                for key in ("event_id","pair_checksum","pixel_ix","pixel_iy","bucket","gas_model")
                    read(sg[key])==read(ag[key]) || error("Stage mismatch: $key")
                end
                q=read(sg["sensor_pe_equiv"]);a=read(ag["readout_pe_equiv"])
                detected=read(sg["detected_pe"]);dark=read(sg["dark_pulses"])
                all(isfinite,a) && all(a.<=limit+1e-8) || error("ASIC range")
                all(0 .<=q .<=detected.+dark.+1e-8) || error("Sensor charge range")
                sum(detected)==read(sg["detected_pe_total"]) || error("Detected PE ledger")
                sum(dark)==read(sg["dark_pulses_total"]) || error("Dark ledger")
                isapprox(sum(q),read(sg["sensor_pe_total"]);rtol=1e-12) || error("Sensor sum")
                isapprox(sum(a),read(ag["readout_pe_total"]);rtol=1e-12) || error("ASIC sum")
                quiet=a[q.==0]
                noise_n+=length(quiet);noise_sum+=sum(quiet);noise_sumsq+=sum(abs2,quiet)
                cells+=length(q);clipped+=count(identity,read(ag["saturated"]))
                repairs+=read(sg["statistical_state_repairs"]);limits+=read(sg["statistical_variance_limits"])
                batches+=read(sg["statistical_batches"]);maxrepair=max(maxrepair,read(sg["maximum_statistical_state_repair"]))
                push!(sensor_totals,sum(q));push!(asic_totals,sum(a));push!(max_pulses,maximum(q;init=0.))
            end
        end
    end
    noise_mean=noise_sum/noise_n;noise_rms=sqrt(noise_sumsq/noise_n-noise_mean^2)
    expected=s.pilot.electronics.noise_rms_pe
    abs(noise_mean)<=6expected/sqrt(noise_n) || error("Quiet-cell noise mean")
    abs(noise_rms-expected)<=6expected/sqrt(2noise_n) || error("Quiet-cell noise RMS")
    (;event_id=ids,sensor_totals,asic_totals,max_pulses,state_repairs=repairs,variance_limits=limits,
      maximum_state_repair=maxrepair,batches,acquired_cells=cells,clipped_cells=clipped,
      clipped_cell_fraction=clipped/cells,quiet_cells=noise_n,noise_mean_pe=noise_mean,noise_rms_pe=noise_rms,
      mean_sensor_total_pe=mean(sensor_totals),mean_asic_total_pe=mean(asic_totals),
      total_charge_retention=sum(asic_totals)/sum(sensor_totals),all_checks_passed=true)
end

"""Time fresh statistical SiPM plus ASIC pilots for both gases at the selected E/P points.

Separate warm-up from sensor simulation/write and ASIC simulation/write. Reuse
completed pilot shards to build the standard production manifest and summary.
The 8,000-event estimate scales measured data processing and finalization; source
resolution is counted once. Benchmark validation and report rendering are separate.
The default tags select E/P=2.5 and 3; output_dir selects the numerical summary directory.
"""
function benchmark_response(;tags=("ep2p5","ep3"),output_dir=OUT)
    mkpath(output_dir);write(joinpath(output_dir,".gitignore"),"*.h5\n")
    data=Dict{String,Any}();manifests=Dict{String,Any}();warmed=false;warm_s=0.
    for tag in tags
        card="config/specs/runs/sipm_bb0nu_statistical_$(tag)_pilot100.toml"
        s=nothing;branches=nothing
        setup_s=@elapsed begin;s=Prod.settings(card);branches=Prod.sources(s);end
        s.pilot.c["response_model"]=="statistical" || error("Statistical card required")
        for b in branches
            isdir(b.sensor_dir) && !isempty(readdir(b.sensor_dir)) && error("Cached response: fresh timing requires a new event selection")
            isdir(b.asic_dir) && !isempty(readdir(b.asic_dir)) && error("Cached ASIC response")
        end
        if !warmed
            scratch=joinpath(ROOT,"work/sipm_production_warmup");mkpath(scratch)
            warm_s=@elapsed mktempdir(scratch) do dir
                b=first(branches)
                warm=merge(b,(;sensor_dir=joinpath(dir,"sensor"),asic_dir=joinpath(dir,"asic")))
                path=Prod.sensor_shard(s,warm,b.ids[1:1],1)
                Prod.asic_shard(s,warm,path,b.ids[1:1],1)
            end
            println("Warm-up: $(round(warm_s;digits=2)) s");flush(stdout);warmed=true
        end
        rows=Dict{String,Any}();processing_s=0.
        for b in branches
            println("$tag / $(b.src.gas): sensor, $(length(b.ids)) events");flush(stdout)
            sensor="";asic=""
            sensor_s=@elapsed sensor=Prod.sensor_shard(s,b,b.ids,1)
            println("$tag / $(b.src.gas): sensor saved in $(round(sensor_s;digits=2)) s; ASIC next");flush(stdout)
            asic_s=@elapsed asic=Prod.asic_shard(s,b,sensor,b.ids,1)
            processing_s+=sensor_s+asic_s
            metrics=validate_response(sensor,asic,b.ids,s)
            rows[b.src.gas]=(;metrics...,events=length(b.ids),sensor_s,asic_s,total_s=sensor_s+asic_s,
                projected_8000_s=(sensor_s+asic_s)*8000/length(b.ids),
                sensor_path=sensor,asic_path=asic,sensor_bytes=filesize(sensor),asic_bytes=filesize(asic),
                projected_8000_bytes=round(Int,(filesize(sensor)+filesize(asic))*8000/length(b.ids)))
            println("$tag / $(b.src.gas): ASIC $(round(asic_s;digits=2)) s; saved-history checks passed");flush(stdout)
        end
        manifest="";finalize_s=@elapsed manifest=Prod.run(card)
        n=length(first(branches).ids)
        electron_card=h5open(first(branches).src.path,"r") do f
            String(read(f["metadata/branch_card_toml"]))
        end
        data[tag]=(;rows,events=n,setup_s,processing_s,finalize_s,
            total_s=setup_s+processing_s+finalize_s,
            projected_8000_s=setup_s+(processing_s+finalize_s)*8000/n,
            projected_8000_minutes=(setup_s+(processing_s+finalize_s)*8000/n)/60,
            run_card=read(s.path,String),sensor_card=read(s.pilot.cp,String),
            electronics_card=read(s.pilot.ep,String),electron_card,
            manifest)
        manifests[tag]=Dict("manifest"=>manifest,"photon_index"=>s.run["photon_index"],"events"=>n)
    end
    path=joinpath(output_dir,"timing.h5")
    h5open(path*".tmp","w") do f
        m=create_group(f,"metadata");m["complete"]=false
        m["schema"]="nausicaa.sipm_production_timing.v1";m["sample"]="bb0nu"
        m["created_utc"]=string(now(UTC));m["warmup_s"]=warm_s;m["benchmark_source"]=read(@__FILE__,String)
        for (tag,row) in data
            g=create_group(f,tag)
            for (k,v) in pairs(row)
                k==:rows && continue
                g[string(k)]=v
            end
            for (gas,metrics) in row.rows
                gg=create_group(g,gas)
                for (k,v) in pairs(metrics);gg[string(k)]=v;end
            end
        end
        m["projected_all_8000_minutes"]=sum(row.projected_8000_minutes for row in values(data))
        write(m["complete"],true)
    end
    mv(path*".tmp",path;force=true)
    open(joinpath(output_dir,"files.toml"),"w") do io
        TOML.print(io,Dict("summary"=>path,"operating_points"=>manifests))
    end
    println("Timing summary: $path")
end
if abspath(PROGRAM_FILE)==@__FILE__;benchmark_response();end
