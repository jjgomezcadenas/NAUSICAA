#!/usr/bin/env julia
include("run_electron_250ns.jl")
const PHOTON_BENCH_OUT=joinpath(ROOT,"results/photon_production_pilot")

"""Read every saved event and check photon sums, identities and the electron ledger."""
function validate_photon_pilot(path,gas,n)
    h5open(path,"r") do f
        m=f["metadata"]
        read(m["complete"]) && read(m["event_count"])==n || error("Incomplete pilot")
        String(read(m["signal_unit"]))=="incident_photons" || error("Photon units")
        length(keys(f["events"]))==n || error("Event count")
        total=Int64(0)
        for id in keys(f["events"])
            h=read_photon_history(f["events/$id"])
            h.placement.event_id==parse(Int64,id) && h.gas_model==Symbol(gas) || error("Event identity")
            all(>=(0),h.cells.count) || error("Negative photon count")
            total+=sum(h.cells.count)
        end
        total==read(m["incident_photons"]) || error("Photon sum")
        read(m["surviving_electrons"])<=read(m["initial_pairs"]) &&
            total<=read(m["emitted_photons"]) || error("Charge/photon ledger")
    end
    true
end

"""Write a TOML index of the measured pilot files and matching full photon files."""
function write_photon_pilot_index(summary::AbstractString)
    entries=Dict{String,Any}()
    h5open(summary,"r") do f
        read(f["metadata/complete"]) || error("Incomplete timing summary")
        for tag in ("ep2p5","ep3")
            e=TOML.parse(String(read(f["$tag/electron_card_toml"])))
            entry=Dict{String,Any}("events"=>read(f["$tag/events"]),
                "reduced_field_kV_cm_bar"=>e["reduced_field_kV_cm_bar"],
                "bucket_width_us"=>e["bucket_width_us"])
            for gas in ("pure_xe","xe_he15")
                g=f["$tag/$gas"]
                paths=Dict{String,Any}("pilot"=>String(read(g["pilot_path"])),
                    "pilot_sha256"=>String(read(g["pilot_sha256"])))
                if read(g["matching_full_exists"])
                    paths["full"]=String(read(g["full_path"]))
                    paths["full_events"]=read(f["$tag/target_events"])
                end
                entry[gas]=paths
            end
            entries[tag]=entry
        end
    end
    path=joinpath(dirname(summary),"files.toml")
    open(path,"w") do io
        TOML.print(io,Dict("schema"=>"nausicaa.photon_production_pilot_index.v1",
            "summary"=>String(summary),"sample"=>"bb0nu","operating_points"=>entries))
    end
    path
end

"""Benchmark fresh bb0nu photon production for both gases at the two current E/P settings.

Warm up Julia separately, then time 500 events per operating point by default.
The existing production runner supplies the physics, seeds, storage and checks.
Save pilot photon files and measured/extrapolated timings under results/.
Completed matching full productions are inspected and linked in the summary.
A cached pilot is rejected because it would not measure generation time; choose
another event limit to make a new measurement.
"""
function benchmark_photons(n::Int=500)
    mkpath(PHOTON_BENCH_OUT);write(joinpath(PHOTON_BENCH_OUT,".gitignore"),"*.h5\n")
    hash=electron_replay_code_hash()
    specs=[load_specs(ROOT,"config/specs/runs/electron_250ns_$tag.toml") for tag in ("ep2p5","ep3")]
    for spec in specs
        1<n<Int(spec.realisation.data["max_events_per_sample"]) || error("Pilot event limit")
        p=electron_paths(spec,"bb0nu",hash,n)
        any(isfile,(p.pure_xe,p.xe_he15)) && error("Pilot already cached; choose another event limit")
    end
    warmn=findfirst(1:30) do count
        p=electron_paths(first(specs),"bb0nu",hash,count)
        !any(isfile,(p.pure_xe,p.xe_he15))
    end
    isnothing(warmn) && error("Warm-up limits already cached")
    println("Warming up on $warmn events, both gases");flush(stdout)
    warm_s=@elapsed run_electron_sample(first(specs),"bb0nu",hash,warmn)
    summary=joinpath(PHOTON_BENCH_OUT,"timing.h5")
    h5open(summary*".tmp","w") do f
        m=create_group(f,"metadata");m["complete"]=false
        m["schema"]="nausicaa.photon_production_timing.v1"
        m["sample"]="bb0nu";m["pilot_events"]=n;m["warmup_events"]=warmn;m["warmup_s"]=warm_s
        m["replay_code_hash"]=hash;m["benchmark_source"]=read(@__FILE__,String)
        m["julia_version"]=string(VERSION);m["created_utc"]=string(now(UTC))
        for (tag,spec) in zip(("ep2p5","ep3"),specs)
            result=nothing
            wall=@elapsed result=run_electron_sample(spec,"bb0nu",hash,n)
            target=Int(spec.realisation.data["max_events_per_sample"])
            full=electron_paths(spec,"bb0nu",hash,target)
            g=create_group(f,tag)
            g["events"]=n;g["target_events"]=target;g["elapsed_s"]=wall
            g["seconds_per_event_both_gases"]=wall/n;g["estimated_full_s"]=wall*target/n
            g["estimated_full_minutes"]=wall*target/n/60
            g["estimated_cold_full_s"]=warm_s+wall*target/n
            g["run_card_toml"]=read(spec.run_path,String)
            g["electron_card_toml"]=read(spec.electrons.path,String)
            g["realisation_card_toml"]=read(spec.realisation.path,String)
            case_dir=joinpath(PHOTON_BENCH_OUT,tag);mkpath(case_dir)
            for gas in ("pure_xe","xe_he15")
                path=getproperty(result.paths,Symbol(gas))
                validate_photon_pilot(path,gas,n)
                saved=joinpath(case_dir,"photons_$(gas)_$(n).h5")
                cp(path,saved;force=false)
                gg=create_group(g,gas);gg["pilot_path"]=saved
                gg["pilot_sha256"]=content_hash(saved;file=true)
                gg["pilot_bytes"]=filesize(saved);gg["estimated_full_bytes"]=round(Int,filesize(saved)*target/n)
                fullpath=getproperty(full,Symbol(gas));gg["matching_full_exists"]=isfile(fullpath)
                if isfile(fullpath)
                    totals=completed_photon_totals(fullpath,full,gas,target)
                    gg["full_path"]=fullpath;gg["full_recorded_elapsed_s"]=totals.elapsed_s
                    gg["full_bytes"]=filesize(fullpath)
                end
            end
            println("$tag: $n events, BOTH gases, $(round(wall;digits=3)) s; estimate $target events $(round(wall*target/n;digits=1)) s")
            flush(stdout)
        end
        write(m["complete"],true)
    end
    mv(summary*".tmp",summary;force=true)
    write_photon_pilot_index(summary)
    println("Timing summary: $summary")
end
if abspath(PROGRAM_FILE)==@__FILE__
    benchmark_photons(isempty(ARGS) ? 500 : parse(Int,only(ARGS)))
end
