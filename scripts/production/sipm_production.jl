"""Sharded, restartable SiPM and ASIC production with PE-per-bucket outputs."""
module SiPMProduction
using HDF5, TOML, SHA, Random, Statistics
include("sipm_pipeline.jl")
const Pilot=SiPMPipeline
const S=Pilot.S
const P=Pilot.P
const E=Pilot.E
const EA=Pilot.EA
const SENSOR_SCHEMA="nausicaa.sipm_pe_file.v1"
const ASIC_SCHEMA="nausicaa.sipm_readout_pe_file.v1"
const ROOT=normpath(joinpath(@__DIR__,"../.."))
sha(path)=open(io->bytes2hex(sha256(io)),path)
str(g,key)=String(read(g[key]))
schema(h)=String(read(attributes(h)["schema"]))
require(ok,message)=ok || throw(ArgumentError(message))
resolve(path)=normpath(isabspath(path) ? path : joinpath(ROOT,path))

"""Read a production card and reuse the validated physical pilot model."""
function settings(card::AbstractString)
    path=resolve(card);run=TOML.parsefile(path)
    require(run["schema"]=="nausicaa.sipm_production_run.v1","Production card schema")
    require(run["event_selection"] in ("all_sorted_source_ids","first_n_sorted_source_ids"),
            "Event selection policy")
    require(run["shard_events"]>0,"Shard size must be positive")
    require(run["sample"] in ("bb0nu","single_electrons"),"Unknown sample")
    require(!isempty(run["gases"]) && allunique(run["gases"]) &&
            all(g->g in ("pure_xe","xe_he15"),run["gases"]),"Unknown gas selection")
    pilot=Pilot.settings(ROOT,run["pilot_card"])
    require(resolve(run["photon_index"])==resolve(pilot.run["photon_index"]),
            "Production and pilot photon indices differ")
    codehash=sha(@__FILE__)
    (;path,run,pilot,codehash)
end

"""Resolve source IDs and preserve the pilot's selection-independent RNG stream."""
function sources(s)
    output=NamedTuple[]
    for src in Pilot.sources(s.pilot)
        src.sample==s.run["sample"] && src.gas in s.run["gases"] || continue
        ids=h5open(src.path,"r") do h
            require(schema(h)=="nausicaa.photon_sipm_file.pilot.v1" &&
                    Bool(read(h["metadata/complete"])),"Incomplete photon source")
            found=sort!(parse.(Int64,collect(keys(h["events"]))))
            require(length(found)==read(h["metadata/event_count"]),"Source event count")
            found
        end
        if s.run["event_selection"]=="first_n_sorted_source_ids"
            n=Int(s.run["max_events"])
            require(0<n<=length(ids),"Invalid benchmark event limit")
            ids=ids[1:n]
        end
        idhash=S.content_hash(S.canonical(ids))
        sensorid=S.content_hash(join((src.stream,idhash,s.codehash,SENSOR_SCHEMA),"|"))
        asicid=S.content_hash(join((sensorid,s.pilot.electronics_model,s.codehash,ASIC_SCHEMA),"|"))
        base=dirname(dirname(dirname(dirname(src.path))))
        parent=joinpath(base,"response","electrons",basename(dirname(src.path)))
        sensor_dir=joinpath(parent,"sensor_pe",sensorid)
        asic_dir=joinpath(parent,"asic_pe",asicid)
        push!(output,(;src,ids,idhash,sensorid,asicid,sensor_dir,asic_dir))
    end
    require(length(output)==length(s.run["gases"]),"Missing gas source")
    output
end

function put(g,values)
    for (key,value) in pairs(values);g[string(key)]=value;end
end
function vector(g,key,value)
    if isempty(value);g[key]=value;else;g[key,compress=4]=value;end
end

"""Complete shards are immutable; interrupted temporary files may be replaced."""
function shard_ready(path,expected_schema,stage_id,ids;upstream_hash=nothing)
    isfile(path) || return false
    h5open(path,"r") do h
        require(schema(h)==expected_schema && Bool(read(h["metadata/complete"])),
                "Incomplete or incompatible shard: $path")
        m=h["metadata"]
        require(str(m,"stage_id")==stage_id && Int64.(read(m["event_ids"]))==ids &&
                length(keys(h["events"]))==length(ids),"Shard identity/event mismatch: $path")
        if upstream_hash!==nothing
            require(str(m,"source_sha256")==upstream_hash,"Upstream shard changed: $path")
        end
    end
    true
end

function atomic_shard(writefile,path,expected_schema,stage_id,ids)
    mkpath(dirname(path));temp=path*".tmp"
    isfile(temp) && rm(temp)
    h5open(temp,"w") do h
        writefile(h)
        write(h["metadata/complete"],true)
    end
    shard_ready(temp,expected_schema,stage_id,ids)
    mv(temp,path;force=true)
    path
end

"""Simulate selected SiPM channels and save integer and charge-equivalent PE."""
function sensor_shard(s,branch,ids,number)
    path=joinpath(branch.sensor_dir,"shard_$(lpad(number,4,'0')).h5")
    shard_ready(path,SENSOR_SCHEMA,branch.sensorid,ids) && return path
    src=branch.src;pilot=s.pilot;q0=P.charge_per_pe_C(pilot.sensor)
    atomic_shard(path,SENSOR_SCHEMA,branch.sensorid,ids) do out
        attributes(out)["schema"]=SENSOR_SCHEMA
        m=create_group(out,"metadata")
        put(m,(;complete=false,stage_id=branch.sensorid,stream_identity=src.stream,
            source_path=src.path,source_sha256=src.sourcehash,sample=src.sample,gas_model=src.gas,
            event_ids=ids,event_count=length(ids),bucket_width_s=src.dt,
            signal_unit="PE",time_unit="s",charge_per_pe_C=q0,
            sensor_card_toml=read(pilot.cp,String),sensor_code_hash=pilot.sensor_code,
            production_code_hash=s.codehash,julia_version=string(VERSION)))
        events=create_group(out,"events")
        h5open(src.path,"r") do input
            for id in ids
                original=input["events/$id"];history=EA.read_photon_history(original)
                require(history.placement.event_id==id && string(history.gas_model)==src.gas,
                        "Photon event identity")
                grid=Pilot.acquisition(history,pilot.c["acquisition"])
                nb=length(grid.buckets);nc=length(grid.channels);rows=nb*nc
                detected=zeros(Int64,rows);dark=zeros(Int64,rows);sensor_pe=zeros(Float64,rows)
                seeds=zeros(UInt64,nc)
                for (j,ch) in enumerate(grid.channels)
                    seed=S.event_seed(src.stream,id,"sensor:$(ch[1]):$(ch[2]):$(pilot.c["seed"])")
                    seeds[j]=seed
                    trace=P.simulate_channel(Xoshiro(seed),view(grid.counts,:,j),src.dt,pilot.sensor)
                    selection=(j-1)*nb+1:j*nb
                    detected[selection]=trace.detected_pe
                    dark[selection]=trace.dark_pulses
                    sensor_pe[selection]=trace.charge_C./q0
                end
                require(all(detected.<=vec(grid.counts)) && all(sensor_pe.>=0) &&
                        all(sensor_pe.<=detected.+dark.+1e-8),"SiPM ledger")
                g=create_group(events,string(id));attributes(g)["schema"]="nausicaa.sipm_pe_history.v1"
                put(g,(;event_id=id,gas_model=src.gas,pair_checksum=str(original,"pair_checksum"),
                    first_bucket=first(grid.buckets),last_bucket=last(grid.buckets),
                    ion_drift_cm=history.placement.ion_drift_cm,
                    electron_drift_cm=history.placement.electron_drift_cm,
                    incident_photons_total=sum(grid.counts),detected_pe_total=sum(detected),
                    dark_pulses_total=sum(dark),sensor_pe_total=sum(sensor_pe),
                    acquired_cells=rows,channel_ix=first.(grid.channels),
                    channel_iy=last.(grid.channels),channel_seed=seeds))
                vector(g,"pixel_ix",repeat(first.(grid.channels);inner=nb))
                vector(g,"pixel_iy",repeat(last.(grid.channels);inner=nb))
                vector(g,"bucket",repeat(grid.buckets,nc))
                vector(g,"detected_pe",detected)
                vector(g,"dark_pulses",dark)
                vector(g,"sensor_pe_equiv",sensor_pe)
            end
        end
    end
end

"""Apply ASIC noise/range to a saved sensor shard; retain PE per bucket."""
function asic_shard(s,branch,sensorpath,ids,number)
    path=joinpath(branch.asic_dir,"shard_$(lpad(number,4,'0')).h5")
    sensorhash=sha(sensorpath)
    shard_ready(path,ASIC_SCHEMA,branch.asicid,ids;upstream_hash=sensorhash) && return path
    pilot=s.pilot;q0=P.charge_per_pe_C(pilot.sensor)
    atomic_shard(path,ASIC_SCHEMA,branch.asicid,ids) do out
        attributes(out)["schema"]=ASIC_SCHEMA
        m=create_group(out,"metadata")
        put(m,(;complete=false,stage_id=branch.asicid,sensor_stage_id=branch.sensorid,
            source_path=sensorpath,source_sha256=sensorhash,sample=branch.src.sample,
            gas_model=branch.src.gas,event_ids=ids,event_count=length(ids),
            bucket_width_s=branch.src.dt,signal_unit="PE_equivalent",time_unit="s",
            charge_per_pe_C=q0,current_limit_A=pilot.electronics.current_limit_A,
            electronics_card_toml=read(pilot.ep,String),electronics_code_hash=pilot.electronics_code,
            production_code_hash=s.codehash,julia_version=string(VERSION)))
        events=create_group(out,"events")
        h5open(sensorpath,"r") do input
            for id in ids
                ig=input["events/$id"];sensor_pe=Float64.(read(ig["sensor_pe_equiv"]))
                ix=Int64.(read(ig["channel_ix"]));iy=Int64.(read(ig["channel_iy"]))
                nb=Int(read(ig["last_bucket"])-read(ig["first_bucket"])+1)
                rows=length(sensor_pe);require(rows==nb*length(ix),"Sensor matrix dimensions")
                measured=zeros(Float64,rows);saturated=falses(rows);seeds=zeros(UInt64,length(ix))
                clip_loss=0.0
                for (j,ch) in enumerate(zip(ix,iy))
                    seed=S.event_seed(branch.src.stream*pilot.electronics_model,id,
                                      "electronics:$(ch[1]):$(ch[2]):$(pilot.e["seed"])")
                    seeds[j]=seed;selection=(j-1)*nb+1:j*nb
                    trace=E.apply_electronics(Xoshiro(seed),view(sensor_pe,selection).*q0,
                                               q0,pilot.electronics)
                    measured[selection]=trace.current_A.*branch.src.dt./q0
                    saturated[selection]=trace.saturated
                    clip_loss+=sum(max.(trace.current_unclipped_A.-trace.current_A,0.))*branch.src.dt/q0
                end
                require(all(isfinite,measured) && all(measured.<=
                        pilot.electronics.current_limit_A*branch.src.dt/q0+1e-8),
                        "ASIC PE range")
                g=create_group(events,string(id));attributes(g)["schema"]="nausicaa.sipm_readout_pe_history.v1"
                put(g,(;event_id=id,gas_model=branch.src.gas,pair_checksum=str(ig,"pair_checksum"),
                    first_bucket=read(ig["first_bucket"]),last_bucket=read(ig["last_bucket"]),
                    ion_drift_cm=read(ig["ion_drift_cm"]),electron_drift_cm=read(ig["electron_drift_cm"]),
                    incident_photons_total=read(ig["incident_photons_total"]),
                    detected_pe_total=read(ig["detected_pe_total"]),
                    dark_pulses_total=read(ig["dark_pulses_total"]),
                    sensor_pe_total=read(ig["sensor_pe_total"]),
                    readout_pe_total=sum(measured),clip_loss_pe_equiv=clip_loss,
                    clipped_cells=count(identity,saturated),acquired_cells=rows,
                    channel_ix=ix,channel_iy=iy,channel_seed=seeds))
                vector(g,"pixel_ix",read(ig["pixel_ix"]))
                vector(g,"pixel_iy",read(ig["pixel_iy"]))
                vector(g,"bucket",read(ig["bucket"]))
                vector(g,"readout_pe_equiv",measured)
                vector(g,"saturated",Vector{Bool}(saturated))
            end
        end
    end
end

"""Run or resume complete shards and save a numerical event summary/manifest."""
function run(card::AbstractString)
    s=settings(card);branches=sources(s)
    runid=S.content_hash(S.canonical(s.run)*s.codehash*join(b.asicid for b in branches))[1:16]
    resultdir=resolve(joinpath(s.run["results_directory"],runid));mkpath(resultdir)
    entries=Dict{String,Any}[]
    for branch in branches
        chunks=Iterators.partition(branch.ids,Int(s.run["shard_events"]))
        for (number,chunk) in enumerate(chunks)
            ids=collect(chunk)
            sensor=sensor_shard(s,branch,ids,number)
            asic=asic_shard(s,branch,sensor,ids,number)
            push!(entries,Dict("sample"=>branch.src.sample,"gas"=>branch.src.gas,
                "event_ids"=>ids,"sensor"=>sensor,"sensor_sha256"=>sha(sensor),
                "asic"=>asic,"asic_sha256"=>sha(asic)))
            println("$(branch.src.sample)/$(branch.src.gas): $(last(ids)) / $(last(branch.ids))")
            flush(stdout)
        end
    end
    summary=joinpath(resultdir,"summary.h5");temp=summary*".tmp"
    h5open(temp,"w") do out
        attributes(out)["schema"]="nausicaa.sipm_production_summary.v1"
        m=create_group(out,"metadata")
        put(m,(;complete=false,run_card_toml=read(s.path,String),run_id=runid,
            production_code_hash=s.codehash,photon_index=resolve(s.run["photon_index"]),
            photon_index_sha256=sha(resolve(s.run["photon_index"]))))
        for branch in branches
            g=create_group(out,"$(branch.src.sample)/$(branch.src.gas)")
            ids=Int64[];incident=Int64[];detected=Int64[];dark=Int64[]
            sensor_pe=Float64[];readout_pe=Float64[];loss=Float64[]
            clipped=Int64[];acquired=Int64[]
            for entry in entries
                entry["gas"]==branch.src.gas && entry["sample"]==branch.src.sample || continue
                h5open(entry["asic"],"r") do h
                    for id in entry["event_ids"]
                        ev=h["events/$id"]
                        push!(ids,id);push!(incident,read(ev["incident_photons_total"]))
                        push!(detected,read(ev["detected_pe_total"]))
                        push!(dark,read(ev["dark_pulses_total"]))
                        push!(sensor_pe,read(ev["sensor_pe_total"]))
                        push!(readout_pe,read(ev["readout_pe_total"]))
                        push!(loss,read(ev["clip_loss_pe_equiv"]))
                        push!(clipped,read(ev["clipped_cells"]))
                        push!(acquired,read(ev["acquired_cells"]))
                    end
                end
            end
            require(ids==branch.ids,"Aggregate event ledger")
            if s.run["event_selection"]=="all_sorted_source_ids"
                expected=h5open(resolve(s.run["photon_index"]),"r") do h
                    Int64(read(h["$(branch.src.sample)/$(branch.src.gas)_incident"]))
                end
                require(sum(incident)==expected,"Aggregate incident-photon ledger")
            end
            put(g,(;event_id=ids,incident_photons=incident,detected_pe=detected,
                dark_pulses=dark,sensor_pe_equiv=sensor_pe,readout_pe_equiv=readout_pe,
                clip_loss_pe_equiv=loss,clipped_cells=clipped,acquired_cells=acquired,
                source_sha256=branch.src.sourcehash,sensor_stage_id=branch.sensorid,
                asic_stage_id=branch.asicid,event_count=length(ids)))
        end
        write(m["complete"],true)
    end
    mv(temp,summary;force=true)
    manifest=joinpath(resultdir,"manifest.toml")
    open(manifest*".tmp","w") do io
        TOML.print(io,Dict("schema"=>"nausicaa.sipm_production_manifest.v1",
            "run_id"=>runid,"complete"=>true,"summary"=>summary,"summary_sha256"=>sha(summary),
            "run_card"=>s.path,"shards"=>entries))
    end
    mv(manifest*".tmp",manifest;force=true)
    println("SiPM production manifest: $manifest")
    manifest
end
end
