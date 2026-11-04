"""Card resolution, photon-history adapter and independently reusable response files."""
module SiPMPipeline
using HDF5, TOML, SHA, Random, Statistics
include("specs.jl")
include("../../src/event_analysis/EventAnalysis.jl")
include("../../src/ion_signal/sipm_response.jl")
include("../../src/ion_signal/sipm_electronics.jl")
const S=ProductionSpecs
const P=SiPMResponse
const E=SiPMElectronics
const EA=EventAnalysis
sha(path)=open(io->bytes2hex(sha256(io)),path)
str(g,k)=String(read(g[k]))
schema(g)=String(read(attributes(g)["schema"]))

function settings(root,runpath)
    runfile=joinpath(root,runpath); run=TOML.parsefile(runfile)
    run["schema"]=="nausicaa.sipm_run.v1" && run["pilot_only"] || error("Pilot run card required")
    ids=Int64.(run["event_ids"])
    0<length(ids)<=16 && allunique(ids) || error("Pilot requires 1–16 unique event IDs")
    sort!(ids)
    cp=joinpath(root,run["sipm_card"]); ep=joinpath(root,run["electronics_card"])
    c=TOML.parsefile(cp); e=TOML.parsefile(ep)
    c["schema"]=="nausicaa.sipm_spec.v1" && e["schema"]=="nausicaa.sipm_electronics_spec.v1" || error("Response card schema")
    for (key,value) in (("arrival_model","uniform_in_bucket"),("illumination","uniform_microcells"),
        ("recovery_pde","constant"),("initial_state","stationary_dark"),("charge_assignment","avalanche_bucket"))
        c[key]==value || error("Unsupported sensor policy: $key")
    end
    !any(c[k] for k in ("crosstalk","afterpulsing","gain_spread")) || error("Additional sensor effects require a model")
    e["transfer"]=="bucket_average_current" && e["noise"]=="independent_gaussian_pe" &&
        e["baseline_pe"]==0 && !e["lower_clipping"] || error("Unsupported ASIC policy")
    all(x->isfinite(x)&&x>0,(c["diameter_m"],c["microcell_pitch_m"],c["dark_reference_area_m2"])) || error("Sensor geometry")
    c["seed"]>=0 && e["seed"]>=0 || error("Seeds must be nonnegative")
    a=c["acquisition"]
    a["channels"]=="illuminated_in_event" && a["quiet_before_buckets"]>=0 &&
        a["quiet_after_buckets"]>=0 && a["empty_last_bucket"]>=a["empty_first_bucket"] || error("Acquisition policy")
    all(x->length(x)==2 && all(v->v isa Integer,x),a["additional_channels"]) || error("Additional channel indices")
    sensor=P.Sensor(c["pde"],c["microcell_count"],c["gain"],c["recovery_s"],
        P.dark_rate_Hz(c["dark_reference_rate_Hz"],c["dark_reference_area_m2"],c["diameter_m"]))
    electronics=E.Electronics(e["bucket_width_s"],e["noise_rms_pe"],e["current_limit_A"])
    common=[@__FILE__,joinpath(@__DIR__,"specs.jl"),joinpath(root,"src/event_analysis/production_readout.jl"),
        joinpath(root,"src/event_analysis/EventAnalysis.jl"),joinpath(root,"Manifest.toml")]
    sensor_files=vcat(common,[joinpath(root,"src/ion_signal/sipm_response.jl")])
    electronics_files=vcat(common,[joinpath(root,"src/ion_signal/sipm_electronics.jl")])
    sensor_code=S.content_hash(join(sha.(sensor_files),"|")*string(VERSION))
    electronics_code=S.content_hash(join(sha.(electronics_files),"|")*string(VERSION))
    sensor_model=S.content_hash(S.canonical(c)*sensor_code)
    electronics_model=S.content_hash(S.canonical(e)*electronics_code)
    (;root,run,runfile,ids,cp,ep,c,e,sensor,electronics,sensor_code,electronics_code,sensor_model,electronics_model,
      sensor_files,electronics_files)
end

"""Owned dense counts on the selected channels and full signed acquisition interval."""
function acquisition(history,a)
    c=history.cells
    channels=sort!(unique(vcat(collect(zip(c.pixel_ix,c.pixel_iy)),
        Tuple{Int64,Int64}[(Int64(x[1]),Int64(x[2])) for x in a["additional_channels"]])))
    lo=isempty(c.bucket) ? Int64(a["empty_first_bucket"]) : minimum(c.bucket)-a["quiet_before_buckets"]
    hi=isempty(c.bucket) ? Int64(a["empty_last_bucket"]) : maximum(c.bucket)+a["quiet_after_buckets"]
    counts=zeros(Int64,hi-lo+1,length(channels))
    positions=Dict(ch=>j for (j,ch) in enumerate(channels))
    for k in eachindex(c.count)
        counts[c.bucket[k]-lo+1,positions[(c.pixel_ix[k],c.pixel_iy[k])]]=c.count[k]
    end
    sum(counts)==sum(c.count) || error("Acquisition lost photons")
    (;channels,counts,buckets=collect(Int64,lo:hi))
end

function sources(s)
    out=NamedTuple[]
    h5open(joinpath(s.root,s.run["photon_index"]),"r") do idx
        for sample in s.run["samples"],gas in s.run["gases"]
            path=str(idx,"$sample/$(gas)_path"); sourcehash=sha(path)
            h5open(path,"r") do h
                schema(h)=="nausicaa.photon_sipm_file.pilot.v1" || error("Photon file schema")
                m=h["metadata"]
                Bool(read(m["complete"])) && str(m,"sample")==sample && str(m,"branch")==gas || error("Photon input identity")
                str(m,"signal_unit")=="incident_photons" && str(m,"time_unit")=="us" || error("Photon input units")
                card=TOML.parse(str(m,"branch_card_toml")); dt=card["bucket_width_us"]*1e-6
                isapprox(dt,s.electronics.bucket_width_s;rtol=1e-12) || error("Photon and ASIC buckets differ")
                isapprox(card["sipm_diameter_mm"]*1e-3,s.c["diameter_m"];rtol=1e-12) || error("Sensor must match photon collection aperture")
                n=read(m["event_count"])
                n==length(keys(h["events"])) && all(id->haskey(h["events"],string(id)),s.ids) || error("Input event count/selection")
                stream=S.content_hash(s.sensor_model*sourcehash) # independent of event selection
                sensorid=S.content_hash(stream*S.canonical(s.ids))
                eid=S.content_hash(sensorid*s.electronics_model)
                base=dirname(dirname(dirname(dirname(path))))
                parent=joinpath(base,"response","electrons",basename(dirname(path)))
                sensorpath=joinpath(parent,sensorid,"sipm_response.h5")
                response=joinpath(parent,eid,"electron_response.h5")
                push!(out,(;sample,gas,path,sourcehash,dt,card,n,stream,sensorid,eid,sensorpath,response))
            end
        end
    end
    out
end

function put(g,values)
    for (key,value) in pairs(values); g[string(key)]=value; end
end

function atomic_h5(writefile,path)
    mkpath(dirname(path)); tmp=path*".tmp"
    h5open(tmp,"w") do h
        writefile(h)
        h["metadata/complete"]=true
    end
    h5open(tmp,"r") do h
        Bool(read(h["metadata/complete"])) || error("Incomplete response")
        length(keys(h["events"]))==read(h["metadata/event_count"]) || error("Response event count")
    end
    mv(tmp,path;force=true)
    path
end

function valid_cache(path,id,expected_schema)
    isfile(path) || return false
    h5open(path,"r") do h
        schema(h)==expected_schema && Bool(read(h["metadata/complete"])) &&
            str(h["metadata"],"response_id")==id || error("Incompatible cached response: $path")
    end
    true
end

function sensor_file(s,src)
    valid_cache(src.sensorpath,src.sensorid,"nausicaa.sipm_file.v1") && return src.sensorpath
    atomic_h5(src.sensorpath) do out
        attributes(out)["schema"]="nausicaa.sipm_file.v1"
        m=create_group(out,"metadata")
        put(m,(;response_id=src.sensorid,stream_identity=src.stream,source_path=src.path,source_sha256=src.sourcehash,
            sample=src.sample,gas_model=src.gas,event_count=length(s.ids),event_ids=s.ids,
            card_toml=read(s.cp,String),card_hash=S.content_hash(S.canonical(s.c)),code_hash=s.sensor_code,
            julia_version=string(VERSION),bucket_width_s=src.dt,charge_per_pe_C=P.charge_per_pe_C(s.sensor),
            dark_rate_Hz=s.sensor.dark_rate_Hz,signal_unit="C",time_unit="s"))
        provenance=create_group(m,"code")
        for (j,path) in enumerate(s.sensor_files)
            g=create_group(provenance,string(j));put(g,(;path=relpath(path,s.root),sha256=sha(path),contents=read(path,String)))
        end
        events=create_group(out,"events")
        h5open(src.path,"r") do input
            up=create_group(m,"photon_metadata")
            for k in keys(input["metadata"]); up[k]=read(input["metadata"][k]); end
            for id in s.ids
                original=input["events/$id"]; history=EA.read_photon_history(original)
                history.placement.event_id==id && string(history.gas_model)==src.gas || error("Event identity mismatch")
                isapprox(Float64(history.grid.bucket_width_mm)/src.card["drift_speed_mm_us"]*1e-6,src.dt;rtol=1e-6) || error("Event bucket width")
                grid=acquisition(history,s.c["acquisition"]); nb=length(grid.buckets); nc=length(grid.channels)
                data=Dict(k=>zeros(T,nb*nc) for (k,T) in (("detected_pe",Int64),("dark_pulses",Int64),
                    ("optical_charge_C",Float64),("dark_charge_C",Float64),("charge_C",Float64)))
                seeds=UInt64[]; elapsed=@elapsed for (j,ch) in enumerate(grid.channels)
                    seed=S.event_seed(src.stream,id,"sensor:$(ch[1]):$(ch[2]):$(s.c["seed"])");push!(seeds,seed)
                    trace=P.simulate_channel(Xoshiro(seed),view(grid.counts,:,j),src.dt,s.sensor)
                    rows=(j-1)*nb+1:j*nb
                    for k in keys(data); data[k][rows]=getproperty(trace,Symbol(k)); end
                end
                all(data["detected_pe"].<=vec(grid.counts)) || error("PDE count bound")
                all(data["charge_C"].>=0) && all(data["charge_C"].<=
                    P.charge_per_pe_C(s.sensor).*(data["detected_pe"]+data["dark_pulses"]).*(1+1e-12)) || error("Charge bound")
                g=create_group(events,string(id));attributes(g)["schema"]="nausicaa.sipm_history.v1"
                EA.write_photon_history(create_group(g,"photon_input"),history)
                put(g,(;event_id=id,gas_model=src.gas,pair_checksum=str(original,"pair_checksum"),
                    pixel_ix=repeat(first.(grid.channels);inner=nb),pixel_iy=repeat(last.(grid.channels);inner=nb),
                    bucket=repeat(grid.buckets,nc),incident_photons=vec(grid.counts),
                    channel_ix=first.(grid.channels),channel_iy=last.(grid.channels),channel_seed=seeds,
                    first_bucket=first(grid.buckets),last_bucket=last(grid.buckets),simulation_elapsed_s=elapsed))
                put(g,data)
                println("SiPM $(src.sample)/$(src.gas) event $id: $(sum(data["detected_pe"])) PE, $(round(elapsed;digits=2)) s")
            end
        end
    end
end

function electronics_file(s,src)
    if valid_cache(src.response,src.eid,"nausicaa.electron_response_file.v1")
        matches=h5open(src.response,"r") do h
            str(h["metadata"],"source_sha256")==sha(src.sensorpath)
        end
        matches && return src.response
        # A reconstructed sensor cache can have different timing metadata. Rebuild
        # electronics from the current saved sensor contents and update provenance.
    end
    atomic_h5(src.response) do out
        attributes(out)["schema"]="nausicaa.electron_response_file.v1"
        m=create_group(out,"metadata"); sensorhash=sha(src.sensorpath)
        put(m,(;response_id=src.eid,sensor_response_id=src.sensorid,source_path=src.sensorpath,
            source_sha256=sensorhash,sample=src.sample,gas_model=src.gas,event_count=length(s.ids),event_ids=s.ids,
            card_toml=read(s.ep,String),card_hash=S.content_hash(S.canonical(s.e)),code_hash=s.electronics_code,
            julia_version=string(VERSION),bucket_width_s=src.dt,charge_per_pe_C=P.charge_per_pe_C(s.sensor),
            signal_unit="A",time_unit="s"))
        code=create_group(m,"code")
        for (j,path) in enumerate(s.electronics_files)
            g=create_group(code,string(j));put(g,(;path=relpath(path,s.root),sha256=sha(path),contents=read(path,String)))
        end
        events=create_group(out,"events")
        h5open(src.sensorpath,"r") do input
            for id in s.ids
                ig=input["events/$id"]; q=read(ig["charge_C"])
                nc=length(read(ig["channel_seed"])); nb=read(ig["last_bucket"])-read(ig["first_bucket"])+1
                noise=zeros(length(q)); raw=similar(noise); current=similar(noise); clipped=falses(length(q));seeds=UInt64[]
                for (j,ch) in enumerate(zip(read(ig["channel_ix"]),read(ig["channel_iy"])))
                    seed=S.event_seed(src.stream*s.electronics_model,id,"electronics:$(ch[1]):$(ch[2]):$(s.e["seed"])")
                    push!(seeds,seed); rows=(j-1)*nb+1:j*nb
                    trace=E.apply_electronics(Xoshiro(seed),view(q,rows),P.charge_per_pe_C(s.sensor),s.electronics)
                    noise[rows]=trace.noise_pe;raw[rows]=trace.current_unclipped_A
                    current[rows]=trace.current_A;clipped[rows]=trace.saturated
                end
                all(current .<= s.electronics.current_limit_A) || error("ASIC range")
                g=create_group(events,string(id));attributes(g)["schema"]="nausicaa.electron_response_history.v1"
                for key in ("event_id","gas_model","pair_checksum","pixel_ix","pixel_iy","bucket","incident_photons",
                            "detected_pe","dark_pulses","optical_charge_C","dark_charge_C","charge_C",
                            "channel_ix","channel_iy","first_bucket","last_bucket")
                    g[key]=read(ig[key])
                end
                put(g,(;channel_seed=seeds,noise_pe=noise,current_unclipped_A=raw,current_A=current,saturated=Vector{Bool}(clipped)))
                EA.write_photon_history(create_group(g,"photon_input"),EA.read_photon_history(ig["photon_input"]))
            end
        end
    end
end

"""Write numerical pilot summaries and physical parameters for Python presentation."""
function run(root,runpath)
    s=settings(root,runpath); srcs=sources(s)
    study=S.content_hash(join(x.eid for x in srcs)*S.canonical(s.run))[1:16]
    dir=joinpath(root,s.run["results_directory"],study);mkpath(dir)
    write(joinpath(dir,".gitignore"),"*.h5\n")
    for src in srcs; sensor_file(s,src);electronics_file(s,src); end
    path=joinpath(dir,"summary.h5")
    h5open(path*".tmp","w") do out
        attributes(out)["schema"]="nausicaa.sipm_pilot_summary.v1"
        p=create_group(out,"parameters")
        put(p,(;pde=s.sensor.pde,microcells=s.sensor.microcells,gain=s.sensor.gain,recovery_s=s.sensor.recovery_s,
            diameter_m=s.c["diameter_m"],dark_rate_Hz=s.sensor.dark_rate_Hz,charge_per_pe_C=P.charge_per_pe_C(s.sensor),
            bucket_width_s=s.electronics.bucket_width_s,noise_rms_pe=s.electronics.noise_rms_pe,
            current_limit_A=s.electronics.current_limit_A,
            current_limit_pe=s.electronics.current_limit_A*s.electronics.bucket_width_s/P.charge_per_pe_C(s.sensor),
            dark_mean_per_bucket=s.sensor.dark_rate_Hz*s.electronics.bucket_width_s,
            noise_rms_A=s.electronics.noise_rms_pe*P.charge_per_pe_C(s.sensor)/s.electronics.bucket_width_s))
        for src in srcs
            g=create_group(out,"$(src.sample)/$(src.gas)")
            put(g,(;sensor_path=src.sensorpath,electronics_path=src.response,source_path=src.path,source_sha256=src.sourcehash))
            h5open(src.response,"r") do h
                for id in s.ids
                    ev=h["events/$id"]; q=read(ev["charge_C"]);pe=read(ev["detected_pe"])
                    dark=read(ev["dark_pulses"]);clip=read(ev["saturated"])
                    eg=create_group(g,string(id))
                    put(eg,(;event_id=id,incident_photons=sum(read(ev["incident_photons"])),detected_pe=sum(pe),
                        dark_pulses=sum(dark),sensor_charge_pe=sum(q)/P.charge_per_pe_C(s.sensor),
                        recovery_loss_fraction=1-sum(q)/(P.charge_per_pe_C(s.sensor)*sum(pe+dark)),
                        max_sensor_pe=maximum(q;init=0.)/P.charge_per_pe_C(s.sensor),
                        saturated_buckets=count(identity,clip),acquired_channel_buckets=length(q),
                        saturation_fraction=isempty(q) ? 0. : mean(clip),
                        noise_mean_pe=isempty(q) ? 0. : mean(read(ev["noise_pe"])),
                        noise_rms_pe=isempty(q) ? 0. : sqrt(mean(abs2,read(ev["noise_pe"])))))
                end
            end
        end
    end
    mv(path*".tmp",path;force=true)
    open(joinpath(dir,"manifest.toml"),"w") do io
        TOML.print(io,Dict("study_id"=>study,"run_card_toml"=>read(s.runfile,String),"summary"=>path,
            "sensor_model_id"=>s.sensor_model,"electronics_model_id"=>s.electronics_model,
            "outputs"=>[Dict("sample"=>x.sample,"gas"=>x.gas,"sensor"=>x.sensorpath,"electronics"=>x.response) for x in srcs]))
    end
    open(joinpath(root,s.run["results_directory"],"latest.toml"),"w") do io
        TOML.print(io,Dict("summary"=>path,"manifest"=>joinpath(dir,"manifest.toml")))
    end
    println("Pilot summary: $path")
    path
end
end
