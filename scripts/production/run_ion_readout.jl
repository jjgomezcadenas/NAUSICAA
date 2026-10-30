#!/usr/bin/env julia
# Replay a saved pair realisation into a new ion grid without running electrons.
include(joinpath(@__DIR__,"run_detector_pilot.jl"))

"""Hash the Julia sources that determine this ion-only replay and its output."""
function ion_replay_code_hash()
    files=("src/event_analysis/EventAnalysis.jl","src/event_analysis/types.jl",
           "src/event_analysis/io.jl","src/event_analysis/production.jl",
           "src/event_analysis/production_readout.jl",
           "src/event_analysis/grid_readout.jl","src/ion_signal/event_io.jl",
           "scripts/production/specs.jl","scripts/production/run_detector_pilot.jl",
           "scripts/production/run_ion_readout.jl")
    return content_hash(join((content_hash(joinpath(ROOT,file);file=true) for file in files),"|"))
end

"""Resolve the saved common realisation and the ion output for one sample/card."""
function ion_paths(spec,sample::String,code_hash::String)
    r=spec.realisation.data;i=spec.ions.data
    source=spec.sample_paths[sample]
    source_sha=content_hash(source;file=true)
    realisation_values=Dict(k=>v for (k,v) in r if k!="samples")
    stream_values=Dict(k=>v for (k,v) in realisation_values if k!="max_events_per_sample")
    stream_identity=content_hash(join((sample,source_sha,canonical(stream_values)),"|"))
    realisation_id=content_hash(join((stream_identity,canonical(realisation_values)),"|"))
    ion_id=content_hash(join((realisation_id,canonical(i),code_hash),"|"))
    base=joinpath(spec.results_dir,sample,realisation_id)
    return (;source,source_sha,stream_identity,realisation_id,ion_id,
            realisation=joinpath(base,"realisation.h5"),
            ions=joinpath(base,"detector","ions",ion_id,"ions.h5"))
end

"""Verify that a completed common realisation matches this source and run card."""
function check_realisation(file::HDF5.File,paths,n::Int,seed::Int)
    String(read(attributes(file)["schema"]))=="nausicaa.event_realisation.pilot.v1" ||
        error("Realisation schema mismatch")
    m=file["metadata"]
    Bool(read(m["complete"])) && Int(read(m["event_count"]))==n ||
        error("Realisation is incomplete or has the wrong event count")
    String(read(m["source_sha256"]))==paths.source_sha || error("Source content mismatch")
    String(read(m["card_hash"]))==paths.realisation_id || error("Realisation card mismatch")
    String(read(m["stream_identity"]))==paths.stream_identity || error("RNG stream mismatch")
    Int(read(m["seed"]))==seed || error("Realisation seed mismatch")
    return nothing
end

"""Replay initial pair RNG draws, verify them against the common realisation,
and collect ions into one new pad/time grid for every event. The ion diffusion
seed remains independent of grid pitch and bucket size. The two electron-gas
outputs and the common realisation are read-only throughout this operation.
"""
function run_ion_sample(spec,sample::String,code_hash::String)
    r=spec.realisation.data;i=spec.ions.data
    paths=ion_paths(spec,sample,code_hash)
    isfile(paths.realisation) || error("Run the common realisation/electron branch first: $(paths.realisation)")
    n=Int(r["max_events_per_sample"]);seed=Int(r["seed"])
    if isfile(paths.ions)
        h5open(paths.ions,"r") do existing
            String(read(attributes(existing)["schema"]))=="nausicaa.ion_pad_file.pilot.v1" ||
                error("Existing ion file has another schema")
            m=existing["metadata"]
            Bool(read(m["complete"])) && Int(read(m["event_count"]))==n &&
                String(read(m["card_hash"]))==paths.ion_id ||
                error("Existing ion file does not match this completed run")
        end
        println("$sample: reused completed ion readout $(paths.ions)")
        return paths.ions
    end
    mkpath(dirname(paths.ions))
    temporary=paths.ions*".tmp"
    isfile(temporary) && rm(temporary)
    total_pairs=Int64(0);total_ions=Int64(0)
    started=time_ns()
    h5open(paths.source,"r") do source
        ids=read(source["events/event_id"])
        length(ids)>=n || error("Source has fewer than $n events")
        h5open(paths.realisation,"r") do common
            check_realisation(common,paths,n,seed)
            h5open(temporary,"w") do output
                write_header(output,"nausicaa.ion_pad_file.pilot.v1",sample,
                    paths.source,paths.source_sha,paths.ion_id,code_hash,seed,
                    "ions",read(spec.realisation.path,String),read(spec.ions.path,String))
                output["metadata/stream_identity"]=paths.stream_identity
                for (index,raw_id) in enumerate(ids[1:n])
                    id=Int64(raw_id);key=string(id)
                    hits,_=read_event(source,id)
                    saved=common["events/$key"]
                    pairs=sample_fano_pairs(Xoshiro(event_seed(paths.stream_identity,id,"pairs")),
                        hits,Float32(r["w_eV"]),Float32(r["fano_factor"]))
                    count=foldl(Base.Checked.checked_add,pairs;init=Int64(0))
                    checksum=pair_checksum(pairs)
                    count==Int64(read(saved["initial_pairs"])) &&
                        checksum==String(read(saved["pair_checksum"])) ||
                        error("Pair replay differs from the saved realisation at event $id")
                    placement=EventPlacement(id,Float32(read(saved["ion_drift_cm"])),
                        Float32(read(saved["electron_drift_cm"])))
                    rx,ry=transverse_reference(hits)
                    phase=sample_grid_offset_mm(
                        Xoshiro(event_seed(paths.stream_identity,id,"ion_phase")),
                        Float32(i["pad_pitch_mm"]))
                    bucket_mm=Float32(i["drift_speed_m_s"]*i["bucket_width_ms"])
                    grid=GridSpec(Float32(i["pad_pitch_mm"]),bucket_mm,
                        (rx,ry,minimum(hits.z_mm)),phase,Int8(1))
                    sigma=sigma_t_ion_mm(Float32(i["temperature_K"]),
                        placement.ion_drift_cm,Float32(i["drift_field_V_cm"]))
                    history=collect_ion_pad_history(
                        Xoshiro(event_seed(paths.stream_identity,id,"ion_diffusion")),
                        hits,pairs,placement,grid;sigma_t_mm=sigma,sigma_l_mm=sigma)
                    group=create_group(output["events"],key)
                    write_ion_history(group,history)
                    group["pair_checksum"]=checksum
                    total_pairs+=count
                    total_ions+=sum(history.cells.count)
                    if index%100==0 || index==n
                        println("$sample ion $(i["pad_pitch_mm"]) mm: $index/$n events")
                        flush(stdout)
                    end
                end
                total_ions==total_pairs || error("Ion charge ledger failed")
                output["metadata/total_pairs"]=total_pairs
                output["metadata/total_ions"]=total_ions
                output["metadata/event_count"]=n
                output["metadata/elapsed_s"]=(time_ns()-started)/1e9
                write(output["metadata/complete"],true)
            end
        end
    end
    mv(temporary,paths.ions;force=true)
    elapsed=(time_ns()-started)/1e9
    println("$sample ion $(i["pad_pitch_mm"]) mm: $n events in $(round(elapsed;digits=2)) s; "*
            "pairs=$total_pairs, ions=$total_ions; $(paths.ions)")
    return paths.ions
end

function main_ions(card::String)
    spec=load_specs(ROOT,card)
    hash=ion_replay_code_hash()
    for sample in get(spec.run_spec,"selected_samples",["bb0nu","single_electrons"])
        run_ion_sample(spec,sample,hash)
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("Usage: julia --project=. scripts/production/run_ion_readout.jl RUN.toml")
    main_ions(ARGS[1])
end
