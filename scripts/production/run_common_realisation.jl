#!/usr/bin/env julia
# Materialize shared event placement and primary-pair ledgers before either arm.
include("run_detector_pilot.jl")

"""Write immutable shared placement and primary ledgers for one sample.

Ion reference drift is uniform within the realization card's min/max bounds;
electron drift is its chamber-length complement. With a primary reference,
each replayed pair count/checksum is verified against the original realization.
No carrier positions are stored. Both collection arms consume this same file.
"""
function run_common(spec,sample)
    r=spec.realisation.data;source=spec.sample_paths[sample]
    source_sha=content_hash(source;file=true)
    identity=realisation_identities(ROOT,sample,source_sha,r)
    path=joinpath(spec.results_dir,sample,identity.realisation_id,"realisation.h5")
    n=Int(r["max_events_per_sample"])
    code=identity.primary_code_hash
    if isfile(path)
        h5open(path,"r") do f
            check_primary_code(ROOT,f["metadata"])
            @assert read(f["metadata/complete"]) && read(f["metadata/event_count"])==n
            @assert read(f["metadata/card_hash"])==identity.realisation_id
            @assert read(f["metadata/source_sha256"])==source_sha && read(f["metadata/code_hash"])==code
        end
        return path
    end
    reference=nothing
    if haskey(r,"primary_reference")
        ref=TOML.parsefile(joinpath(ROOT,r["primary_reference"]))
        refids=realisation_identities(ROOT,sample,source_sha,ref)
        reference=joinpath(spec.results_dir,sample,refids.realisation_id,"realisation.h5")
        if !isfile(reference)
            reference=joinpath(spec.results_dir,sample,refids.legacy_realisation_id,"realisation.h5")
        end
        isfile(reference) || error("Primary reference ledger missing")
    end
    # A same-card legacy ledger is replay-verified, never overwritten.
    legacy=joinpath(spec.results_dir,sample,identity.legacy_realisation_id,"realisation.h5")
    previous=isfile(legacy) ? h5open(legacy,"r") : nothing
    if previous!==nothing
        m=previous["metadata"]
        @assert read(m["complete"]) && read(m["event_count"])==n
        @assert read(m["source_sha256"])==source_sha && read(m["card_hash"])==identity.legacy_realisation_id
        @assert read(m["stream_identity"])==identity.stream_identity
    end
    mkpath(dirname(path));started=time()
    h5open(source,"r") do src
      ref=isnothing(reference) ? nothing : h5open(reference,"r")
      try
        ids=Int64.(read(src["events/event_id"]));length(ids)>=n || error("Insufficient events")
        h5open(path*".tmp","w") do f
            write_header(f,"nausicaa.event_realisation.v2",sample,source,source_sha,
                identity.realisation_id,code,Int(r["seed"]),"realisation",read(spec.realisation.path,String),"")
            m=f["metadata"];m["primary_code_hash"]=code;m["stream_identity"]=identity.stream_identity
            m["placement_mode"]="uniform_ion_reference"
            m["minimum_ion_drift_cm"]=r["min_ion_drift_cm"]
            m["maximum_ion_drift_cm"]=get(r,"max_ion_drift_cm",r["chamber_length_cm"])
            m["chamber_length_cm"]=r["chamber_length_cm"]
            if previous!==nothing
                m["verified_legacy_path"]=legacy
                m["verified_legacy_sha256"]=content_hash(legacy;file=true)
            end
            if reference!==nothing
                m["primary_reference_path"]=reference
                m["primary_reference_sha256"]=content_hash(reference;file=true)
            end
            for (j,id) in enumerate(ids[1:n])
                hits,_=read_event(src,id)
                pairs=sample_fano_pairs(Xoshiro(event_seed(identity.stream_identity,id,"pairs")),hits,
                    Float32(r["w_eV"]),Float32(r["fano_factor"]))
                count=foldl(Base.Checked.checked_add,pairs;init=Int64(0));checksum=pair_checksum(pairs)
                if ref!==nothing
                    @assert read(ref["events/$id/initial_pairs"])==count
                    @assert read(ref["events/$id/pair_checksum"])==checksum
                end
                li=sample_drift_cm(Xoshiro(event_seed(identity.stream_identity,id,"placement")),
                    Float32(r["min_ion_drift_cm"]),Float32(get(r,"max_ion_drift_cm",r["chamber_length_cm"])))
                le=Float32(r["chamber_length_cm"])-li
                if previous!==nothing
                    old=previous["events/$id"]
                    @assert read(old["initial_pairs"])==count && read(old["pair_checksum"])==checksum
                    @assert read(old["ion_drift_cm"])==li && read(old["electron_drift_cm"])==le
                end
                g=create_group(f["events"],string(id));g["event_id"]=id
                g["ion_drift_cm"]=li;g["electron_drift_cm"]=le
                g["initial_pairs"]=count;g["pair_checksum"]=checksum
                if j%1000==0 || j==n;println("Common realization: $j/$n");flush(stdout);end
            end
            m["event_count"]=n;m["elapsed_s"]=time()-started;write(m["complete"],true)
        end
      finally
        ref===nothing || close(ref)
        previous===nothing || close(previous)
      end
    end
    mv(path*".tmp",path);println("Common realization: $path");path
end
if abspath(PROGRAM_FILE)==@__FILE__
    spec=load_specs(ROOT,only(ARGS))
    for sample in get(spec.run_spec,"selected_samples",["bb0nu"]);run_common(spec,sample);end
end
