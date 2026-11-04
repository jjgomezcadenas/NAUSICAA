"""Production orchestration and owned readers for reconstructed ion-pad histories."""
module IonCSAProduction
using HDF5, TOML, SHA, LinearAlgebra
include("../../scripts/production/specs.jl")
const S = ProductionSpecs
export ReconstructedIonEvent, read_event, read_covariance, validate_response, settings, run
const CHANNEL_POLICY = "pixels_with_positive_input_charge"
const ABSENT_POLICY = "no_input_charge_noise_not_simulated"
const PRODUCTION_SCHEMA = "nausicaa.ion_csa_response.v1"
sha(path) = open(io -> bytes2hex(sha256(io)), path)
str(g,k) = String(read(g[k]))
schema(h) = String(read(attributes(h)["schema"]))
require(ok,message) = ok || throw(ArgumentError(message))
resolve(root,path) = normpath(isabspath(path) ? path : joinpath(root,path))

"""Owned frame-by-pad matrices in ion-charge units, with explicit pixel and time axes.

Each pixel received input charge somewhere in this event. Every saved bucket,
including quiet buckets, is retained. Time edges are seconds. Coordinates and
counts preserve their exact signed Int64 keys and integer counts; signal and
noise matrices preserve the validated Float32 storage values.
"""
struct ReconstructedIonEvent
    event_id::Int64
    pixel_ix::Vector{Int64}
    pixel_iy::Vector{Int64}
    first_bucket::Int64
    edge_time_s::Vector{Float64}
    frame_s::Float64
    pitch_mm::Float64
    reference_mm::Vector{Float32}
    phase_mm::Vector{Float32}
    bucket_width_mm::Float32
    drift_direction::Int8
    ion_drift_cm::Float32
    electron_drift_cm::Float32
    input_ions::Matrix{Int64}
    noiseless_ions::Matrix{Float32}
    noise_ions::Matrix{Float32}
    reconstructed_ions::Matrix{Float32}
    noise_sigma_ions::Vector{Float64}
    model_used::Vector{Int8}
    noise_seed::Vector{UInt64}
    pair_checksum::String
end

function check_header(h)
    require(schema(h) in ("nausicaa.ion_csa_response.pilot.v2",PRODUCTION_SCHEMA),
            "Requires a corrected ion CSA response")
    m=h["metadata"]
    require(Bool(read(m["complete"])),"Incomplete ion CSA response")
    require(str(m,"signal_unit")=="ion_charges" && str(m,"time_unit")=="s",
            "Response units must be ion charges and seconds")
    require(str(m,"matrix_axes_julia")=="frame,pad" && str(m,"matrix_axes_hdf5")=="pad,frame",
            "Response matrix axes mismatch")
    require(str(m,"reconstruction_method")=="three_frame_inverse","Unsupported reconstruction")
    if schema(h)==PRODUCTION_SCHEMA
        require(str(m,"channel_selection")==CHANNEL_POLICY && str(m,"absent_channel_policy")==ABSENT_POLICY,
                "Unsupported channel selection")
        require(!Bool(read(m["measured_zero_suppression"])) && !Bool(read(m["negative_estimate_clipping"])),
                "Requires signed, unsuppressed estimates on input-charged pixels")
    end
    m
end

function read_event(h::HDF5.File,id::Integer)
    m=check_header(h);require(haskey(h,"events/$id"),"Event ID is absent: $id")
    g=h["events/$id"]
    require(Int64(read(g["event_id"]))==id,"Event ID/key mismatch")
    px=Int64.(read(g["pixel_ix"]));py=Int64.(read(g["pixel_iy"]))
    n=Int64.(read(g["counts"]));q=Float32.(read(g["noiseless_estimated_ions"]))
    z=Float32.(read(g["noise_ions"]));r=Float32.(read(g["recovered_ions"]))
    sigma=Float64.(read(g["noise_sigma_ions"]));edges=Float64.(read(g["edge_time_s"]))
    first=Int64(read(g["first_bucket"]));dt=Float64(read(m["frame_s"]))
    np=length(px);nf=size(n,1);keys=collect(zip(px,py))
    require(length(py)==np && np>0 && nf>0 && size(n,2)==np,"Pad/count dimensions mismatch")
    require(issorted(keys) && allunique(keys),"Pixel keys must be sorted and unique")
    require(all(size(x)==size(n) for x in (q,z,r)),"Signal/noise matrix dimensions mismatch")
    require(all(>=(0),n) && all(>(0),vec(sum(n;dims=1))),"Rows must select pixels with input charge")
    require(sum(n)==Int64(read(g["initial_pairs"])),"Input ion ledger mismatch")
    require(length(sigma)==nf && all(x->isfinite(x)&&x>0,sigma),"Frame noise sigma mismatch")
    require(length(edges)==nf+1 && all(isfinite,edges) && isfinite(dt) && dt>0,"Time axis dimensions/units")
    require(all(isapprox.(diff(edges),dt;rtol=1e-10,atol=1e-12)) &&
            isapprox(edges[1],first*dt;rtol=1e-10,atol=1e-12),"Bucket keys and time edges disagree")
    require(all(x->all(isfinite,x),(q,z,r)),"Nonfinite reconstructed signal")
    require(maximum(abs.(Float64.(r).-Float64.(q).-Float64.(z))./sigma)<0.001,
            "Reconstruction is inconsistent with mean plus noise")
    model=Int8.(read(g["model_used"]));seed=UInt64.(read(g["noise_seed"]))
    require(length(model)==np && all(x->x in (0,1),model) && length(seed)==np,"Pad model/seed dimensions")
    reference=Float32.(read(g["reference_mm"]));phase=Float32.(read(g["phase_mm"]))
    width=Float32(read(g["bucket_width_mm"]));direction=Int8(read(g["drift_direction"]))
    require(length(reference)==3 && length(phase)==2 && all(isfinite,reference) && all(isfinite,phase) && width>0 && direction==1,
            "Invalid ion grid geometry")
    ReconstructedIonEvent(Int64(id),px,py,first,edges,dt,Float64(read(m["pitch_mm"])),
        reference,phase,width,direction,Float32(read(g["ion_drift_cm"])),Float32(read(g["electron_drift_cm"])),
        n,q,z,r,sigma,model,seed,str(g,"pair_checksum"))
end

"""Read one event into independent arrays; subsequent mutation leaves the HDF5 input intact."""
read_event(path::AbstractString,id::Integer) = h5open(h->read_event(h,id),path,"r")

"""Read the exact reconstructed temporal noise covariance, in elementary-charge squared."""
function read_covariance(h::HDF5.File,nframes::Integer)
    check_header(h)
    path="covariance_by_frames/$nframes/covariance_e2"
    require(haskey(h,path),"No covariance for $nframes frames")
    cv=Float64.(read(h[path]))
    require(size(cv)==(nframes,nframes) && all(isfinite,cv) && issymmetric(cv),"Invalid frame covariance")
    cv
end
read_covariance(path::AbstractString,nframes::Integer) = h5open(h->read_covariance(h,nframes),path,"r")

"""Replay all source cells against a corrected response and return a compact validation ledger.

Checks the upstream content hash, event/pair identities, exact pad and bucket
assignment, zero buckets, noise covariance and input charge. An optional source
path supports relocation while preserving its content identity.
"""
function validate_response(path::AbstractString;source_path=nothing)
    h5open(path,"r") do h
        m=check_header(h);source=source_path===nothing ? str(m,"source_path") : String(source_path)
        require(isfile(source) && sha(source)==str(m,"source_sha256"),"Upstream ion content changed")
        ids=Int64.(read(m["selected_event_ids"]));actual=sort!(parse.(Int64,collect(keys(h["events"]))))
        require(ids==sort!(unique(ids)) && actual==ids && length(ids)==read(m["event_count"]),"Event selection mismatch")
        total=Int64(0);pads=0;negatives=0;empty=0;variances=Dict{Int,Vector{Float64}}();pairs=String[]
        h5open(source,"r") do upstream
            require(schema(upstream)=="nausicaa.ion_pad_file.pilot.v1" && Bool(read(upstream["metadata/complete"])),"Incomplete/unsupported ion input")
            require(str(upstream["metadata"],"sample")==str(m,"sample"),"Sample identity mismatch")
            require(read(upstream["metadata/event_count"])==read(m["source_event_count"]),"Source event count mismatch")
            card=TOML.parse(str(upstream["metadata"],"branch_card_toml"))
            require(isapprox(card["pad_pitch_mm"],read(m["pitch_mm"])) &&
                    isapprox(card["bucket_width_ms"]*1e-3,read(m["frame_s"])),"Detector/response geometry mismatch")
            for id in ids
                e=read_event(h,id);u=upstream["events/$id"]
                require(e.pair_checksum==str(u,"pair_checksum"),"Primary-pair checksum mismatch")
                for key in ("reference_mm","phase_mm","ion_drift_cm","electron_drift_cm","drift_direction")
                    require(read(h["events/$id/$key"])==read(u[key]),"Event placement changed: $key")
                end
                x=Int64.(read(u["pixel_ix"]));y=Int64.(read(u["pixel_iy"]));b=Int64.(read(u["bucket"]));c=Int64.(read(u["count"]))
                positions=Dict(k=>i for (i,k) in enumerate(zip(e.pixel_ix,e.pixel_iy)))
                require(Set(zip(x,y))==Set(keys(positions)),"Input-charged pixel selection changed")
                require(length(x)==length(y)==length(b)==length(c) && all(>(0),c),"Invalid sparse input columns")
                require(allunique(collect(zip(x,y,b))),"Duplicate sparse input cells")
                for k in eachindex(c)
                    row=b[k]-e.first_bucket+1
                    require(1<=row<=size(e.input_ions,1) && e.input_ions[row,positions[(x[k],y[k])]]==c[k],"Sparse-to-dense cell mismatch")
                end
                require(count(>(0),e.input_ions)==length(c) && sum(c)==sum(e.input_ions)==read(u["initial_pairs"]),"Missing/extra input cells")
                nf=size(e.input_ions,1)
                diagonal=get!(variances,nf) do
                    cv=read_covariance(h,nf)
                    require(isposdef(Symmetric(cv)),"Nonpositive frame covariance")
                    diag(cv)
                end
                require(isapprox(diagonal,e.noise_sigma_ions.^2;rtol=1e-10),"Noise sigma/covariance mismatch")
                total+=sum(c);pads+=length(e.pixel_ix);negatives+=count(<(0),e.reconstructed_ions);empty+=count(==(0),e.input_ions)
                push!(pairs,"$id:$(e.pair_checksum)")
            end
        end
        require(total==read(m["input_ion_total"]),"Response total ion ledger mismatch")
        (;sample=str(m,"sample"),pitch_mm=Int(read(m["pitch_mm"])),frame_s=Float64(read(m["frame_s"])),
          selected_event_ids=ids,event_count=length(ids),input_ion_total=total,pads,negative_estimates=negatives,
          zero_input_buckets=empty,pair_identity=S.content_hash(join(pairs,"|")),
          source_path=source,source_sha256=str(m,"source_sha256"),response_id=str(m,"response_id"))
    end
end

"""Resolve production cards and saved-study branches; numerical acquisition remains in the validated backend."""
function settings(root::AbstractString,card::AbstractString)
    root=abspath(root);path=resolve(root,card);r=TOML.parsefile(path)
    require(r["schema"]=="nausicaa.ion_csa_pipeline_run.v1","Production run schema mismatch")
    require(r["channel_selection"]==CHANNEL_POLICY && r["absent_channel_policy"]==ABSENT_POLICY,"Channel selection policy mismatch")
    require(!r["measured_zero_suppression"] && !r["negative_estimate_clipping"],"Signed input-selected output required")
    cp=resolve(root,r["response_card"]);rp=resolve(root,r["reconstruction_card"])
    c=TOML.parsefile(cp);rc=TOML.parsefile(rp)
    require(c["schema"]=="nausicaa.ion_csa_spec.v1" && rc["schema"]=="nausicaa.ion_csa_reconstruction.v1","Response/reconstruction schema mismatch")
    require(rc["method"]=="three_frame_inverse" && !rc["signal_to_noise_cut"] && !rc["nonnegative_constraint"],"Signed three-frame reconstruction required")
    require(r["selection"]=="evenly_spaced_sorted_event_ids" && r["events_per_sample"]>0,"Event selection mismatch")
    require(!isempty(r["samples"]) && allunique(r["samples"]) && all(x->x in ("bb0nu","single_electrons"),r["samples"]),"Sample selection mismatch")
    branches=r["branches"];pitches=[b["pitch_mm"] for b in branches]
    require(!isempty(pitches) && allunique(pitches) && all(x->x in (2,3,4),pitches),"Select unique available pitches")
    require(isfile(resolve(root,r["ion_index"])),"Ion index is absent")
    files=[@__FILE__,joinpath(root,"Manifest.toml")]
    append!(files,[joinpath(root,p) for p in (
        "scripts/production/specs.jl","scripts/production/run_ion_csa.jl",
        "scripts/production/run_ion_csa_pilot.py","scripts/production/ion_csa_pilot.jl",
        "scripts/production/run_ion_csa_correction.py","scripts/production/ion_csa_correction.jl",
        "src/ion_signal/ion_csa_response.jl","src/ion_signal/ion_csa_correction.jl",
        "src/ion_signal/ion_csa_summary.jl","src/ion_signal/pad_history_readout.jl",
        "src/ion_signal/bucket_study.jl","src/ion_signal/averaging_scan.jl",
        "src/ion_signal/discharge_calibration.jl","src/ion_signal/bucket_noise_covariance.jl",
        "src/ion_signal/spice_acquisition.py","src/ion_signal/spice.py",
        "src/event_analysis/production_readout.jl")])
    codehash=S.content_hash(join(sha.(files),"|")*string(VERSION))
    calibration_file=resolve(root,joinpath(c["calibration_directory"],"analysis.h5"))
    require(isfile(calibration_file),"Selected circuit calibration is absent; generate it with the calibration card before running production")
    calibration_hash=sha(calibration_file)
    identity=S.content_hash(S.canonical(r)*S.canonical(c)*S.canonical(rc)*codehash*sha(resolve(root,r["ion_index"]))*calibration_hash)
    dir=resolve(root,joinpath(r["results_directory"],identity[1:16]))
    (;root,path,r,cp,rp,c,rc,branches,files,codehash,calibration_file,calibration_hash,identity,dir)
end

function write_toml(path,values)
    mkpath(dirname(path));temp=path*".tmp"
    open(temp,"w") do io;TOML.print(io,values);end
    mv(temp,path;force=true);path
end

function study_responses(s,study)
    path=resolve(s.root,study);require(isfile(joinpath(path,"summary.h5")),"Corrected study summary is absent: $path")
    h5open(joinpath(path,"summary.h5"),"r") do h
        paths=String[]
        for name in keys(h)
            g=h[name]
            g isa HDF5.Group && haskey(attributes(g),"response_path") || continue
            a=attributes(g);p=String(read(a["response_path"]))
            require(haskey(a,"response_sha256") && isfile(p) && sha(p)==String(read(a["response_sha256"])),
                    "Saved corrected response differs from its study summary: $p")
            push!(paths,p)
        end
        paths
    end
end

"""Run one selected geometry using the existing LTspice acquisition and Julia reconstruction stages."""
function simulate_branch(s,b;python,simulator=nothing,workers=4)
    pitch=b["pitch_mm"];base=joinpath(s.dir,"execution","pitch$pitch");mkpath(base)
    raw=Dict("schema"=>"nausicaa.ion_csa_run.v1","response_card"=>s.r["response_card"],"ion_index"=>s.r["ion_index"],
        "samples"=>s.r["samples"],"pitches_mm"=>[pitch],"events_per_sample"=>s.r["events_per_sample"],"selection"=>s.r["selection"],
        "results_directory"=>relpath(joinpath(base,"filtered"),s.root),"report_directory"=>relpath(joinpath(base,"filtered_report"),s.root))
    rawcard=write_toml(joinpath(base,"filtered_run.toml"),raw)
    args=[python,joinpath(s.root,"scripts/production/run_ion_csa_pilot.py"),"--config",rawcard,"--workers",string(workers)]
    simulator===nothing || append!(args,["--simulator",simulator])
    Base.run(Cmd(Cmd(args);dir=s.root))
    choices=filter(p->isfile(joinpath(p,"summary.h5")),readdir(joinpath(base,"filtered");join=true))
    require(length(choices)==1,"Expected one completed filtered study")
    reconstruction=Dict("schema"=>"nausicaa.ion_csa_reconstruction_run.v1","parent_study"=>relpath(only(choices),s.root),
        "reconstruction_card"=>s.r["reconstruction_card"],"results_directory"=>relpath(joinpath(base,"corrected"),s.root),
        "report_directory"=>relpath(joinpath(base,"corrected_report"),s.root))
    rcard=write_toml(joinpath(base,"reconstruction_run.toml"),reconstruction)
    Base.run(Cmd(Cmd([python,joinpath(s.root,"scripts/production/run_ion_csa_correction.py"),"--config",rcard]);dir=s.root))
    choices=filter(p->isfile(joinpath(p,"summary.h5")),readdir(joinpath(base,"corrected");join=true))
    require(length(choices)==1,"Expected one completed corrected study")
    only(choices)
end

function check_selection(s,path,ledger,pitch)
    require(ledger.pitch_mm==pitch && ledger.sample in s.r["samples"],"Saved study sample/pitch mismatch")
    h5open(resolve(s.root,s.r["ion_index"]),"r") do h
        require(sha(str(h,"$(ledger.sample)/ions/$pitch/path"))==ledger.source_sha256,"Saved study is from a different ion index")
    end
    h5open(ledger.source_path,"r") do h
        ids=sort!(parse.(Int64,collect(keys(h["events"]))));n=s.r["events_per_sample"]
        require(n<=length(ids),"Event selection exceeds source")
        expected=n==1 ? ids[1:1] : ids[round.(Int,range(1,length(ids);length=n))]
        require(expected==ledger.selected_event_ids,"Saved response event selection differs")
    end
    h5open(path,"r") do h
        require(S.canonical(TOML.parse(str(h["metadata"],"response_card_toml")))==S.canonical(s.c),"Saved circuit/noise card differs")
        require(S.canonical(TOML.parse(str(h["metadata"],"reconstruction_card_toml")))==S.canonical(s.rc),"Saved reconstruction card differs")
        m=h["metadata"]
        expected=S.content_hash(join((S.canonical(s.c),s.calibration_hash,str(m,"parent_code_hash")),"|"))
        require(expected==str(m,"physical_model_id"),"Saved response used a different circuit calibration")
    end
end

function publish(s,parent,ledger)
    parenthash=sha(parent)
    nodeid=S.content_hash(parenthash*s.codehash*CHANNEL_POLICY*ABSENT_POLICY*PRODUCTION_SCHEMA)
    # Corrected responses are response/ions/<ion_output_id>/<response_id>/ion_response.h5.
    target=joinpath(dirname(dirname(parent)),nodeid,"ion_response.h5")
    mkpath(dirname(target));write(joinpath(dirname(target),".gitignore"),"*.h5\n*.h5.tmp\n")
    if isfile(target)
        manifest=joinpath(dirname(target),"production_manifest.toml")
        require(isfile(manifest) && get(TOML.parsefile(manifest),"sha256","")==sha(target),"Published response checksum mismatch")
        h5open(target,"r") do h
            check_header(h)
            require(str(h["metadata"],"parent_response_sha256")==parenthash && str(h["metadata"],"response_id")==nodeid,"Published cache identity mismatch")
        end
        validate_response(target)
        return target
    end
    tmp=target*".tmp";cp(parent,tmp;force=true)
    h5open(tmp,"r+") do h
        HDF5.delete_attribute(h,"schema");attributes(h)["schema"]=PRODUCTION_SCHEMA;m=h["metadata"];write(m["complete"],false)
        function put(k,v)
            haskey(m,k) && HDF5.delete_object(m,k)
            m[k]=v
        end
        put("parent_response_path",parent);put("parent_response_sha256",parenthash)
        put("response_id",nodeid);put("production_code_hash",s.codehash)
        put("channel_selection",CHANNEL_POLICY);put("absent_channel_policy",ABSENT_POLICY)
        put("bucket_selection","complete_event_interval_with_quiet_and_recovery_frames")
        put("measured_zero_suppression",false);put("negative_estimate_clipping",false)
        put("numerical_payload","unchanged_from_corrected_parent")
        put("production_run_card_toml",read(s.path,String))
        write(m["complete"],true)
    end
    validate_response(tmp)
    mv(tmp,target)
    write_toml(joinpath(dirname(target),"production_manifest.toml"),Dict("schema"=>"nausicaa.ion_response_node.v1",
        "response_id"=>nodeid,"response_file"=>target,"sha256"=>sha(target),"parent_response"=>parent,
        "parent_sha256"=>parenthash,"channel_selection"=>CHANNEL_POLICY,"absent_channel_policy"=>ABSENT_POLICY,"complete"=>true))
    target
end

"""Produce or register all selected ion nodes and return the complete run manifest path.

Saved-study branches publish existing corrected arrays after full validation.
Branches without corrected_study run acquisition and reconstruction first.
The output manifest exposes reusable ion nodes for a later combined detector run.
"""
function run(root,card;python=joinpath(root,"work/event-viewer-venv/bin/python"),simulator=nothing,workers=4)
    s=settings(root,card);mkpath(s.dir);write(joinpath(s.dir,".gitignore"),"execution/\n*.h5\n*.tmp\n")
    responses=Dict{String,Any}[];shared=Dict{String,Tuple{Vector{Int64},String}}()
    for b in s.branches
        study=get(b,"corrected_study",nothing)
        study===nothing && (study=simulate_branch(s,b;python,simulator,workers))
        paths=study_responses(s,study)
        require(length(paths)==length(s.r["samples"]),"Study sample count differs")
        seen=Set{String}()
        for parent in paths
            ledger=validate_response(parent);check_selection(s,parent,ledger,b["pitch_mm"])
            require(!(ledger.sample in seen),"Duplicate sample in study");push!(seen,ledger.sample)
            signature=(ledger.selected_event_ids,ledger.pair_identity)
            require(get(shared,ledger.sample,signature)==signature,"Event IDs/primary pairs differ across pitches")
            shared[ledger.sample]=signature
            output=publish(s,parent,ledger)
            values=Dict{String,Any}(string(k)=>v for (k,v) in pairs(ledger))
            values["parent_corrected_response"]=parent;values["path"]=output;values["sha256"]=sha(output)
            h5open(output,"r") do h;values["response_id"]=str(h["metadata"],"response_id");end
            push!(responses,values)
            println("Published $(ledger.sample), $(ledger.pitch_mm) mm: $output")
        end
        require(seen==Set(s.r["samples"]),"Study sample set differs")
    end
    graphid=S.content_hash(s.identity*S.canonical([Dict("path"=>r["path"],"sha256"=>r["sha256"]) for r in responses]))
    manifest=Dict("schema"=>"nausicaa.ion_response_run.v1","graph_id"=>graphid,"run_identity"=>s.identity,
        "complete"=>true,"run_card"=>relpath(s.path,s.root),"run_card_toml"=>read(s.path,String),
        "production_code_hash"=>s.codehash,"calibration_file"=>s.calibration_file,"calibration_sha256"=>s.calibration_hash,
        "channel_selection"=>CHANNEL_POLICY,"absent_channel_policy"=>ABSENT_POLICY,
        "matrix_axes_hdf5"=>"pad,frame","matrix_axes_julia"=>"frame,pad","responses"=>responses)
    # Each source realisation exposes the same ion-stage manifest for final graph assembly.
    for base in unique(dirname(dirname(dirname(dirname(r["source_path"])))) for r in responses)
        write_toml(joinpath(base,"runs",graphid,"ion_manifest.toml"),manifest)
    end
    for p in s.files
        dest=joinpath(s.dir,"source_snapshot",sha(p),basename(p));mkpath(dirname(dest));cp(p,dest;force=true)
    end
    path=write_toml(joinpath(s.dir,"manifest.toml"),manifest)
    write_toml(joinpath(dirname(s.dir),"latest.toml"),Dict("manifest"=>path,"run_identity"=>s.identity,"graph_id"=>graphid))
    println("Ion response manifest: $path")
    path
end
end
