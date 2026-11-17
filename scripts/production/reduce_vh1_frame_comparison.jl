#!/usr/bin/env julia
# VH1 frame comparison: the same accepted events read with two run cards that
# differ in the ion drift velocity, the ion bucket and the VH1 frame. For each
# card the ion detector response and the VH1 readout of every accepted event of
# the given batches are recomputed with the production seeds (the electron arm
# is not needed), and per-event diagnostics, the response-matrix conditioning
# and histograms are saved. Event energies and hits come from the saved hit
# files of both runs.
#   julia --project=. scripts/production/reduce_vh1_frame_comparison.jl config/runs/<run>.toml config/runs/<reference>.toml FIRST LAST
# Output: results/readout/<run_label>/vh1_frame_comparison/comparison.h5
include("run_hit_production.jl")
using Statistics,LinearAlgebra

const SNR_EDGES=collect(0.:.5:40.)            # signal-to-noise per pad-bucket, overflow in the last bin
const RESIDUAL_EDGES=collect(-150.:3.:150.)   # reconstructed minus input ions per pad-bucket

"""Histogram counts of x on fixed edges; values outside land in the end bins."""
function fixed_histogram!(counts,edges,x)
    for v in x;counts[clamp(searchsortedlast(edges,v),1,length(counts))]+=1;end
    counts
end

"""Selection of one batch: accepted event numbers, placements and frames, after the MC source check."""
function batch_selection(rc)
    hp=input_path(ROOT,rc["hits"]);sha=bytes2hex(open(sha256,hp))
    source=joinpath(ROOT,"results/readout",rc["run_label"],"detector_response/controls.h5")
    h5open(source,"r") do f
        read(f["source_sha256"],String)==sha||error("MC source changed: $source")
        g=f["events"];ids=read(g["event_number"]);n=rc["hit_production"]["accepted_events"];n>0&&(ids=ids[1:min(n,length(ids))])
        ps=[EventPlacement(id,read(g["ion_drift_cm"])[j],read(g["electron_drift_cm"])[j]) for (j,id) in enumerate(ids)]
        fs=[EventFrame(mc_vertex_frame,ntuple(k->read(g["frame_$("xyz"[k])_min"])[j],3),ntuple(k->read(g["frame_$("xyz"[k])_max"])[j],3)) for j in eachindex(ids)]
        (hits=hp,ids=ids,placements=ps,frames=fs)
    end
end

"""Ion arm and VH1 readout of the accepted events of one run card, batch by batch, with the production seeds."""
function vh1_diagnostics(cardpath,batches)
    r0=resolve_run(cardpath;root=ROOT);response=r0.run["response"]
    vh=load_vh1_settings(ROOT,joinpath(ROOT,response["vh1_card"]);minimum_samples=get(response,"vh1_minimum_samples",32768))
    rows=NamedTuple[];plans=Dict{Int,Any}()
    snr=zeros(Int64,length(SNR_EDGES)-1);empty=zeros(Int64,length(RESIDUAL_EDGES)-1);signal=zeros(Int64,length(RESIDUAL_EDGES)-1)
    for b in batches
        r=resolve_run(cardpath;root=ROOT);rc=r.run;select_batch!(rc,b);ec=electronics_card(r)
        sel=batch_selection(rc)
        h5open(sel.hits,"r") do input
            rowmap=Dict(id=>j for (j,id) in enumerate(read(input["events/event_id"])))
            for (j,id) in enumerate(sel.ids)
                p=sel.placements[j];frame=sel.frames[j]
                det=@timed begin
                    h,info=read_mc_event(input,id;event_row=rowmap[id]);h=hits_in_frame(center_hits(h,info),frame);ds=event_settings(r,p)
                    rng(tag)=Xoshiro(seedfor(rc["seed"],id,tag))
                    pairs=sample_fano_pairs(rng("pairs"),h,Float32(r.detector["energy"]["W_eV"]),Float32(r.detector["energy"]["fano_factor"]))
                    sim_ion_detector_response(rng("ions"),h,pairs,p,frame,ds;max_vector_bytes=rc["max_vector_bytes"])
                end
                ions=det.value;erng(tag)=Xoshiro(seedfor(ec["seed"],id,tag))
                fit=@timed begin
                    adc=sim_vh1_electronics(erng("vh1_noise"),erng("vh1_adc"),ions,vh)
                    nb=Int(ions.ions.time.bucket_count)
                    plan=get!(plans,nb) do;prepare_vh1_reconstruction(adc.time,ions.ions.time,vh);end
                    reco_vh1_ions(adc,plan),adc,plan
                end
                ri,adc,plan=fit.value
                counts=ions.ions.counts;valid=ri.valid;noise=sqrt.(diag(plan.ion_covariance))
                has_input=(counts.>0).&valid;is_empty=(counts.==0).&valid
                q=Float64.(ri.ions)
                res_signal=[q[i,k]-counts[i,k] for k in axes(q,2),i in axes(q,1) if has_input[i,k]]
                res_empty=[q[i,k] for k in axes(q,2),i in axes(q,1) if is_empty[i,k]]
                snr_cells=[q[i,k]/noise[i] for k in axes(q,2),i in axes(q,1) if has_input[i,k]]
                fixed_histogram!(snr,SNR_EDGES,snr_cells);fixed_histogram!(empty,RESIDUAL_EDGES,res_empty);fixed_histogram!(signal,RESIDUAL_EDGES,res_signal)
                push!(rows,(event_number=id,batch=b,ion_buckets=nb,ion_pads=size(counts,2),ions_input=sum(counts),
                    signal_pad_buckets=count(has_input),signal_input_ions=sum(counts[has_input]),signal_reconstructed_ions=sum(q[has_input];init=0.),signal_snr_sum=sum(snr_cells;init=0.),
                    saturated_pads=count(adc.saturated),adc_limit_pads=count(ri.checks.adc_at_limit),invalid_pads=count(.!ri.checks.converged),
                    ions_on_invalid_pads=sum(counts[.!valid];init=0),reconstructed_ions=sum(q[valid];init=0.),
                    bucket_residual_rms=sqrt(mean(abs2,q[valid].-counts[valid])),signal_residual_rms=isempty(res_signal) ? 0. : sqrt(mean(abs2,res_signal)),
                    empty_cells=length(res_empty),empty_sum=sum(res_empty;init=0.),empty_sumsq=sum(abs2,res_empty;init=0.),
                    predicted_bucket_noise=sqrt(mean(diag(plan.ion_covariance))),detector_s=det.time,vh1_s=fit.time))
                println("$(rc["run_label"]) event $id ($j/$(length(sel.ids))): buckets $nb, saturated $(count(adc.saturated)), invalid $(count(.!ri.checks.converged))");flush(stdout)
            end
        end
    end
    conditioning=Dict{String,Any}()
    for (nb,plan) in sort(collect(plans);by=first)
        W=plan.lower\plan.response;g=vh.gain_V_per_ion
        conditioning["buckets_$nb"]=Dict("buckets"=>nb,"readings"=>size(plan.response,1),"condition_response"=>cond(plan.response),"condition_whitened"=>cond(W),
            "predicted_bucket_noise_ions"=>sqrt.(diag(plan.ion_covariance)),"reading_drop_noise_ions"=>sqrt.(diag(plan.covariance))./g,
            "response_matrix_ions"=>plan.response./g)
    end
    card=TOML.parsefile(joinpath(ROOT,response["vh1_card"]))
    (rows=rows,snr=snr,empty=empty,signal=signal,conditioning=conditioning,run_label=r0.run["run_label"],
     frame_ms=card["frame_ms"],conversion_window_ms=card["conversion_window_ms"],ion_drift_velocity_cm_s=r0.detector["drift"]["ion_drift_velocity_cm_s"],
     ion_bucket_ms=r0.detector["ions"]["ion_bucket_ms"],run_card=read(cardpath,String),vh1_card=read(joinpath(ROOT,response["vh1_card"]),String),
     overrides=Dict(string(k)=>v for (k,v) in r0.overrides))
end

"""Event energies and hits of a run from the hit files of the given batches, and the itrack hits at or above `cut_keV`."""
function hit_energies(label,batches,cut_keV)
    out=Dict{String,Any}("event_number"=>Int64[],"ion_energy_keV"=>Float32[],"ions"=>Float32[],"invalid_ion_buckets"=>Int64[],"itrack_hits"=>Int64[],"itrack_hits_above_cut"=>Int64[],"hit_cut_keV"=>cut_keV)
    for b in batches
        h5open(joinpath(ROOT,"results/readout",label,"batches",batch_tag(b),"vh1_hits/hits.h5"),"r") do f
            g=f["events"];ids=read(g["event_number"]);e=read(f["itrack/energy_keV"]);off=read(f["itrack/offsets"])
            append!(out["event_number"],ids);append!(out["ion_energy_keV"],read(g["ion_energy_keV"]));append!(out["ions"],read(g["ions"]))
            append!(out["invalid_ion_buckets"],read(g["invalid_ion_buckets"]))
            append!(out["itrack_hits"],[off[j+1]-off[j] for j in eachindex(ids)])
            append!(out["itrack_hits_above_cut"],[count(>=(cut_keV),view(e,off[j]+1:off[j+1])) for j in eachindex(ids)])
            cut=read(f["configuration/ion_hit_cut_keV"]);haskey(out,"output_cut_keV")&&out["output_cut_keV"]!=cut&&error("Output cut differs between batches")
            out["output_cut_keV"]=cut
        end
    end
    out
end

function write_run!(f,name,d,energies)
    g=create_group(f,name)
    for k in (:run_label,:frame_ms,:conversion_window_ms,:ion_drift_velocity_cm_s,:ion_bucket_ms,:run_card,:vh1_card);g[string(k)]=getfield(d,k);end
    write_control_tree(create_group(g,"overrides"),d.overrides)
    columns(create_group(g,"events"),d.rows)
    h=create_group(g,"histograms");h["snr_edges"]=SNR_EDGES;h["snr_counts"]=d.snr
    h["residual_edges"]=RESIDUAL_EDGES;h["empty_residual_counts"]=d.empty;h["signal_residual_counts"]=d.signal
    write_control_tree(create_group(g,"conditioning"),d.conditioning)
    write_control_tree(create_group(g,"hits"),energies)
end

function main(run_card,reference_card,first,last;cut_keV=5.0)
    BLAS.set_num_threads(1)
    batches=first:last
    new=vh1_diagnostics(abspath(run_card),batches);ref=vh1_diagnostics(abspath(reference_card),batches)
    ids_new=[r.event_number for r in new.rows];ids_ref=[r.event_number for r in ref.rows]
    ids_new==ids_ref||error("The two runs do not hold the same accepted events")
    e_new=hit_energies(new.run_label,batches,cut_keV);e_ref=hit_energies(ref.run_label,batches,cut_keV)
    e_new["event_number"]==e_ref["event_number"]||error("The chained hit files do not hold the same events")
    out=joinpath(ROOT,"results/readout",new.run_label,"vh1_frame_comparison");mkpath(out);target=joinpath(out,"comparison.h5")
    isfile(target)&&error("Comparison exists: $target")
    h5open(target,"w") do f
        f["schema"]="itaca.vh1_frame_comparison.v1";f["batches"]=collect(batches);f["events"]=length(new.rows)
        write_run!(f,"run",new,e_new);write_run!(f,"reference",ref,e_ref)
    end
    println("FINISHED $target: $(length(new.rows)) events, batches $(first)-$(last)")
    target
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==4||error("Usage: reduce_vh1_frame_comparison.jl run_card.toml reference_card.toml FIRST LAST")
    main(ARGS[1],ARGS[2],parse(Int,ARGS[3]),parse(Int,ARGS[4]))
end
