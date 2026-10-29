# Paired numerical/selection scans, reusing explicit event RNG streams.
using HDF5,DataFrames,TOML,SHA,Dates,Random,Statistics
module Step3
include(joinpath(@__DIR__,"validate_event_step3.jl"))
end
using .Step3.EventAnalysis
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const CARD=isempty(ARGS) ? joinpath(ROOT,"config/event_analysis.toml") : abspath(ARGS[1])
const OUT=joinpath(ROOT,"results/event_analysis/step4_scans.h5")

"""Generate baseline plus one-parameter scan variants from the TOML card."""
function scan_variants(config)
    baseline=(id="baseline",nbins=Int(config["diffusion"]["nbins"]),nsigma=Float32(config["diffusion"]["nsigma"]),
        k=Int(config["graph"]["k"]),threshold=Int64(config["graph"]["min_ions_per_bin"]),scale=1f0)
    result=[baseline]
    for (key,field,prefix,T) in (("histogram_bins",:nbins,"bins",Int),("histogram_padding_sigma",:nsigma,"padding",Float32),
            ("knn_k",:k,"k",Int),("min_ions_per_bin",:threshold,"ions",Int64),("graph_scale",:scale,"scale",Float32))
        for raw in sort(config["scans"][key])
            value=T(raw)
            value==getproperty(baseline,field) && continue
            patch=NamedTuple{(field,)}((value,))
            push!(result,merge(baseline,patch,(id=prefix*"_"*string(value),)))
        end
    end
    allunique([v.id for v in result]) || error("Repeated scan values")
    return result
end

"""
Diagnostic collection with the baseline accepted track's physical grid reference.
Extend axes only by whole cells; express keys in that fixed lattice so changing
coverage/time origin cannot shift the comparison. The third index is a signed
geometric z-cell label, used only for this paired diagnostic, not ensemble time.
"""
function fixed_grid_counts(base,diffused,config,detector)
    base.selection.accepted || return nothing
    a=pad_axes(diffused.bins,only(base.selection.tracks);pitch_mm=Float32(detector["geometry"]["pitch_m"])*1000f0,
        z_bin_width_mm=Float32(config["pads"]["z_bin_width_mm"]),grid_offset_mm=base.offset,
        minimum_xy_span_mm=Float32(config["pads"]["minimum_xy_span_mm"]))
    h=compute_occupancy_analysis(diffused.bins;x_edges_mm=a.x_edges_mm,y_edges_mm=a.y_edges_mm,
        z_edges_mm=a.z_edges_mm,drift_speed_mm_ms=Float32(detector["gas"]["bulk_speed_m_s"]),drift_direction=1)
    return DataFrame(pad_dx=h.counts.pad_ix.-1 .+a.first_edge_indices[1],
        pad_dy=h.counts.pad_iy.-1 .+a.first_edge_indices[2],
        arrival_bin=h.counts.slice_index.-1 .+a.first_edge_indices[3],ions=copy(h.counts.ions))
end

"""Reconstruct and collect a cached diffusion table under one graph variant."""
function collect_variant(base,diffused,variant,config,detector)
    b,d = base.b*variant.scale,base.distance*variant.scale
    params=DiffusionParams(base.L,base.sigma,base.sigma_z,b,d,
        Float32(variant.threshold)*Float32(config["ionization"]["ionisation_energy_eV"])/1000f0,variant.nbins,variant.nsigma)
    selection=select_diffused_event(diffused.bins;bin_width_mm=b,max_distance_mm=d,
        min_ions_per_bin=variant.threshold,k=variant.k,method=config["graph"]["method"],diffusion=params)
    axes,history,aligned=nothing,nothing,nothing
    if selection.accepted
        axes=pad_axes(diffused.bins,only(selection.tracks);pitch_mm=Float32(detector["geometry"]["pitch_m"])*1000f0,
            z_bin_width_mm=Float32(config["pads"]["z_bin_width_mm"]),grid_offset_mm=base.offset,
            minimum_xy_span_mm=Float32(config["pads"]["minimum_xy_span_mm"]))
        history=compute_occupancy_analysis(diffused.bins;x_edges_mm=axes.x_edges_mm,y_edges_mm=axes.y_edges_mm,
            z_edges_mm=axes.z_edges_mm,drift_speed_mm_ms=Float32(detector["gas"]["bulk_speed_m_s"]),
            drift_direction=Int(config["pads"]["drift_direction"]))
        aligned=align_pad_history(history,axes.reference_mm)
        @assert sum(aligned.ions)==diffused.retained_ions
    end
    return (;selection,axes,history,aligned)
end

"""
Normalised L1 distance between paired aligned pad/slice count tables:
sum(abs(A-B))/(sum(A)+sum(B)). Zero means identical; one means disjoint charge.
This measures the full pipeline response, including changes in its grid reference.
"""
function paired_distance(a,b)
    difference=Dict{NTuple{3,Int64},Int64}()
    for (sign,table) in ((1,a),(-1,b)),r in eachrow(table)
        key=(r.pad_dx,r.pad_dy,r.arrival_bin)
        difference[key]=get(difference,key,Int64(0))+sign*r.ions
    end
    denominator=sum(a.ions)+sum(b.ions)
    return denominator>0 ? Float32(sum(abs,values(difference)))/Float32(denominator) : 0f0
end

function run_scan()
    config=TOML.parsefile(CARD)
    config["ensemble"]["alignment"]=="central_pad" || error("This checkpoint implements central_pad alignment")
    detector_path=joinpath(ROOT,config["inputs"]["detector_card"])
    detector=TOML.parsefile(detector_path)
    tests=TOML.parsefile(joinpath(ROOT,"work/event_analysis/step4_validation.toml"))
    tests["status"]=="passed" || error("Run step-four tests first")
    for (p,h) in tests["files"];bytes2hex(sha256(read(joinpath(ROOT,p))))==h || error("Changed after test: $p");end
    variants=scan_variants(config)
    mkpath(dirname(OUT))
    h5open(OUT,"w") do f
        attributes(f)["schema"]="nausicaa.event.scans.v1"
        attributes(f)["status"]="running"
        attributes(f)["created_utc"]=string(now(UTC))
        attributes(f)["julia_version"]=string(VERSION)
        attributes(f)["pairing"]="Same events, drift, stage seeds and grid offsets; identical ion draws across histogram variants; cached diffusion for graph-only variants"
        f["configuration/card"]=read(CARD,String)
        f["configuration/detector"]=read(detector_path,String)
        f["tests/passed"]=tests["counts"]["passed"]
        paths=sort(collect(keys(tests["files"])))
        append!(paths,[relpath(CARD,ROOT),relpath(detector_path,ROOT),"scripts/validation/scan_event_analysis.jl","scripts/validation/validate_event_step3.jl"])
        f["provenance/paths"]=paths
        f["provenance/sha256"]=[bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in paths]
        Step3.save_table(f,"variants",DataFrame(variants))
        for input_relative in config["inputs"]["samples"]
            input=joinpath(ROOT,input_relative)
            sample=replace(basename(input),"_hitdf.h5"=>"")
            sg=create_group(f,"samples/"*sample)
            sg["input_sha256"]=bytes2hex(open(sha256,input))
            accepted_tables=Dict(v.id=>DataFrame[] for v in variants)
            measurements=Dict(v.id=>NamedTuple[] for v in variants)
            # Save every selected count with sample, variant and event identity.
            histories=Dict(v.id=>DataFrame(event_id=Int64[],pad_dx=Int64[],pad_dy=Int64[],arrival_bin=Int64[],ions=Int64[]) for v in variants)
            h5open(input,"r") do src
                ids=read(src["events/event_id"])
                n=min(Int(config["scans"]["events_per_sample"]),length(ids))
                n>0 || error("scans.events_per_sample must be positive")
                ids=ids[unique(round.(Int,range(1,length(ids);length=n)))]
                for (index,id) in enumerate(ids)
                    hits,_=read_event(src,id)
                    base=Step3.process_event(hits,config,detector,sample)
                    base_aligned=base.selection.accepted ? align_pad_history(base.history,base.axes.reference_mm) : nothing
                    cache=Dict((variants[1].nbins,variants[1].nsigma)=>base.diffused)
                    fixed_cache=Dict((variants[1].nbins,variants[1].nsigma)=>fixed_grid_counts(base,base.diffused,config,detector))
                    event_results=Dict{String,Any}()
                    for v in variants
                        key=(v.nbins,v.nsigma)
                        if !haskey(cache,key)
                            cache[key]=diffuse_xyz_image_mc(Xoshiro(base.seeds.diffusion),hits,base.counts.sampled;
                                sigma_t_mm=base.sigma,sigma_l_mm=base.sigma_z,nbins=v.nbins,nsigma=v.nsigma,
                                ionisation_energy_eV=Float32(config["ionization"]["ionisation_energy_eV"]))
                        end
                        d=cache[key]
                        if !haskey(fixed_cache,key);fixed_cache[key]=fixed_grid_counts(base,d,config,detector);end
                        result=collect_variant(base,d,v,config,detector)
                        s=result.selection
                        if v.id=="baseline"
                            @assert s.accepted==base.selection.accepted
                            s.accepted && @assert result.aligned==base_aligned
                        end
                        # Changing only bin count preserves every retained ion.
                        v.nsigma==variants[1].nsigma && @assert d.retained_ions==base.diffused.retained_ions
                        shared=s.accepted && base.selection.accepted
                        distance=shared ? paired_distance(base_aligned,result.aligned) : NaN32
                        fixed_distance=shared ? paired_distance(fixed_cache[(variants[1].nbins,variants[1].nsigma)],fixed_cache[key]) : NaN32
                        reference_id="baseline"
                        if startswith(v.id,"bins_") && v.nbins>variants[1].nbins
                            previous=maximum(filter(n->n<v.nbins,Int.(config["scans"]["histogram_bins"])))
                            reference_id=previous==variants[1].nbins ? "baseline" : "bins_"*string(previous)
                        end
                        event_results[v.id]=(;accepted=s.accepted,fixed=fixed_cache[key])
                        prior=event_results[reference_id]
                        resolution_shared=shared && prior.accepted
                        resolution_distance=resolution_shared ? paired_distance(prior.fixed,fixed_cache[key]) : NaN32
                        push!(measurements[v.id],(event_id=id,accepted=Int8(s.accepted),baseline_accepted=Int8(base.selection.accepted),
                            n_tracks=length(s.tracks),drift_cm=base.L,seed_diffusion=base.seeds.diffusion,
                            seed_drift=base.seeds.drift,seed_grid=base.seeds.grid,
                            offset_x_mm=base.offset[1],offset_y_mm=base.offset[2],
                            sampled_ions=d.sampled_ions,retained_ions=d.retained_ions,outside_ions=d.outside_ions,
                            graph_retained_ions=s.graph_retained_ions,graph_below_threshold_ions=s.graph_below_threshold_ions,
                            histogram_dx_mm=maximum(diff(d.x_edges_mm)),histogram_dy_mm=maximum(diff(d.y_edges_mm)),
                            histogram_dz_mm=maximum(diff(d.z_edges_mm)),paired_pad_distance=distance,
                            fixed_grid_pad_distance=fixed_distance,resolution_reference=reference_id,
                            resolution_pair_shared=Int8(resolution_shared),resolution_fixed_grid_distance=resolution_distance))
                        if s.accepted
                            push!(accepted_tables[v.id],result.aligned)
                            rows=copy(result.aligned);insertcols!(rows,1,:event_id=>fill(id,nrow(rows)))
                            append!(histories[v.id],rows)
                        end
                    end
                    if index%10==0 || index==length(ids)
                        println(sample,": ",index,"/",length(ids)," paired events across ",length(variants)," settings");flush(stdout)
                    end
                end
            end
            summaries=NamedTuple[]
            for v in variants
                g=create_group(sg,v.id)
                table=DataFrame(measurements[v.id])
                Step3.save_table(g,"events",table)
                Step3.save_table(g,"aligned_histories",histories[v.id])
                ensemble=ensemble_pad_maps(accepted_tables[v.id];processed_events=nrow(table),
                    pitch_mm=Float32(detector["geometry"]["pitch_m"])*1000f0,
                    z_bin_width_mm=Float32(config["pads"]["z_bin_width_mm"]),drift_speed_mm_ms=Float32(detector["gas"]["bulk_speed_m_s"]),
                    minimum_half_cells=Int(config["ensemble"]["minimum_half_cells"]))
                write_ensemble(create_group(g,"ensemble"),ensemble)
                keep=table.accepted.==1;basekeep=table.baseline_accepted.==1
                shared=keep.&basekeep
                @assert ensemble.total_ions==sum(table.retained_ions[keep])
                push!(summaries,(variant=v.id,processed=nrow(table),accepted=count(keep),
                    lost=count(basekeep.&.!keep),gained=count(.!basekeep.&keep),shared=count(shared),
                    acceptance_fraction=Float32(count(keep))/Float32(nrow(table)),
                    mean_shared_pad_distance=any(shared) ? mean(table.paired_pad_distance[shared]) : NaN32,
                    mean_fixed_grid_pad_distance=any(shared) ? mean(table.fixed_grid_pad_distance[shared]) : NaN32,
                    resolution_reference=first(table.resolution_reference),resolution_shared=count(table.resolution_pair_shared.==1),
                    mean_resolution_fixed_grid_distance=any(table.resolution_pair_shared.==1) ? mean(table.resolution_fixed_grid_distance[table.resolution_pair_shared.==1]) : NaN32,
                    accepted_ions=ensemble.total_ions,mean_ions_accepted=any(keep) ? Float32(ensemble.total_ions)/Float32(count(keep)) : NaN32,
                    mean_selected_ions_processed=Float32(ensemble.total_ions)/Float32(nrow(table)),
                    retained_ions_all=sum(table.retained_ions),outside_ions_all=sum(table.outside_ions)))
            end
            Step3.save_table(sg,"summary",DataFrame(summaries))
        end
        HDF5.delete_attribute(f,"status");attributes(f)["status"]="complete"
    end
    h5open(OUT,"r") do f
        for sample in keys(f["samples"]),v in variants
            g=f["samples/"*sample*"/"*v.id]
            e=read_ensemble(g["ensemble"])
            @assert sum(read(g["aligned_histories/ions"]))==e.total_ions
        end
    end
    println(OUT)
end
run_scan()
