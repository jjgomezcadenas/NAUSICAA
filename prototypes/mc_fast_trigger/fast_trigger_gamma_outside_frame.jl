using HDF5, TOML, SHA, Statistics, Random
const RUN_START=time()
include("fast_trigger_gamma_outside_frame_functions.jl")
using .MCFastTrigger
const EA=MCFastTrigger.EA
const ROOT=normpath(joinpath(@__DIR__,"../.."))
card=TOML.parsefile(joinpath(@__DIR__,"fast_trigger_gamma_outside_frame.toml"))
out=joinpath(@__DIR__,"results","fast_trigger_gamma_outside_frame.h5")
settings=EA.FastTriggerSettings(Float32(card["detector_length_mm"]),Float32(card["vertex_fiducial_min_mm"]),Float32(card["vertex_fiducial_max_mm"]),ntuple(_->Float32(card["illustrative_margin_mm"]),3),Float32(card["gamma_energy_cut_keV"]))
function columns(g,rows)
    isempty(rows) && return
    for k in keys(first(rows));g[string(k)]=[getproperty(r,k) for r in rows];end
end
function distributions(g,rows)
    isempty(rows) && return
    for k in keys(first(rows))
        x=[getproperty(r,k) for r in rows]
        eltype(x)<:Number || continue
        q=create_group(g,string(k));s=MCFastTrigger.summary(x)
        for a in keys(s);q[string(a)]=getproperty(s,a);end
        low,high=extrema(Float64.(x));high==low && (high=low+1)
        edges=collect(range(low,high;length=card["histogram_bins"]+1))
        if low>=0 && (occursin("outside",string(k)) || occursin("nearest_hit",string(k)) || string(k)=="energy_keV")
            edges=expm1.(range(log1p(low),log1p(high);length=card["histogram_bins"]+1));edges[1]=low;edges[end]=high
        end
        counts=zeros(Int64,length(edges)-1)
        for v in x;counts[clamp(searchsortedlast(edges,v),1,length(counts))]+=1;end
        q["edges"]=edges;q["counts"]=counts
    end
end
h5open(out*".tmp","w") do output
    output["card"]=read(joinpath(@__DIR__,"fast_trigger_gamma_outside_frame.toml"),String)
    output["schema"]="itaca.mc_fast_trigger.fiducial.v4"
    codefiles=[joinpath(@__DIR__,"fast_trigger_gamma_outside_frame_functions.jl"),@__FILE__,joinpath(ROOT,"src/event_analysis/fast_trigger_types.jl"),joinpath(ROOT,"src/event_analysis/mc_fast_trigger.jl"),joinpath(ROOT,"src/event_analysis/io.jl")]
    output["implementation_files"]=relpath.(codefiles,ROOT)
    output["implementation_file_sha256"]=[bytes2hex(sha256(read(p))) for p in codefiles]
    output["implementation_sha256"]=bytes2hex(sha256(join([bytes2hex(sha256(read(p))) for p in codefiles],"|")))
    for species in card["species"]
        rng=Xoshiro(card["placement_seed"])
        t=time();cfg=card[species];hp=joinpath(ROOT,cfg["hits"]);pp=expanduser(cfg["particles"])
        g=create_group(output,species);g["hits_path"]=hp;g["particles_path"]=pp
        rawhash=open(sha256,pp)|>bytes2hex;g["particles_sha256"]=rawhash
        h5open(hp,"r") do h
            rawhash==read(h["metadata/source_sha256"]) || error("Source checksum mismatch")
            ids=read(h["events/event_id"])
            card["events_per_species"]=="all" || (ids=ids[1:card["events_per_species"]])
            starts=read(h["events/first_row"]);sizes=read(h["events/n_hits"])
            parts=Dict(id=>Any[] for id in ids)
            h5open(pp,"r") do raw
                ds=raw["MC/particles"];step=card["particle_block_rows"]
                for start in 1:step:length(ds)
                    for p in ds[start:min(start+step-1,length(ds))]
                        haskey(parts,p.event_id) && push!(parts[p.event_id],p)
                    end
                end
            end
            outputs=Dict{Int64,EA.FastTriggerOutput}();trigger_stats=EA.FastTriggerStatistics();evs=NamedTuple[];all_evs=NamedTuple[];gammas=NamedTuple[];outside_ids=Int64[];details=create_group(g,"event_details");calc=0.0
            for (i,id) in enumerate(ids)
                ix=starts[i]:starts[i]+sizes[i]-1
                names=(:x_mm,:y_mm,:z_mm,:energy_MeV,:particle_id)
                hits=NamedTuple{names}(Tuple(h["hits/"*string(k)][ix] for k in names))
                t0=time();r=mc_fast_trigger(hits,parts[id];margin_mm=Float32(card["illustrative_margin_mm"]),expected_primaries=cfg["expected_primary_electrons"]);calc+=time()-t0
                placement=EA.place_event(rng,id,settings.chamber_length_mm)
                vertex_z=10f0*placement.ion_drift_cm
                mh,mi=EA.read_mc_event(h,id;event_row=i)
                mp=EA.MCParticle[EA.MCParticle(p.particle_id,p.mother_id,String(p.particle_name)) for p in parts[id]]
                decision=EA.mc_fast_trigger(EA.center_hits(mh,mi),mi,mp,placement,settings)
                outputs[id]=decision;trigger_stats=EA.update_statistics(trigger_stats,decision)
                selection=(contained=decision.decision!=EA.outside_detector,accepted=decision.decision in (EA.trigger_accepted,EA.outside_frame_gamma))
                row=merge((event_number=id,event_in_file=i,trigger_decision=UInt8(decision.decision),vertex_z_mm=vertex_z,ion_drift_cm=placement.ion_drift_cm,electron_drift_cm=placement.electron_drift_cm,frame_x_min_mm=decision.frame.lower_mm[1],frame_y_min_mm=decision.frame.lower_mm[2],frame_z_min_mm=decision.frame.lower_mm[3],frame_x_max_mm=decision.frame.upper_mm[1],frame_y_max_mm=decision.frame.upper_mm[2],frame_z_max_mm=decision.frame.upper_mm[3],contained=UInt8(selection.contained),accepted=UInt8(selection.accepted)),r.event)
                push!(all_evs,row)
                selection.accepted || continue
                push!(evs,row)
                any(v->max(v.outside_frame_x_mm,v.outside_frame_y_mm,v.outside_frame_z_mm)>0,r.gammas) && push!(outside_ids,id)
                append!(gammas,[merge((event_number=id,),v) for v in r.gammas])

            end
            columns(create_group(g,"events"),evs);columns(create_group(g,"gammas"),gammas)
            distributions(create_group(g,"event_statistics"),evs);distributions(create_group(g,"gamma_statistics"),gammas)
            distributions(create_group(g,"inside_frame_gamma_statistics"),filter(v->v.inside_frame==1,gammas))
            outside_set=Set(outside_ids)
            inside_event_gammas=filter(v->!(v.event_number in outside_set),gammas)
            columns(create_group(g,"all_gammas_inside_frame"),inside_event_gammas)
            distributions(create_group(g,"all_gammas_inside_frame_statistics"),inside_event_gammas)
            g["all_gammas_inside_frame_event_count"]=length(evs)-length(outside_ids)
            g["all_gammas_inside_frame_events_with_deposits"]=length(unique(v.event_number for v in inside_event_gammas))
            columns(create_group(g,"all_events"),all_evs)
            sel=create_group(g,"selection")
            total=length(all_evs);contained=count(e->e.contained==1,all_evs);accepted=length(evs)
            for (k,v) in [("input",total),("contained",contained),("rejected_containment",total-contained),("accepted",accepted),("rejected_vertex",contained-accepted)]
                sel[k]=v
            end
            sel["containment_percent"]=100contained/total;sel["fiducial_percent_of_contained"]=100accepted/contained
            sel["normalization_events"]=contained
            sel["normalization_condition"]="Only events with all deposited hits inside detector Z = 0–1500 mm enter any efficiency numerator or denominator."
            threshold=card["gamma_energy_cut_keV"]
            deposits_by_event=Dict(e.event_number=>NamedTuple[] for e in evs)
            for v in gammas;push!(deposits_by_event[v.event_number],v);end
            survivors=filter(e->outputs[e.event_number].decision==EA.trigger_accepted,evs)
            survivor_ids=Set(e.event_number for e in survivors)
            surviving_outside=filter(v->v.inside_frame==0 && v.event_number in survivor_ids,gammas)
            columns(create_group(g,"after_gamma_cut_events"),survivors)
            columns(create_group(g,"surviving_outside_gamma_deposits"),surviving_outside)
            distributions(create_group(g,"surviving_outside_gamma_statistics"),surviving_outside)
            survivor_outside_ids=sort(unique(Int64[v.event_number for v in surviving_outside]))
            g["surviving_outside_example_numbers"]=survivor_outside_ids[1:min(5,length(survivor_outside_ids))]
            @assert length(survivors)==trigger_stats.accepted_events
            @assert contained==trigger_stats.input_events-trigger_stats.outside_detector_events
            for field in fieldnames(EA.FastTriggerStatistics);sel["count_"*string(field)]=getfield(trigger_stats,field);end
            sel["gamma_energy_cut_keV"]=threshold
            sel["accepted_gamma_cut"]=length(survivors)
            sel["rejected_gamma_cut"]=length(evs)-length(survivors)
            sel["gamma_cut_percent_of_fiducial"]=100length(survivors)/length(evs)
            sel["final_percent"]=100length(survivors)/contained
            sel["surviving_events_with_outside_gamma"]=length(survivor_outside_ids)
            sel["surviving_outside_gamma_deposits"]=length(surviving_outside)
            ordered=sort(evs;by=e->e.gamma_free_dz_mm)
            indices=unique(round.(Int,range(1,length(ordered);length=min(10,length(ordered)))))
            g["example_event_numbers"]=[ordered[i].event_number for i in indices[1:min(5,length(indices))]]
            g["outside_frame_event_numbers"]=outside_ids[1:min(10,length(outside_ids))]
            sel["events_with_gamma_outside_frame"]=length(outside_ids)
            # Save hit arrays only for displayed examples. Tables retain every event.
            example_ids=unique(vcat([ordered[i].event_number for i in indices[1:min(5,length(indices))]],survivor_outside_ids[1:min(5,length(survivor_outside_ids))]))
            for id in example_ids
                i=findfirst(==(id),ids);ix=starts[i]:starts[i]+sizes[i]-1
                names=(:x_mm,:y_mm,:z_mm,:energy_MeV,:particle_id)
                hits=NamedTuple{names}(Tuple(h["hits/"*string(k)][ix] for k in names))
                r=mc_fast_trigger(hits,parts[id];margin_mm=Float32(card["illustrative_margin_mm"]),expected_primaries=cfg["expected_primary_electrons"])
                vertex_z=all_evs[i].vertex_z_mm
                d=create_group(details,string(id))
                for k in keys(r.hits);d[string(k)]=(k==:z_mm ? r.hits.z_mm.+vertex_z : getproperty(r.hits,k));end
                d["primary_low_mm"]=r.lo.+Float32[0,0,vertex_z];d["primary_high_mm"]=r.hi.+Float32[0,0,vertex_z];d["vertex_z_mm"]=vertex_z
            end
            g["calculation_seconds"]=calc
        end
        g["elapsed_seconds"]=time()-t
        println(species," completed in ",round(time()-t;digits=2)," s");flush(stdout)
    end
end
filesize(out*".tmp")<=card["maximum_output_bytes"] || error("Prototype output exceeds size limit")
mv(out*".tmp",out;force=true)
println(out)
println("Julia total seconds including compilation: ",round(time()-RUN_START;digits=3))
