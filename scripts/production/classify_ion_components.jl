#!/usr/bin/env julia
include("prepare_hit_presentation.jl")

"""Compare detached graph components with the original MC gamma deposits."""
function classify(path, graph_card, mc_path, raw_path, output; arm="itrack")
    c=TOML.parsefile(graph_card)["track_analysis"]
    prefix=arm=="itrack" ? "ion" : "electron"
    cut=Float32(c["$(prefix)_hit_cut_keV"])
    width=Float32(c["$(prefix)_voxel_width_mm"])
    radius=Float32(c["$(prefix)_connection_mm"])
    r=open_hit_input(path)
    ids=Int64.(read(r.file["events/event_number"]))
    particle_tables=EventAnalysis.read_mc_particles(raw_path,ids)
    records=NamedTuple[]
    try
        h5open(mc_path,"r") do mc
            rows=Dict(id=>j for (j,id) in enumerate(read(mc["events/event_id"])))
            for id in ids
                event=read_hit_event(r,id)
                h=getfield(event,Symbol(arm))
                keep=h.energy_keV .> cut
                df=EventAnalysis.DataFrame(x=h.x_mm[keep],y=h.y_mm[keep],z=h.z_mm[keep],energy=h.energy_keV[keep]./1000f0)
                tracks=make_tracks(voxelize_event(df,width);method="KDT",max_distance_mm=radius,energy_threshold_kev=0f0)
                length(tracks)<=1 && continue

                source,info=EventAnalysis.read_mc_event(mc,id;event_row=rows[id])
                source=EventAnalysis.center_hits(source,info)
                roots=EventAnalysis.gamma_map(particle_tables[id],info)
                primary=findall(i->roots[source.particle_id[i]]==0,eachindex(source.particle_id))
                gamma=findall(i->roots[source.particle_id[i]]!=0,eachindex(source.particle_id))
                # Classify every secondary component by its proximity to MC hits.
                fractions=Float64[]
                for track in tracks[2:end]
                    v=track.voxels
                    assigned=0.0;total=0.0
                    for row in eachrow(v)
                        e=Float64(row.energy);total+=e
                        dg=isempty(gamma) ? Inf : minimum((row.x-source.x_mm[i])^2+(row.y-source.y_mm[i])^2+(row.z-source.z_mm[i])^2 for i in gamma)
                        dp=minimum((row.x-source.x_mm[i])^2+(row.y-source.y_mm[i])^2+(row.z-source.z_mm[i])^2 for i in primary)
                        if dg<=radius^2 && dg<dp
                            assigned+=e
                        end
                    end
                    push!(fractions,total==0 ? 0.0 : assigned/total)
                end
                energies=Float64[1000sum(t.voxels.energy) for t in tracks]
                push!(records,(event_number=id,component_count=length(tracks),second_keV=energies[2],
                    secondary_gamma_fraction=maximum(fractions),all_secondary_gamma=all(>=(0.8),fractions),
                    any_secondary_gamma=any(>=(0.8),fractions)))
            end
        end
    finally
        close_hit_input(r)
    end
    h5open(output,"w") do f
        f["schema"]="itaca.component_mc_classification.v1"
        f["arm"]=arm
        f["hit_file"]=path
        f["graph_card"]=read(graph_card,String)
        f["rule"]="A secondary graph component is matched to a gamma deposit when at least 80% of its energy lies within the graph radius of MC gamma hits and is nearer to those hits than to MC non-gamma hits."
        for key in keys(first(records))
            f[string(key)]=[getfield(x,key) for x in records]
        end
    end
    println("MULTI $(length(records)) ALL_GAMMA $(count(x->x.all_secondary_gamma,records)) ANY_GAMMA $(count(x->x.any_secondary_gamma,records)) NONE_GAMMA $(count(x->!x.any_secondary_gamma,records))")
end

if abspath(PROGRAM_FILE)==@__FILE__
    classify(abspath(ARGS[1]),abspath(ARGS[2]),abspath(ARGS[3]),abspath(ARGS[4]),abspath(ARGS[5]);arm=length(ARGS)>5 ? ARGS[6] : "itrack")
end
