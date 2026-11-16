#!/usr/bin/env julia
# Second pass on the itrack with shorter connection radii.
#   julia --project=. scripts/production/prepare_second_pass.jl results/readout/<run_label> [RADII_MM]
# For every accepted event of every batch (hits.h5), the itrack hits above the
# display cut are connected, at their positions on the readout lattice, with
# each radius of RADII_MM (default 3.5,4.6,5.5,6.7,8.0). Writes
# results/readout/<run_label>/<output_subdir>/second_pass/components_lattice.h5: for
# each radius, the hits, energy and largest extent of every component of
# every event. The events are in the order of the run statistics (events.h5).
include("prepare_batch_features.jl")

"""Components of the itrack of one event for one connection radius: hits, energy (keV) and largest extent (mm) of each, largest first.

Hits above the display cut are connected when their centres are closer than
the radius. The hits are used at their own positions on the readout lattice
(pad pitch in X and Y, one frame in Z), without voxelising them again.
"""
function itrack_components(h::SavedHits,tc::AbstractDict,radius::Real)
    keep=findall(>(Float32(tc["ion_hit_cut_keV"])),h.energy_keV);n=length(keep)
    n==0&&return (voxels=Int32[],energy_keV=Float32[],extent_mm=Float32[])
    x=Float64.(h.x_mm[keep]);y=Float64.(h.y_mm[keep]);z=Float64.(h.z_mm[keep]);e=h.energy_keV[keep]
    groups=Dict{Int,Vector{Int}}();for (i,l) in enumerate(point_components(x,y,z,radius));push!(get!(groups,l,Int[]),i);end
    members=sort(collect(values(groups));by=g->(-length(g),first(g)))
    extent(g)=Float32(max(maximum(x[g])-minimum(x[g]),maximum(y[g])-minimum(y[g]),maximum(z[g])-minimum(z[g])))
    (voxels=Int32[length(g) for g in members],energy_keV=Float32[sum(e[g]) for g in members],extent_mm=Float32[extent(g) for g in members])
end

function prepare_second_pass(run_dir::AbstractString,radii::Vector{Float64}=[3.5,4.6,5.5,6.7,8.0])
    run_dir=abspath(run_dir);batch_root=joinpath(run_dir,"batches");tags=sort(filter(t->occursin(r"^b\d{5}$",t),readdir(batch_root)))
    subdir=batch_output_subdir(joinpath(batch_root,first(tags)));out=joinpath(run_dir,subdir,"second_pass");mkpath(out);target=joinpath(out,"components_lattice.h5");isfile(target)&&error("Preserve existing second pass: $target")
    started=time();ids=Int64[];data=[(offsets=Int64[0],voxels=Int32[],energy=Float32[],extent=Float32[]) for _ in radii]
    for t in tags
        r=open_hit_input(joinpath(batch_root,t,subdir,"hits.h5"))
        try
            tc=TOML.parse(read(r.file["configuration/cards/resolved"],String))["run"]["track_analysis"]
            for id in read(r.file["events/event_number"])
                h=read_hit_event(r,id).itrack;push!(ids,id)
                for (k,radius) in enumerate(radii)
                    c=itrack_components(h,tc,radius);d=data[k];append!(d.voxels,c.voxels);append!(d.energy,c.energy_keV);append!(d.extent,c.extent_mm);push!(d.offsets,length(d.voxels))
                end
            end
        finally
            close_hit_input(r)
        end
        println("batch $t: $(length(ids)) events, $(round(time()-started;digits=1)) s");flush(stdout)
    end
    tmp=target*".tmp"
    h5open(tmp,"w") do f
        f["schema"]="itaca.second_pass_components.v1";f["event_number"]=ids;f["radii_mm"]=radii
        for (k,d) in enumerate(data);g="radius_$(k)";f["$g/offsets"]=d.offsets;f["$g/voxels"]=d.voxels;f["$g/energy_keV"]=d.energy;f["$g/extent_mm"]=d.extent;end
    end
    Base.rename(tmp,target);println("SECOND PASS $target events=$(length(ids)) time_s=$(round(time()-started;digits=1))");target
end

if abspath(PROGRAM_FILE)==@__FILE__
    1<=length(ARGS)<=2||error("Usage: prepare_second_pass.jl results/readout/<run_label> [RADII_MM]")
    length(ARGS)==2 ? prepare_second_pass(ARGS[1],parse.(Float64,split(ARGS[2],","))) : prepare_second_pass(ARGS[1])
end
