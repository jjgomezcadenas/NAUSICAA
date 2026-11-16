#!/usr/bin/env julia
# Track selection on the readout lattice for every accepted event of a run.
#   julia --project=. scripts/production/prepare_track_selection.jl results/readout/<run_label>
# Reads the hit file of each batch and applies `itrack_selection`
# (prepare_batch_features.jl) with the settings of the run card, or their
# defaults when the card predates them. Writes
# results/readout/<run_label>/<output_subdir>/track_selection/events.h5 with one
# row per event, in the order of the run statistics. Feature tables written
# since 2026-11-16 already hold these columns.
include("prepare_batch_features.jl")

function prepare_track_selection(run_dir::AbstractString)
    run_dir=abspath(run_dir);batch_root=joinpath(run_dir,"batches");tags=sort(filter(t->occursin(r"^b\d{5}$",t),readdir(batch_root)))
    subdir=batch_output_subdir(joinpath(batch_root,first(tags)));out=joinpath(run_dir,subdir,"track_selection");mkpath(out);target=joinpath(out,"events.h5");isfile(target)&&error("Preserve existing selection: $target")
    ids=Int64[];rows=NamedTuple[];settings=nothing
    for t in tags
        r=open_hit_input(joinpath(batch_root,t,subdir,"hits.h5"))
        try
            tc=TOML.parse(read(r.file["configuration/cards/resolved"],String))["run"]["track_analysis"]
            settings=Dict("ion_segment_connection_mm"=>Float64(get(tc,"ion_segment_connection_mm",5.5)),"ion_segment_min_hits"=>Int(get(tc,"ion_segment_min_hits",6)),
                "ion_deposit_min_keV"=>Float64(get(tc,"ion_deposit_min_keV",10.0)),"ion_max_segments"=>Int(get(tc,"ion_max_segments",3)),"ion_connection_mm"=>Float64(tc["ion_connection_mm"]),"ion_hit_cut_keV"=>Float64(tc["ion_hit_cut_keV"]))
            for id in read(r.file["events/event_number"]);push!(ids,id);push!(rows,itrack_selection(read_hit_event(r,id).itrack,tc));end
        finally
            close_hit_input(r)
        end
    end
    tmp=target*".tmp"
    h5open(tmp,"w") do f
        f["schema"]="itaca.track_selection.v1";f["event_number"]=ids
        for (k,v) in settings;f["settings/$k"]=v;end
        f["segments"]=[r.segments for r in rows];f["deposits"]=[r.deposits for r in rows];f["deposit_keV"]=[r.deposit_keV for r in rows]
        f["segments_joined"]=Int8[r.segments_joined for r in rows];f["selected"]=Int8[r.selected for r in rows]
    end
    Base.rename(tmp,target);println("TRACK SELECTION $target events=$(length(ids)) selected=$(count(r->r.selected,rows))");target
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1||error("Usage: prepare_track_selection.jl results/readout/<run_label>")
    prepare_track_selection(ARGS[1])
end
