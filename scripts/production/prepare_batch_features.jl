#!/usr/bin/env julia
# Track features of the accepted events of one batch.
#   julia --project=. scripts/production/prepare_batch_features.jl results/readout/<run_label>/batches/b0000N
# Reads <output_subdir>/hits.h5 and detector_response/controls.h5 of the batch;
# writes <output_subdir>/features.h5 with one row per accepted event, the
# settings and the selection counts of the batch. The table is enough to
# select events and to build the statistics of a run without its hit files.
include("prepare_hit_presentation.jl")

const FEATURE_SCHEMA="itaca.batch_features.v1"

"""Graph features of one track of one event, from its retained hits.

Hits above the display cut are voxelised and connected as in the event
presentation (`[track_analysis]` of the run card). The main track is the
component with most voxels. Its length is the length of the tree path between
its two endpoints; its extensions are the ranges of its voxel centres plus one
voxel width. E1 >= E2 are the energies in spheres around the two endpoints.
"""
function track_features(h::SavedHits,tc::AbstractDict,arm::AbstractString)
    prefix=arm=="itrack" ? "ion" : "electron"
    cut=Float32(tc["$(prefix)_hit_cut_keV"]);width=Float32(tc["$(prefix)_voxel_width_mm"])
    radius=Float32(tc["$(prefix)_connection_mm"]);blob=Float64(tc["$(prefix)_blob_radius_mm"])
    keep=h.energy_keV .> cut
    nothing_found=(hits=count(keep),components=0,main_keV=0f0,second_keV=0f0,main_voxels=0,length_mm=NaN32,
        extent_x_mm=NaN32,extent_y_mm=NaN32,extent_z_mm=NaN32,E1_keV=NaN32,E2_keV=NaN32)
    any(keep)||return nothing_found
    df=EventAnalysis.DataFrame(x=h.x_mm[keep],y=h.y_mm[keep],z=h.z_mm[keep],energy=h.energy_keV[keep]./1000f0)
    tracks=make_tracks(voxelize_event(df,width);method="KDT",max_distance_mm=radius,energy_threshold_kev=0f0)
    isempty(tracks)&&return nothing_found
    energies=sort(Float32[1000sum(t.voxels.energy) for t in tracks];rev=true)
    track=tracks[argmax([size(t.voxels,1) for t in tracks])];v=track.voxels
    u,w,path,_,_=EventAnalysis.find_extremes_mst_diameter(track,EventAnalysis.extract_coords(track))
    len=sum(sqrt((v.x[path[i+1]]-v.x[path[i]])^2+(v.y[path[i+1]]-v.y[path[i]])^2+(v.z[path[i+1]]-v.z[path[i]])^2) for i in 1:length(path)-1;init=0.0)
    b=EventAnalysis.energy_in_spheres_around_extremes(track,(extremes=(v[u,:],v[w,:]),),blob)
    e1,e2=1000 .*(max(b.blob1_energy,b.blob2_energy),min(b.blob1_energy,b.blob2_energy))
    extent(c)=Float32(maximum(c)-minimum(c)+width)
    (hits=count(keep),components=length(tracks),main_keV=1000f0*Float32(sum(v.energy)),second_keV=length(energies)<2 ? 0f0 : energies[2],
        main_voxels=size(v,1),length_mm=Float32(len),extent_x_mm=extent(v.x),extent_y_mm=extent(v.y),extent_z_mm=extent(v.z),
        E1_keV=Float32(e1),E2_keV=Float32(e2))
end

batch_number(tag::AbstractString)=parse(Int,tag[2:end])

"""Folder of a batch that holds its hit file or its feature table."""
function batch_output_subdir(batch_dir::AbstractString)
    found=filter(d->isdir(joinpath(batch_dir,d))&&any(f->isfile(joinpath(batch_dir,d,f)),("hits.h5","features.h5")),readdir(batch_dir))
    length(found)==1||error("Expected one hit output folder in $batch_dir, found $found")
    only(found)
end

"""Feature table of one batch, computed from its hit file and selection controls.

Returns a dictionary: `events` (column name => vector, one entry per accepted
event), `selection` (counts), `settings` (selection settings), the settings as
TOML text (`track_analysis`, `detector`), the two saved hit cuts, `species`,
`run_label`, `batch` and the SHA-256 of the hit file.
"""
function compute_batch_features(batch_dir::AbstractString)
    batch_dir=abspath(batch_dir);tag=basename(batch_dir);b=batch_number(tag);path=joinpath(batch_dir,batch_output_subdir(batch_dir),"hits.h5")
    rows=NamedTuple[];r=open_hit_input(path);out=Dict{String,Any}("batch"=>tag,"hit_file_sha256"=>bytes2hex(open(sha256,path)))
    try
        read(r.file["schema"],String)==SignalReconstruction.HIT_SCHEMA_LOCATED||error("Expected hit files with event locations: $path")
        resolved=TOML.parse(read(r.file["configuration/cards/resolved"],String));d=resolved["detector"];tc=resolved["run"]["track_analysis"]
        dc=Dict{String,Any}("chamber_length_mm"=>10d["geometry"]["fiducial_length_cm"],"detector_radius_mm"=>5d["geometry"]["fiducial_diameter_cm"],
            "ion_pitch_mm"=>d["ions"]["pad_pitch_mm"],"EL_pitch_mm"=>d["el"]["hole_pitch_mm"],"ion_bucket_ms"=>d["ions"]["ion_bucket_ms"],
            "electron_bucket_us"=>d["sipm"]["electron_bucket_us"]*resolved["run"]["electrons"]["rebin_factor"])
        out["track_analysis"]=sprint(io->TOML.print(io,tc;sorted=true));out["detector"]=sprint(io->TOML.print(io,dc;sorted=true))
        out["ion_hit_cut_keV"]=read(r.file["configuration/ion_hit_cut_keV"]);out["electron_hit_cut_keV"]=read(r.file["configuration/electron_hit_cut_keV"])
        out["species"]=read(r.file["configuration/species"],String);out["run_label"]=dirname(dirname(read(r.file["configuration/run_label"],String)))
        for id in read(r.file["events/event_number"])
            e=read_hit_event(r,id);where=read_hit_location(r,id);f=e.info.frame;c=where.fixed_frame
            ion=track_features(e.itrack,tc,"itrack");ele=track_features(e.etrack,tc,"etrack")
            push!(rows,merge((event_number=id,batch=b,ion_energy_keV=e.energy.ion_energy_keV,electron_energy_keV=e.energy.electron_energy_keV,
                invalid_ion_buckets=e.energy.invalid_ion_buckets,invalid_electron_buckets=e.energy.invalid_electron_buckets,
                ion_drift_cm=e.info.ion_drift_cm,electron_drift_cm=e.info.electron_drift_cm,
                vertex_x_mm=where.vertex_mm[1],vertex_y_mm=where.vertex_mm[2],vertex_z_mm=where.vertex_mm[3],
                frame_side_mm=c.upper_mm[1]-c.lower_mm[1],response_box_side_mm=maximum(f.upper_mm[i]-f.lower_mm[i] for i in 1:3),
                frame_shift_mm=maximum(abs((f.upper_mm[i]+f.lower_mm[i])/2-(c.upper_mm[i]+c.lower_mm[i])/2) for i in 1:3),
                itrack_saved_hits=length(e.itrack.x_mm),etrack_saved_hits=length(e.etrack.x_mm)),
                NamedTuple{Tuple(Symbol("itrack_",k) for k in keys(ion))}(values(ion)),NamedTuple{Tuple(Symbol("etrack_",k) for k in keys(ele))}(values(ele))))
        end
    finally
        close_hit_input(r)
    end
    isempty(rows)&&error("No accepted events in $path")
    out["events"]=Dict{String,Any}(String(k)=>[getfield(row,k) for row in rows] for k in keys(first(rows)))
    selection=Dict{String,Int64}();settings=Dict{String,Any}()
    h5open(joinpath(batch_dir,"detector_response/controls.h5"),"r") do f
        for k in keys(f["selection"])
            f["selection"][k] isa HDF5.Dataset||continue
            v=read(f["selection"][k]);v isa Integer ? (selection[k]=v) : (settings[k]=v)
        end
    end
    out["selection"]=selection;out["settings"]=settings
    out
end

"""Write the feature table of a batch beside its hit file; an existing table is preserved."""
function write_batch_features(batch_dir::AbstractString)
    batch_dir=abspath(batch_dir);target=joinpath(batch_dir,batch_output_subdir(batch_dir),"features.h5");isfile(target)&&error("Preserve existing features: $target")
    started=time();t=compute_batch_features(batch_dir);budget=16*1024^2+2Base.summarysize(t["events"]);GuardedFiles.space_check(dirname(target),budget,1024^3);tmp=target*".tmp"
    try
        h5open(tmp,"w") do f
            f["schema"]=FEATURE_SCHEMA
            for k in ("batch","hit_file_sha256","track_analysis","detector","ion_hit_cut_keV","electron_hit_cut_keV","species","run_label");f[k]=t[k];end
            for (k,v) in t["events"];f["events/$k"]=v;end
            for (k,v) in t["selection"];f["selection/counts/$k"]=v;end
            for (k,v) in t["settings"];f["selection/settings/$k"]=v;end
            f["time_s"]=time()-started;flush(f);GuardedFiles.check_size(tmp,budget)
        end
        Base.rename(tmp,target)
    catch
        isfile(tmp)&&rm(tmp);rethrow()
    end
    println("FEATURES $target events=$(length(t["events"]["event_number"])) time_s=$(round(time()-started;digits=1))")
    target
end

"""Read a feature table written by `write_batch_features`, in the form of `compute_batch_features`."""
function read_batch_features(path::AbstractString)
    h5open(path,"r") do f
        read(f["schema"])==FEATURE_SCHEMA||error("Unexpected schema in $path")
        t=Dict{String,Any}(k=>read(f[k]) for k in ("batch","hit_file_sha256","track_analysis","detector","ion_hit_cut_keV","electron_hit_cut_keV","species","run_label"))
        t["events"]=Dict{String,Any}(k=>read(f["events"][k]) for k in keys(f["events"]))
        t["selection"]=Dict{String,Int64}(k=>read(f["selection/counts"][k]) for k in keys(f["selection/counts"]))
        t["settings"]=Dict{String,Any}(k=>read(f["selection/settings"][k]) for k in keys(f["selection/settings"]))
        t
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1||error("Usage: prepare_batch_features.jl results/readout/<run_label>/batches/b0000N")
    write_batch_features(ARGS[1])
end
