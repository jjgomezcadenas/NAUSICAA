#!/usr/bin/env julia
# Statistics of all accepted events of a run kept as per-batch hit files.
#   julia --project=. scripts/production/prepare_batch_statistics.jl results/readout/<run_label>
# Reads results/readout/<run_label>/batches/*/<output_subdir>/hits.h5 and the
# selection counts of each batch; writes
# results/readout/<run_label>/<output_subdir>/statistics/events.h5 with one row
# per accepted event, the energy fits and the selection counts. No detector or
# electronics calculation is repeated. present_batch_statistics.py makes the page.
include("prepare_hit_presentation.jl")

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

"""Gaussian peak fits of a quantity in slices of another; NaN where a slice has too few events."""
function sliced_fits(values,by,edges;minimum_events=200)
    n=length(edges)-1;out=Dict("edges"=>collect(Float64.(edges)),"entries"=>zeros(Int64,n),"mu"=>fill(NaN,n),"sigma"=>fill(NaN,n),"mu_error"=>fill(NaN,n),"sigma_error"=>fill(NaN,n))
    for k in 1:n
        x=[v for (v,b) in zip(values,by) if isfinite(v)&&edges[k]<=b<edges[k+1]];out["entries"][k]=length(x)
        length(x)<minimum_events&&continue
        fit=fit_peak(x);out["mu"][k]=fit.mu;out["sigma"][k]=fit.sigma;out["mu_error"][k]=fit.mu_error;out["sigma_error"][k]=fit.sigma_error
    end
    out
end

function prepare_batch_statistics(run_dir::AbstractString)
    run_dir=abspath(run_dir);batch_root=joinpath(run_dir,"batches");isdir(batch_root)||error("No batches under $run_dir")
    tags=sort(filter(t->occursin(r"^b\d{5}$",t),readdir(batch_root)));isempty(tags)&&error("No batch folders under $batch_root")
    subdirs=unique(vcat([filter(d->isfile(joinpath(batch_root,t,d,"hits.h5")),readdir(joinpath(batch_root,t))) for t in tags]...))
    length(subdirs)==1||error("Expected one hit output folder per batch, found $subdirs")
    subdir=only(subdirs);files=[joinpath(batch_root,t,subdir,"hits.h5") for t in tags]
    all(isfile,files)||error("Missing hit files: $(tags[.!isfile.(files)])")
    out=joinpath(run_dir,subdir,"statistics");mkpath(out);target=joinpath(out,"events.h5");isfile(target)&&error("Preserve existing statistics: $target")
    started=time();rows=NamedTuple[];species="";label="";tc=nothing;dc=nothing;cuts=nothing
    for (b,path) in enumerate(files)
        r=open_hit_input(path)
        try
            read(r.file["schema"],String)==SignalReconstruction.HIT_SCHEMA_LOCATED||error("Expected hit files with event locations: $path")
            resolved=TOML.parse(read(r.file["configuration/cards/resolved"],String));d=resolved["detector"]
            thistc=resolved["run"]["track_analysis"];thisdc=Dict{String,Any}("chamber_length_mm"=>10d["geometry"]["fiducial_length_cm"],"detector_radius_mm"=>5d["geometry"]["fiducial_diameter_cm"],
                "ion_pitch_mm"=>d["ions"]["pad_pitch_mm"],"EL_pitch_mm"=>d["el"]["hole_pitch_mm"],"ion_bucket_ms"=>d["ions"]["ion_bucket_ms"],
                "electron_bucket_us"=>d["sipm"]["electron_bucket_us"]*resolved["run"]["electrons"]["rebin_factor"])
            thiscuts=(read(r.file["configuration/ion_hit_cut_keV"]),read(r.file["configuration/electron_hit_cut_keV"]))
            b==1 ? (tc=thistc;dc=thisdc;cuts=thiscuts;species=read(r.file["configuration/species"],String);label=dirname(dirname(read(r.file["configuration/run_label"],String)))) :
                (tc==thistc&&dc==thisdc&&cuts==thiscuts||error("Settings differ in $path"))
            for id in read(r.file["events/event_number"])
                e=read_hit_event(r,id);where=read_hit_location(r,id);f=e.info.frame;c=where.fixed_frame
                ion=track_features(e.itrack,tc,"itrack");ele=track_features(e.etrack,tc,"etrack")
                push!(rows,merge((event_number=id,batch=b-1,ion_energy_keV=e.energy.ion_energy_keV,electron_energy_keV=e.energy.electron_energy_keV,
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
        println("batch $(tags[b]): $(length(rows)) events so far, $(round(time()-started;digits=1)) s");flush(stdout)
    end
    allunique(r.event_number for r in rows)||error("Duplicate event numbers across batches")
    # Selection counts of each batch, from its detector-response controls
    selection=Dict{String,Vector{Int64}}();settings=Dict{String,Any}()
    for t in tags
        h5open(joinpath(batch_root,t,"detector_response/controls.h5"),"r") do f
            for k in keys(f["selection"])
                f["selection"][k] isa HDF5.Dataset||continue
                v=read(f["selection"][k])
                v isa Integer ? push!(get!(selection,k,Int64[]),v) : (get!(settings,k,v)==v||error("Selection setting $k differs in batch $t"))
            end
        end
    end
    budget=16*1024^2+Base.summarysize(rows)*2;GuardedFiles.space_check(out,budget,1024^3);tmp=target*".tmp"
    try
        h5open(tmp,"w") do f
            f["schema"]="itaca.batch_statistics.v1";f["species"]=species;f["run_label"]=label;f["batches"]=tags
            f["hit_files"]=relpath.(files,out);f["hit_file_sha256"]=[bytes2hex(open(sha256,p)) for p in files]
            f["track_analysis"]=sprint(io->TOML.print(io,tc;sorted=true));f["detector"]=sprint(io->TOML.print(io,dc;sorted=true))
            f["ion_hit_cut_keV"]=cuts[1];f["electron_hit_cut_keV"]=cuts[2]
            for k in keys(first(rows));f["events/$k"]=[getfield(r,k) for r in rows];end
            for (k,v) in selection;f["selection/batch/$k"]=v;f["selection/total/$k"]=sum(v);end
            for (k,v) in settings;f["selection/settings/$k"]=v;end
            drift_edges=range(0,dc["chamber_length_mm"]/10;length=7)
            for (key,drift) in (("ion_energy_keV","ion_drift_cm"),("electron_energy_keV","electron_drift_cm"))
                values=Float64[getfield(r,Symbol(key)) for r in rows];fit=fit_peak(filter(isfinite,values))
                for (k,v) in Base.pairs(fit);f["energy_fits/$key/$k"]=v;end
                for (k,v) in sliced_fits(values,[getfield(r,Symbol(drift)) for r in rows],drift_edges);f["energy_fits_by_drift/$key/$k"]=v;end
            end
            f["time_s"]=time()-started;flush(f);GuardedFiles.check_size(tmp,budget)
        end
        Base.rename(tmp,target)
    catch
        isfile(tmp)&&rm(tmp);rethrow()
    end
    println("STATISTICS $target events=$(length(rows)) batches=$(length(tags)) time_s=$(round(time()-started;digits=1))")
    target
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1||error("Usage: prepare_batch_statistics.jl results/readout/<run_label>")
    prepare_batch_statistics(ARGS[1])
end
