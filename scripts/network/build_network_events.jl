#!/usr/bin/env julia
# Training events of one run for the classification network.
#   julia --project=. scripts/network/build_network_events.jl results/readout/<run_label> [OUTPUT.h5]
# For each batch it reads <output_subdir>/features.h5 and hits.h5, keeps the
# events with one component in the itrack and one in the etrack, and places the
# saved hits of each track on the readout lattice inside the fixed frame.
# Writes results/network/<run_label>/events.h5: for each track the voxel
# indices and energies of every event (sparse, one hit per cell), and one row
# per event with its
# features. train_network.py fills the dense grids from this file.
include("../production/prepare_batch_features.jl")

const NETWORK_EVENTS_SCHEMA="itaca.network_events.v1"

"""Number of lattice cells along one axis: cells whose centre lies within half a pitch of the frame."""
lattice_cells(side_mm::Real,pitch_mm::Real)=floor(Int,side_mm/pitch_mm)+2

"""Cell index (from 0) of hit coordinates `x` on a lattice of the given pitch.

`lower` is the lower face of the frame. Cell 0 is centred within half a pitch
of it. The lattice of an event is the one most of its hits lie on. An error is
raised when a coordinate is off that lattice by more than `tolerance` (in
units of the pitch), or when a cell falls outside 0:cells-1. The electron Z
axis needs a tolerance of half a pitch: the last time group of an event is
shorter than the others and its hits are centred on its actual span.
"""
function lattice_indices(x::AbstractVector{<:Real},lower::Real,pitch::Real,cells::Int;tolerance::Real=0.02)
    out=Vector{UInt8}(undef,length(x));isempty(x)&&return out
    off(a,b)=(s=(Float64(a)-Float64(b))/pitch;abs(s-round(s)))
    candidates=unique(round.(mod.(Float64.(x),pitch)./pitch;digits=2));anchors=[x[findfirst(v->round(mod(Float64(v),pitch)/pitch;digits=2)==c,x)] for c in candidates]
    x0=Float64(anchors[argmax([count(v->off(v,a)<0.02,x) for a in anchors])]);first_cell=floor(Int,(x0-lower+pitch/2)/pitch)
    for i in eachindex(x)
        steps=(Float64(x[i])-x0)/pitch;k=round(Int,steps)
        abs(steps-k)<tolerance||error("Hit at $(x[i]) mm is off the lattice of pitch $pitch mm")
        cell=first_cell+k;0<=cell<cells||error("Cell $cell outside 0:$(cells-1)")
        out[i]=cell
    end
    out
end

"""Positions of the hits to keep so that each cell holds one hit: of several hits in a cell, the one with most energy.

Two hits share a cell when a noise hit was drawn for a pad and frame that also holds a signal hit.
"""
function one_hit_per_cell(ix,iy,iz,energy)
    best=Dict{NTuple{3,UInt8},Int}()
    for i in eachindex(energy)
        k=(ix[i],iy[i],iz[i]);j=get(best,k,0);(j==0||energy[i]>energy[j])&&(best[k]=i)
    end
    sort!(collect(values(best)))
end

"""Pitches (x, y, z in mm) of the two readout lattices, from the resolved card of a hit file."""
function lattice_pitches(resolved::AbstractDict)
    d=resolved["detector"]
    (itrack=(Float64(d["ions"]["pad_pitch_mm"]),Float64(d["ions"]["pad_pitch_mm"]),d["ions"]["ion_bucket_ms"]*d["drift"]["ion_drift_velocity_cm_s"]/100),
     etrack=(Float64(d["el"]["hole_pitch_mm"]),Float64(d["el"]["hole_pitch_mm"]),d["sipm"]["electron_bucket_us"]*resolved["run"]["electrons"]["rebin_factor"]*d["drift"]["electron_drift_velocity_mm_us"]))
end

function build_network_events(run_dir::AbstractString,target::AbstractString="")
    run_dir=abspath(run_dir);batch_root=joinpath(run_dir,"batches")
    tags=sort(filter(t->occursin(r"^b\d{5}$",t),readdir(batch_root)));isempty(tags)&&error("No batch folders under $batch_root")
    isempty(target)&&(target=joinpath(ROOT,"results/network",basename(run_dir),"events.h5"));target=abspath(target)
    isfile(target)&&error("Preserve existing events: $target");mkpath(dirname(target))
    started=time();columns=Dict{String,Vector}();accepted=0;pitches=nothing;side=0.0;first_table=nothing;shas=String[]
    arms=Dict(a=>(offsets=Int64[0],ix=UInt8[],iy=UInt8[],iz=UInt8[],energy=Float32[]) for a in ("itrack","etrack"));cells=Dict{String,NTuple{3,Int}}();dropped=Dict("itrack"=>0,"etrack"=>0)
    for t in tags
        dir=joinpath(batch_root,t,batch_output_subdir(joinpath(batch_root,t)));feat=read_batch_features(joinpath(dir,"features.h5"));hits=joinpath(dir,"hits.h5")
        bytes2hex(open(sha256,hits))==feat["hit_file_sha256"]||error("Feature table and hit file of $t do not match")
        first_table===nothing ? (first_table=feat) : for k in ("track_analysis","detector","ion_hit_cut_keV","electron_hit_cut_keV","species","run_label","settings")
            feat[k]==first_table[k]||error("Setting $k differs in batch $t")
        end
        e=feat["events"];keep=findall((e["itrack_components"].==1).&(e["etrack_components"].==1));accepted+=length(e["event_number"]);push!(shas,feat["hit_file_sha256"])
        h5open(hits,"r") do f
            read(f["events/event_number"])==e["event_number"]||error("Event order differs between features and hits in $t")
            thispitch=lattice_pitches(TOML.parse(read(f["configuration/cards/resolved"],String)))
            lower=[read(f["events/fixed_frame_$(c)_min_mm"]) for c in "xyz"];upper=read(f["events/fixed_frame_x_max_mm"])
            thisside=Float64(upper[1]-lower[1][1]);all(i->isapprox(upper[i]-lower[1][i],thisside;atol=1f-2),eachindex(upper))||error("Frame side varies in $t")
            pitches===nothing ? (pitches=thispitch;side=thisside) : (pitches==thispitch&&side==thisside||error("Lattice or frame differs in $t"))
            for arm in ("itrack","etrack")
                p=getfield(pitches,Symbol(arm));n=ntuple(i->lattice_cells(side,p[i]),3);all(<=(256),n)||error("Grid too large for 8-bit indices: $n");cells[arm]=n
                off=read(f["$arm/offsets"]);xyz=[read(f["$arm/$(c)_mm"]) for c in "xyz"];energy=read(f["$arm/energy_keV"]);a=arms[arm]
                for row in keep
                    r=off[row]+1:off[row+1];idx=[lattice_indices(view(xyz[i],r),lower[i][row],p[i],n[i];tolerance=(arm=="etrack"&&i==3) ? 0.5 : 0.02) for i in 1:3]
                    en=view(energy,r);kept=one_hit_per_cell(idx[1],idx[2],idx[3],en);dropped[arm]+=length(en)-length(kept)
                    append!(a.ix,idx[1][kept]);append!(a.iy,idx[2][kept]);append!(a.iz,idx[3][kept]);append!(a.energy,en[kept]);push!(a.offsets,length(a.energy))
                end
            end
        end
        for (k,v) in e;append!(get!(columns,k,similar(v,0)),v[keep]);end
        println("batch $t: $(length(columns["event_number"])) of $accepted events kept, $(round(time()-started;digits=1)) s");flush(stdout)
    end
    allunique(columns["event_number"])||error("Duplicate event numbers across batches")
    budget=16*1024^2+2*(Base.summarysize(columns)+sum(Base.summarysize(a) for a in values(arms)));GuardedFiles.space_check(dirname(target),budget,1024^3);tmp=target*".tmp"
    try
        h5open(tmp,"w") do f
            f["schema"]=NETWORK_EVENTS_SCHEMA;f["species"]=first_table["species"];f["run_label"]=first_table["run_label"];f["batches"]=tags;f["hit_file_sha256"]=shas
            f["selection"]="one component in the itrack and one in the etrack";f["accepted_events"]=accepted;f["frame_side_mm"]=side
            f["track_analysis"]=first_table["track_analysis"];f["detector"]=first_table["detector"]
            f["ion_hit_cut_keV"]=first_table["ion_hit_cut_keV"];f["electron_hit_cut_keV"]=first_table["electron_hit_cut_keV"]
            for (k,v) in columns;f["events/$k"]=v;end
            for (arm,a) in arms
                f["$arm/pitch_mm"]=collect(getfield(pitches,Symbol(arm)));f["$arm/cells"]=collect(cells[arm]);f["$arm/offsets"]=a.offsets;f["$arm/hits_dropped_from_shared_cells"]=dropped[arm]
                for (k,v) in (("ix",a.ix),("iy",a.iy),("iz",a.iz),("energy_keV",a.energy));f["$arm/$k",chunk=(min(length(v),65536),),shuffle=true,deflate=4]=v;end
            end
            f["time_s"]=time()-started;flush(f);GuardedFiles.check_size(tmp,budget)
        end
        Base.rename(tmp,target)
    catch
        isfile(tmp)&&rm(tmp);rethrow()
    end
    println("NETWORK EVENTS $target hits dropped from shared cells: $dropped events=$(length(columns["event_number"])) of $accepted bytes=$(filesize(target)) time_s=$(round(time()-started;digits=1))")
    target
end

if abspath(PROGRAM_FILE)==@__FILE__
    1<=length(ARGS)<=2||error("Usage: build_network_events.jl results/readout/<run_label> [OUTPUT.h5]")
    build_network_events(ARGS...)
end
