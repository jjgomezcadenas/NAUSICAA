using HDF5
using ..EventAnalysis: EventFrame, mc_vertex_frame
include("../ion_signal/guarded_files.jl")

"""Four retained-hit columns in vertex-relative mm and keV; one entry per hit."""
struct SavedHits
    x_mm::Vector{Float32}
    y_mm::Vector{Float32}
    z_mm::Vector{Float32}
    energy_keV::Vector{Float32}
end
"""Signed event totals before hit cuts, stored at Float32 precision.
`corrected_pe` includes mean attachment correction. Invalid-bucket counts
refer to the reconstructed vectors used for the respective totals.
"""
struct EventEnergy
    ions::Float32
    corrected_pe::Float32
    ion_energy_keV::Float32
    electron_energy_keV::Float32
    invalid_ion_buckets::Int64
    invalid_electron_buckets::Int64
end
"""MC event number, complementary drift distances, acquisition frame and
transverse diffusion widths. Widths reproduce subsequent graph distance cuts."""
struct HitEventInfo
    event_number::Int64
    ion_drift_cm::Float32
    electron_drift_cm::Float32
    frame::EventFrame
    ion_sigma_t_mm::Float32
    electron_sigma_t_mm::Float32
end
"""One event's information, pre-cut totals and both retained hit sets."""
struct HitEvent
    info::HitEventInfo
    energy::EventEnergy
    itrack::SavedHits
    etrack::SavedHits
end
"""Hit thresholds (keV), deflate level, hit chunk length and hard file-size cap.
The minimum permitted threshold is 1 keV. Retained energies exceed the threshold.
Compression levels 0–9 are accepted; level zero is useful for validation.
"""
struct HitOutputSettings
    ion_hit_cut_keV::Float32
    electron_hit_cut_keV::Float32
    compression_level::Int64
    chunk_hits::Int64
    max_output_bytes::Int64
    function HitOutputSettings(i,e,c,n,b)
        all(x->isfinite(x)&&x>=1,(i,e))||throw(ArgumentError("Hit cuts must be at least 1 keV"))
        0<=c<=9&&n>0&&b>0||throw(ArgumentError("Invalid storage settings"))
        new(i,e,c,n,b)
    end
end
"""Run label, species, gas, complete TOML card text and labelled SHA-256 checksums."""
struct HitFileMetadata
    run_label::String
    species::String
    gas_mixture::GasMixture
    cards::Dict{String,String}
    checksums::Dict{String,String}
end
"""Append-only temporary HDF5 output. Publish with `finish_hit_output`;
on exceptions call `abort_hit_output`. Completed files are never overwritten."""
mutable struct HitWriter
    file::HDF5.File
    temporary_path::String
    output_path::String
    settings::HitOutputSettings
    event_count::Int64
    ion_hit_count::Int64
    electron_hit_count::Int64
    event_numbers::Set{Int64}
end
"""Open hit file with a lookup from MC event number to one-based event row."""
struct HitReader
    file::HDF5.File
    event_rows::Dict{Int64,Int64}
end
const HIT_SCHEMA="itaca.reconstructed_hits.v1"
const HIT_COLUMNS=("x_mm","y_mm","z_mm","energy_keV")
const EVENT_FLOAT_COLUMNS=("ion_drift_cm","electron_drift_cm","frame_x_min_mm","frame_x_max_mm",
 "frame_y_min_mm","frame_y_max_mm","frame_z_min_mm","frame_z_max_mm","ion_sigma_t_mm",
 "electron_sigma_t_mm","ions","corrected_pe","ion_energy_keV","electron_energy_keV")
const EVENT_INT_COLUMNS=("event_number","invalid_ion_buckets","invalid_electron_buckets")

"""Extract the four saved columns without copying; the writer consumes them synchronously."""
saved_hits(h::ReconstructedHits)=SavedHits(h.x_mm,h.y_mm,h.z_mm,h.energy_keV)

function create_column(g,name,T,chunk,level)
    d=create_dataset(g,name,T,((0,),(-1,));chunk=(chunk,),shuffle=(),deflate=level)
    close(d)
end
function append_column(f,path,values)
    d=f[path]
    try
        n=length(d);m=length(values)
        m==0&&return
        HDF5.set_extent_dims(d,(n+m,));d[n+1:n+m]=values
    finally
        close(d)
    end
end
"""Create an extendible compressed hit file beside its final path.
Checks the output cap and a 1 GiB free-space reserve before creating a unique
run's `.tmp` file. Refuses an existing final file or unfinished temporary.
"""
function open_hit_output(path::AbstractString,metadata::HitFileMetadata,settings::HitOutputSettings)::HitWriter
    target=abspath(path);tmp=target*".tmp";mkpath(dirname(target))
    (ispath(target)||ispath(tmp))&&throw(ArgumentError("Output or temporary already exists: $target"))
    GuardedFiles.space_check(dirname(target),settings.max_output_bytes,1024^3)
    f=h5open(tmp,"w")
    try
        f["schema"]=HIT_SCHEMA
        g=create_group(f,"events")
        for key in EVENT_FLOAT_COLUMNS;create_column(g,key,Float32,256,settings.compression_level);end
        for key in EVENT_INT_COLUMNS;create_column(g,key,Int64,256,settings.compression_level);end
        close(g)
        for arm in ("itrack","etrack")
            g=create_group(f,arm)
            for key in HIT_COLUMNS;create_column(g,key,Float32,settings.chunk_hits,settings.compression_level);end
            create_column(g,"offsets",Int64,256,settings.compression_level);close(g)
            append_column(f,"$arm/offsets",Int64[0])
        end
        g=create_group(f,"configuration")
        for key in (:run_label,:species);g[string(key)]=getfield(metadata,key);end
        g["gas_mixture"]=string(metadata.gas_mixture);g["coordinate_frame"]="mc_vertex_frame"
        g["position_unit"]="mm";g["energy_unit"]="keV"
        for key in fieldnames(HitOutputSettings);g[string(key)]=getfield(settings,key);end
        for key in (:cards,:checksums)
            sub=create_group(g,string(key))
            for (name,value) in getfield(metadata,key);sub[name]=value;end
            close(sub)
        end
        close(g);flush(f);GuardedFiles.check_size(tmp,settings.max_output_bytes)
        HitWriter(f,tmp,target,settings,0,0,0,Set{Int64}())
    catch
        close(f);rm(tmp;force=true);rethrow()
    end
end

"""Append one event's totals and hits and update its zero-based cumulative offsets.
Validate columns, coordinate frame, threshold and duplicate event number before
writing. Check free space and the file budget after flushing each event.
On a write failure the caller must abort the temporary output.
"""
function write_hit_event(w::HitWriter,event::HitEvent)::Nothing
    info=event.info;energy=event.energy
    info.event_number in w.event_numbers&&throw(ArgumentError("Duplicate MC event number"))
    info.frame.frame==mc_vertex_frame||throw(ArgumentError("Expected vertex-relative frame"))
    for (h,cut) in ((event.itrack,w.settings.ion_hit_cut_keV),(event.etrack,w.settings.electron_hit_cut_keV))
        n=length(h.x_mm)
        all(k->length(getfield(h,Symbol(k)))==n,HIT_COLUMNS)||throw(ArgumentError("Unequal hit columns"))
        all(k->all(isfinite,getfield(h,Symbol(k))),HIT_COLUMNS)||throw(ArgumentError("Nonfinite hit"))
        all(>(cut),h.energy_keV)||throw(ArgumentError("Hit at or below saved cut"))
    end
    GuardedFiles.space_check(dirname(w.output_path),max(0,w.settings.max_output_bytes-filesize(w.temporary_path)),1024^3)
    values=Dict{String,Any}("event_number"=>info.event_number,"ion_drift_cm"=>info.ion_drift_cm,
      "electron_drift_cm"=>info.electron_drift_cm,"ion_sigma_t_mm"=>info.ion_sigma_t_mm,
      "electron_sigma_t_mm"=>info.electron_sigma_t_mm)
    for (i,axis) in enumerate(("x","y","z"))
        values["frame_$(axis)_min_mm"]=info.frame.lower_mm[i]
        values["frame_$(axis)_max_mm"]=info.frame.upper_mm[i]
    end
    for k in fieldnames(EventEnergy);values[string(k)]=getfield(energy,k);end
    for (k,v) in values;append_column(w.file,"events/$k",[v]);end
    for (arm,h) in (("itrack",event.itrack),("etrack",event.etrack))
        for key in HIT_COLUMNS;append_column(w.file,"$arm/$key",getfield(h,Symbol(key)));end
        count=(arm=="itrack" ? w.ion_hit_count : w.electron_hit_count)+length(h.x_mm)
        append_column(w.file,"$arm/offsets",[count])
    end
    w.event_count+=1;w.ion_hit_count+=length(event.itrack.x_mm);w.electron_hit_count+=length(event.etrack.x_mm)
    push!(w.event_numbers,info.event_number);flush(w.file)
    GuardedFiles.check_size(w.temporary_path,w.settings.max_output_bytes)
    nothing
end
"""Check schema lengths and offsets, close, and atomically publish the completed file."""
function finish_hit_output(w::HitWriter)::String
    for key in (EVENT_FLOAT_COLUMNS...,EVENT_INT_COLUMNS...)
        d=w.file["events/$key"];try length(d)==w.event_count||error("Incomplete event column") finally close(d) end
    end
    for (arm,count) in (("itrack",w.ion_hit_count),("etrack",w.electron_hit_count))
        off=read(w.file["$arm/offsets"])
        length(off)==w.event_count+1&&first(off)==0&&last(off)==count&&issorted(off)||error("Invalid offsets")
        for key in HIT_COLUMNS
            d=w.file["$arm/$key"];try length(d)==count||error("Incomplete hit column") finally close(d) end
        end
    end
    flush(w.file);GuardedFiles.check_size(w.temporary_path,w.settings.max_output_bytes);close(w.file)
    ispath(w.output_path)&&error("Output appeared during publication")
    Base.rename(w.temporary_path,w.output_path)
    w.output_path
end
"""Close the writer and remove its unpublished temporary file; preserve completed output."""
function abort_hit_output(w::HitWriter)::Nothing
    isopen(w.file)&&close(w.file)
    isfile(w.temporary_path)&&rm(w.temporary_path)
    nothing
end
"""Open the supported schema and index MC event numbers. Close with `close_hit_input`."""
function open_hit_input(path::AbstractString)::HitReader
    f=h5open(path,"r")
    try
        read(f["schema"],String)==HIT_SCHEMA||throw(ArgumentError("Unsupported hit schema"))
        ids=read(f["events/event_number"]);length(unique(ids))==length(ids)||error("Duplicate event numbers")
        HitReader(f,Dict(id=>Int64(i) for (i,id) in enumerate(ids)))
    catch;close(f);rethrow();end
end
"""Read pre-cut totals and retained hits by MC event number, not file ordinal."""
function read_hit_event(r::HitReader,event_number::Int64)::HitEvent
    row=get(r.event_rows,event_number,0);row>0||throw(KeyError(event_number))
    value(key)=begin d=r.file["events/$key"];try d[row] finally close(d) end end
    frame=EventFrame(mc_vertex_frame,ntuple(i->value("frame_$("xyz"[i])_min_mm"),3),ntuple(i->value("frame_$("xyz"[i])_max_mm"),3))
    info=HitEventInfo(event_number,value("ion_drift_cm"),value("electron_drift_cm"),frame,value("ion_sigma_t_mm"),value("electron_sigma_t_mm"))
    energy=EventEnergy((value(string(k)) for k in fieldnames(EventEnergy))...)
    function hits(arm)
        d=r.file["$arm/offsets"];off=try d[row:row+1] finally close(d) end
        columns=map(HIT_COLUMNS) do key
            d=r.file["$arm/$key"];try d[off[1]+1:off[2]] finally close(d) end
        end
        SavedHits(columns...)
    end
    HitEvent(info,energy,hits("itrack"),hits("etrack"))
end
"""Release an open hit reader."""
close_hit_input(r::HitReader)=close(r.file)
