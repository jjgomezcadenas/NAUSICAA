# Chain the hit files of the batches of a run into one file.
using HDF5,SHA,TOML
include("../../src/ion_signal/guarded_files.jl")

"""Chain per-batch hit files into `target`.

Event and hit columns are concatenated in batch order and the hit offsets
continued. All batches must share the schema, the hit cuts and every checksum
except those of their own inputs (`source`, `selection`); event numbers must be
unique. The configuration of the first batch is kept, and a group `batches`
records each file's run label, input checksums and event count. The output is
bounded by the summed input sizes plus an allowance, written to a temporary
file and renamed.
"""
function merge_hit_batches(files::Vector{String},target::AbstractString;compression_level::Int=4,chunk_hits::Int=65536)::String
    isempty(files)&&throw(ArgumentError("No hit files to chain"))
    ispath(target)||ispath(target*".tmp") ? throw(ArgumentError("Output or temporary already exists: $target")) : mkpath(dirname(target))
    own=("source","selection");budget=sum(filesize,files)+16*1024^2
    GuardedFiles.space_check(dirname(target),budget,1024^3)
    tmp=target*".tmp";out=h5open(tmp,"w")
    try
        first_file=h5open(files[1],"r")
        schema=read(first_file["schema"],String);reference=Dict(k=>read(first_file["configuration/checksums"][k]) for k in keys(first_file["configuration/checksums"]) if !(k in own))
        cuts=(read(first_file["configuration/ion_hit_cut_keV"]),read(first_file["configuration/electron_hit_cut_keV"]))
        event_keys=collect(keys(first_file["events"]));hit_keys=[k for k in keys(first_file["itrack"]) if k!="offsets"]
        out["schema"]=schema;HDF5.copy_object(first_file["configuration"],out,"configuration");close(first_file)
        columns=Dict{String,Vector}();offsets=Dict("itrack"=>Int64[0],"etrack"=>Int64[0]);labels=String[];counts=Int64[];shas=Dict(k=>String[] for k in own)
        for path in files
            h5open(path,"r") do f
                read(f["schema"],String)==schema||error("Schema differs in $path")
                (read(f["configuration/ion_hit_cut_keV"]),read(f["configuration/electron_hit_cut_keV"]))==cuts||error("Hit cuts differ in $path")
                sums=f["configuration/checksums"]
                Set(k for k in keys(sums) if !(k in own))==Set(keys(reference))||error("Checksum names differ in $path")
                for (k,v) in reference;read(sums[k])==v||error("Checksum $k differs in $path");end
                Set(keys(f["events"]))==Set(event_keys)||error("Event columns differ in $path")
                for k in event_keys;append!(get!(columns,"events/$k",similar(read(f["events/$k"]),0)),read(f["events/$k"]));end
                for arm in ("itrack","etrack")
                    off=read(f["$arm/offsets"]);first(off)==0||error("Invalid offsets in $path")
                    append!(offsets[arm],off[2:end].+last(offsets[arm]))
                    for k in hit_keys;append!(get!(columns,"$arm/$k",similar(read(f["$arm/$k"]),0)),read(f["$arm/$k"]));end
                end
                push!(labels,read(f["configuration/run_label"],String));push!(counts,length(read(f["events/event_number"])))
                for k in own;push!(shas[k],read(sums[k]));end
            end
        end
        ids=columns["events/event_number"];allunique(ids)||error("Duplicate event numbers across batches")
        for arm in ("itrack","etrack")
            length(offsets[arm])==length(ids)+1&&issorted(offsets[arm])||error("Invalid chained offsets")
            all(k->length(columns["$arm/$k"])==last(offsets[arm]),hit_keys)||error("Hit columns and offsets disagree")
            columns["$arm/offsets"]=offsets[arm]
        end
        for (name,values) in columns
            chunk=(max(1,min(length(values),startswith(name,"events/")||endswith(name,"offsets") ? 256 : chunk_hits)),)
            isempty(values) ? (out[name]=values) : write(create_dataset(out,name,eltype(values),(length(values),);chunk,shuffle=(),deflate=compression_level),values)
        end
        g=create_group(out,"batches");g["run_label"]=labels;g["events"]=counts;g["file_sha256"]=[bytes2hex(open(sha256,p)) for p in files]
        for k in own;g["$(k)_sha256"]=shas[k];end
        flush(out);GuardedFiles.check_size(tmp,budget);close(out)
        Base.rename(tmp,target)
    catch
        isopen(out)&&close(out);rm(tmp;force=true);rethrow()
    end
    String(target)
end

"""Add the selection counters of the batches' detector-response controls; write them as TOML."""
function merge_selection_counts(controls::Vector{String},target::AbstractString)::Dict{String,Any}
    total=Dict{String,Any}()
    for path in controls
        h5open(path,"r") do f
            g=f["selection"]
            for k in keys(g)
                g[k] isa HDF5.Dataset||continue
                v=read(g[k])
                v isa Integer ? (total[k]=get(total,k,0)+v) : (get!(total,k,v)==v||error("Selection setting $k differs in $path"))
            end
        end
    end
    total["batches"]=length(controls)
    open(io->TOML.print(io,total;sorted=true),target,"w")
    total
end
