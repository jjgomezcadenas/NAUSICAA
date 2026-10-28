# Persistent application worker. stdin/stdout carry request IDs; logs use stderr.
using HDF5, TOML, SHA
include(joinpath(@__DIR__,"../event_analysis/EventAnalysis.jl"))
using .EventAnalysis

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const WORK=abspath(ARGS[1])
const CARD=abspath(ARGS[2])
const CONFIG=TOML.parsefile(CARD)
const CODE_PATHS=sort(vcat(["Project.toml","Manifest.toml","src/event_viewer/worker.jl","src/ion_signal/event_io.jl"],
    [relpath(joinpath(dir,file),ROOT) for (dir,_,files) in walkdir(joinpath(ROOT,"src/event_analysis"))
     for file in files if endswith(file,".jl")]))
const CODE_HASHES=[bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in CODE_PATHS]
const SOURCE_HASH=bytes2hex(sha256(join(CODE_PATHS.*":".*CODE_HASHES,"\n")))

"""Write owned table columns to an HDF5 group with their Julia datatypes."""
function save_columns(g,table)
    for column in propertynames(table);g[string(column)]=table[!,column];end
end

"""Persist one completed MC display snapshot in the agreed viewer v1 schema.
The worker closes the temporary file before replacing current.h5. Float32
physical arrays and Int64 indices/counts retain their Julia types. HDF5 axis
attributes state the Julia and Python matrix conventions. All input arrays are
preserved. The reader owns and retains its original input-file handle.
"""
function write_snapshot(hits,info,maps,source,sample,position,event_count,request_id,nbins)
    target=joinpath(WORK,"current.h5");temporary=joinpath(WORK,"current.tmp.h5")
    h5open(temporary,"w") do f
        attributes(f)["schema"]="nausicaa.event.viewer.v1"
        attributes(f)["status"]="complete"
        for (k,v) in ("request_id"=>request_id,"sample"=>sample,"source_path"=>source,
                "source_schema"=>"nausicaa.hitdf.v3",
                "source_size_bytes"=>Int64(filesize(source)),"source_mtime_ns"=>round(Int64,mtime(source)*1e9),
                "code_sha256"=>SOURCE_HASH,"julia_version"=>string(VERSION))
            f["context/"*k]=v
        end
        f["provenance/code_paths"]=CODE_PATHS;f["provenance/code_sha256"]=CODE_HASHES
        f["configuration/card_toml"]=read(CARD,String)
        f["configuration/display_nbins"]=Int64(nbins)
        f["event/event_id"]=hits.event_id
        f["event/position"]=Int64(position);f["event/count"]=Int64(event_count)
        f["event/frame"]=string(hits.frame)
        f["event/info/vertex_mm"]=collect(info.vertex_mm)
        f["event/info/primary_particle_ids"]=info.primary_particle_ids
        for s in (:x_mm,:y_mm,:z_mm,:energy_MeV,:time_ns,:particle_id,:hit_id)
            f["event/hits/"*string(s)]=getproperty(hits,s)
        end
        f["event/summary/hit_count"]=maps.hit_count
        f["event/summary/total_energy_MeV"]=maps.total_energy_MeV
        for s in (:x_edges_mm,:y_edges_mm,:z_edges_mm);f["display/axes/"*string(s)]=getproperty(maps,s);end
        for (prefix,axes) in (("xy","x,y"),("xz","x,z"),("yz","y,z")),suffix in ("energy_MeV","hit_count")
            name=prefix*"_"*suffix;f["display/projections/"*name]=getproperty(maps,Symbol(name))
            attributes(f["display/projections/"*name])["julia_axes"]=axes
            attributes(f["display/projections/"*name])["python_axes"]=join(reverse(split(axes,",")),",")
        end
        save_columns(create_group(f,"display/occupied"),maps.occupied)
        attributes(f["display/occupied"])["index_base"]=1
        f["reconstruction/status"]="pending"
    end
    mv(temporary,target;force=true)
    return target
end

"""Run ordered TOML requests, retaining the sample handle and latest event.
An unsuccessful request returns its ID and error without replacing the last
completed display snapshot. `open_sample` writes a catalogue of original IDs;
`load_event` loads the requested ID and builds maps; `close` releases the handle.
"""
function main()
    source_file=nothing;sample="";source="";ids=Int64[];hits=nothing;info=nothing
    println("READY");flush(stdout)
    for request_id in eachline(stdin)
        occursin(r"^[a-zA-Z0-9_-]+$",request_id) || continue
        response=Dict{String,Any}("request_id"=>request_id)
        finish=false
        try
            r=TOML.parsefile(joinpath(WORK,request_id*".request.toml"))
            op=r["op"]
            if op=="close"
                finish=true
            elseif op in ("open_sample","load_event")
                requested=r["sample"]
                if requested!=sample
                    path=normpath(joinpath(ROOT,CONFIG["samples"][requested]))
                    new_file=h5open(path,"r")
                    try
                        new_ids=read(new_file["events/event_id"])
                        read(attributes(new_file)["schema"])=="nausicaa.hitdf.v3" || error("Viewer requires hitdf v3")
                        source_file!==nothing && close(source_file)
                        source_file=new_file;ids=new_ids;sample=requested;source=path;hits=nothing;info=nothing
                    catch
                        close(new_file);rethrow()
                    end
                end
                if op=="open_sample"
                    path=joinpath(WORK,"catalogue.h5")
                    h5open(path,"w") do f
                        f["event_id"]=ids;f["sample"]=sample
                    end
                    response["path"]=path
                else
                    id=parse(Int64,string(r["event_id"]));position=findfirst(==(id),ids)
                    position===nothing && error("Event ID $id is absent from $sample")
                    if hits===nothing || hits.event_id!=id;hits,info=read_event(source_file,id);end
                    nbins=Int(r["nbins"])
                    1<=nbins<=500 || error("Display bin target must be between 1 and 500")
                    maps=bin_mc_energy(hits;nbins)
                    response["path"]=write_snapshot(hits,info,maps,source,sample,position,length(ids),request_id,nbins)
                end
            else
                error("Unknown viewer operation: $op")
            end
            response["status"]="ok"
        catch e
            response["status"]="error";response["message"]=sprint(showerror,e)
            showerror(stderr,e,catch_backtrace());println(stderr)
        end
        open(io->TOML.print(io,response),joinpath(WORK,request_id*".response.toml"),"w")
        println(request_id);flush(stdout)
        finish && break
    end
    source_file!==nothing && close(source_file)
end
main()
