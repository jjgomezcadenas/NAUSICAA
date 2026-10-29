# Persistent application worker. stdin/stdout carry request IDs; logs use stderr.
using HDF5, TOML, SHA
include(joinpath(@__DIR__,"../event_analysis/EventAnalysis.jl"))
using .EventAnalysis

const ROOT=normpath(joinpath(@__DIR__,"../.."))
const WORK=abspath(ARGS[1])
const CARD=abspath(ARGS[2])
const CONFIG=TOML.parsefile(CARD)
include(joinpath(@__DIR__,"transport.jl"))
include(joinpath(@__DIR__,"pads.jl"))
const CODE_PATHS=sort(vcat(["Project.toml","Manifest.toml","src/event_viewer/worker.jl","src/event_viewer/transport.jl","src/event_viewer/pads.jl","src/ion_signal/event_io.jl"],
    [relpath(joinpath(dir,file),ROOT) for (dir,_,files) in walkdir(joinpath(ROOT,"src/event_analysis"))
     for file in files if endswith(file,".jl")]))
const CODE_HASHES=[bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in CODE_PATHS]
const SOURCE_HASH=bytes2hex(sha256(join(CODE_PATHS.*":".*CODE_HASHES,"\n")))

"""Write owned table columns to an HDF5 group with their Julia datatypes."""
function save_columns(g,table)
    for column in propertynames(table);g[string(column)]=table[!,column];end
end

"""Build the current event's MC graph with the TOML settings; own all outputs.
Energy fractions use the full original event energy as denominator. Julia saves
the energy retained by the voxel threshold and the energy below that threshold.
"""
function reconstruct(hits,request_id)
    s=MCGraphSettings(; (Symbol(k)=>v for (k,v) in CONFIG["mc_graph"])...)
    distance=resolved_mc_distance(s)
    voxels=voxelize_event(to_petit_dataframe(hits),s.voxel_width_mm)
    components=build_mc_components(voxels;method=s.method,max_distance_mm=distance,
        energy_threshold_keV=s.energy_threshold_keV,k=s.k)
    keep=voxels.energy .>= s.energy_threshold_keV*1f-3
    total=sum(hits.energy_MeV)
    fractions=Float32[total>0f0 ? c.energy_MeV/total : 0f0 for c in components]
    return (;settings=s,distance,voxels,components,fractions,request_id,
        retained_energy_MeV=sum(voxels.energy[keep]),below_threshold_energy_MeV=sum(voxels.energy[.!keep]))
end

"""Write graph-local one-based indices, applied settings and component summaries.
All physical arrays are Julia Float32; counts and persisted indices are Int64.
The build request ID identifies the cached reconstruction across display rebins.
"""
function save_graph(f,g)
    f["reconstruction/status"]="complete"
    f["reconstruction/build_request_id"]=g.request_id
    for name in fieldnames(MCGraphSettings)
        f["configuration/mc_graph/"*string(name)]=getfield(g.settings,name)
    end
    f["configuration/mc_graph/resolved_distance_mm"]=g.distance
    v=g.voxels
    for (name,col) in (("event_id",v.event_id),("x_mm",v.x),("y_mm",v.y),("z_mm",v.z),("energy_MeV",v.energy))
        f["reconstruction/voxels/"*name]=col
    end
    f["reconstruction/summary/component_count"]=Int64(length(g.components))
    f["reconstruction/summary/retained_energy_MeV"]=g.retained_energy_MeV
    f["reconstruction/summary/below_threshold_energy_MeV"]=g.below_threshold_energy_MeV
    create_group(f,"reconstruction/components")
    for c in g.components
        base="reconstruction/components/$(c.component_id)/"
        save_columns(create_group(f,base*"nodes"),c.nodes)
        for (name,indices) in (("edges/src",c.edge_src),("edges/dst",c.edge_dst),
                ("backbone/path_indices",c.path_indices),
                ("backbone/endpoint_indices",c.endpoint_indices===nothing ? Int64[] : collect(c.endpoint_indices)))
            f[base*name]=indices;attributes(f[base*name])["index_base"]=1
        end
        f[base*"backbone/s_mm"]=c.path_s_mm
        for (name,value) in (("energy_MeV",c.energy_MeV),("energy_fraction",g.fractions[c.component_id]),
                ("node_count",Int64(length(c.nodes.x_mm))),("edge_count",Int64(length(c.edge_src))),
                ("length_mm",c.length_mm),("confidence",c.confidence))
            f[base*"summary/"*name]=value
        end
    end
end

"""Persist one completed MC display snapshot in the agreed viewer v1 schema.
The worker closes the temporary file before replacing current.h5. Float32
physical arrays and Int64 indices/counts retain their Julia types. HDF5 axis
attributes state the Julia and Python matrix conventions. All input arrays are
preserved. The reader owns and retains its original input-file handle.
"""
function write_snapshot(hits,info,maps,source,sample,position,event_count,request_id,nbins;graph=nothing,transport=nothing,active_transport=false,pads=nothing,active_pads=false)
    target=joinpath(WORK,"current.h5");temporary=joinpath(WORK,"current.tmp.h5")
    h5open(temporary,"w") do f
        attributes(f)["schema"]="nausicaa.event.viewer.v2"
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
        if graph===nothing;f["reconstruction/status"]="pending";else;save_graph(f,graph);end
        f["context/active_transport"]=active_transport
        f["context/active_pads"]=active_pads
        pads===nothing || save_pads(f,pads)
        transport===nothing || save_transport(f,transport)
    end
    mv(temporary,target;force=true)
    return target
end

"""Run ordered TOML requests, retaining the sample handle and latest event.
An unsuccessful request returns its ID and error without replacing the last
completed display snapshot. `open_sample` writes a catalogue of original IDs;
requests may supply source_path for a user-selected reduced HDF5 file.
`prepare_file` accepts a reduced v3 file or a NEXUS file, creates and verifies
the _hitdf.h5 through EventIO when needed, and returns its path. Raw archive
inputs map to project events/; other inputs map beside their source.
An existing reduced companion is reused after checking its schema and source
SHA-256. A mismatch returns regeneration guidance and preserves both files.
`load_event` loads the requested ID and builds maps; `reconstruct` also builds
the graph on first request. Maps and graph are cached for the current event;
display rebinning preserves the graph. `close` releases the handle.
"""
function main()
    source_file=nothing;sample="";source="";ids=Int64[];hits=nothing;info=nothing
    graph=nothing;maps=nothing;map_nbins=0
    transported=nothing;transport_key=nothing;collected=nothing;pad_key=nothing
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
            elseif op=="prepare_file"
                path=abspath(r["path"])
                kind=h5open(path,"r") do f
                    haskey(f,"MC/hits") ? :nexus :
                        (haskey(attributes(f),"schema") && read(attributes(f)["schema"])=="nausicaa.hitdf.v3" ? :hitdf : :unknown)
                end
                kind==:unknown && error("Select a NEXUS event file or a hitdf v3 file")
                if kind==:nexus
                    target=EventAnalysis.EventIO.default_output(path)
                    if isfile(target)
                        h5open(target,"r") do f
                            read(attributes(f)["schema"])=="nausicaa.hitdf.v3" || error("Existing reduced file needs regeneration: $target")
                            read(f["metadata/source_sha256"])==EventAnalysis.EventIO.file_sha256(path) ||
                                error("Existing reduced file belongs to different source contents: $target. Regenerate it with scripts/events/extract_hits.jl --overwrite.")
                        end
                    else
                        EventAnalysis.EventIO.convert_hits(path)
                    end
                    path=target
                end
                response["path"]=path
            elseif op in ("open_sample","load_event","reconstruct","transport","pads")
                requested=r["sample"]
                if requested!=sample
                    path=get(r,"source_path",normpath(joinpath(ROOT,get(CONFIG["samples"],requested,""))))
                    new_file=h5open(path,"r")
                    try
                        new_ids=read(new_file["events/event_id"])
                        read(attributes(new_file)["schema"])=="nausicaa.hitdf.v3" || error("Viewer requires hitdf v3")
                        source_file!==nothing && close(source_file)
                        source_file=new_file;ids=new_ids;sample=requested;source=path;hits=nothing;info=nothing
                        graph=nothing;maps=nothing;map_nbins=0
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
                    if hits===nothing || hits.event_id!=id
                        hits,info=read_event(source_file,id);graph=nothing;maps=nothing;map_nbins=0
                    end
                    nbins=Int(r["nbins"])
                    1<=nbins<=500 || error("Display bin target must be between 1 and 500")
                    if maps===nothing || map_nbins!=nbins;maps=bin_mc_energy(hits;nbins);map_nbins=nbins;end
                    if op=="reconstruct" && graph===nothing;graph=reconstruct(hits,request_id);end
                    if op in ("transport","pads")
                        key=(source,id,r["drift_mode"],Float32(r["fixed_cm"]),Int(r["realisation"]))
                        if key!=transport_key
                            transported=transport_event(hits,source,key[3],key[4],key[5]);transport_key=key
                        end
                    end
                    if op=="pads" && pad_key!=transport_key
                        collected=collect_pads(transported);pad_key=transport_key
                    end
                    current_transport=transport_key!==nothing && transport_key[1:2]==(source,id) ? transported : nothing
                    response["path"]=write_snapshot(hits,info,maps,source,sample,position,length(ids),request_id,nbins;graph,
                        transport=current_transport,active_transport=op in ("transport","pads"),
                        pads=pad_key==transport_key && current_transport!==nothing ? collected : nothing,active_pads=op=="pads")
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
