module PixelArea
using TOML,HDF5,SHA
export read_devices, pack_devices, validate_boxes, analyze

function spice_value(token,parameters)
    startswith(token,"{") && return parameters[token[2:end-1]]
    m=match(r"^([+-]?[0-9.]+(?:[eE][+-]?\d+)?)([a-zA-Z]*)$",token)
    m===nothing && error("Unsupported SPICE value $token")
    factors=Dict(""=>1.,"u"=>1e-6,"n"=>1e-9,"p"=>1e-12,"f"=>1e-15,"m"=>1e-3)
    parse(Float64,m[1])*factors[lowercase(m[2])]
end
function read_devices(path)
    lines=readlines(path);parameters=Dict{String,Float64}()
    for line in lines
        m=match(r"^\.param(?:eter)?\s+(\w+)=(\S+)",line)
        m===nothing || (parameters[m[1]]=spice_value(m[2],parameters))
    end
    devices=NamedTuple[]
    for line in lines
        startswith(line,"M") || continue
        tokens=split(line);values=Dict(split(x,'=') for x in tokens[7:end])
        w=spice_value(values["w"],parameters);l=spice_value(values["l"],parameters)
        multiplier=spice_value(values["m"],parameters)
        @assert multiplier==round(multiplier) && multiplier>0 && w>0 && l>0
        for j in 1:Int(multiplier)
            name=multiplier==1 ? String(tokens[1]) : "$(tokens[1])_$j"
            push!(devices,(name=name,instance=String(tokens[1]),W_m=w,L_m=l))
        end
    end
    isempty(devices) && error("No MOS devices found in $path")
    devices
end

"""Constructive shelf packing of independent channel envelopes, rotated as needed.
The shorter side faces the narrow strip. Margins surround every device; the
passive reserve has its own fixed dimensions. No folding or shared diffusion.
"""
function pack_devices(devices; width, margin, gap, outer, passive_side)
    @assert width>0 && margin>=0 && gap>=0 && outer>=0 && passive_side>0
    inner=width-2outer;inner>0 || error("Margins exhaust CSA strip width")
    boxes=[(name=d.name,w=min(d.W_m,d.L_m)+2margin,h=max(d.W_m,d.L_m)+2margin,kind="MOS envelope") for d in devices]
    push!(boxes,(name="Cf reserve",w=passive_side,h=passive_side,kind="Passive allocation"))
    sort!(boxes;by=b->(-b.h,-b.w,b.name))
    placed=NamedTuple[];x=outer;y=outer;rowheight=0.
    for box in boxes
        box.w<=inner+1e-15 || error("$(box.name) exceeds strip width")
        if x>outer && x+box.w>width-outer+1e-15
            x=outer;y+=rowheight+gap;rowheight=0.
        end
        push!(placed,(name=box.name,x=x,y=y,w=box.w,h=box.h,kind=box.kind))
        x+=box.w+gap;rowheight=max(rowheight,box.h)
    end
    height=y+rowheight+outer
    validate_boxes(placed,width,height)
    placed,height
end

function validate_boxes(boxes,width,height)
    for box in boxes
        @assert min(box.x,box.y)>=-1e-15 && min(box.w,box.h)>0
        @assert box.x+box.w<=width+1e-15 && box.y+box.h<=height+1e-15
    end
    for i in eachindex(boxes),j in 1:i-1
        a=boxes[i];b=boxes[j]
        overlap=min(a.x+a.w,b.x+b.w)-max(a.x,b.x)>1e-15 && min(a.y+a.h,b.y+b.h)-max(a.y,b.y)>1e-15
        @assert !overlap
    end
    true
end
function write_boxes(group,boxes)
    group["names"]=[b.name for b in boxes];group["kinds"]=[b.kind for b in boxes]
    for field in (:x,:y,:w,:h);group[string(field)*"_m"]=[getproperty(b,field) for b in boxes];end
    group["area_m2"]=[b.w*b.h for b in boxes]
end

function analyze(root)
    configpath=joinpath(root,"config/pixel_area.toml");c=TOML.parsefile(configpath)
    netlist=joinpath(root,c["inputs"]["netlist"]);detector=joinpath(root,c["inputs"]["detector"])
    geometry=TOML.parsefile(detector)["geometry"];budget=c["budget"];refs=c["references"]
    devices=read_devices(netlist)
    pitch=geometry["pitch_m"];pad=geometry["pad_width_m"];side=budget["envelope_side_m"]
    clear=budget["pad_clearance_m"];keep=pad+2clear
    left=budget["shared_strip_width_m"];right=budget["csa_strip_width_m"]
    middle=side-left-right;edge=(side-keep)/2
    @assert side<=pitch && middle>=keep && edge>0
    regions=[(name="Pad and clearance",x=left,y=edge,w=middle,h=keep,kind="Geometric reserve"),
        (name="CSA strip",x=side-right,y=0.,w=right,h=side,kind="Circuit packing capacity"),
        (name="Shared services",x=0.,y=0.,w=left,h=side,kind="Unimplemented reserve"),
        (name="Buffer and selection",x=left,y=side-edge,w=middle,h=edge,kind="Unimplemented reserve"),
        (name="Control and routing",x=left,y=0.,w=middle,h=edge,kind="Unimplemented reserve")]
    validate_boxes(regions,side,side)
    @assert isapprox(sum(r.w*r.h for r in regions),side^2;rtol=1e-12)
    totalgate=sum(d.W_m*d.L_m for d in devices)
    yang=refs["yang_chip_width_m"]*refs["yang_chip_height_m"]
    liang=refs["liang_chip_width_m"]*refs["liang_chip_height_m"]
    array_side=round(Int,sqrt(budget["array_pixels"]))
    @assert array_side^2==budget["array_pixels"]
    tile_side=array_side*pitch
    shared_area=left*side*budget["array_pixels"]
    edge_band=shared_area/(2tile_side)
    @assert edge_band<=(pitch-side)/2
    outdir=joinpath(root,"results/pixel_area");mkpath(outdir)
    h5open(joinpath(outdir,"analysis.h5"),"w") do o
        a=attributes(o)
        for (k,v) in ("pitch_m"=>pitch,"pad_width_m"=>pad,"cell_area_m2"=>pitch^2,"pad_area_m2"=>pad^2,"envelope_side_m"=>side,"envelope_area_m2"=>side^2,"gate_area_m2"=>totalgate,"mos_instances"=>length(devices),"pad_clearance_m"=>clear,"yang_chip_area_m2"=>yang,"liang_chip_area_m2"=>liang,"array_pixels"=>budget["array_pixels"],"yang_area_per_array_pixel_m2"=>yang/budget["array_pixels"],"shared_reserve_total_m2"=>left*side*budget["array_pixels"],"minimum_pixels_for_yang_area_proxy"=>ceil(Int,yang/(left*side)),"spare_cell_area_m2"=>pitch^2-side^2,"array_side_pixels"=>array_side,"tile_side_m"=>tile_side,"shared_edge_band_m"=>edge_band,"edge_free_margin_m"=>(pitch-side)/2)
            a[k]=v
        end
        a["architecture"]=c["inputs"]["readout_architecture"]
        a["limits"]="Constructive rectangle packing with stated margins. Readout/control/shared blocks are allocations, not transistor-level layouts. External/shared ADC. No circuit under collection-pad keepout. No DRC, routing, parasitic extraction, package, die edge or multi-pixel tiling validation."
        files=[configpath,netlist,detector,@__FILE__]
        o["source_paths"]=[relpath(p,root) for p in files];o["source_sha256"]=[bytes2hex(sha256(read(p))) for p in files]
        r=create_group(o,"references");for (k,v) in refs;attributes(r)[k]=v;end
        d=create_group(o,"devices");d["names"]=[d.name for d in devices];d["W_m"]=[d.W_m for d in devices];d["L_m"]=[d.L_m for d in devices];d["channel_area_m2"]=[d.W_m*d.L_m for d in devices]
        write_boxes(create_group(o,"regions"),regions)
        for s in c["scenarios"]
            placed,height=pack_devices(devices;width=right,margin=s["device_margin_m"],gap=s["row_gap_m"],outer=s["outer_margin_m"],passive_side=budget["passive_reserve_side_m"])
            g=create_group(o,s["name"]);b=attributes(g)
            for (k,v) in s;b[k]=v;end
            footprint=right*height;remaining=side^2-right*side+footprint
            for (k,v) in ("csa_width_m"=>right,"csa_height_m"=>height,"csa_area_m2"=>footprint,"padded_device_area_m2"=>sum(x.w*x.h for x in placed),"fits_strip"=>height<=side,"remaining_strip_area_m2"=>right*(side-height),"budget_used_m2"=>remaining)
                b[k]=v
            end
            write_boxes(g,placed)
            @assert height<=side
            # Independent dimensional scaling test, preserving every margin.
            doubled=[merge(d,(W_m=2d.W_m,L_m=2d.L_m)) for d in devices]
            p2,h2=pack_devices(doubled;width=2right,margin=2s["device_margin_m"],gap=2s["row_gap_m"],outer=2s["outer_margin_m"],passive_side=2budget["passive_reserve_side_m"])
            @assert isapprox(h2,2height;rtol=1e-12)
            @assert all(isapprox(p2[i].x,2placed[i].x;atol=1e-15) && isapprox(p2[i].y,2placed[i].y;atol=1e-15) for i in eachindex(placed))
            println(s["name"],": CSA ",right*1e6," x ",height*1e6," um; area=",footprint*1e6," mm2; budget used=",remaining*1e6," mm2")
        end
        println("Channels: ",totalgate*1e12," um2; cell ",pitch^2*1e6," mm2; Topmetal-S ",yang*1e6," and ",liang*1e6," mm2")
    end
end
end
