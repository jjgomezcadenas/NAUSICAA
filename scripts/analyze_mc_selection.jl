#!/usr/bin/env julia
# Julia batch reduction and process benchmark. Run from any working directory.
using HDF5, TOML, SHA, Statistics
const ROOT=normpath(joinpath(@__DIR__,".."))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis
if "--pixel-charge" in ARGS
    include(joinpath(@__DIR__,"run_pixel_charge.jl"))
    include(joinpath(@__DIR__,"reduce_pixel_charge.jl"))
    run_pixel_charge()
    pixel_card=TOML.parsefile(get(ENV,"PIXEL_CHARGE_CARD",joinpath(ROOT,"config/pixel_charge.toml")))
    reduce_pixel_charge(joinpath(ROOT,pixel_card["run"]["output_directory"]))
    summarize_pixel_charge(joinpath(ROOT,pixel_card["run"]["output_directory"]))
    exit()
end
const CARD=abspath(get(ENV,"MC_SELECTION_CARD",joinpath(ROOT,"config/mc_selection.toml")))
const CFG=TOML.parsefile(CARD)
const SETTINGS=MCGraphSettings(;(Symbol(k)=>v for (k,v) in CFG["mc_graph"])...)

"""Read ordered IDs and reduce a contiguous catalogue subset with one owned HDF5 reader.
Return per-event results and measured read/analysis seconds. Compilation is warmed
by the caller. Physical calculations remain in analyze_mc_event.
"""
function batch(path,ids)
    results=Any[];read_seconds=0.0;analysis_seconds=0.0
    h5open(path) do f
        for id in ids
            t=time_ns();hits,_=read_event(f,id);read_seconds+=(time_ns()-t)/1e9
            t=time_ns();push!(results,analyze_mc_event(hits,SETTINGS));analysis_seconds+=(time_ns()-t)/1e9
        end
    end
    return results,read_seconds,analysis_seconds
end
"""Write event-summary columns from owned results; preserve input order and types."""
function write_events(f,results)
    for name in keys(first(results).event)
        f["events/"*string(name)]=[getproperty(r.event,name) for r in results]
    end
end
"""Count half-open bins, including the final upper edge, and explicit under/overflow."""
function histogram(values,edges)
    counts=zeros(Int64,length(edges)-1);under=Int64(0);over=Int64(0)
    for v in values
        if v<first(edges);under+=1
        elseif v>last(edges);over+=1
        else;counts[min(searchsortedlast(edges,v),length(counts))]+=1
        end
    end
    return counts,under,over
end
"""Save common histogram edges and before/after counts; voxel entries are pooled over events."""
function write_histograms(f,results)
    for (name,values,units) in (
        ("voxel_energy_MeV",r->r.voxel_energy_MeV,"MeV"),
        ("retained_voxel_energy_MeV",r->r.retained_voxel_energy_MeV,"MeV"),
        ("n_voxels",r->[r.event.n_voxels],"1"),
        ("n_retained_voxels",r->[r.event.n_retained_voxels],"1"),
        ("track_length_mm",r->[r.event.track_length_mm],"mm"),
        ("component_length_mm",r->r.component_length_mm,"mm"),
        ("voxel_energy_sum_MeV",r->[r.event.voxel_energy_MeV],"MeV"),
        ("retained_energy_sum_MeV",r->[r.event.retained_energy_MeV],"MeV"),
        ("n_components",r->[r.event.n_components],"1"))
        before=reduce(vcat,(values(r) for r in results);init=Float64[])
        after=reduce(vcat,(values(r) for r in results if r.event.accepted);init=Float64[])
        # Shared, fixed edges across species. Overflow is stored and reported.
        edges= name in ("voxel_energy_MeV","retained_voxel_energy_MeV") ? collect(0.0:.002:1.0) :
            occursin("energy_sum",name) ? collect(0.0:.025:3.0) :
            name=="n_components" ? collect(-.5:1:49.5) :
            occursin("length",name) ? collect(0.0:5:500.0) : collect(-.5:20:3999.5)
        base="histograms/"*name;f[base*"/edges"]=edges
        attributes(f[base])["units"]=units
        for (label,data) in (("before",before),("after",after))
            c,u,o=histogram(data,edges)
            f[base*"/"*label*"/counts"]=c;f[base*"/"*label*"/underflow"]=u;f[base*"/"*label*"/overflow"]=o
        end
    end
end

if !isempty(ARGS) && ARGS[1]=="--worker"
    path=ARGS[2];ids=parse.(Int64,split(ARGS[3],','));output=ARGS[4];gate=ARGS[5]
    warm=@elapsed batch(path,ids[1:1])
    touch(output*".ready")
    while !isfile(gate);sleep(.02);end
    elapsed=@elapsed results,reading,analysis=batch(path,ids)
    h5open(output,"w") do f
        write_events(f,results);f["reading_seconds"]=reading;f["analysis_seconds"]=analysis
        f["batch_seconds"]=elapsed;f["warmup_seconds"]=warm
    end
    exit()
end

out=joinpath(ROOT,CFG["run"]["output_directory"]);mkpath(out)
work=mktempdir(joinpath(ROOT,"work");prefix="mc-selection-")
for sample in ("bb0nu","electrons")
    path=joinpath(ROOT,CFG["samples"][sample])
    ids=h5open(f->read(f["events/event_id"]),path)
    ids=ids[1:min(length(ids),CFG["run"]["n_events"])]
    println(sample,": warming serial calculation");flush(stdout)
    warm=@elapsed batch(path,ids[1:1])
    serial=@elapsed results,reading,analysis=batch(path,ids)
    println(sample,": serial ",round(serial;digits=3)," s; starting process benchmark");flush(stdout)
    nw=min(CFG["run"]["workers"],length(ids));gate=joinpath(work,sample*".start")
    processes=[];outputs=String[];logs=IO[];cold_start=time_ns()
    for i in 1:nw
        subset=ids[i:nw:end];output=joinpath(work,"$(sample)_$(i).h5");push!(outputs,output)
        log=open(output*".log","w");push!(logs,log)
        cmd=`$(Base.julia_cmd()) --startup-file=no --project=$ROOT $(@__FILE__) --worker $path $(join(subset,',')) $output $gate`
        push!(processes,run(pipeline(cmd;stdout=log,stderr=log);wait=false))
    end
    while !all(isfile(p*".ready") for p in outputs)
        any(process_exited(p) for p in processes) && error("Worker exited during warm-up; see $work")
        sleep(.05)
    end
    startup=(time_ns()-cold_start)/1e9
    t=time_ns();touch(gate);foreach(wait,processes);parallel=(time_ns()-t)/1e9
    all(success,processes) || error("Worker failed; see $work")
    foreach(close,logs)
    # Compare every saved event scalar against the serial reference by event ID.
    reference=Dict(r.event.event_id=>r.event for r in results)
    for output in outputs
        h5open(output) do f
            for (i,id) in enumerate(read(f["events/event_id"])),name in keys(reference[id])
                @assert read(f["events/"*string(name)])[i]==getproperty(reference[id],name)
            end
        end
    end
    target=joinpath(out,sample*".h5")
    h5open(target,"w") do f
        attributes(f)["schema"]="nausicaa.mc_selection.v1"
        write_events(f,results);write_histograms(f,results)
        for (k,v) in CFG["mc_graph"];f["settings/"*k]=v;end
        f["provenance/source_path"]=path
        f["provenance/source_sha256"]=open(io->bytes2hex(sha256(io)),path)
        f["provenance/card_toml"]=read(CARD,String)
        f["provenance/code_sha256"]=bytes2hex(sha256(read(joinpath(ROOT,"src/event_analysis/mc_selection.jl"))))
        f["provenance/julia_version"]=string(VERSION)
        accepted=count(r->r.event.accepted,results);n=length(results);eff=accepted/n
        f["selection/processed"]=n;f["selection/accepted"]=accepted;f["selection/efficiency"]=eff
        f["selection/binomial_standard_error"]=sqrt(eff*(1-eff)/n)
        # Wilson score interval, z=1 (approximately 68% coverage).
        centre=(eff+1/(2n))/(1+1/n);half=sqrt(eff*(1-eff)/n+1/(4n^2))/(1+1/n)
        f["selection/wilson_68_low"]=centre-half;f["selection/wilson_68_high"]=centre+half
        for (key,value) in (("serial_warmup_seconds",warm),("serial_batch_seconds",serial),
            ("serial_read_seconds",reading),("serial_analysis_seconds",analysis),
            ("parallel_startup_warmup_seconds",startup),("parallel_batch_seconds",parallel),
            ("parallel_total_seconds",startup+parallel),("workers",nw),("serial_parallel_identical",true))
            f["timing/"*key]=value
        end
        f["validation/max_energy_closure_MeV"]=maximum(abs(r.event.energy_closure_MeV) for r in results)
        f["validation/max_voxel_error_vs_Float64_MeV"]=maximum(r.event.voxel_reduction_error_MeV for r in results)
        @assert maximum(r.event.voxel_reduction_error_MeV for r in results)<1.25e-4 # 1% of 12.5 keV resolution
    end
    println(sample,": ",count(r->r.event.accepted,results),"/",length(ids)," accepted; serial=",serial,
        " s; parallel warm=",parallel," s; startup+warmup=",startup," s; results=",target);flush(stdout)
end
println("Worker logs: ",work)
