# Batch orchestration invoked by analyze_mc_selection.jl --pixel-charge.
using Random, DataFrames
include(joinpath(@__DIR__,"pixel_diffusion_cache.jl"))

"""Run full-cloud transport and pixel/frame reduction for both 100 μm samples.

Use separate deterministic drift/diffusion/grid-phase streams matching the viewer
at realisation zero. Stream events through one input reader and retain sparse
pixel/frame counts and per-event axes in HDF5. Retain every positive pixel-frame cell for subsequent histogram reductions.
Align pixel indices to the pixel containing the cloud reference. All events
reach collection. Build shared histogram edges after processing both species.
"""
function run_pixel_charge()
    card=abspath(get(ENV,"PIXEL_CHARGE_CARD",joinpath(ROOT,"config/pixel_charge.toml")))
    c=TOML.parsefile(card);tc=c["transport"];pc=c["pads"];detpath=joinpath(ROOT,tc["detector_card"])
    det=TOML.parsefile(detpath);pitch=Float32(get(pc,"pitch_mm",Float32(det["geometry"]["pitch_m"])*1000f0))
    speed=Float32(det["gas"]["bulk_speed_m_s"]);w=Float32(tc["ionisation_energy_eV"])
    output=joinpath(ROOT,c["run"]["output_directory"]);mkpath(output)
    n_min=Int64(c["statistics"]["n_min"]);n_min>=0 || error("n_min must be nonnegative")
    for sample in ("bb0nu","electrons")
        input=realpath(joinpath(ROOT,c["samples"][sample]));target=joinpath(output,sample*".h5")
        temporary=tempname(joinpath(ROOT,"work"))*".h5"
        events=DataFrame(event_id=Int64[],drift_cm=Float32[],sigma_mm=Float32[],seed=String[],
            produced_ions=Int64[],retained_ions=Int64[],outside_ions=Int64[],collected_ions=Int64[],
            cells_with_charge=Int64[],pulse_count=Int64[],has_pulse=Bool[],
            maximum_ions=Int64[],minimum_ions=Int64[],grid_offset_x_mm=Float32[],grid_offset_y_mm=Float32[])
        cells=DataFrame(event_id=Int64[],pixel_x=Int64[],pixel_y=Int64[],time_slice=Int64[],n_ions=Int64[])
        source=h5open(input);ids=read(source["events/event_id"]);ids=ids[1:min(length(ids),c["run"]["n_events"])]
        cache,cache_path=open_pixel_diffusion_cache(c,input,sample,det,ROOT)
        println(sample,": diffusion cache ",cache_path);flush(stdout)
        cache_hits=0;cache_misses=0;cache_read_seconds=0.0;cache_write_seconds=0.0
        started=time_ns();reading=0.0;transport_time=0.0;collection_time=0.0;writing=0.0
        h5open(temporary,"w") do f
            attributes(f)["schema"]="nausicaa.pixel_charge.v1"
            for (i,id) in enumerate(ids)
                t=time_ns();cached=read_cached_diffusion(cache,id);cache_read_seconds+=(time_ns()-t)/1e9
                if cached===nothing
                    "--repixelize-only" in ARGS && error("No cached diffusion for $sample event $id. Run without --repixelize-only to build missing clouds.")
                    cache_misses+=1
                    t=time_ns();hits,_=read_event(source,id);reading+=(time_ns()-t)/1e9
                    token=bytes2hex(sha256(join((input,string(id),string(c["run"]["seed"]),"0"),"|")))
                    seed=parse(UInt64,token[1:16];base=16)
                    t=time_ns();drift=sample_drift_cm(Xoshiro(seed),Float32(tc["min_cm"]),Float32(tc["max_cm"]))
                    sigma=sigma_t_ion_mm(Float32(det["gas"]["temperature_K"]),drift,Float32(det["bias"]["bulk_field_V_m"])/100f0)
                    n=ion_counts(hits,w)
                    d=diffuse_xyz_image_mc(Xoshiro(seed⊻0x9e3779b97f4a7c15),hits,n.sampled;
                        sigma_t_mm=sigma,sigma_l_mm=sigma,nbins=Int(tc["nbins"]),nsigma=Float32(tc["nsigma"]),ionisation_energy_eV=w)
                    transport_time+=(time_ns()-t)/1e9
                    t=time_ns();save_cached_diffusion(cache,id,d,drift,sigma,seed);cache_write_seconds+=(time_ns()-t)/1e9
                else
                    cache_hits+=1;d=cached.d;drift=cached.drift;sigma=cached.sigma;seed=cached.seed
                end
                t=time_ns();phase=sample_grid_offset_mm(Xoshiro(seed⊻0xd1b54a32d192ed03),pitch)
                axes,history=collect_ion_cloud(d.bins;pitch_mm=pitch,z_bin_width_mm=Float32(pc["z_bin_width_mm"]),
                    grid_offset_mm=phase,drift_speed_mm_ms=speed,drift_direction=Int(pc["drift_direction"]),
                    minimum_xy_span_mm=Float32(pc["minimum_xy_span_mm"]))
                aligned=align_pad_history(history,axes.reference_mm)
                stats=pixel_frame_statistics(history;n_min)
                @assert history.total_ions==d.retained_ions && d.sampled_ions==d.retained_ions+d.outside_ions
                for row in eachrow(aligned)
                    push!(cells,(id,row.pad_dx,row.pad_dy,row.arrival_bin,row.ions))
                end
                push!(events,(id,drift,sigma,string(seed),d.sampled_ions,d.retained_ions,d.outside_ions,
                    history.total_ions,nrow(aligned),length(stats.counts),!isempty(stats.counts),
                    something(stats.maximum_ions,0),something(stats.minimum_ions,0),phase...))
                collection_time+=(time_ns()-t)/1e9
                t=time_ns();base="event_axes/"*string(id)*"/"
                for (key,value) in (("x_edges_mm",history.x_edges_mm),("y_edges_mm",history.y_edges_mm),
                    ("z_edges_mm",history.z_edges_mm),("time_edges_ms",history.time_edges_ms),
                    ("reference_mm",collect(axes.reference_mm)),("grid_offset_mm",collect(phase)))
                    f[base*key]=value
                end
                writing+=(time_ns()-t)/1e9
                if i%100==0;println(sample,": ",i,"/",length(ids)," events; ",round((time_ns()-started)/1e9;digits=1)," s");flush(stdout);end
            end
            close(source);close(cache)
            for (group,table) in (("events",events),("cells",cells)),name in propertynames(table)
                f[group*"/"*string(name)]=table[!,name]
            end
            attributes(f["cells"])["index_convention"]="pixel_x/y relative to reference pixel; time_slice zero-based relative arrival frame"
            f["summary/processed_events"]=length(ids)
            f["summary/events_with_charge"]=count(events.has_pulse)
            f["summary/pulse_count"]=sum(events.pulse_count)
            f["summary/produced_ions"]=sum(events.produced_ions)
            f["summary/outside_ions"]=sum(events.outside_ions)
            f["summary/collected_ions"]=sum(events.collected_ions)
            f["provenance/diffusion_cache_path"]=cache_path
            f["timing/diffusion_cache_hits"]=cache_hits;f["timing/diffusion_cache_misses"]=cache_misses
            f["timing/cache_read_seconds"]=cache_read_seconds;f["timing/cache_write_seconds"]=cache_write_seconds
            f["settings/n_min"]=n_min;f["settings/pitch_mm"]=pitch
            f["settings/frame_ms"]=Float32(pc["z_bin_width_mm"])/speed
            f["settings/drift_direction"]=pc["drift_direction"]
            f["provenance/card_toml"]=read(card,String);f["provenance/detector_toml"]=read(detpath,String)
            f["provenance/source_path"]=input;f["provenance/source_sha256"]=open(io->bytes2hex(sha256(io)),input)
            paths=sort(vcat([relpath(@__FILE__,ROOT),"scripts/analyze_mc_selection.jl","scripts/pixel_diffusion_cache.jl","Manifest.toml"],
                [relpath(joinpath(dir,file),ROOT) for (dir,_,files) in walkdir(joinpath(ROOT,"src/event_analysis")) for file in files if endswith(file,".jl")]))
            f["provenance/code_paths"]=paths;f["provenance/code_sha256"]=[bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in paths]
            f["provenance/julia_version"]=string(VERSION)
            f["timing/read_seconds"]=reading;f["timing/transport_seconds"]=transport_time
            f["timing/collection_reduction_seconds"]=collection_time;f["timing/event_axes_write_seconds"]=writing
            f["timing/elapsed_seconds"]=Float64((time_ns()-started)/1e9)
        end
        mv(temporary,target;force=true)
        println(sample,": saved ",nrow(events)," events, ",nrow(cells)," pixel-frame cells to ",target);flush(stdout)
    end
end
