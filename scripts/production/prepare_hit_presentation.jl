#!/usr/bin/env julia
include("run_hit_production.jl")

"""Calculate histograms, projections and graph backbones from a completed hit file.
No detector or electronics calculation is repeated. Retained hits are read from
the source file; this output saves only derived presentation data.
"""
function prepare_hit_presentation(path;graph_card=nothing,output_dir=nothing)
    r=open_hit_input(path);base=dirname(abspath(path));out=output_dir===nothing ? joinpath(base,"presentation") : abspath(output_dir);mkpath(out)
    target=joinpath(out,"controls.h5");isfile(target)&&error("Preserve existing presentation data")
    rc=TOML.parse(read(r.file["configuration/cards/run"],String));dc=TOML.parse(read(r.file["configuration/cards/detector_card"],String))
    graph_card===nothing || (rc["track_analysis"]=TOML.parsefile(graph_card)["track_analysis"])
    ids=read(r.file["events/event_number"]);examples=Dict{String,Any}();histograms=Dict{String,Any}()
    tc=rc["track_analysis"];started=time()
    try
      for (j,id) in enumerate(ids)
        e=read_hit_event(r,id);ex=Dict{String,Any}("ordinal"=>j);hh=Dict{String,Any}()
        for (arm,key,h,pitch,dz,sigma,k,mult) in (
          ("itrack","ions",e.itrack,dc["ion_pitch_mm"],dc["ion_bucket_ms"]*dc["ion_drift_velocity_mm_ms"],e.info.ion_sigma_t_mm,tc["ion_graph_k"],tc["ion_edge_sigma_scale"]),
          ("etrack","pe",e.etrack,dc["EL_pitch_mm"],dc["electron_bucket_us"]*rc["electrons"]["rebin_factor"]*dc["electron_drift_velocity_mm_us"],e.info.electron_sigma_t_mm,tc["electron_graph_k"],tc["edge_sigma_scale"]))
            cut=Float32(tc[arm=="itrack" ? "ion_hit_cut_keV" : "electron_hit_cut_keV"])
            keep=h.energy_keV .> cut
            h=SavedHits(h.x_mm[keep],h.y_mm[keep],h.z_mm[keep],h.energy_keV[keep])
            h0=ReconstructedHits(h.x_mm,h.y_mm,h.z_mm,h.energy_keV,h.energy_keV)
            placement=EventPlacement(id,e.info.ion_drift_cm,e.info.electron_drift_cm)
            z0=10Float64(placement.ion_drift_cm)
            if arm=="itrack"
                t=ion_readout_time(e.info.frame,placement,Float32(dc["chamber_length_mm"]),Float32(dc["ion_drift_velocity_mm_ms"]),Float32(dc["ion_bucket_ms"]))
                zs=Float32[(Float64(t.start_ms)+(i-.5)*t.bucket_ms)*Float32(dc["ion_drift_velocity_mm_ms"])-z0 for i in 1:t.bucket_count]
            else
                t=electron_readout_time(e.info.frame,placement,Float32(dc["chamber_length_mm"]),Float32(dc["electron_drift_velocity_mm_us"]),Float32(dc["electron_bucket_us"]))
                width=Float32(t.bucket_us*rc["electrons"]["rebin_factor"])
                last_us=Float32(Float64(t.start_us)+Float64(t.bucket_count)*t.bucket_us)
                zs=Float32[]
                for i in 1:cld(t.bucket_count,rc["electrons"]["rebin_factor"])
                    lo=Float64(t.start_us)+(i-1)*width;hi=min(lo+width,Float64(last_us))
                    push!(zs,Float32(dc["chamber_length_mm"])-(lo+hi)/2*Float32(dc["electron_drift_velocity_mm_us"])-z0)
                end
            end
            ex[arm]=hit_projections(h0,pitch,dz;z_centres=zs)
            df=EventAnalysis.DataFrame(x=h.x_mm,y=h.y_mm,z=h.z_mm,energy=h.energy_keV./1000f0)
            tracks=make_tracks(df;method="kNN",k=Int(k),max_distance_mm=Float32(mult*sigma),energy_threshold_kev=0f0)
            paths=Dict{String,Any}()
            for (n,track) in enumerate(tracks)
                u,w,indices,_,_=EventAnalysis.find_extremes_mst_diameter(track,EventAnalysis.extract_coords(track));v=track.voxels
                paths[string(n)]=Dict("x_mm"=>v.x[indices],"y_mm"=>v.y[indices],"z_mm"=>v.z[indices],"energy_MeV"=>sum(v.energy),"hit_count"=>size(v,1))
            end
            ex[arm]["backbones"]=paths
            blobs=Dict{String,Any}()
            if !isempty(tracks)
                track=tracks[argmax([size(t.voxels,1) for t in tracks])]
                u,w,_,_,_=EventAnalysis.find_extremes_mst_diameter(track,EventAnalysis.extract_coords(track))
                v=track.voxels;walk=(extremes=(v[u,:],v[w,:]),)
                for scale in (2,3)
                    radius=Float64(scale*sigma)
                    b=EventAnalysis.energy_in_spheres_around_extremes(track,walk,radius)
                    overlap=sum(v.energy[i] for i in 1:size(v,1) if
                        sum((v[i,axis]-b.start_center[j])^2 for (j,axis) in enumerate((:x,:y,:z)))<=radius^2 &&
                        sum((v[i,axis]-b.end_center[j])^2 for (j,axis) in enumerate((:x,:y,:z)))<=radius^2;init=0.0)
                    blobs[string(scale)]=Dict("radius_mm"=>radius,"centres_mm"=>hcat(collect(b.blob1_center),collect(b.blob2_center)),
                        "energy_keV"=>1000 .* [b.blob1_energy,b.blob2_energy],"overlap_keV"=>1000overlap)
                end
            end
            ex[arm]["blobs"]=blobs
            hh[key]=control_histogram(h.energy_keV,collect(0.:.25:200.))
        end
        examples[string(id)]=ex;histograms[string(id)]=hh
      end
      budget=16*1024^2+3Base.summarysize(examples)+3Base.summarysize(histograms)
      GuardedFiles.space_check(out,budget,1024^3);tmp=target*".tmp"
      try
        h5open(tmp,"w") do f
          f["schema"]="itaca.hit_presentation.v2";f["hit_file"]=relpath(path,out);f["hit_file_sha256"]=bytes2hex(open(sha256,path))
          f["card"]=sprint(io->TOML.print(io,rc));f["cards/detector_card"]=read(r.file["configuration/cards/detector_card"],String)
          f["events/event_number"]=ids
          for (targetkey,sourcekey) in (("ion_energy_keV","ion_energy_keV"),("electron_energy_keV","electron_energy_keV"),("invalid_ion_buckets","invalid_ion_buckets"),("invalid_pe_buckets","invalid_electron_buckets"))
            f["events/$targetkey"]=read(r.file["events/$sourcekey"])
          end
          summary=Dict{String,Any}()
          for key in ("ion_energy_keV","electron_energy_keV")
            values=filter(isfinite,read(r.file["events/$key"]));summary[key]=control_histogram(values,histogram_edges(values,80))
          end
          write_control_tree(create_group(f,"track_summary"),summary)
          write_control_tree(create_group(f,"track_examples"),examples)
          write_control_tree(create_group(f,"hit_histograms"),histograms)
          f["postprocessing_s"]=time()-started
          flush(f);GuardedFiles.check_size(tmp,budget)
        end
        Base.rename(tmp,target)
      catch;isfile(tmp)&&rm(tmp);rethrow();end
    finally;close_hit_input(r);end
    println("PRESENTATION $target time_s=$(time()-started)");target
end
if abspath(PROGRAM_FILE)==@__FILE__
 prepare_hit_presentation(abspath(ARGS[1]);graph_card=length(ARGS)>1 ? abspath(ARGS[2]) : nothing,output_dir=length(ARGS)>2 ? abspath(ARGS[3]) : nothing)
end
