#!/usr/bin/env julia
# Run bb0nu framed detector controls: summaries and example projections.
# Optional second CLI argument limits input events.
using HDF5,TOML,Random,SHA,Statistics
include("../../src/event_analysis/EventAnalysis.jl")
using .EventAnalysis
include("../../src/detector_response/DetectorResponse.jl")
using .DetectorResponse
include("../../src/ion_signal/guarded_files.jl")
const ROOT=normpath(joinpath(@__DIR__,"../.."))
seedfor(seed,id,part)=parse(UInt64,bytes2hex(sha256("framed-v1|$seed|$id|$part"))[1:16];base=16)

"""Resolve event diffusion and optics from card values; grids have independent fixed random offsets."""
function event_settings(c,p)
    F(k)=Float32(c[k]);rng(tag)=Xoshiro(seedfor(c["seed"],p.event_id,tag))
    iong=ReadoutGrid(F("ion_pitch_mm"),sample_grid_offset_mm(rng("ion_grid"),F("ion_pitch_mm")))
    elg=ReadoutGrid(F("EL_pitch_mm"),sample_grid_offset_mm(rng("electron_grid"),F("EL_pitch_mm")))
    si=sigma_t_ion_mm(F("temperature_K"),p.ion_drift_cm,F("drift_field_V_cm"))
    scale=sqrt(p.electron_drift_cm/F("pressure_bar"))
    a=elm_light_acceptance(F("hole_diameter_mm")/2,F("hole_depth_mm");
        tpb_efficiency=F("tpb_efficiency"),forward_fraction=F("forward_fraction"),
        guide_index=F("guide_index"),gap_index=F("gap_index"),
        propagation_transmission=F("propagation_transmission"),coupling_transmission=F("coupling_transmission"),electrode_transmission=F("electrode_transmission"))
    y=el_photon_yield(F("EL_gain_photons_e_kV"),F("EL_threshold_kV_cm_bar"),F("EL_reduced_field_kV_cm_bar"),F("pressure_bar"),F("hole_depth_mm"))
    c["gas_mixture"]=="Xe_85_He_15" || error("This run card specifies Xe–15% He")
    DetectorResponseSettings(F("chamber_length_mm"),F("W_eV"),F("fano_factor"),Xe_85_He_15,
      iong,elg,F("ion_bucket_ms"),F("electron_bucket_us"),IonTransportSettings(F("ion_drift_velocity_mm_ms"),si,si),
      ElectronTransportSettings(F("electron_drift_velocity_mm_us"),F("electron_dT_mm_sqrtbar_sqrtcm")*scale,F("electron_dL_mm_sqrtbar_sqrtcm")*scale),
      F("electron_lifetime_ms"),y,F("EL_variance_ratio"),a,Symbol(c["photon_wavelength"]))
end
fraction(n,d)=d==0 ? 0. : n/d
"""Per-event controls, reduced in one pass without sorting events or pulses."""
function metrics(r,v)
    a=vector_statistics(v);l=r.losses
    dt=hasproperty(v.time,:bucket_ms) ? v.time.bucket_ms : v.time.bucket_us
    (;collected=a.total,mean=a.mean,maximum=a.maximum,positive_cells=a.positive_cells,
      acquired_channels=length(v.channels.ix),active_channels=a.channels_with_signal,
      active_fraction=fraction(a.channels_with_signal,length(v.channels.ix)),
      total=l.total,loss_xy=l.excluded_xy,loss_time=l.excluded_time,loss=l.excluded_total,
      loss_xy_fraction=fraction(l.excluded_xy,l.total),loss_time_fraction=fraction(l.excluded_time,l.total),loss_fraction=fraction(l.excluded_total,l.total),
      outside_channels=l.outside_channels_with_signal,signal_channels=l.channels_with_signal,
      outside_channel_fraction=fraction(l.outside_channels_with_signal,l.channels_with_signal),
      buckets=Int(v.time.bucket_count),duration=Float64(dt)*v.time.bucket_count,
      dense_bytes=sizeof(v.counts))
end
"""Small XY and XZ count projections for an example; also project excluded cells."""
function projections(r,v,g,velocity,ion;rebin::Int=1)
    rebin>0 || throw(ArgumentError("Positive display rebin required"))
    display_bucket(k)=fld(k-1,rebin)+1
    a,dt,nb=DetectorResponse._time(v.time)
    first=round(Int64,Float64(a)/dt)
    xy=Dict{Tuple{Int64,Int64},Int64}();xz=Dict{Tuple{Int64,Int64},Int64}()
    for j in axes(v.counts,2),k in axes(v.counts,1)
        c=v.counts[k,j];c==0&&continue
        x=Int64(v.channels.ix[j]);y=Int64(v.channels.iy[j])
        xy[(x,y)]=get(xy,(x,y),0)+c;xz[(x,display_bucket(k))]=get(xz,(x,display_bucket(k)),0)+c
    end
    outxy=Dict{Tuple{Int64,Int64},Int64}();outxz=Dict{Tuple{Int64,Int64},Int64}()
    for i in eachindex(r.excluded.count)
        x=r.excluded.pixel_ix[i];y=r.excluded.pixel_iy[i];k=r.excluded.bucket[i];c=r.excluded.count[i]
        outxy[(x,y)]=get(outxy,(x,y),0)+c;outxz[(x,display_bucket(k))]=get(outxz,(x,display_bucket(k)),0)+c
    end
    result=Dict{String,Any}()
    for (name,d,isxz) in (("xy",xy,false),("xz",xz,true),("outside_xy",outxy,false),("outside_xz",outxz,true))
        keys_=sort!(collect(keys(d)))
        result[name*"_x"]=Float32[g.edge_mm[1]+(q[1]+.5f0)*g.pitch_mm for q in keys_]
        if isxz
            midpoint(k)=begin
                low=(k-1)*rebin
                high=0<=low<nb ? min(k*rebin,nb) : k*rebin
                first+(low+high)/2
            end
            result[name*"_y"]=Float32[(ion ? (midpoint(q[2])*dt*velocity-10r.placement.ion_drift_cm) :
                (10r.placement.electron_drift_cm-midpoint(q[2])*dt*velocity)) for q in keys_]
        else
            result[name*"_y"]=Float32[g.edge_mm[2]+(q[2]+.5f0)*g.pitch_mm for q in keys_]
        end
        result[name*"_count"]=Int64[d[q] for q in keys_]
    end
    @assert sum(values(xy);init=0)==sum(values(xz);init=0)==sum(v.counts)
    @assert sum(values(outxy);init=0)==sum(values(outxz);init=0)==sum(r.excluded.count)
    result["display_rebin"]=rebin
    result["display_z_bin_mm"]=Float32(rebin*dt*velocity)
    result["acquired_xy_bounds"]=Float32[minimum(v.channels.x_mm)-g.pitch_mm/2,maximum(v.channels.x_mm)+g.pitch_mm/2,
        minimum(v.channels.y_mm)-g.pitch_mm/2,maximum(v.channels.y_mm)+g.pitch_mm/2]
    result
end
function columns(g,rows)
    isempty(rows)&&return
    for k in keys(first(rows));g[string(k)]=[r[k] for r in rows];end
end
"""Save histogram bins and numerics in Julia; presentation reads them unchanged."""
function histograms(g,rows,nb)
    isempty(rows)&&return
    for k in keys(first(rows))
        x=Float64[r[k] for r in rows];all(isfinite,x)||error("Nonfinite metric $k")
        low,high=extrema(x);high==low&&(high=low+max(1.,abs(low)*.01))
        edges=collect(range(low,high;length=nb+1));counts=zeros(Int64,nb)
        for v in x;counts[clamp(searchsortedlast(edges,v),1,nb)]+=1;end
        q=create_group(g,string(k));q["edges"]=edges;q["counts"]=counts;q["entries"]=length(x)
        q["mean"]=mean(x);q["rms"]=sqrt(mean(abs2,x));q["std"]=std(x;corrected=false)
    end
end
function main(card,limit=0)
    started=time();c=TOML.parsefile(card);limit>0&&(c["max_events"]=limit;c["run_label"]*="_first$(limit)")
    out=joinpath(ROOT,"results/readout",c["run_label"]);mkpath(joinpath(out,"detector_response"))
    target=joinpath(out,"detector_response/controls.h5");isfile(target)&&error("Output exists: $target")
    hp=joinpath(ROOT,c["hits"]);pp=expanduser(c["particles"])
    source_sha=open(sha256,hp)|>bytes2hex;particle_sha=open(sha256,pp)|>bytes2hex
    ids=h5open(hp,"r") do f
        read(f["metadata/source_sha256"])==particle_sha||error("Particle source mismatch")
        read(f["events/event_id"])
    end
    c["max_events"]>0&&(ids=ids[1:min(length(ids),c["max_events"])])
    parts=read_mc_particles(pp,ids)
    trigger=FastTriggerSettings(Float32(c["chamber_length_mm"]),Float32(c["vertex_min_mm"]),Float32(c["vertex_max_mm"]),ntuple(_->Float32(c["frame_margin_mm"]),3),Float32(c["gamma_energy_cut_keV"]))
    rng=Xoshiro(c["seed"]);stats=FastTriggerStatistics();ev=NamedTuple[];im=NamedTuple[];em=NamedTuple[];allrows=NamedTuple[]
    examples=Dict{Int64,Any}();worst=Dict{Int64,Any}();checks=Dict{Int64,Any}();selection_s=0.;response_s=0.;per_event=Float64[]
    h5open(hp,"r") do f
        for (row,id) in enumerate(ids)
            t=time();h,info=read_mc_event(f,id;event_row=row);h=center_hits(h,info)
            p=place_event(rng,id,trigger.chamber_length_mm);o=mc_fast_trigger(h,info,parts[id],p,trigger)
            stats=update_statistics(stats,o);push!(allrows,(event_number=id,decision=Int(o.decision),ion_drift_cm=p.ion_drift_cm,electron_drift_cm=p.electron_drift_cm))
            selection_s+=time()-t
            o.decision==trigger_accepted||continue
            t=time();frame=acquisition_frame(o.frame,p,trigger.chamber_length_mm);selected=hits_in_frame(h,frame)
            gm=gamma_map(parts[id],info)
            main_outside=sum(h.energy_keV[j] for j in eachindex(h.x_mm) if gm[h.particle_id[j]]==0 && !contains(frame,(h.x_mm[j],h.y_mm[j],h.z_mm[j]));init=0f0)
            iszero(main_outside)||error("Main-track energy outside frame in event $id")
            pairs=sample_fano_pairs(Xoshiro(seedfor(c["seed"],id,"pairs")),selected,Float32(c["W_eV"]),Float32(c["fano_factor"]))
            s=event_settings(c,p)
            ri=sim_ion_detector_response(Xoshiro(seedfor(c["seed"],id,"ions")),selected,pairs,p,frame,s;max_vector_bytes=c["max_vector_bytes"])
            re=sim_electron_detector_response(Xoshiro(seedfor(c["seed"],id,"electrons")),selected,pairs,p,frame,s;max_vector_bytes=c["max_vector_bytes"])
            elapsed=time()-t;response_s+=elapsed;push!(per_event,elapsed)
            push!(im,metrics(ri,ri.ions));push!(em,metrics(re,re.photons))
            inside=sum(selected.energy_keV);outside=sum(h.energy_keV[i] for i in eachindex(h.x_mm) if !contains(frame,(h.x_mm[i],h.y_mm[i],h.z_mm[i]));init=0f0)
            push!(ev,(event_number=id,energy_inframe_keV=inside,energy_outside_frame_keV=outside,main_track_energy_outside_frame_keV=main_outside,e_ion_pairs=sum(pairs),
                e_lost_attachment=re.e_lost_attachment,EL_photons=re.EL_photons,uncollected_photons=re.uncollected_photons,
                ion_drift_cm=p.ion_drift_cm,electron_drift_cm=p.electron_drift_cm,response_s=elapsed,
                frame_x_min=frame.lower_mm[1],frame_x_max=frame.upper_mm[1],frame_y_min=frame.lower_mm[2],frame_y_max=frame.upper_mm[2],frame_z_min=frame.lower_mm[3],frame_z_max=frame.upper_mm[3],
                ion_phase_x=s.ion_grid.edge_mm[1],ion_phase_y=s.ion_grid.edge_mm[2],electron_phase_x=s.electron_grid.edge_mm[1],electron_phase_y=s.electron_grid.edge_mm[2]))
            score=max(last(im).loss_fraction,last(em).loss_fraction)
            check_event=get(c,"sample","")=="single_electrons" && id==604
            if check_event || length(examples)<c["example_events"] || length(worst)<3 || score>minimum(v["score"] for v in values(worst))
                d=Dict("frame_low"=>collect(frame.lower_mm),"frame_high"=>collect(frame.upper_mm),"score"=>score,
                    "ions"=>projections(ri,ri.ions,s.ion_grid,s.ion_transport.ion_drift_velocity_mm_ms,true),
                    "photons"=>projections(re,re.photons,s.electron_grid,s.electron_transport.electron_drift_velocity_mm_us,false;rebin=Int(get(c,"electron_display_rebin",1))))
                length(examples)<c["example_events"]&&(examples[id]=d)
                check_event&&(checks[id]=d)
                worst[id]=d
                length(worst)>3&&delete!(worst,first(sort!(collect(keys(worst));by=k->worst[k]["score"])))
            end
            if length(ev)%100==0
                println("Accepted $(length(ev)); input $row/$(length(ids)); response $(round(response_s;digits=2)) s");flush(stdout)
                GC.gc(false)
            end
        end
    end
    # Budget: scalar columns/histograms plus projected examples, with HDF5 metadata allowance.
    budget=16*1024^2+length(ids)*4096
    GuardedFiles.space_check(out,budget,1024^3)
    tmp=target*".tmp"
    try
      h5open(tmp,"w") do f
        f["schema"]="itaca.framed_detector_controls.v1";f["card"]=sprint(io->TOML.print(io,c))
        f["source_sha256"]=source_sha;f["particles_sha256"]=particle_sha
        files=[@__FILE__]
        for dir in ("src/detector_response","src/event_analysis")
          for (root,_,names) in walkdir(joinpath(ROOT,dir)),name in names
            endswith(name,".jl")&&push!(files,joinpath(root,name))
          end
        end
        f["code_files"]=relpath.(files,ROOT);f["code_sha256"]=[bytes2hex(open(sha256,p)) for p in files]
        columns(create_group(f,"selection/events"),allrows)
        for k in fieldnames(FastTriggerStatistics);f["selection/"*string(k)]=getfield(stats,k);end
        columns(create_group(f,"events"),ev);columns(create_group(f,"ions/events"),im);columns(create_group(f,"photons/events"),em)
        histograms(create_group(f,"event_histograms"),ev,c["histogram_bins"])
        histograms(create_group(f,"ions/histograms"),im,c["histogram_bins"]);histograms(create_group(f,"photons/histograms"),em,c["histogram_bins"])
        f["timing/selection_s"]=selection_s;f["timing/response_s"]=response_s
        f["timing/response_mean_s"]=mean(per_event);f["timing/response_mean_excluding_first_s"]=length(per_event)>1 ? mean(per_event[2:end]) : mean(per_event)
        f["timing/response_median_s"]=median(per_event);f["timing/total_s"]=time()-started
        for (label,rows) in (("ions",im),("photons",em))
            for (name,field) in (("xy","loss_xy"),("time","loss_time"),("any","loss"))
                f["$label/summary/affected_$name"]=count(r->r[Symbol(field)]>0,rows)
                f["$label/summary/fraction_$name"]=sum(r[Symbol(field)] for r in rows)/sum(r.total for r in rows)
            end
        end
        for (label,ex) in (("examples",examples),("largest_losses",worst),("frame_checks",checks))
          for (id,d) in ex
            g=create_group(f,"$label/$id")
            for k in ("frame_low","frame_high","score");g[k]=d[k];end
            for carrier in ("ions","photons"), (k,v) in d[carrier];g["$carrier/$k"]=v;end
          end
        end
        flush(f);GuardedFiles.check_size(tmp,budget)
      end
      mv(tmp,target)
    catch
      isfile(tmp)&&rm(tmp);rethrow()
    end
    open(joinpath(out,"run.toml"),"w") do io;TOML.print(io,c);end
    println("FINISHED: $(stats.input_events) input; $(stats.accepted_events) accepted; $(round(response_s;digits=2)) s response; $(round(time()-started;digits=2)) s total; $target")
end
if abspath(PROGRAM_FILE)==@__FILE__
 main(isempty(ARGS) ? joinpath(ROOT,"config/readout/bb0nu_detector_response_3mm.toml") : abspath(ARGS[1]),length(ARGS)>1 ? parse(Int,ARGS[2]) : 0)
end
