#!/usr/bin/env julia
include("run_framed_reconstruction.jl")
include("../../src/vh1/VH1Readout.jl")
using .VH1Readout
include("frame_noise_hits.jl")

"""Parse the explicit diagnostic selection. No diagnostics are constructed in `none` mode."""
function hit_control_settings(card)
    c=get(card,"controls",Dict("mode"=>"none"));mode=c["mode"]
    mode in ("none","debug")||throw(ArgumentError("controls.mode must be none or debug"))
    mode=="none"&&return (mode=mode,limit=0,diagnostics=String[])
    names=String.(get(c,"diagnostics",String[]));allowed=["noise","hit_cut","graph_scan","electronics_examples"]
    all(in(allowed),names)||throw(ArgumentError("Unknown diagnostic"))
    limit=Int(get(c,"event_limit",0));limit>0||throw(ArgumentError("Debug mode requires positive event_limit"))
    (mode=mode,limit=limit,diagnostics=names)
end

"""Calculate only requested debug results, using the existing diagnostic functions."""
function hit_debug(c;ri,pe,ions,photons,adc,asic,p,ds,totals,s,card)
    result=Dict{String,Any}()
    for name in c.diagnostics
        if name=="electronics_examples"
            result[name]=Dict("adc_counts"=>adc.counts[:,1:min(3,size(adc.counts,2))],
              "asic_pe"=>asic.response_pe[:,1:min(3,size(asic.response_pe,2))])
            continue
        end
        ei=grouped_input(photons.photons.counts,card["electrons"]["rebin_factor"]).*s.sipm.sipm.pde
        pairs=(("ions",ri.ions,ri.valid,ions.ions.counts,ri.channels,totals.ion_keV_per_count,ri,ds.ion_transport.sigma_t_mm),
          ("pe",pe.pe,pe.valid,ei,pe.channels,totals.attachment_correction*totals.pe_keV_per_count,pe,ds.electron_transport.sigma_t_mm))
        selected=Dict{String,Any}()
        for (arm,v,valid,input,channels,scale,signal,sigma) in pairs
            cut=card["hit_output"][arm=="ions" ? "ion_hit_cut_keV" : "electron_hit_cut_keV"]
            if name=="noise"
                selected[arm]=channel_histograms(v,valid,input,100)
            elseif name=="hit_cut"
                selected[arm]=Dict(string(k)=>value for (k,value) in Base.pairs(hit_cut_statistics(v,valid,input,scale,cut)))
            elseif name=="graph_scan"
                tc=card["track_analysis"];hits=reconstruct_hits(signal,p,ds;charge_correction=arm=="ions" ? 1. : totals.attachment_correction,
                  keV_per_count=arm=="ions" ? totals.ion_keV_per_count : totals.pe_keV_per_count)
                ks=get(tc,"graph_k_scan",[tc[arm=="ions" ? "ion_graph_k" : "electron_graph_k"]])
                selected[arm]=Dict("k_$k"=>hit_cut_topology(v,valid,input,channels,hits,scale,cut,sigma,k,
                  arm=="ions" ? tc["ion_edge_sigma_scale"] : tc["edge_sigma_scale"];compute_endpoints=false) for k in ks)
            end
        end
        result[name]=selected
    end
    result
end

"""Run accepted MC events through detector, electronics and reconstruction.
Write compact event totals and retained hits. In controls.mode=none no control
histogram, graph or example routine is called. Optional on_event sees one event
for validation and must not alter it. Timings are printed to the run log.
"""
function run_hit_production(cardpath;on_event=event->nothing,batch=nothing)
    BLAS.set_num_threads(1)
    r=resolve_run(cardpath;root=ROOT);rc=r.run;batch===nothing||select_batch!(rc,batch);stage=rc["hit_production"];response=rc["response"];ec=electronics_card(r)
    debug=hit_control_settings(rc)
    use_vh1=haskey(response,"vh1_card")
    use_vh1&&debug.mode!="none"&&throw(ArgumentError("VH1 hit production requires controls.mode=none"))
    hp=input_path(ROOT,rc["hits"]);sha=bytes2hex(open(sha256,hp))
    source=joinpath(ROOT,get(stage,"selection_controls",joinpath("results/readout",rc["run_label"],"detector_response/controls.h5")))
    out=joinpath(ROOT,"results/readout",rc["run_label"],stage["output_subdir"])
    write_configuration(r,ROOT;stage="hits")
    c=rc["hit_output"];settings=HitOutputSettings(c["ion_hit_cut_keV"],c["electron_hit_cut_keV"],c["compression_level"],c["chunk_hits"],c["max_output_bytes"])
    ids,positions,frames,locations=h5open(source,"r") do f
        read(f["source_sha256"],String)==sha||error("MC source changed")
        g=f["events"];ids=read(g["event_number"]);stage["accepted_events"]>0&&(ids=ids[1:stage["accepted_events"]])   # 0: all accepted events
        ps=[EventPlacement(id,read(g["ion_drift_cm"])[j],read(g["electron_drift_cm"])[j]) for (j,id) in enumerate(ids)]
        fs=[EventFrame(mc_vertex_frame,ntuple(k->read(g["frame_$("xyz"[k])_min"])[j],3),ntuple(k->read(g["frame_$("xyz"[k])_max"])[j],3)) for j in eachindex(ids)]
        # Fixed-frame selection: the vertex in detector coordinates and the fixed frame of each event
        ls=haskey(g,"fixed_frame_x_min") ? [HitEventLocation(ntuple(k->read(g["vertex_$("xyz"[k])_mm"])[j],3),
            EventFrame(mc_vertex_frame,ntuple(k->read(g["fixed_frame_$("xyz"[k])_min"])[j],3),ntuple(k->read(g["fixed_frame_$("xyz"[k])_max"])[j],3))) for j in eachindex(ids)] : nothing
        ids,ps,fs,ls
    end
    cards=configuration_texts(r,ROOT)
    hashes=Dict("source"=>sha,"selection"=>bytes2hex(open(sha256,source)),"calibration"=>bytes2hex(open(sha256,joinpath(ROOT,ec["nonlinear_calibration"]))))
    for (name,h) in r.set.sha256;hashes["detector/"*name]=h;end
    if use_vh1
        hashes["vh1_card"]=bytes2hex(open(sha256,joinpath(ROOT,response["vh1_card"])))
        hashes["vh1_code"]=bytes2hex(open(sha256,joinpath(ROOT,"src/vh1/VH1Readout.jl")))
    end
    files=vcat([@__FILE__,joinpath(@__DIR__,"framed_electronics_settings.jl"),joinpath(@__DIR__,"run_framed_detector_response.jl")],
      [joinpath(d,n) for sub in ("src/event_analysis","src/detector_response","src/electronics","src/reconstruction","src/ion_signal") for (d,_,ns) in walkdir(joinpath(ROOT,sub)) for n in ns if endswith(n,".jl")])
    for file in files;hashes[replace(relpath(file,ROOT),"/"=>"__")]=bytes2hex(open(sha256,file));end
    metadata=HitFileMetadata(rc["run_label"],get(rc,"sample","bb0nu"),Xe_85_He_15,cards,hashes)
    setup=@timed prepare_electronics_run(ec);prepared=setup.value;workspace=ER.CSAWorkspace()
    vh=use_vh1 ? load_vh1_settings(ROOT,joinpath(ROOT,response["vh1_card"]);minimum_samples=get(response,"vh1_minimum_samples",32768)) : nothing
    c=rc["ions"];is=IonReconstructionSettings(c["baseline_ms"],c["average_ms"],c["nonlinear_below_V"],c["max_iterations"],c["convergence_ions"])
    es=ElectronReconstructionSettings(prepared.sipm.sipm,prepared.sipm.illuminated_microcells,rc["electrons"]["subtract_mean_dark"])
    located=locations!==nothing
    w=open_hit_output(joinpath(out,"hits.h5"),metadata,settings;located)
    timing=zeros(6);debug_results=debug.mode=="debug" ? Dict{String,Any}() : nothing
    started=time()
    try
        h5open(hp,"r") do input
            rowmap=Dict(id=>j for (j,id) in enumerate(read(input["events/event_id"])))
            for (j,id) in enumerate(ids)
                t=@elapsed begin
                    h,info=read_mc_event(input,id;event_row=rowmap[id]);p=positions[j];frame=frames[j]
                    h=hits_in_frame(center_hits(h,info),frame);ds=event_settings(r,p)
                    rng(tag)=Xoshiro(seedfor(rc["seed"],id,tag));erng(tag)=Xoshiro(seedfor(ec["seed"],id,tag))
                    pairs=sample_fano_pairs(rng("pairs"),h,Float32(r.detector["energy"]["W_eV"]),Float32(r.detector["energy"]["fano_factor"]))
                    ions=sim_ion_detector_response(rng("ions"),h,pairs,p,frame,ds;max_vector_bytes=rc["max_vector_bytes"])
                    photons=sim_electron_detector_response(rng("electrons"),h,pairs,p,frame,ds;max_vector_bytes=rc["max_vector_bytes"])
                end;timing[1]+=t
                t=@elapsed begin
                    s=electronics_settings(prepared,ions.ions.time)
                    if use_vh1
                        adc=sim_vh1_electronics(erng("vh1_noise"),erng("vh1_adc"),ions,vh)
                    else
                        adc,_=ER.sim_ion_electronics_prepared(erng("csa_noise"),erng("adc_noise"),ions,s.ion,workspace)
                    end
                end;timing[2]+=t
                t=@elapsed asic,_=sim_sipm_readout(erng("sipm"),erng("asic"),photons,s.sipm);timing[3]+=t
                t=@elapsed begin
                    if use_vh1
                        plan=prepare_vh1_reconstruction(adc.time,ions.ions.time,vh)
                        ri=reco_vh1_ions(adc,plan)
                    else
                        cal=IonReconstructionCalibration(s.ion.adc,s.ion.csa.current_pulse,s.ion.csa.linear_response,s.ion.csa.nonlinear_response,s.ion.csa.noise)
                        plan=SR.prepare_ion_reconstruction(adc.time,ions.ions.time,cal,is)
                        ri=SR.reconstruct_ions(adc,plan)
                    end
                    fine=reco_electron_signal(asic,es)
                    totals=event_charge_estimates(ri,fine,p,ds,s.sipm.sipm.pde)
                    pe=rebin_electron_signal(fine,Int64(rc["electrons"]["rebin_factor"]))
                end;timing[4]+=t
                t=@elapsed begin
                    ih=reconstruct_hits(ri,p,ds;keV_per_count=totals.ion_keV_per_count,minimum_energy_keV=settings.ion_hit_cut_keV)
                    if located&&use_vh1   # noise hits of the frame cells outside the simulated pads and buckets
                        sigma=sqrt(sum(plan.ion_covariance[k,k] for k in axes(plan.ion_covariance,1))/size(plan.ion_covariance,1))
                        ih=join_hits(ih,frame_noise_hits(erng("frame_noise"),locations[j].fixed_frame,ri,p,ds,sigma,totals.ion_keV_per_count,settings.ion_hit_cut_keV))
                    end
                    eh=reconstruct_hits(pe,p,ds;charge_correction=totals.attachment_correction,keV_per_count=totals.pe_keV_per_count,minimum_energy_keV=settings.electron_hit_cut_keV)
                    energy=EventEnergy(totals.ions_total,totals.pe_total,totals.ion_energy_keV,totals.electron_energy_keV,totals.invalid_ion_buckets,totals.invalid_pe_buckets)
                    info=HitEventInfo(id,p.ion_drift_cm,p.electron_drift_cm,frame,ds.ion_transport.sigma_t_mm,ds.electron_transport.sigma_t_mm)
                    event=HitEvent(info,energy,saved_hits(ih),saved_hits(eh))
                end;timing[5]+=t
                timing[6]+=@elapsed write_hit_event(w,event;location=located ? locations[j] : nothing)
                on_event(event)
                if debug.mode=="debug"&&j<=debug.limit
                    debug_results[string(id)]=hit_debug(debug;ri,pe,ions,photons,adc,asic,p,ds,totals,s,card=rc)
                end
                println("Event $id ($j/$(length(ids))): itrack $(length(ih.x_mm)), etrack $(length(eh.x_mm)) hits");flush(stdout)
            end
        end
        target=finish_hit_output(w)
        if debug.mode=="debug"
            dir=joinpath(out,"debug");mkpath(dir);path=joinpath(dir,"controls.h5");ispath(path)&&error("Debug output exists")
            budget=1024^2+3Base.summarysize(debug_results);GuardedFiles.space_check(dir,budget,1024^3)
            try
                h5open(path*".tmp","w") do f
                    write_control_tree(f,debug_results);f["card"]=read(cardpath,String);flush(f);GuardedFiles.check_size(path*".tmp",budget)
                end
                Base.rename(path*".tmp",path)
            catch;isfile(path*".tmp")&&rm(path*".tmp");rethrow();end
        end
        println("TIMING setup_s=$(setup.time) event_loop_s=$(time()-started) events=$(length(ids))")
        for (name,t) in zip(("detector",use_vh1 ? "VH1_ADC" : "CSA_ADC","SiPM_ASIC","reconstruction","hits","write"),timing)
            println("TIMING $name total_s=$t per_event_s=$(t/length(ids))")
        end
        println("OUTPUT $target bytes=$(filesize(target))")
        target
    catch
        abort_hit_output(w);rethrow()
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    batch,args=batch_argument(ARGS)
    length(args)==1||error("Usage: run_hit_production.jl run_card.toml [--batch=N]")
    run_hit_production(abspath(args[1]);batch)
end
