module IonCSASummary
using HDF5, Statistics, TOML, Test, SHA

"""Count a saved diagnostic distribution with explicit under/overflow ledgers."""
function histogram(values,edges)
 bins=zeros(Int64,length(edges)-1);under=0;over=0
 for x in values
  j=searchsortedlast(edges,x)
  if j<1;under+=1;elseif j>=length(edges);over+=1;else;bins[j]+=1;end
 end
 @assert sum(bins)+under+over==length(values)
 (;bins,under,over)
end
function savehist(g,name,values,edges)
 h=create_group(g,name);r=histogram(values,edges);h["edges"]=edges;h["count"]=r.bins
 attributes(h)["underflow"]=r.under;attributes(h)["overflow"]=r.over
end

"""Reduce complete ion-response files to per-event quality, charge and compression diagnostics."""
function run(C,root,card,dir)
 s=C.settings(root,card);@assert s.dir==dir
 run_context(C,s,dir,card)
end

"""Reduce a resolved response context, including versioned reconstructions of saved data."""
function run_context(C,s,dir,card)
 destination=joinpath(dir,"summary.h5")
 @testset "Ion CSA pilot summary and event accounting" begin
 h5open(destination*".tmp","w") do out
  attributes(out)["study_id"]=s.studyid;attributes(out)["model_id"]=s.modelid
  attributes(out)["summary_source_sha256"]=C.sha(@__FILE__)
  attributes(out)["run_card"]=card
  create_group(out,"calibration");create_group(out,"validation")
  h5open(joinpath(dir,"calibration.h5"),"r") do cal
   for name in keys(cal)
    C.HDF5.copy_object(cal[name],out,"calibration/"*name)
    if cal[name] isa HDF5.Group && !haskey(out["calibration/"*name],"linear_reconstruction_kernel")
     cg=out["calibration/"*name];a=C.attrs(cg);v=read(cg["average_V_per_ion"])
     taps=[v[1],v[2]-a["alpha"]*v[1],v[3]-a["alpha"]*v[2]]./a["gain_V_per_ion"]
     cg["linear_reconstruction_kernel"]=taps
     attributes(cg)["predicted_linear_total_bias_fraction"]=sum(taps)-1
    end
   end
  end
  h5open(joinpath(dir,"solver_validation.h5"),"r") do h
   C.HDF5.copy_object(h["maximum_error_ENC"],out,"validation/refinement_error_ENC")
  end
  v=TOML.parsefile(joinpath(dir,"fast_path_validation.toml"))
  for (key,value) in v;attributes(out["validation"])[key]=value;end
  names=String[]
  for source in s.sources
   name=source.sample*"_p"*string(source.pitch);push!(names,name);g=create_group(out,name)
   h5open(source.response,"r") do h
    @test Bool(read(h["metadata/complete"]))
    @test read(h["metadata/selected_event_ids"])==source.selected
    for key in ("sample","frame_s","pitch_mm","pads","direct_spice_pads","negative_estimates","response_seconds","input_ion_total","maximum_storage_error")
     attributes(g)[key]=read(h["metadata/"*key])
    end
    attributes(g)["response_path"]=source.response;attributes(g)["response_sha256"]=C.sha(source.response)
    enc=read(attributes(h["calibration"])["ENC_e_rms"]);attributes(g)["ENC_e_rms"]=enc
    expected=Float64[];truthsnr=Float64[];eventsnr=Float64[];linearbias=Float64[];meanbias=Float64[];noisybias=Float64[];deficitfraction=Float64[]
    initial=Int64[];cells=Int64[];padcounts=Int64[];spicecounts=Int64[];negative=Int64[]
    levels=length(s.c["compression_fractions"]);stories=zeros(Int64,length(source.selected),levels);affectedcharge=zeros(length(source.selected),levels)
    examples=NamedTuple[];ledger=Int64(0)
    for (i,id) in enumerate(source.selected)
     e=h["events/$id"];n=read(e["counts"]);linear=Float64.(read(e["linear_estimated_ions"]));mean=Float64.(read(e["noiseless_estimated_ions"]));noisy=Float64.(read(e["recovered_ions"]));flags=read(e["compression_level"])
     ni=Int64(read(e["initial_pairs"]));@test sum(n)==ni;ledger+=ni;pos=n.>0
     sigma=haskey(e,"noise_sigma_ions") ? repeat(read(e["noise_sigma_ions"]),1,size(n,2)) : fill(enc,size(n))
     snr=mean[pos]./sigma[pos];append!(expected,snr);append!(truthsnr,n[pos]./sigma[pos]);push!(eventsnr,Statistics.mean(snr))
     push!(initial,ni);push!(cells,count(pos));push!(padcounts,size(n,2));push!(spicecounts,count(==(1),read(e["model_used"])));push!(negative,count(<(0),noisy))
     push!(linearbias,(sum(linear)-ni)/ni);push!(meanbias,(sum(mean)-ni)/ni);push!(noisybias,(sum(noisy)-ni)/ni)
     deficit=max.(linear.-mean,0.);push!(deficitfraction,sum(deficit[pos])/ni)
     for j in 1:levels
      mask=flags.>=j;stories[i,j]=count(vec(any(mask;dims=1)));affectedcharge[i,j]=sum(n[mask])/ni
     end
     for col in axes(n,2)
      positive=n[:,col].>0
      push!(examples,(;event=id,column=col,mean_count=Statistics.mean(n[positive,col]),maximum_count=maximum(n[:,col]),maximum_flag=maximum(flags[:,col]),deficit=sum(deficit[positive,col])))
     end
    end
    @test ledger==read(h["metadata/input_ion_total"])
    for (key,value) in ("event_id"=>source.selected,"initial_ions"=>initial,"positive_cells"=>cells,"pads_per_event"=>padcounts,"spice_pads_per_event"=>spicecounts,"negative_estimates_per_event"=>negative,"event_mean_snr"=>eventsnr,"linear_charge_bias_fraction"=>linearbias,"noiseless_charge_bias_fraction"=>meanbias,"noisy_charge_bias_fraction"=>noisybias,"compression_deficit_fraction"=>deficitfraction,"compressed_histories_per_event"=>stories,"affected_input_charge_fraction"=>affectedcharge,"compression_fractions"=>s.c["compression_fractions"])
     g[key]=value
    end
    for (key,value) in ("inclusive_mean_snr"=>Statistics.mean(expected),"inclusive_median_snr"=>median(expected),"event_weighted_mean_snr"=>Statistics.mean(eventsnr),"mean_ions_per_positive_cell"=>ledger/length(expected),"positive_cells"=>length(expected),"nonpositive_snr_cells"=>count(<=(0),expected),"mean_noiseless_charge_bias_fraction"=>Statistics.mean(meanbias),"mean_noisy_charge_bias_fraction"=>Statistics.mean(noisybias),"mean_compression_deficit_fraction"=>Statistics.mean(deficitfraction))
     attributes(g)[key]=value
    end
    g["mean_compressed_histories_per_event"]=vec(Statistics.mean(stories;dims=1));g["fraction_events_with_compression"]=vec(Statistics.mean(stories.>0;dims=1));g["mean_affected_input_charge_fraction"]=vec(Statistics.mean(affectedcharge;dims=1))
    savehist(g,"snr",expected,10.0.^range(-3,3.5;length=101));savehist(g,"input_snr",truthsnr,10.0.^range(-3,3.5;length=101))
    savehist(g,"snr_low",expected,collect(-5.:.25:15.));savehist(g,"event_snr",eventsnr,10.0.^range(-3,3.;length=81))
    for (key,values) in ("linear_bias"=>linearbias,"mean_bias"=>meanbias,"noisy_bias"=>noisybias,"compression_deficit"=>deficitfraction)
     # Symmetric binning preserves signed charge bias and noise fluctuations.
     width=max(.005,maximum(abs.(values))*1.02);savehist(g,key,values,collect(range(-width,width;length=81)))
    end
    for j in 1:levels
     savehist(g,"stories_$j",stories[:,j],collect(-.5:1.:maximum(stories[:,j])+.5))
    end
    order=sortperm(examples;by=x->x.deficit,rev=true);chosen=order[1:min(6,length(order))]
    target=ledger/length(expected);typical=argmin([abs(x.mean_count-target) for x in examples]);pushfirst!(chosen,typical)
    quiet=findall(x->x.maximum_flag==0,examples)
    !isempty(quiet) && push!(chosen,quiet[argmax([examples[j].maximum_count for j in quiet])])
    weak=argmin([abs(x.mean_count-5) for x in examples]);push!(chosen,weak);unique!(chosen)
    for (j,k) in enumerate(chosen)
     x=examples[k];e=h["events/$(x.event)"];eg=create_group(g,"examples/$j");col=x.column
     attributes(eg)["event_id"]=x.event;attributes(eg)["pixel_ix"]=read(e["pixel_ix"])[col];attributes(eg)["pixel_iy"]=read(e["pixel_iy"])[col]
     attributes(eg)["label"]=k==typical ? "Typical signal" : k==weak ? "Weak signal" : x.maximum_flag>0 ? "Compressed history" : "Strong uncompressed signal"
     for field in ("counts","linear_estimated_ions","noiseless_estimated_ions","recovered_ions","noise_ions","output_at_edges_V","average_at_edges_V","compression_level")
      eg[field]=e[field][:,col]
     end
     eg["edge_time_s"]=read(e["edge_time_s"])
     n=read(e["counts"]);nz=findall(vec(any(n.>0;dims=2)));times=read(e["edge_time_s"])
     eg["event_window_s"]=[times[max(1,first(nz)-2)],times[min(length(times),last(nz)+3)]]
     eg["model_used"]=read(e["model_used"])[col]
    end
    println(name," positive cells=",length(expected)," mean S/N=",Statistics.mean(expected)," event-mean S/N=",Statistics.mean(eventsnr)," compression50 events=",100Statistics.mean(stories[:,3].>0),"% mean deficit=",100Statistics.mean(deficitfraction),"%");flush(stdout)
   end
  end
  out["groups"]=names
 end
 mv(destination*".tmp",destination;force=true)
 end
 println(destination)
end
end
