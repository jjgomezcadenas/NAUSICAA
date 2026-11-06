module IonCSAResponse
using HDF5, TOML, SHA, Random, LinearAlgebra, Statistics, Test
include("pad_history_readout.jl")
const P=PadHistoryReadout
const D=P.D
const B=D.B
const A=D.A
include("bucket_noise_covariance.jl")
include("../event_analysis/EventAnalysis.jl")
include("../../scripts/production/specs.jl")
include("ion_csa_nonlinear_io.jl")
include("ion_csa_correction.jl")
const NL=IonCSANonlinearIO
const NC=IonCSACorrection
const EA=EventAnalysis
const S=ProductionSpecs
sha(path)=open(io->bytes2hex(sha256(io)),path)
attrs(g)=Dict(k=>read(attributes(g)[k]) for k in keys(attributes(g)))
str(g,key)=String(read(g[key]))

"""Restore dense Int64 pad/frame counts, signed bucket origin and a complete ion ledger."""
function dense_history(history,frame,card)
 cells=history.cells
 keys=collect(zip(cells.pixel_ix,cells.pixel_iy));pads=unique(keys)
 isempty(pads) && error("Pilot requires an event with ion counts")
 firstbucket=minimum(cells.bucket)-card["quiet_before_frames"]
 lastbucket=maximum(cells.bucket)+ceil(Int,card["recovery_tail_s"]/frame)
 n=lastbucket-firstbucket+1;counts=zeros(Int64,n,length(pads));positions=Dict(p=>j for (j,p) in enumerate(pads))
 for i in eachindex(cells.count)
  counts[cells.bucket[i]-firstbucket+1,positions[keys[i]]]=cells.count[i]
 end
 @assert sum(counts)==history.initial_pairs
 (;counts,pads,firstbucket,edges=collect(firstbucket:lastbucket+1).*frame)
end

"""Integrate filtered one-sided PSD covariance at all requested integer frame lags.

Cosine recurrence shares spectral quadrature across lags. The caller bounds the
omitted high-frequency covariance independently through its integrated variance.
"""
function covariance(f,asd,width,frame,alpha,maxlag;order=8)
 b=[k/sqrt(4k^2-1) for k in 1:order-1];eg=eigen(SymTridiagonal(zeros(order),b))
 edges=sort!(unique(vcat(f,collect(first(f):(.5/(maxlag*frame+width+frame)):last(f)))))
 lf=log.(f);lp=2log.(asd);out=zeros(maxlag+1)
 for j in 1:length(edges)-1
  half=(edges[j+1]-edges[j])/2;mid=(edges[j+1]+edges[j])/2
  for q in eachindex(eg.values)
   x=mid+half*eg.values[q];k=clamp(searchsortedlast(f,x),1,length(f)-1)
   psd=exp(lp[k]+(lp[k+1]-lp[k])*(log(x)-lf[k])/(lf[k+1]-lf[k]))
   w=half*2eg.vectors[1,q]^2*psd*sinc(x*width)^2*((1-alpha)^2+4alpha*sinpi(x*frame)^2)
   out[1]+=w
   if maxlag>0
    cosine=cospi(2x*frame);old=1.;current=cosine;out[2]+=w*current
    for lag in 2:maxlag
     next=2cosine*current-old;out[lag+1]+=w*next;old=current;current=next
    end
   end
  end
 end
 out
end

"""Read validated cards and derive response identities from content and implementation."""
function settings(root,runpath)
 runfile=joinpath(root,runpath);run=TOML.parsefile(runfile)
 get(run,"schema","")=="nausicaa.ion_csa_run.v1" || error("Ion CSA run schema")
 cp=joinpath(root,run["response_card"]);c=NL.resolve_card(root,TOML.parsefile(cp))
 get(c,"schema","")=="nausicaa.ion_csa_spec.v1" || error("Ion CSA card schema")
 c["input_timing"]=="uniform_mean_arrivals_with_finite_gap_current" || error("Unknown arrival model")
 c["channel_noise"]=="independent_pads" || error("Unknown pad covariance model")
 all(>(0),[c[k] for k in ("averaging_s","pad_side_m","gap_m","gap_speed_m_s","current_step_s","solver_step_s")]) || error("Nonpositive CSA parameter")
 c["quiet_before_frames"]>=1 && c["recovery_tail_s"]>=0 || error("Invalid quiet interval")
 c["noise_low_Hz"]>0 && c["noise_high_Hz"]>c["noise_low_Hz"] || error("Noise band")
 c["noise_seed"]>=0 || error("Negative noise seed")
 run["selection"]=="evenly_spaced_sorted_event_ids" || error("Unknown event selection")
 all(p->p in (2,3,4),run["pitches_mm"]) && run["events_per_sample"]>0 || error("Pilot grid/event count")
 caldir=joinpath(root,c["calibration_directory"]);calfile=joinpath(caldir,"analysis.h5")
 files=[@__FILE__,joinpath(@__DIR__,"pad_history_readout.jl"),joinpath(@__DIR__,"bucket_study.jl"),joinpath(@__DIR__,"averaging_scan.jl"),joinpath(@__DIR__,"discharge_calibration.jl"),joinpath(@__DIR__,"bucket_noise_covariance.jl"),joinpath(root,"src/event_analysis/production_readout.jl"),joinpath(root,"scripts/production/specs.jl"),joinpath(root,"scripts/production/run_ion_csa_pilot.py"),joinpath(root,"scripts/production/ion_csa_pilot.jl"),joinpath(@__DIR__,"spice_acquisition.py"),joinpath(@__DIR__,"spice.py")]
 append!(files,[joinpath(@__DIR__,f) for f in ("ion_csa_nonlinear.jl","ion_csa_nonlinear_io.jl","ion_csa_nonlinear_cases.jl","ion_csa_correction.jl")])
 codehash=bytes2hex(sha256(join(sha.(files),"|")))
 modelid=bytes2hex(sha256(join((S.canonical(c),sha(calfile),codehash),"|")))
 sources=NamedTuple[]
 h5open(joinpath(root,run["ion_index"]),"r") do idx
  for sample in run["samples"],pitch in run["pitches_mm"]
   path=str(idx,"$sample/ions/$pitch/path")
   h5open(path,"r") do h
    String(read(attributes(h)["schema"]))=="nausicaa.ion_pad_file.pilot.v1" || error("Ion source schema")
    m=h["metadata"];Bool(read(m["complete"])) || error("Incomplete ion input")
    str(m,"sample")==sample && str(m,"signal_unit")=="ions" && str(m,"time_unit")=="ms" || error("Ion source units/sample")
    upstream=TOML.parse(str(m,"branch_card_toml"));frame=upstream["bucket_width_ms"]*.001
    isapprox(upstream["pad_pitch_mm"],pitch) || error("Ion pitch differs")
    frame>c["averaging_s"] || error("Aperture exceeds frame")
    ids=sort!(parse.(Int64,collect(keys(h["events"]))));n=run["events_per_sample"]
    length(ids)==read(m["event_count"]) && n<=length(ids) || error("Event count mismatch")
    selected=n==1 ? ids[1:1] : ids[round.(Int,range(1,length(ids);length=n))];@assert allunique(selected)
    sourcehash=sha(path)
    nodeid=bytes2hex(sha256(join((modelid,sourcehash,S.canonical(selected)),"|")))
    base=dirname(dirname(dirname(dirname(path))))
    response=joinpath(base,"response","ions",basename(dirname(path)),nodeid,"ion_response.h5")
    push!(sources,(;sample,pitch,path,frame,selected,sourcehash,nodeid,response,source_events=length(ids)))
   end
  end
 end
 studyid=bytes2hex(sha256(join((modelid,S.canonical(run),join(s.nodeid for s in sources)),"|")))[1:16]
 dir=joinpath(root,run["results_directory"],studyid);mkpath(dir)
 write(joinpath(dir,".gitignore"),"*.h5\n*.h5.tmp\nspice_cases/\n")
 (;root,run,c,cp,runfile,caldir,calfile,codehash,modelid,sources,studyid,dir,files)
end

"""Save the v3 frame covariance, finite-path charge and independent numerical controls."""
function calibration(s,maxframes)
 path=joinpath(s.dir,"calibration.h5");isfile(path) && return path
 cfg=Dict("circuit_analysis"=>relpath(s.calfile,s.root),"circuit_results"=>relpath(joinpath(s.caldir,"variants",s.c["variant"]*".h5"),s.root),"circuit_name"=>s.c["variant"],"prehistory_s"=>.01,"averaging_s"=>s.c["averaging_s"])
 @testset "Pilot calibration and covariance" begin
 h5open(path*".tmp","w") do out
  attributes(out)["model_id"]=s.modelid
  h5open(s.calfile,"r") do h
   v=h[s.c["variant"]];f=read(v["frequency_Hz"]);asd=read(v["output_asd_V_sqrtHz"])
   @test first(f)≈s.c["noise_low_Hz"]
   @test last(f)≈s.c["noise_high_Hz"]
   q=B.E*(1-B.phi(s.c["gap_m"],s.c["pad_side_m"]/2));out["finite_path_charge_C"]=q
   t,i=B.pwl(.03,s.c["pad_side_m"]/2,s.c["gap_m"],s.c["gap_speed_m_s"],.01,.15,2.5e-6)
   @test B.integral(t,i)≈q rtol=1e-9
   for frame in sort!(unique(vcat([.02,.04],[x.frame for x in s.sources])))
    key="$(round(Int,1000frame))ms";source=v[key];k=P.kernel(s.root,cfg,frame)
    sigma=A.noise_rms(f,asd,s.c["averaging_s"],frame,k.alpha)
    @test sigma/only(read(source["noise_rms_V"]))≈1 rtol=1e-8
    @test k.gain/only(read(source["gain_V_per_ion"]))≈1 rtol=1e-8
    keep=f.<=s.c["covariance_max_frequency_Hz"]
    cv=covariance(f[keep],asd[keep],s.c["averaging_s"],frame,k.alpha,maxframes)
    omitted=max(0.,sigma^2-cv[1]);fraction=omitted/sigma^2
    @test fraction<s.c["covariance_tail_fraction_tolerance"]
    refined=covariance(f[keep],asd[keep],s.c["averaging_s"],frame,k.alpha,maxframes;order=12)
    @test maximum(abs.(refined.-cv))/sigma^2<1e-6
    reference=BucketNoiseCovariance.covariance(f,asd,s.c["averaging_s"],frame,k.alpha,3)
    @test maximum(abs.(reference.-cv[1:4]))/sigma^2<fraction+1e-6
    cv[1]=sigma^2
    matrix=[cv[abs(i-j)+1]/k.gain^2 for i in 1:64,j in 1:64]
    @test minimum(eigvals(Symmetric(matrix)))>0
    L=cholesky(Symmetric(matrix)).L
    draws=L*randn(Xoshiro(s.c["noise_seed"]+round(Int,frame*1e6)),64,12000)
    variance=mean(abs2,draws);empirical=[mean(draws[1:end-l,:].*draws[1+l:end,:]) for l in 0:5]
    @test abs(variance/(sigma/k.gain)^2-1)<.025
    @test maximum(abs.(empirical.-cv[1:6]./k.gain^2))/(sigma/k.gain)^2<.025
    @test abs(mean(draws))/(sigma/k.gain)<.01
    g=create_group(out,key)
    for (name,value) in ("frame_s"=>frame,"tau_s"=>k.tau,"alpha"=>k.alpha,"gain_V_per_ion"=>k.gain,"baseline_V"=>k.baseline,"noise_rms_V"=>sigma,"ENC_e_rms"=>sigma/k.gain,"omitted_covariance_bound_fraction"=>fraction,"noise_control_mean_e"=>mean(draws),"noise_control_variance_e2"=>variance)
     attributes(g)[name]=value
    end
    g["raw_V_per_ion"]=k.raw;g["average_V_per_ion"]=k.avg;g["covariance_V2"]=cv;g["covariance_e2"]=cv./k.gain^2
    g["noise_control_empirical_covariance_e2"]=empirical
    edges=collect(-5.:.1:5.).*(sigma/k.gain);g["noise_control_hist_edges_e"]=edges
    bins=zeros(Int64,length(edges)-1);under=0;over=0
    for value in draws
     j=searchsortedlast(edges,value)
     if j<1;under+=1;elseif j>=length(edges);over+=1;else;bins[j]+=1;end
    end
    g["noise_control_hist_count"]=bins;attributes(g)["noise_control_underflow"]=under;attributes(g)["noise_control_overflow"]=over
    println(key," gain=",k.gain," ENC=",sigma/k.gain," correlations=",cv[1:6]./cv[1]," omitted bound=",fraction);flush(stdout)
   end
  end
 end
 mv(path*".tmp",path)
 end
 path
end

"""Create a deterministic pilot plan and linear pad predictions from owned ion histories."""
function prepare(root,runpath)
 s=settings(root,runpath);plan=joinpath(s.dir,"plan.h5")
 if isfile(plan);println(s.dir);return s.dir;end
 maxframes=64
 for source in s.sources
  h5open(source.path,"r") do h
   for id in source.selected
    b=read(h["events/$id/bucket"])
    maxframes=max(maxframes,maximum(b)-minimum(b)+s.c["quiet_before_frames"]+ceil(Int,s.c["recovery_tail_s"]/source.frame)+1)
   end
  end
 end
 calibration(s,maxframes)
 names=String[];events=Int64[];xs=Int64[];ys=Int64[];columns=Int64[];minimums=Float64[];maxima=Int64[];groups=String[];candidate=Bool[]
 h5open(plan*".tmp","w") do out
  attributes(out)["study_id"]=s.studyid;attributes(out)["model_id"]=s.modelid;attributes(out)["code_hash"]=s.codehash
  out["run_card_toml"]=read(s.runfile,String);out["response_card_toml"]=read(s.cp,String)
  h5open(joinpath(s.dir,"calibration.h5"),"r") do cal
   for source in s.sources
    name=source.sample*"_p"*string(source.pitch);g=create_group(out,name);at=attributes(g)
    for field in (:sample,:pitch,:path,:frame,:sourcehash,:nodeid,:response,:source_events);at[string(field)]=getproperty(source,field);end
    g["selected_event_ids"]=source.selected
    cg=cal["$(round(Int,1000source.frame))ms"];a=attrs(cg)
    k=(;frame=source.frame,tau=a["tau_s"],alpha=a["alpha"],gain=a["gain_V_per_ion"],baseline=a["baseline_V"],raw=read(cg["raw_V_per_ion"]),avg=read(cg["average_V_per_ion"]))
    h5open(source.path,"r") do input
     for id in source.selected
      history=EA.read_ion_history(input["events/$id"]);dense=dense_history(history,source.frame,s.c)
      e=create_group(g,"events/$id");e["counts",compress=6]=dense.counts;e["edge_time_s"]=dense.edges
      e["pixel_ix"]=first.(dense.pads);e["pixel_iy"]=last.(dense.pads);attributes(e)["first_bucket"]=dense.firstbucket;attributes(e)["initial_pairs"]=history.initial_pairs
      linear=zeros(size(dense.counts));voltage=zeros(size(dense.counts,1)+1,size(dense.counts,2));average=similar(voltage)
      for col in axes(dense.counts,2)
       rv,av,di=P.response(view(dense.counts,:,col),k)
       linear[:,col]=di[1:size(linear,1)]./k.gain;voltage[:,col]=rv[1:size(voltage,1)];average[:,col]=av[1:size(voltage,1)]
       push!(names,name);push!(events,id);push!(xs,dense.pads[col][1]);push!(ys,dense.pads[col][2]);push!(columns,col);push!(minimums,minimum(rv));push!(maxima,maximum(view(dense.counts,:,col)));push!(groups,"$name/events/$id");push!(candidate,minimum(rv)<s.c["initial_spice_guard_V"])
      end
      e["linear_estimated_ions",compress=6]=linear;e["linear_output_at_edges_V",compress=6]=voltage;e["linear_average_V",compress=6]=average
     end
    end
   end
  end
  for (key,value) in ("name"=>names,"event_id"=>events,"pixel_ix"=>xs,"pixel_iy"=>ys,"column"=>columns,"minimum_V"=>minimums,"maximum_count"=>maxima,"event_group"=>groups,"candidate"=>candidate)
   out["pads/$key"]=value
  end
  controls=fill(false,length(names));edges=s.c["control_voltage_edges_V"]
  for name in unique(names),j in 1:length(edges)-1
   get(s.c,"mode","validation")=="production" && continue
   possible=findall((names.==name).&(minimums.>=edges[j]).&(minimums.<edges[j+1]))
   isempty(possible) && continue
   sort!(possible;by=i->minimums[i])
   chosen=possible[round.(Int,range(1,length(possible);length=min(length(possible),s.c["controls_per_stratum"])))]
   controls[chosen].=true
  end
  out["pads/control"]=controls
 end
 mv(plan*".tmp",plan)
 for sourcefile in s.files
  dest=joinpath(s.dir,"source_snapshot",sha(sourcefile),basename(sourcefile));mkpath(dirname(dest));cp(sourcefile,dest;force=true)
 end
 open(joinpath(s.dir,"manifest.toml"),"w") do io
  TOML.print(io,Dict("study_id"=>s.studyid,"model_id"=>s.modelid,"code_hash"=>s.codehash,"run_card"=>runpath,"response_card"=>s.run["response_card"],"calibration_sha256"=>sha(joinpath(s.dir,"calibration.h5")),"source_paths"=>[x.path for x in s.sources],"source_hashes"=>[x.sourcehash for x in s.sources]))
 end
 jobs(root,runpath,s.dir)
 println(s.dir);s.dir
end

const PLAN_COLUMN_CACHE=IdDict{Any,Any}()
"""Look up one pad request using owned, read-only columns cached per open plan."""
function caseinfo(plan,i)
 data=get!(PLAN_COLUMN_CACHE,plan) do
  (;names=read(plan["pads/name"]),events=read(plan["pads/event_id"]),xs=read(plan["pads/pixel_ix"]),ys=read(plan["pads/pixel_iy"]),columns=read(plan["pads/column"]),groups=read(plan["pads/event_group"]))
 end
 name=data.names[i];event=data.events[i];x=data.xs[i];y=data.ys[i]
 (;name,event,x,y,key="$(name)_$(event)_$(x)_$(y)",column=data.columns[i],group=data.groups[i])
end

"""Write finite-path currents for the selected transistor simulations."""
function jobs(root,runpath,dir;guard=nothing)
 s=settings(root,runpath);@assert s.dir==dir
 planfile=joinpath(dir,"plan.h5");records=Dict{String,Any}[]
 h5open(planfile,"r") do plan
  minv=read(plan["pads/minimum_V"]);controls=read(plan["pads/control"])
  cutoff=isnothing(guard) ? s.c["initial_spice_guard_V"] : guard
  selected=findall((minv.<cutoff).|controls)
  for i in selected
   info=caseinfo(plan,i);folder=joinpath(dir,"spice_cases",info.key);mkpath(folder)
   if !isfile(joinpath(folder,"input.toml"))
    e=plan[info.group];n=e["counts"][:,info.column];edges=read(e["edge_time_s"]);frame=read(attributes(plan[info.name])["frame"])
    active=findall(>(0),n);firstbin=first(active);pre=s.c["simulation_prehistory_s"];offset=edges[firstbin]
    u=collect(range(-s.c["gap_m"]/s.c["gap_speed_m_s"],0.;length=ceil(Int,s.c["gap_m"]/s.c["gap_speed_m_s"]/s.c["current_step_s"])+1))
    ut=vcat(u,frame.+u);ui=B.current.(ut,frame,s.c["pad_side_m"]/2,s.c["gap_m"],s.c["gap_speed_m_s"])
    stop=last(edges)-offset+pre+s.c["averaging_s"]
    times=sort!(unique(round.(vcat(0.,[x+pre+(j-firstbin)*frame for j in active for x in ut],stop);digits=12)))
    current=[sum(n[j]*((first(ut)<=x-pre-(j-firstbin)*frame<=last(ut)) ? B.interp(ut,ui,x-pre-(j-firstbin)*frame) : 0.) for j in active) for x in times]
    charge=sum(n)*B.E*(1-B.phi(s.c["gap_m"],s.c["pad_side_m"]/2));error=abs(B.integral(times,current)/charge-1)
    @assert error<1e-8
    open(joinpath(folder,"current.pwl"),"w") do io
     for j in eachindex(times);println(io,times[j]," ",current[j]);end
    end
    open(joinpath(folder,"input.toml"),"w") do io
     TOML.print(io,Dict("pad_row"=>i,"offset_s"=>offset,"prehistory_s"=>pre,"stoptime_s"=>stop,"expected_charge_C"=>charge,"pwl_charge_error"=>error,"frame_s"=>frame))
    end
   end
   push!(records,Dict("row"=>i,"folder"=>folder,"key"=>info.key))
  end
  open(joinpath(dir,"jobs.toml"),"w") do io;TOML.print(io,Dict("guard_V"=>cutoff,"jobs"=>records));end
 end
 println("Requested ",length(records)," transistor histories, guard ",isnothing(guard) ? s.c["initial_spice_guard_V"] : guard);flush(stdout)
end

"""Reduce transistor waveforms to the same sampled mean voltages and ion estimator."""
function reduce_cases(root,runpath,dir;refined=false)
 s=settings(root,runpath);jobcard=TOML.parsefile(joinpath(dir,"jobs.toml"));failures=String[]
 h5open(joinpath(dir,"plan.h5"),"r") do plan
 h5open(joinpath(dir,"calibration.h5"),"r") do cal
  for job in jobcard["jobs"]
   folder=job["folder"];destination=joinpath(folder,refined ? "refined_response.h5" : "response.h5");isfile(destination) && continue
   rawpath=joinpath(folder,refined ? "refined_raw.h5" : "raw.h5");
   if refined && !isfile(rawpath);continue;end
   isfile(rawpath) || error("Missing acquisition $folder")
   info=caseinfo(plan,job["row"]);input=TOML.parsefile(joinpath(folder,"input.toml"));e=plan[info.group]
   frame=read(attributes(plan[info.name])["frame"]);k=attrs(cal["$(round(Int,1000frame))ms"])
   h5open(rawpath,"r") do raw
    t=read(raw["time"]);v=read(raw["V(out)"]);ii=read(raw["I(Is1)"]);charge_error=abs(B.integral(t,ii)/input["expected_charge_C"]-1)
    if charge_error>=1e-7
     open(joinpath(dir,"charge_retry.toml"),"w") do io
      TOML.print(io,Dict("row"=>job["row"],"folder"=>folder,"refined"=>refined,"charge_error"=>charge_error))
     end
     error("Charge quadrature failed: $folder, error $charge_error")
    end
    baseline=k["baseline_V"];cum=vcat(0.,cumsum(diff(t).*((v[1:end-1].-baseline).+(v[2:end].-baseline))./2))
    function primitive(x)
     x<=first(t) && return 0.
     @assert x<=last(t)+1e-9
     j=clamp(searchsortedlast(t,x),1,length(t)-1);dx=x-t[j]
     cum[j]+(v[j]-baseline)*dx+(v[j+1]-v[j])/(t[j+1]-t[j])*dx^2/2
    end
    edges=read(e["edge_time_s"]);times=edges.-input["offset_s"].+input["prehistory_s"];width=s.c["averaging_s"]
    average=[(primitive(x+width/2)-primitive(x-width/2))/width for x in times]
    voltage=[x<first(t) ? baseline : B.interp(t,v,x) for x in times]
    signal=-(average[2:end].-k["alpha"].*average[1:end-1]);ions=signal./k["gain_V_per_ion"]
    linear=e["linear_estimated_ions"][:,info.column]
    h5open(destination*".tmp","w") do out
     for (key,value) in ("charge_error"=>charge_error,"maximum_error_ENC"=>maximum(abs.(ions.-linear))/k["ENC_e_rms"],"raw_sha256"=>sha(rawpath),"minimum_output_V"=>minimum(v))
      attributes(out)[key]=value
     end
     for (key,value) in attrs(raw);attributes(out)[key]=value;end
     out["estimated_ions"]=ions;out["output_at_edges_V"]=voltage;out["average_at_edges_V"]=average.+baseline
    end
    mv(destination*".tmp",destination)
   end
  end
 end
 end
end

"""Expand the direct-SPICE range until control errors are below the saved noise tolerance."""
function assess(root,runpath,dir)
 s=settings(root,runpath);current=TOML.parsefile(joinpath(dir,"jobs.toml"))["guard_V"]
 guard=current;errors=Float64[];bins=s.c["control_voltage_edges_V"]
 h5open(joinpath(dir,"plan.h5"),"r") do plan
  controls=read(plan["pads/control"]);minv=read(plan["pads/minimum_V"])
  for i in findall(controls)
   info=caseinfo(plan,i);p=joinpath(dir,"spice_cases",info.key,"response.h5")
   err=h5open(h->read(attributes(h)["maximum_error_ENC"]),p,"r")
   if minv[i]>=current
    push!(errors,err)
    if err>s.c["linear_error_tolerance_ENC"]
     j=clamp(searchsortedlast(bins,minv[i]),1,length(bins)-1);guard=max(guard,bins[j+1])
    end
   end
  end
 end
 open(joinpath(dir,"fast_path_validation.toml"),"w") do io
  production=get(s.c,"mode","validation")=="production"
  TOML.print(io,Dict("guard_V"=>current,"proposed_guard_V"=>guard,"maximum_remaining_control_error_ENC"=>isempty(errors) ? NaN : maximum(errors),"control_errors_ENC"=>errors,"passed"=>guard==current,
   "mode"=>get(s.c,"mode","validation"),"validation_source"=>production ? "saved_calibration" : "run_controls",
   "nonlinear_calibration_sha256"=>get(s.c,"nonlinear_calibration_sha256","")))
 end
 if guard>current;jobs(root,runpath,dir;guard);println("EXPANDED ",guard);else;println("VALIDATED ",guard);end
end

"""Sample stationary correlated Gaussian noise with a stable per-pad seed."""
function noise_draw(covariance_e2,n,seed;factor_cache=Dict{Int,Any}())
 L=get!(factor_cache,n) do
  matrix=[covariance_e2[abs(i-j)+1] for i in 1:n,j in 1:n]
  cholesky(Symmetric(matrix)).L
 end
 L*randn(Xoshiro(seed),n)
end

"""Write complete versioned pilot responses atomically, preserving all signed estimates."""
function finish(root,runpath,dir)
 s=settings(root,runpath);validation=TOML.parsefile(joinpath(dir,"fast_path_validation.toml"));@assert validation["passed"]
 @testset "Pilot response ledgers and storage" begin
 h5open(joinpath(dir,"plan.h5"),"r") do plan
 h5open(joinpath(dir,"calibration.h5"),"r") do calibration
  names=read(plan["pads/name"]);ev=read(plan["pads/event_id"]);cols=read(plan["pads/column"])
  for source in s.sources
   name=source.sample*"_p"*string(source.pitch);g=plan[name]
   output=source.response;mkpath(dirname(output));write(joinpath(dirname(output),".gitignore"),"*.h5\n*.h5.tmp\n")
   if isfile(output)
    h5open(output,"r") do existing;@assert Bool(read(existing["metadata/complete"]));end
    continue
   end
   factor_cache=Dict{Int,Any}();key="$(round(Int,1000source.frame))ms";ca=attrs(calibration[key]);cv=read(calibration[key*"/covariance_e2"])
   count_ledger=Int64(0);maximum_storage_error=0.;spicepads=0;nonlinearpads=0;allpads=0;negative=0
   started=time_ns()
   h5open(source.path,"r") do upstream
   h5open(output*".tmp","w") do out
    attributes(out)["schema"]="nausicaa.ion_csa_response.pilot.v1"
    m=create_group(out,"metadata");m["complete"]=false;m["source_path"]=source.path;m["source_sha256"]=source.sourcehash
    m["source_event_count"]=source.source_events;m["selected_event_ids"]=source.selected;m["event_count"]=length(source.selected)
    m["sample"]=source.sample;m["model_id"]=s.modelid;m["response_id"]=source.nodeid;m["study_id"]=s.studyid;m["code_hash"]=s.codehash
    m["response_model"]=get(s.c,"response_model","hybrid_spice")
    m["mode"]=get(s.c,"mode","validation")
    m["model_used_codes"]="0=linear,1=SPICE,2=calibrated_nonlinear"
    if haskey(s.c,"nonlinear_calibration_sha256")
     m["nonlinear_calibration_sha256"]=s.c["nonlinear_calibration_sha256"]
     m["nonlinear_calibration"]=s.c["nonlinear_calibration"]
    end
    m["response_card_toml"]=read(s.cp,String);m["run_card_toml"]=read(s.runfile,String);m["signal_unit"]="ion_charges";m["time_unit"]="s"
    m["matrix_axes_hdf5"]="pad,frame";m["matrix_axes_julia"]="frame,pad";m["nominal_discharge_time_s"]=s.c["nominal_discharge_time_s"]
    m["frame_s"]=source.frame;m["pitch_mm"]=source.pitch;m["noise_seed"]=s.c["noise_seed"];m["noise_assumption"]="Stationary quiescent-point Gaussian frame noise; independent pads; full saved temporal covariance."
    HDF5.copy_object(calibration[key],out,"calibration")
    for id in source.selected
     e=plan["$name/events/$id"];counts=read(e["counts"]);linear=read(e["linear_estimated_ions"]);mean=copy(linear)
     voltage=read(e["linear_output_at_edges_V"]);average=read(e["linear_average_V"])
     model=zeros(Int8,size(counts,2));noise=zeros(size(mean));seeds=zeros(UInt64,size(counts,2));px=read(e["pixel_ix"]);py=read(e["pixel_iy"])
     rows=findall((names.==name).&(ev.==id))
     for row in rows
      info=caseinfo(plan,row);p=joinpath(dir,"spice_cases",info.key,"response.h5");col=cols[row]
      if isfile(p)
       h5open(p,"r") do h
        mean[:,col]=read(h["estimated_ions"]);voltage[:,col]=read(h["output_at_edges_V"]);average[:,col]=read(h["average_at_edges_V"])
        model[col]=haskey(attributes(h),"response_model") && read(attributes(h)["response_model"])=="calibrated_nonlinear" ? 2 : 1
       end
       spicepads+=model[col]==1;nonlinearpads+=model[col]==2
      end
      # Model identity excludes the pilot event selection, giving stable replay.
      seed=S.event_seed(s.modelid*"|"*source.sourcehash,id,"pad:$((px[col],py[col])):$(s.c["noise_seed"])")
      seeds[col]=seed;noise[:,col]=noise_draw(cv,size(counts,1),seed;factor_cache)
     end
     recovered=mean.+noise;negative+=count(<(0.),recovered);allpads+=length(model)
     original=upstream["events/$id"];initial=Int64(read(original["initial_pairs"]));@test sum(counts)==initial
     count_ledger+=initial;event=create_group(out,"events/$id");attributes(event)["schema"]="nausicaa.ion_csa_event.v1"
     for field in ("event_id","initial_pairs","pair_checksum","ion_drift_cm","electron_drift_cm","reference_mm","phase_mm","pitch_mm","bucket_width_mm","drift_direction")
      event[field]=read(original[field])
     end
     event["pixel_ix"]=px;event["pixel_iy"]=py;event["first_bucket"]=read(attributes(e)["first_bucket"]);event["edge_time_s"]=read(e["edge_time_s"])
     event["counts",compress=6]=counts;event["model_used"]=model;event["noise_seed"]=seeds
     for (field,array) in ("linear_estimated_ions"=>linear,"noiseless_estimated_ions"=>mean,"noise_ions"=>noise,"recovered_ions"=>recovered,"output_at_edges_V"=>voltage,"average_at_edges_V"=>average)
      saved=Float32.(array);error=maximum(abs.(Float64.(saved).-array));maximum_storage_error=max(maximum_storage_error,error)
      @test error<(occursin("_V",field) ? 1e-7 : .001ca["ENC_e_rms"])
      event[field,compress=6]=saved
     end
     loss=max.(0.,linear.-mean);relative=loss./max.(linear,eps());positive=counts.>0
     flags=zeros(Int8,size(counts))
     for (j,fraction) in enumerate(s.c["compression_fractions"])
      mask=positive .&(relative.>=fraction).&(loss.>=s.c["compression_minimum_deficit_ENC"]*ca["ENC_e_rms"])
      flags[mask].=j
     end
     event["compression_level",compress=6]=flags
     attributes(event)["input_ion_ledger"]=sum(counts);attributes(event)["storage_max_absolute_error"]=maximum_storage_error
    end
    m["input_ion_total"]=count_ledger;m["pads"]=allpads;m["direct_spice_pads"]=spicepads;m["nonlinear_pads"]=nonlinearpads;m["negative_estimates"]=negative
    m["maximum_storage_error"]=maximum_storage_error;m["response_seconds"]=(time_ns()-started)/1e9
    write(m["complete"],true)
   end
   end
   mv(output*".tmp",output)
   open(joinpath(dirname(output),"manifest.toml"),"w") do io
    TOML.print(io,Dict("response_id"=>source.nodeid,"model_id"=>s.modelid,"response_file"=>output,"source_path"=>source.path,"source_sha256"=>source.sourcehash,"code_hash"=>s.codehash,"selected_event_ids"=>source.selected,"complete"=>true))
   end
   println(source.sample," response ",output,"; pads=",allpads," SPICE=",spicepads," negatives=",negative);flush(stdout)
  end
 end
 end
 end
end
"""Check independent solver refinements of the nonlinear histories."""
function validate_cases(root,runpath,dir)
 reduce_cases(root,runpath,dir;refined=true)
 s=settings(root,runpath);rows=Int64[];errors=Float64[];charge=Float64[]
 @testset "Nonlinear solver refinement" begin
 h5open(joinpath(dir,"plan.h5"),"r") do plan
 h5open(joinpath(dir,"calibration.h5"),"r") do cal
  for job in TOML.parsefile(joinpath(dir,"jobs.toml"))["jobs"]
   folder=job["folder"];refined=joinpath(folder,"refined_response.h5");isfile(refined) || continue
   info=caseinfo(plan,job["row"]);frame=read(attributes(plan[info.name])["frame"])
   enc=read(attributes(cal["$(round(Int,1000frame))ms"])["ENC_e_rms"])
   h5open(joinpath(folder,"response.h5"),"r") do a
   h5open(refined,"r") do b
    err=maximum(abs.(read(a["estimated_ions"]).-read(b["estimated_ions"])))/enc
    @test err<.01
    push!(rows,job["row"]);push!(errors,err);push!(charge,read(attributes(b)["charge_error"]))
   end
   end
  end
 end
 end
 @test get(s.c,"mode","validation")=="production" || !isempty(rows)
 end
 h5open(joinpath(dir,"solver_validation.h5"),"w") do out
  out["pad_row"]=rows;out["maximum_error_ENC"]=errors;out["charge_error"]=charge
 end
end

include("ion_csa_nonlinear_cases.jl")

end
