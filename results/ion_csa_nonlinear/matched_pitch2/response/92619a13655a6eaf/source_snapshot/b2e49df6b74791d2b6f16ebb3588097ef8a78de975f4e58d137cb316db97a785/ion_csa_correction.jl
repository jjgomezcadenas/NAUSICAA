module IonCSACorrection
using HDF5,LinearAlgebra,Statistics,Random,SHA,TOML,Test
include("../../scripts/production/specs.jl")
include("ion_csa_summary.jl")
include("ion_csa_nonlinear_io.jl")
include("guarded_files.jl")
const S=ProductionSpecs
sha(p)=open(io->bytes2hex(sha256(io)),p)
attrs(g)=Dict(k=>read(attributes(g)[k]) for k in keys(attributes(g)))
str(g,k)=String(read(g[k]))

"""Construct y[n] = late*x[n-1] + central*x[n] + early*x[n+1].

`taps` is ordered `[early, central, late]`, following the packet response at
the frame before, during and after its collection.
"""
function response_matrix(taps,n)
 length(taps)==3 && n>0 || throw(ArgumentError("Three taps and a positive frame count are required"))
 all(isfinite,taps) && taps[2]>abs(taps[1])+abs(taps[3]) || throw(ArgumentError("Response must be strictly diagonally dominant"))
 Tridiagonal(fill(taps[3],n-1),fill(taps[2],n),fill(taps[1],n-1))
end

"""Build a banded solve and transform the complete input frame covariance."""
function operator(taps,cv,n)
 length(cv)>=n || throw(ArgumentError("Covariance does not cover this history"))
 M=response_matrix(taps,n);F=lu(M)
 old=[cv[abs(i-j)+1] for i in 1:n,j in 1:n]
 intermediate=F\old
 corrected=Matrix(transpose(F\Matrix(transpose(intermediate))))
 corrected=(corrected.+transpose(corrected))./2
 sigma=sqrt.(diag(corrected))
 @assert all(>(0),sigma)
 @assert norm(Matrix(M)*corrected*Matrix(transpose(M))-old)/norm(old)<1e-10
 cholesky(Symmetric(corrected))
 (;M,F,old,covariance=corrected,sigma)
end

"""Variances of normalised-noise sum, sum of squares and adjacent products.

For independent Gaussian pad channels with the given temporal covariance,
these give finite-sample uncertainty of the diagnostic noise moments.
"""
function noise_statistic_variances(covariance,npads)
 n=size(covariance,1)
 n>1 && size(covariance,2)==n && npads>0 || throw(ArgumentError("Noise moments require at least two frames and one pad"))
 sigma=sqrt.(diag(covariance));R=covariance./(sigma*transpose(sigma))
 AR=zeros(n,n)
 AR[1:end-1,:].+=R[2:end,:]./2
 AR[2:end,:].+=R[1:end-1,:]./2
 (;sum=npads*sum(R),squares=2npads*sum(abs2,R),lag=2npads*sum(AR.*transpose(AR)))
end

"""Replace one dataset in a new output, retaining compressed matrix storage."""
function replace(g,key,value)
 haskey(g,key) && HDF5.delete_object(g,key)
 if value isa AbstractArray && ndims(value)>1
  g[key,compress=6]=value
 else
  g[key]=value
 end
end

"""Derive the three reconstructed-ion taps from the saved voltage calibration."""
function taps(g)
 a=attrs(g);v=read(g["average_V_per_ion"])
 [v[1],v[2]-a["alpha"]*v[1],v[3]-a["alpha"]*v[2]]./a["gain_V_per_ion"]
end

function setup(root,runpath)
 path=abspath(joinpath(root,runpath));r=TOML.parsefile(path)
 r["schema"]=="nausicaa.ion_csa_reconstruction_run.v1" || error("Reconstruction run schema")
 cp=joinpath(root,r["reconstruction_card"]);c=TOML.parsefile(cp)
 c["schema"]=="nausicaa.ion_csa_reconstruction.v1" || error("Reconstruction card schema")
 c["method"]=="three_frame_inverse" && c["noise"]=="exact_covariance_transform" || error("Unknown reconstruction")
 !c["nonnegative_constraint"] && !c["signal_to_noise_cut"] || error("This calculation preserves signed, uncut estimates")
 parent=joinpath(root,r["parent_study"]);summary=joinpath(parent,"summary.h5")
 files=[@__FILE__,joinpath(@__DIR__,"guarded_files.jl"),joinpath(@__DIR__,"ion_csa_summary.jl"),joinpath(root,"scripts/production/specs.jl"),joinpath(root,"scripts/production/ion_csa_correction.jl")]
 codehash=bytes2hex(sha256(join(sha.(files),"|")))
 sources=NamedTuple[]
 h5open(summary,"r") do h
  for name in read(h["groups"])
   g=h[name];a=attrs(g);p=String(a["response_path"])
   @assert sha(p)==a["response_sha256"]
   h5open(p,"r") do source
    String(read(attributes(source)["schema"]))=="nausicaa.ion_csa_response.pilot.v1" || error("Requires the original filtered response")
    @assert Bool(read(source["metadata/complete"]))
    parentid=str(source["metadata"],"model_id")
    modelid=bytes2hex(sha256(join((parentid,S.canonical(c),codehash),"|")))
    sourcehash=sha(p);nodeid=bytes2hex(sha256(modelid*"|"*sourcehash))
    target=joinpath(dirname(dirname(p)),nodeid,"ion_response.h5")
    push!(sources,(;name=String(name),sample=String(a["sample"]),pitch=Int(a["pitch_mm"]),frame=Float64(a["frame_s"]),parentpath=p,parenthash=sourcehash,parentid,modelid,nodeid,response=target,selected=read(source["metadata/selected_event_ids"])))
   end
  end
 end
 ids=join(x.nodeid for x in sources)
 studyid=bytes2hex(sha256(join((S.canonical(r),codehash,ids),"|")))[1:16]
 dir=joinpath(root,r["results_directory"],studyid);mkpath(dir)
 write(joinpath(dir,".gitignore"),"*.h5\n*.h5.tmp\nion_csa.html\n")
 (;root,r,c,path,cp,parent,dir,studyid,modelid=first(sources).modelid,codehash,sources,files)
end

"""Save corrected interior calibration and independent zero-input noise controls."""
function calibration(s)
 destination=joinpath(s.dir,"calibration.h5")
 if isfile(destination);return destination;end
 @testset "Three-frame calibration and covariance" begin
 h5open(joinpath(s.parent,"calibration.h5"),"r") do parent
 h5open(destination*".tmp","w") do out
  for key in keys(parent)
   HDF5.copy_object(parent[key],out,key)
   parent[key] isa HDF5.Group || continue
   g=out[key];old=parent[key];tt=taps(old);cv=read(old["covariance_e2"]);n=64
   op=operator(tt,cv,n);mid=n÷2
   n0=zeros(n);n0[mid]=100.;recovered=op.F\(op.M*n0)
   @test maximum(abs.(recovered.-n0))<1e-10
   @test sum(recovered)≈100 rtol=1e-12
   @test norm(op.M*(op.F\Matrix{Float64}(I,n,n))-I)<1e-10
   @test isposdef(Symmetric(op.covariance))
   sigma=op.sigma[mid]
   draws=op.F\(cholesky(Symmetric(op.old)).L*randn(Xoshiro(20261101+parse(Int,key[1:end-2])),n,12000))
   z=draws./op.sigma
   observed_var=mean(abs2,z);observed_mean=mean(z)
   expected_lag=mean(diag(op.covariance,1)./(op.sigma[1:end-1].*op.sigma[2:end]))
   observed_lag=mean(z[1:end-1,:].*z[2:end,:])
   @test abs(observed_var-1)<s.c["sampling_tolerance"]
   @test abs(observed_mean)<s.c["sampling_tolerance"]
   @test abs(observed_lag-expected_lag)<s.c["sampling_tolerance"]
   a=HDF5.attrs(g);a["precorrection_ENC_e_rms"]=read(attributes(old)["ENC_e_rms"])
   a["ENC_e_rms"]=sigma;a["reconstruction_method"]="three_frame_inverse"
   a["covariance_representation"]="Interior display lags; exact history matrices in covariance_by_frames"
   a["predicted_linear_total_bias_fraction"]=0.
   a["precorrection_linear_total_bias_fraction"]=sum(tt)-1
   a["condition_bound"]=sum(abs,tt)/(tt[2]-abs(tt[1])-abs(tt[3]))
   a["noise_control_mean_e"]=mean(draws);a["noise_control_variance_e2"]=mean(abs2,draws)
   a["noise_control_normalised_mean"]=observed_mean;a["noise_control_normalised_variance"]=observed_var
   a["noise_control_lag1_correlation"]=observed_lag;a["expected_control_lag1_correlation"]=expected_lag
   replace(g,"precorrection_covariance_e2",cv);replace(g,"precorrection_linear_reconstruction_kernel",tt)
   replace(g,"linear_reconstruction_kernel",[0.,1.,0.]);replace(g,"covariance_e2",[op.covariance[mid,mid+j] for j in 0:12])
   replace(g,"nominal_covariance_e2",op.covariance)
   edges=collect(-5.:.1:5.).*sigma;bins=zeros(Int64,length(edges)-1);under=0;over=0
   for x in draws
    j=searchsortedlast(edges,x)
    if j<1;under+=1;elseif j>=length(edges);over+=1;else;bins[j]+=1;end
   end
   replace(g,"noise_control_hist_edges_e",edges);replace(g,"noise_control_hist_count",bins)
   a["noise_control_underflow"]=under;a["noise_control_overflow"]=over
   println(key," ENC before=",read(attributes(old)["ENC_e_rms"])," corrected=",sigma," lag1=",op.covariance[mid,mid+1]/sigma^2);flush(stdout)
  end
 end
 end
 end
 mv(destination*".tmp",destination)
 destination
end

"""Re-evaluate saved fast-path and solver controls after the same inverse operation."""
function validate_controls(s)
 original=TOML.parsefile(joinpath(s.parent,"fast_path_validation.toml"));errors=Float64[];refine_errors=Float64[];rows=Int64[]
 if get(original,"mode","validation")=="production"
  card=h5open(h->TOML.parse(String(read(h["response_card_toml"]))),joinpath(s.parent,"plan.h5"),"r")
  resolved=IonCSANonlinearIO.resolve_card(s.root,card)
  @assert resolved["mode"]=="production" && original["passed"]
  @assert original["validation_source"]=="saved_calibration"
  @assert resolved["nonlinear_calibration_sha256"]==original["nonlinear_calibration_sha256"]
  @assert original["guard_V"]==resolved["initial_spice_guard_V"]
  open(joinpath(s.dir,"fast_path_validation.toml"),"w") do io;TOML.print(io,original);end
  h5open(joinpath(s.dir,"solver_validation.h5"),"w") do h
   h["pad_row"]=rows;h["maximum_error_ENC"]=refine_errors
  end
  println("Production: using validated calibration; no repeated SPICE controls")
  return
 end
 tolerance=h5open(joinpath(s.parent,"plan.h5"),"r") do plan
  TOML.parse(String(read(plan["response_card_toml"])))["linear_error_tolerance_ENC"]
 end
 cache=Dict{Tuple{Int,Int},Any}()
 h5open(joinpath(s.parent,"plan.h5"),"r") do plan
 h5open(joinpath(s.parent,"calibration.h5"),"r") do cal
  names=read(plan["pads/name"]);ids=read(plan["pads/event_id"]);xs=read(plan["pads/pixel_ix"]);ys=read(plan["pads/pixel_iy"]);cols=read(plan["pads/column"]);groups=read(plan["pads/event_group"]);controls=read(plan["pads/control"]);minv=read(plan["pads/minimum_V"])
  for row in eachindex(names)
   folder=joinpath(s.parent,"spice_cases","$(names[row])_$(ids[row])_$(xs[row])_$(ys[row])")
   fine=joinpath(folder,"refined_response.h5")
   control=controls[row] && minv[row]>=original["guard_V"]
   (!control && !isfile(fine)) && continue
   e=plan[groups[row]];frame=read(attributes(plan[names[row]])["frame"]);key="$(round(Int,1000frame))ms"
   n=size(e["counts"],1);op=get!(cache,(round(Int,1000frame),n)) do
    operator(taps(cal[key]),read(cal[key*"/covariance_e2"]),n)
   end
   h5open(joinpath(folder,"response.h5"),"r") do coarse
    y=read(coarse["estimated_ions"])
    if control
     reference=e["linear_estimated_ions"][:,cols[row]]
     push!(errors,maximum(abs.(op.F\(y.-reference))./op.sigma))
    end
    if isfile(fine)
     h5open(fine,"r") do h
      push!(refine_errors,maximum(abs.(op.F\(y.-read(h["estimated_ions"])))./op.sigma));push!(rows,row)
     end
    end
   end
  end
 end
 end
 @assert maximum(errors)<tolerance && maximum(refine_errors)<.01
 open(joinpath(s.dir,"fast_path_validation.toml"),"w") do io
  TOML.print(io,Dict("passed"=>true,"guard_V"=>original["guard_V"],"tolerance_ENC"=>tolerance,"maximum_remaining_control_error_ENC"=>maximum(errors),"control_errors_ENC"=>errors,"reconstruction"=>"three_frame_inverse"))
 end
 h5open(joinpath(s.dir,"solver_validation.h5"),"w") do out
  out["pad_row"]=rows;out["maximum_error_ENC"]=refine_errors
 end
 println("Corrected controls: maximum fast error ",maximum(errors)," ENC; refinement ",maximum(refine_errors)," ENC")
end

"""Apply the inverse to saved means and the same noise, retaining voltages and signed counts."""
function correct(s)
 calibration(s);validate_controls(s)
 cache=Dict{Tuple{Int,Int},Any}();output_records=Dict{String,Any}[]
 @testset "Saved pilot reconstruction and noise propagation" begin
 h5open(joinpath(s.parent,"calibration.h5"),"r") do oldcal
 h5open(joinpath(s.dir,"calibration.h5"),"r") do newcal
  for source in s.sources
   key="$(round(Int,1000source.frame))ms";tt=taps(oldcal[key]);cv=read(oldcal[key*"/covariance_e2"])
   output=source.response;mkpath(dirname(output));write(joinpath(dirname(output),".gitignore"),"*.h5\n*.h5.tmp\n")
   if isfile(output)
    h5open(output,"r") do h;@assert Bool(read(h["metadata/complete"]));end
    continue
   end
   started=time_ns();max_storage_ions=0.;max_linear_error=0.;max_additivity_error=0.;negative=0;input_total=Int64(0)
   zsum=0.;zsquare=0.;zcount=Int64(0);lag_sum=0.;lag_expected=0.;lag_count=Int64(0)
   mean_stat_variance=0.;square_stat_variance=0.;lag_stat_variance=0.
   h5open(source.parentpath,"r") do parent
   h5open(output*".tmp","w") do out
    attributes(out)["schema"]="nausicaa.ion_csa_response.pilot.v2"
    HDF5.copy_object(parent["metadata"],out,"metadata");m=out["metadata"];write(m["complete"],false)
    replace(m,"parent_response_path",source.parentpath);replace(m,"parent_response_sha256",source.parenthash)
    replace(m,"physical_model_id",source.parentid);replace(m,"parent_code_hash",str(parent["metadata"],"code_hash"))
    replace(m,"model_id",source.modelid);replace(m,"response_id",source.nodeid);replace(m,"study_id",s.studyid);replace(m,"code_hash",s.codehash)
    replace(m,"reconstruction_card_toml",read(s.cp,String));replace(m,"reconstruction_method","three_frame_inverse")
    replace(m,"covariance_lookup","/covariance_by_frames/<n_frames>/covariance_e2")
    HDF5.copy_object(newcal[key],out,"calibration");create_group(out,"events");create_group(out,"covariance_by_frames")
    lengths=Set{Int}()
    for id in source.selected
     e=parent["events/$id"];counts=read(e["counts"]);n,npads=size(counts)
     op=get!(cache,(round(Int,1000source.frame),n)) do;operator(tt,cv,n);end
     if !(n in lengths)
      g=create_group(out,"covariance_by_frames/$n");g["covariance_e2",compress=6]=op.covariance;g["sigma_ions"]=op.sigma
      push!(lengths,n)
     end
     linear=op.F\Float64.(read(e["linear_estimated_ions"]))
     mean=op.F\Float64.(read(e["noiseless_estimated_ions"]))
     noise=op.F\Float64.(read(e["noise_ions"]))
     measured=op.F\Float64.(read(e["recovered_ions"]));recovered=mean.+noise
     err=maximum(abs.(linear.-counts)./op.sigma)
     additive=maximum(abs.(recovered.-measured)./op.sigma)
     @assert err<s.c["numerical_tolerance_ENC"] && additive<s.c["numerical_tolerance_ENC"]
     @test err<s.c["numerical_tolerance_ENC"]
     @test additive<s.c["numerical_tolerance_ENC"]
     @test abs(sum(linear)/sum(counts)-1)<1e-6
     @assert sum(counts)==read(e["initial_pairs"])
     max_linear_error=max(max_linear_error,err);max_additivity_error=max(max_additivity_error,additive);input_total+=sum(counts)
     z=noise./op.sigma;zsum+=sum(z);zsquare+=sum(abs2,z);zcount+=length(z)
     uncertainty=noise_statistic_variances(op.covariance,npads)
     mean_stat_variance+=uncertainty.sum;square_stat_variance+=uncertainty.squares;lag_stat_variance+=uncertainty.lag
     lag_sum+=sum(z[1:end-1,:].*z[2:end,:]);lag_count+=(n-1)*npads
     lag_expected+=npads*sum(diag(op.covariance,1)./(op.sigma[1:end-1].*op.sigma[2:end]))
     HDF5.copy_object(e,out,"events/$id");g=out["events/$id"];HDF5.attrs(g)["schema"]="nausicaa.ion_csa_event.v2"
     replace(g,"noise_sigma_ions",op.sigma)
     for (field,array) in ("linear_estimated_ions"=>linear,"noiseless_estimated_ions"=>mean,"noise_ions"=>noise,"recovered_ions"=>recovered)
      stored=Float32.(array);rounding=maximum(abs.(Float64.(stored).-array)./op.sigma)
      @assert rounding<s.c["numerical_tolerance_ENC"]
      @test rounding<s.c["numerical_tolerance_ENC"]
      max_storage_ions=max(max_storage_ions,maximum(abs.(Float64.(stored).-array)))
      replace(g,field,stored)
     end
     negative+=count(<(0),recovered)
     loss=max.(linear.-mean,0.);relative=loss./max.(linear,eps());flags=zeros(Int8,size(counts))
     for (j,fraction) in enumerate(s.c["compression_fractions"])
      mask=(counts.>0).&(relative.>=fraction).&(loss.>=s.c["compression_minimum_deficit_ENC"].*op.sigma)
      flags[mask].=j
     end
     replace(g,"compression_level",flags)
    end
    @assert input_total==read(parent["metadata/input_ion_total"])
    normalized_mean=zsum/zcount;normalized_variance=zsquare/zcount
    measured_lag=lag_sum/lag_count;expected_lag=lag_expected/lag_count
    # Small event subsets need their covariance-derived statistical uncertainty.
    mean_limit=max(s.c["sampling_tolerance"],6sqrt(mean_stat_variance)/zcount)
    variance_limit=max(s.c["sampling_tolerance"],6sqrt(square_stat_variance)/zcount)
    lag_limit=max(s.c["sampling_tolerance"],6sqrt(lag_stat_variance)/lag_count)
    @assert abs(normalized_mean)<mean_limit
    @assert abs(normalized_variance-1)<variance_limit
    @assert abs(measured_lag-expected_lag)<lag_limit
    replace(m,"noise_mean_validation_limit",mean_limit);replace(m,"noise_variance_validation_limit",variance_limit);replace(m,"noise_lag1_validation_limit",lag_limit)
    replace(m,"parent_response_seconds",read(parent["metadata/response_seconds"]))
    for (k,v) in ("negative_estimates"=>negative,"maximum_storage_error"=>max_storage_ions,"maximum_storage_error_ions"=>max_storage_ions,"maximum_linear_recovery_error_ENC"=>max_linear_error,"maximum_transform_additivity_error_ENC"=>max_additivity_error,"normalised_noise_mean"=>normalized_mean,"normalised_noise_variance"=>normalized_variance,"normalised_noise_lag1_measured"=>measured_lag,"normalised_noise_lag1_expected"=>expected_lag,"response_seconds"=>(time_ns()-started)/1e9)
     replace(m,k,v)
    end
    write(m["complete"],true)
    println(source.sample," corrected; noise variance=",normalized_variance," lag1=",measured_lag," expected=",expected_lag);flush(stdout)
   end
   end
   mv(output*".tmp",output)
   open(joinpath(dirname(output),"manifest.toml"),"w") do io
    TOML.print(io,Dict("response_id"=>source.nodeid,"model_id"=>source.modelid,"parent_response"=>source.parentpath,"parent_sha256"=>source.parenthash,"response_file"=>output,"reconstruction_card"=>relpath(s.cp,s.root),"selected_event_ids"=>source.selected,"complete"=>true))
   end
  end
 end
 end
 end
 for p in s.files
  target=joinpath(s.dir,"source_snapshot",sha(p),basename(p));mkpath(dirname(target));GuardedFiles.checked_copy(p,target;force=true)
 end
 open(joinpath(s.dir,"manifest.toml"),"w") do io
  TOML.print(io,Dict("study_id"=>s.studyid,"model_id"=>s.modelid,"code_hash"=>s.codehash,"run_card"=>relpath(s.path,s.root),"parent_study_directory"=>s.r["parent_study"],"reconstruction_method"=>"three_frame_inverse","reconstruction_card"=>relpath(s.cp,s.root),"parent_response_paths"=>[x.parentpath for x in s.sources],"parent_response_hashes"=>[x.parenthash for x in s.sources]))
 end
 context=(;studyid=s.studyid,modelid=s.modelid,c=s.c,sources=s.sources)
 IonCSASummary.run_context(IonCSACorrection,context,s.dir,relpath(s.path,s.root))
 println(s.dir)
end
end
