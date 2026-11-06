# Included inside IonCSAResponse; orchestration uses the pure nonlinear module.

"""Read the finite-path piecewise-linear current saved for one pad, in s and A."""
function read_case_current(path)
 t=Float64[];i=Float64[]
 for line in eachline(path)
  values=split(line);push!(t,parse(Float64,values[1]));push!(i,parse(Float64,values[2]))
 end
 t,i
end

"""
    nonlinear_cases(root, runpath, directory)

Evaluate selected pad histories with the validated stateful CSA. Below-guard
histories inside the current/voltage calibration domain become `response.h5`.
Control histories receive an independent prediction in `nonlinear_control.h5`
and continue through SPICE. Frame widths outside validation and histories outside
the numerical domain remain SPICE jobs, with a recorded fallback reason.

Use the exact saved finite-path current, baseline, 5 ms window and filtered
estimator. Cache identity includes the calibration checksum and source current.
The later three-frame inverse and full temporal noise covariance are unchanged.
"""
function nonlinear_cases(root,runpath,dir)
 s=settings(root,runpath)
 get(s.c,"response_model","hybrid_spice")=="hybrid_nonlinear" || return
 k=NL.load_calibration(joinpath(root,s.c["nonlinear_calibration"]);expected_sha256=s.c["nonlinear_calibration_sha256"])
 accepted=0;controls_count=0;fallbacks=Dict{String,String}();started=time_ns()
 h5open(joinpath(dir,"plan.h5"),"r") do plan
 h5open(joinpath(dir,"calibration.h5"),"r") do cal
  controls=read(plan["pads/control"])
  for job in TOML.parsefile(joinpath(dir,"jobs.toml"))["jobs"]
   folder=job["folder"];control=controls[job["row"]]
   dest=joinpath(folder,control ? "nonlinear_control.h5" : "response.h5")
   currentpath=joinpath(folder,"current.pwl")
   if isfile(dest)
    h5open(dest,"r") do h
     a=attributes(h)
     if haskey(a,"nonlinear_calibration_sha256")
      read(a["nonlinear_calibration_sha256"])==k.sha256 && read(a["current_sha256"])==sha(currentpath) || error("Nonlinear case cache mismatch")
      control ? (controls_count+=1) : (accepted+=1)
     end
    end
    continue
   end
   input=TOML.parsefile(joinpath(folder,"input.toml"));frame=input["frame_s"];ms=round(Int,1000frame)
   if !(ms in k.frame_ms) || !isapprox(frame,ms/1000;atol=1e-12)
    fallbacks[job["key"]]="frame_outside_calibration";continue
   end
   ca=attrs(cal["$(ms)ms"])
   isapprox(ca["baseline_V"],k.model.baseline_V;atol=1e-10) || error("Nonlinear baseline mismatch")
   t,i=read_case_current(currentpath)
   if maximum(i)>k.maximum_current_A*(1+1e-12)
    fallbacks[job["key"]]="current_outside_calibration";continue
   end
   info=caseinfo(plan,job["row"]);e=plan[info.group]
   edges=read(e["edge_time_s"]).-input["offset_s"].+input["prehistory_s"]
   r=NL.N.sampled_response(k.model,t,i,edges;averaging_s=s.c["averaging_s"],alpha=ca["alpha"],gain_V_per_ion=ca["gain_V_per_ion"])
   if !r.valid;fallbacks[job["key"]]=r.wave.reason;continue;end
   charge_error=abs(B.integral(t,i)/input["expected_charge_C"]-1)
   charge_error<1e-8 || error("Nonlinear input charge ledger failed")
   linear=e["linear_estimated_ions"][:,info.column]
   h5open(dest*".tmp","w") do out
    for (key,value) in ("response_model"=>"calibrated_nonlinear","nonlinear_calibration_sha256"=>k.sha256,
       "current_sha256"=>sha(currentpath),"charge_error"=>charge_error,"model_sha256"=>k.circuit_sha256,
       "minimum_output_V"=>k.model.baseline_V-r.wave.maximum_excursion_V,
       "maximum_error_ENC"=>maximum(abs.(r.estimated_ions.-linear))/ca["ENC_e_rms"])
     attributes(out)[key]=value
    end
    out["estimated_ions"]=r.estimated_ions;out["output_at_edges_V"]=r.voltage_V;out["average_at_edges_V"]=r.average_V
   end
   mv(dest*".tmp",dest)
   control ? (controls_count+=1) : (accepted+=1)
  end
 end
 end
 open(joinpath(dir,"nonlinear_summary.toml"),"w") do io
  TOML.print(io,Dict("nonlinear_histories"=>accepted,"control_predictions"=>controls_count,
     "fallback_histories"=>length(fallbacks),"fallback_reasons"=>fallbacks,
     "calibration_sha256"=>k.sha256,"elapsed_s"=>(time_ns()-started)/1e9))
 end
 println("Nonlinear histories=",accepted,"; SPICE controls=",controls_count,"; domain fallbacks=",length(fallbacks));flush(stdout)
end

"""
    validate_nonlinear_controls(root, runpath, directory)

Compare independent SPICE controls with their saved nonlinear predictions after
the three-frame inverse. Normalize every bucket by its corrected noise sigma.
Require the configured 0.5 ENC tolerance and save the measured residuals, charge
biases and calibration identity. SPICE remains the output for these controls.
"""
function validate_nonlinear_controls(root,runpath,dir)
 s=settings(root,runpath)
 get(s.c,"mode","validation")=="production" && return
 get(s.c,"response_model","hybrid_spice")=="hybrid_nonlinear" || return
 rows=Int[];errors=Float64[];bias=Float64[];operators=Dict{Tuple{String,Int},Any}()
 h5open(joinpath(dir,"plan.h5"),"r") do plan
 h5open(joinpath(dir,"calibration.h5"),"r") do cal
  for job in TOML.parsefile(joinpath(dir,"jobs.toml"))["jobs"]
   p=joinpath(job["folder"],"nonlinear_control.h5");isfile(p) || continue
   info=caseinfo(plan,job["row"]);frame=read(attributes(plan[info.name])["frame"]);key="$(round(Int,1000frame))ms"
   pred=h5open(h->read(h["estimated_ions"]),p,"r")
   ref=h5open(h->read(h["estimated_ions"]),joinpath(job["folder"],"response.h5"),"r")
   op=get!(operators,(key,length(pred))) do;NC.operator(NC.taps(cal[key]),read(cal[key*"/covariance_e2"]),length(pred));end
   delta=op.F\(pred.-ref);push!(rows,job["row"]);push!(errors,maximum(abs.(delta)./op.sigma));push!(bias,sum(delta))
  end
 end
 end
 isempty(rows) && error("No nonlinear controls inside the calibration domain")
 h5open(joinpath(dir,"nonlinear_control_validation.h5"),"w") do h
  h["pad_row"]=rows;h["maximum_error_ENC"]=errors;h["charge_difference_ions"]=bias
  h["calibration_sha256"]=s.c["nonlinear_calibration_sha256"]
  h["passed"]=maximum(errors)<s.c["linear_error_tolerance_ENC"]
 end
 maximum(errors)<s.c["linear_error_tolerance_ENC"] || error("Nonlinear production controls exceed the requested accuracy")
 println("Nonlinear corrected controls: ",length(rows),"; maximum ",maximum(errors)," ENC");flush(stdout)
end
