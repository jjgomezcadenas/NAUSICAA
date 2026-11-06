# Fit and independently validate the one-state CSA using immutable acquisitions.
using HDF5,TOML,LinearAlgebra,Statistics,SHA,Printf
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/ion_signal/ion_csa_nonlinear.jl"))
include(joinpath(ROOT,"src/ion_signal/ion_csa_correction.jl"))
using .IonCSANonlinear
const C=IonCSACorrection
sha(p)=open(io->bytes2hex(sha256(io)),p)
card=isempty(ARGS) ? joinpath(ROOT,"config/specs/response/ion_csa/nonlinear_calibration_v1.toml") : abspath(only(ARGS))
cfg=TOML.parsefile(card);outdir=joinpath(ROOT,cfg["results_directory"])
mkpath(outdir);write(joinpath(outdir,".gitignore"),"*.h5\n*.tmp\n")
isfile(joinpath(outdir,"calibration.h5")) && error("Choose a fresh results_directory; preserve existing calibration")
function rows(study)
    jobs=TOML.parsefile(joinpath(study,"jobs.toml"))["jobs"]
    h5open(joinpath(study,"plan.h5"),"r") do h
        events=read(h["pads/event_id"]);minimums=read(h["pads/minimum_V"]);groups=read(h["pads/event_group"]);names=read(h["pads/name"]);cols=read(h["pads/column"])
        [(;job,event=events[job["row"]],minv=minimums[job["row"]],group=groups[job["row"]],name=names[job["row"]],column=cols[job["row"]]) for job in jobs]
    end
end
select(rows,n)=begin
    sorted=sort(rows;by=x->x.minv)
    n==0 ? sorted : sorted[unique(round.(Int,range(1,length(sorted);length=min(n,length(sorted)))))]
end
heldout(row)=mod(row.event,cfg["event_split_modulus"])==cfg["validation_remainder"]
function waveform(row)
    h5open(joinpath(row.job["folder"],"raw.h5"),"r") do h
        (;time_s=read(h["time"]),current_A=read(h["I(Is1)"]),voltage_V=read(h["V(out)"]),
          circuit_hash=String(read(attributes(h)["model_sha256"])),simulator_hash=String(read(attributes(h)["simulator_sha256"])))
    end
end
function current_pwl(path)
    times=Float64[];values=Float64[]
    for line in eachline(path)
        a=split(line);push!(times,parse(Float64,a[1]));push!(values,parse(Float64,a[2]))
    end
    times,values
end
study=joinpath(ROOT,cfg["training_study"])
training=select(filter(x->!heldout(x),rows(study)),cfg["training_histories"])
println("Loading ",length(training)," calibration histories; validation events have remainder ",cfg["validation_remainder"]);flush(stdout)
train=waveform.(training)
source=first(train)
length(unique(x.circuit_hash for x in train))==length(unique(x.simulator_hash for x in train))==1 || error("Training circuit/solver mismatch")
baseline=h5open(h->read(attributes(h["40ms"])["baseline_V"]),joinpath(study,"calibration.h5"),"r")
maximum_excursion=maximum(baseline-minimum(w.voltage_V) for w in train)
knots=collect(0.:cfg["knot_spacing_V"]:maximum_excursion)
last(knots)<maximum_excursion && push!(knots,maximum_excursion)
fit=fit_nonlinear_csa(train,knots;baseline_V=baseline,step_s=cfg["maximum_step_s"],ridge=cfg["ridge"])
model=fit.model;maxcurrent=maximum(maximum(w.current_A) for w in train)
fine=NonlinearCSA(knots,fit.gain_V_per_e,fit.recovery_V_s;baseline_V=baseline,step_s=model.step_s/2)
println("Fitted ",fit.points," intervals; excursion range ",last(knots)," V, current bound ",maxcurrent," A");flush(stdout)
train=nothing;GC.gc()
function calibration(cal,key)
    g=cal[key];a=C.attrs(g)
    (;a,taps=C.taps(g),cv=read(g["covariance_e2"]))
end
records=NamedTuple[];casekeys=String[];sourcepaths=String[];sourcesha=String[];eventstats=Dict{Tuple{Int,String,Int},Vector{Float64}}();domains=String[]
for (studyindex,relative) in enumerate(cfg["validation_studies"])
    directory=joinpath(ROOT,relative);selected=select(filter(heldout,rows(directory)),cfg["validation_limit_per_study"])
    isempty(intersect(Set(x.event for x in training),Set(x.event for x in selected))) || error("Event leakage")
    operators=Dict{Tuple{String,Int},Any}();calibrations=Dict{String,Any}()
    # Validate every held-out SPICE candidate in the selected study; unchanged linear pads cancel in event differences.
    h5open(joinpath(directory,"plan.h5"),"r") do plan
    h5open(joinpath(directory,"calibration.h5"),"r") do cal
        for (j,row) in enumerate(selected)
            folder=row.job["folder"];input=TOML.parsefile(joinpath(folder,"input.toml"));frame=input["frame_s"];key="$(round(Int,1000frame))ms"
            k=get!(calibrations,key) do;calibration(cal,key);end
            isapprox(k.a["baseline_V"],baseline;atol=1e-10) || error("Calibration baseline mismatch")
            rawpath=joinpath(folder,"raw.h5")
            h5open(rawpath,"r") do h
                attrs=attributes(h)
                String(read(attrs["model_sha256"]))==source.circuit_hash || error("Validation circuit differs")
                String(read(attrs["simulator_sha256"]))==source.simulator_hash || error("Validation simulator differs")
            end
            t,i=current_pwl(joinpath(folder,"current.pwl"))
            edges=read(plan[row.group*"/edge_time_s"]).-input["offset_s"].+input["prehistory_s"]
            before=time_ns()
            response=maximum(i)<=maxcurrent*(1+1e-12) ? sampled_response(model,t,i,edges;averaging_s=cfg["averaging_s"],alpha=k.a["alpha"],gain_V_per_ion=k.a["gain_V_per_ion"]) : (;valid=false)
            seconds=(time_ns()-before)/1e9;err=NaN;bias=0.;maxpulseerror=NaN;refinement=NaN;reference_sum=NaN
            if response.valid
                ref=h5open(h->read(h["estimated_ions"]),joinpath(folder,"response.h5"),"r")
                n=length(ref)
                op=get!(operators,(key,n)) do;C.operator(k.taps,k.cv,n);end
                truth=op.F\ref;predicted=op.F\response.estimated_ions;delta=predicted.-truth
                err=maximum(abs.(delta)./op.sigma);bias=sum(delta);reference_sum=sum(truth)
                maxpulseerror=maximum(predicted)-maximum(truth)
                if j<=20 || j%100==0
                    other=sampled_response(fine,t,i,edges;averaging_s=cfg["averaging_s"],alpha=k.a["alpha"],gain_V_per_ion=k.a["gain_V_per_ion"])
                    other.valid || error("Fine integration leaves domain")
                    refinement=maximum(abs.(op.F\(response.estimated_ions.-other.estimated_ions))./op.sigma)
                end
            end
            # Model is deployed only below the same guard. Above it, preserve the linear path.
            deployed=response.valid && row.minv<cfg["guard_V"]
            stats=get!(eventstats,(studyindex,String(row.name),Int(row.event)),zeros(2));stats[1]+=deployed ? bias : 0.
            stats[2]=max(stats[2],deployed ? abs(maxpulseerror) : 0.)
            push!(records,(;study=studyindex,event=row.event,frame_ms=round(Int,1000frame),minv=row.minv,
                valid=response.valid,deployed,seconds,error_ENC=err,bias_ions=bias,maximum_pulse_error_ions=maxpulseerror,
                reference_sum_ions=reference_sum,refinement_ENC=refinement))
            push!(casekeys,row.job["key"]);push!(sourcepaths,rawpath);push!(sourcesha,sha(rawpath))
            key in domains || push!(domains,key)
            if j%250==0 || j==length(selected)
                println("Study ",studyindex," ",j,"/",length(selected)," latest error ",err," ENC; valid ",response.valid);flush(stdout)
            end
        end
    end
    end
end
valid=filter(x->x.valid,records);refined=filter(x->isfinite(x.refinement_ENC),records)
passed=!isempty(valid)&&!isempty(refined)&&maximum(x.error_ENC for x in valid)<cfg["maximum_error_ENC"]&&maximum(x.refinement_ENC for x in refined)<cfg["maximum_refinement_error_ENC"]
h5open(joinpath(outdir,"calibration.h5.tmp"),"w") do h
    attributes(h)["schema"]="nausicaa.ion_csa_nonlinear.v1"
    for (key,value) in ("knots_V"=>knots,"gain_V_per_e"=>fit.gain_V_per_e,"recovery_V_s"=>fit.recovery_V_s,"support_time_s"=>fit.support_time_s,
        "baseline_V"=>baseline,"step_s"=>model.step_s,"maximum_current_A"=>maxcurrent,"averaging_s"=>cfg["averaging_s"],"guard_V"=>cfg["guard_V"],
        "maximum_error_ENC"=>cfg["maximum_error_ENC"],"validated_frame_ms"=>sort(unique(x.frame_ms for x in records)),
        "circuit_sha256"=>source.circuit_hash,"simulator_sha256"=>source.simulator_hash,
        "model_code_sha256"=>sha(joinpath(ROOT,"src/ion_signal/ion_csa_nonlinear.jl")),"script_sha256"=>sha(@__FILE__),
        "card_toml"=>read(card,String),"passed"=>passed,"complete"=>true)
        h[key]=value
    end
    h["training/event_ids"]=[x.event for x in training];h["training/case_keys"]=[x.job["key"] for x in training]
    h["training/raw_sha256"]=[sha(joinpath(x.job["folder"],"raw.h5")) for x in training]
    for key in propertynames(first(records));h["validation/"*String(key)]=[getproperty(x,key) for x in records];end
    h["validation/case_keys"]=casekeys;h["validation/raw_paths"]=sourcepaths;h["validation/raw_sha256"]=sourcesha
    entries=sort!(collect(eventstats);by=first)
    h["events/study"]=[p.first[1] for p in entries];h["events/sample_group"]=[p.first[2] for p in entries];h["events/event_id"]=[p.first[3] for p in entries]
    h["events/charge_difference_ions"]=[p.second[1] for p in entries];h["events/maximum_pad_peak_difference_bound_ions"]=[p.second[2] for p in entries]
end
mv(joinpath(outdir,"calibration.h5.tmp"),joinpath(outdir,"calibration.h5"))
summary=Dict("passed"=>passed,"validation_histories"=>length(records),"valid_histories"=>length(valid),"fallback_histories"=>length(records)-length(valid),
    "maximum_error_ENC"=>maximum(x.error_ENC for x in valid),"median_error_ENC"=>median(x.error_ENC for x in valid),
    "maximum_refinement_ENC"=>maximum(x.refinement_ENC for x in refined),"median_model_seconds"=>median(x.seconds for x in valid),
    "maximum_event_charge_difference_ions"=>maximum(abs(x[1]) for x in values(eventstats)),
    "model_sha256"=>sha(joinpath(outdir,"calibration.h5")))
open(joinpath(outdir,"validation.toml"),"w") do io;TOML.print(io,summary);end
println(summary)
passed || error("Calibration has not met production acceptance")
