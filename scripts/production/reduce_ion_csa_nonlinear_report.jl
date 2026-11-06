# Reduce held-out validation and the end-to-end pilot for presentation.
using HDF5,TOML,Statistics,SHA
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/ion_signal/ion_csa_pipeline.jl"))
const P=IonCSAProduction
sha(p)=open(io->bytes2hex(sha256(io)),p)
directory=joinpath(ROOT,"results/ion_csa_nonlinear/v1")
calpath=joinpath(directory,"calibration.h5")
latest=TOML.parsefile(joinpath(ROOT,"results/ion_csa_nonlinear/integration/latest.toml"))["manifest"]
pilot=TOML.parsefile(latest)
reference_manifest=joinpath(ROOT,"results/elm_matched_production/ion_csa/359613d487dfeee8/manifest.toml")
reference=TOML.parsefile(reference_manifest)
newpath=only(pilot["responses"])["path"];oldpath=only(reference["responses"])["path"]
newrows=NamedTuple[];padcounts=zeros(Int64,3)
h5open(newpath,"r") do new
h5open(oldpath,"r") do old
 ids=read(new["metadata/selected_event_ids"])
 for id in ids
    a=P.read_event(new,id);b=P.read_event(old,id)
    @assert a.input_ions==b.input_ions && a.pixel_ix==b.pixel_ix && a.pixel_iy==b.pixel_iy
    @assert a.edge_time_s==b.edge_time_s && a.pair_checksum==b.pair_checksum
    @assert a.noise_sigma_ions≈b.noise_sigma_ions
    delta=Float64.(a.noiseless_ions).-Float64.(b.noiseless_ions)
    nonlinear=a.model_used.==2
    @assert all(b.model_used[nonlinear].==1)
    nl=delta[:,nonlinear]
    positive=a.input_ions.>0
    # Reuse exactly the same noise draw when comparing final reconstructed signals.
    common_noise=Float64.(b.noise_ions)
    noisy_a=Float64.(a.noiseless_ions).+common_noise
    noisy_b=Float64.(b.noiseless_ions).+common_noise
    @assert maximum(abs.((noisy_a.-noisy_b).-delta))<1e-8
    push!(newrows,(;event_id=id,error_ENC=maximum(abs.(delta)./a.noise_sigma_ions),
        charge_difference_ions=sum(delta),reference_charge_ions=sum(Float64.(b.noiseless_ions)),
        nonlinear_charge_difference_ions=sum(nl),
        nonlinear_error_ENC=isempty(nl) ? 0. : maximum(abs.(nl)./a.noise_sigma_ions),
        maximum_pulse_difference_ions=maximum(a.noiseless_ions[positive])-maximum(b.noiseless_ions[positive]),
        noisy_maximum_difference_ions=maximum(noisy_a[positive])-maximum(noisy_b[positive])))
    for j in 0:2;padcounts[j+1]+=count(==(j),a.model_used);end
 end
end
end
@assert maximum(x.error_ENC for x in newrows)<.5
# A changed set of SPICE controls can change a few above-guard pads as well as the nonlinear replacement.
summary=Dict{String,Any}("pilot_events"=>length(newrows),"pilot_linear_pads"=>padcounts[1],
    "pilot_spice_pads"=>padcounts[2],"pilot_nonlinear_pads"=>padcounts[3],
    "pilot_maximum_error_ENC"=>maximum(x.error_ENC for x in newrows),
    "pilot_maximum_charge_difference_ions"=>maximum(abs(x.charge_difference_ions) for x in newrows),
    "pilot_mean_charge_difference_ions"=>mean(x.charge_difference_ions for x in newrows),
    "pilot_maximum_pulse_difference_ions"=>maximum(abs(x.maximum_pulse_difference_ions) for x in newrows))
summary["pilot_nonlinear_maximum_error_ENC"]=maximum(x.nonlinear_error_ENC for x in newrows)
summary["pilot_nonlinear_maximum_charge_difference_ions"]=maximum(abs(x.nonlinear_charge_difference_ions) for x in newrows)
frames=NamedTuple[]
h5open(calpath,"r") do h
    v=h["validation"];valid=read(v["valid"]);width=read(v["frame_ms"]);errors=read(v["error_ENC"])
    seconds=read(v["seconds"]);deployed=read(v["deployed"]);refine=read(v["refinement_ENC"])
    for ms in sort(unique(width))
        mask=width.==ms;accepted=mask.&valid
        push!(frames,(;frame_ms=ms,histories=count(mask),valid=count(accepted),fallbacks=count(mask.&.!valid),
            maximum_error_ENC=maximum(errors[accepted]),median_error_ENC=median(errors[accepted]),
            median_seconds=median(seconds[accepted]),deployed=count(mask.&deployed)))
    end
    summary["training_histories"]=length(read(h["training/event_ids"]))
    summary["validation_events"]=length(unique(zip(read(v["study"]),read(v["event"]))))
    # Gas/sample groups are included separately in the saved event aggregates.
    summary["validation_event_groups"]=length(read(h["events/event_id"]))
    summary["maximum_excursion_V"]=last(read(h["knots_V"]))
    summary["maximum_current_A"]=read(h["maximum_current_A"])
    summary["guard_V"]=read(h["guard_V"]);summary["step_s"]=read(h["step_s"])
    summary["maximum_pulse_difference_ions"]=maximum(abs.(read(v["maximum_pulse_error_ions"])[valid]))
    summary["maximum_event_charge_difference_ions"]=maximum(abs.(read(h["events/charge_difference_ions"])))
    summary["maximum_refinement_ENC"]=maximum(filter(isfinite,refine))
end
h5open(joinpath(directory,"report.h5.tmp"),"w") do h
    h["calibration_path"]=calpath;h["calibration_sha256"]=sha(calpath)
    h["pilot_manifest"]=latest;h["pilot_response_sha256"]=sha(newpath)
    h["reference_manifest"]=reference_manifest;h["reference_response_sha256"]=sha(oldpath)
    h["script_sha256"]=sha(@__FILE__)
    for (k,v) in summary;h["summary/"*k]=v;end
    for k in propertynames(first(frames));h["frames/"*String(k)]=getproperty.(frames,k);end
    for k in propertynames(first(newrows));h["pilot/"*String(k)]=getproperty.(newrows,k);end
    h["passed"]=true
end
mv(joinpath(directory,"report.h5.tmp"),joinpath(directory,"report.h5");force=true)
open(joinpath(directory,"report_summary.toml"),"w") do io;TOML.print(io,summary);end
println(summary)
