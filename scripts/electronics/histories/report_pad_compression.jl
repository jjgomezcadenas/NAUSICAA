using HDF5,Statistics,Test,SHA
length(ARGS)==1 || error("Usage: julia --project=. scripts/electronics/histories/report_pad_compression.jl RESULTS_DIRECTORY")
dir=abspath(ARGS[1]);source=joinpath(dir,"statistics.h5")
h5open(source,"r") do input
 h5open(joinpath(dir,"presentation.h5"),"w") do out
  attributes(out)["statistics_sha256"]=bytes2hex(sha256(read(source)))
  attributes(out)["code_sha256"]=bytes2hex(sha256(read(@__FILE__)))
  losses=read(input["criterion_loss_fraction"]);noise=read(input["criterion_noise_multiple"])
  out["criterion_loss_fraction"]=losses;out["criterion_noise_multiple"]=noise
  @testset "Event compression accounting" begin
   for name in keys(input)
    name in ("criterion_loss_fraction","criterion_noise_multiple") && continue
    g=input[name];a=attributes(g)
    @test read(a["completed"])==read(a["planned"])
    @test read(a["candidate_completed"])==read(a["candidate_planned"])
    @test read(a["control_completed"])==read(a["control_planned"])
    h=create_group(out,name)
    for key in keys(a);attributes(h)[key]=read(a[key]);end
    relative=read(g["pad_maximum_relative_loss"]);stratum=read(g["stratum"])
    h["control_stratum"]=collect(1:4)
    h["control_maximum_relative_loss"]=[maximum(relative[stratum.==s]) for s in 1:4]
    h["control_pads_with_flag"]=read(g["control_pads_with_flag"])
    controlflags=read(g["control_pads_with_flag"]);oneenc=findall(noise.==1.)
    attributes(h)["control_pads_with_10percent_and_1ENC_loss"]=controlflags[first(oneenc)]
    errors=read(g["refinement_error"])
    attributes(h)["refined_histories"]=length(errors)
    attributes(h)["maximum_refinement_error"]=isempty(errors) ? NaN : maximum(errors)
    @test isempty(errors) || maximum(errors)<.005
    for j in eachindex(losses)
     cg=g["criterion_$j"];q=read(cg["affected_ion_fraction"]);p=read(cg["affected_pads"]);b=read(cg["affected_buckets"])
     @test all(0 .<=q.<=1)
     @test all(b.>=p)
     @test count(>(0),p)==read(attributes(cg)["event_count"])
     dest=create_group(h,"criterion_$j")
     for key in keys(attributes(cg));attributes(dest)[key]=read(attributes(cg)[key]);end
     affected=q[q.>0]
     dest["affected_ion_fraction_quantiles"]=isempty(affected) ? zeros(3) : quantile(affected,[.5,.9,1.])
     thresholds=[0.,.01,.05,.1,.25,.5]
     dest["event_charge_threshold"]=thresholds
     dest["event_fraction_above_charge_threshold"]=[count(>(v),q)/length(q) for v in thresholds]
     for v in unique(noise)
      indices=findall(noise.==v)
      vals=[read(attributes(g["criterion_$i"])["event_count"]) for i in indices]
      @test issorted(vals;rev=true)
     end
    end
   end
  end
 end
end
println(joinpath(dir,"presentation.h5"))
