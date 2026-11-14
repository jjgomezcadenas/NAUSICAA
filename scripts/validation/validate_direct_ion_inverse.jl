using Test
include("../production/run_framed_reconstruction.jl")
"""Compare direct and iterative inverses using noiseless calibrated waveforms.
Five saved pad shapes are tested at several amplitudes. ADC quantisation is
checked separately, without CSA or ADC random noise. No full event simulation
or noise-bias study is performed."""
function validate_direct()
 BLAS.set_num_threads(1)
 card=TOML.parsefile(joinpath(ROOT,"config/readout/bb0nu_electronics_pilot_100.toml"))
 prepared=prepare_electronics_run(card)
 dc=TOML.parsefile(joinpath(ROOT,card["detector_card"]))
 base=joinpath(ROOT,"results/readout",dc["run_label"])
 source=joinpath(base,"electronics_100events_efficient/controls.h5")
 outdir=joinpath(base,"reconstruction/direct_inverse");mkpath(outdir)
 target=joinpath(outdir,"controls.h5");isfile(target)&&error("Preserve existing direct-inverse controls")
 rows=NamedTuple[];examples=Dict{String,Any}()
 @testset "Calibrated nonlinear shapes without added noise" begin
 h5open(source,"r") do f
  ids=read(f["events/event_number"]);starts=read(f["events/ion_start_ms"])
  for event in keys(f["examples"]),name in keys(f["examples/$event"])
   startswith(name,"ions")||continue
   g=f["examples/$event/$name"];q=Float64.(read(g["input"]));sum(q)>0||continue
   ts=Float64.(read(g["time_ms"])).*.001
   start=starts[findfirst(==(parse(Int,event)),ids)]
   it=IonReadoutTime(Float32(start+first(read(g["input_time_ms"]))-10),20f0,Int32(length(q)))
   at=IonADCTime(start,100f0,Int64(length(ts)))
   setup=electronics_settings(prepared,it)
   cs=setup.ion.csa
   cal=IonReconstructionCalibration(setup.ion.adc,cs.current_pulse,cs.linear_response,cs.nonlinear_response,cs.noise)
   settings=IonReconstructionSettings(10f0,5f0,.63f0,20,.005f0)
   plan_time=@timed SR.prepare_ion_reconstruction(at,it,cal,settings);p=plan_time.value
   direct_setup=@timed SR.prepare_direct_inverse!(p)
   physical=ion_current(Int64.(q),it,cal.current_pulse,at)
   current_t=vcat(physical.time_s,max(last(physical.time_s),last(p.times))+1e-12)
   current=vcat(physical.current_A,0.)
   for scale in (.1,.25,.5,1.,1.25)
    truth=q.*scale
    forward=nonlinear_waveform(cs.nonlinear_response,current_t,current.*scale;sample_time_s=p.times)
    @test forward.valid
    forward.valid||continue
    x=cs.nonlinear_response.baseline_V.-forward.voltage_V
    direct=SR.direct_ion_inverse(p,x)
    @test direct.valid
    @test maximum(abs.(direct.ions-truth))<.5
    @test abs(sum(direct.ions)-sum(truth))<.5
    y=[mean(view(x,g))-mean(view(x,p.groups[1])) for g in p.groups[2:end]]
    old=@timed SR.fit_nonlinear(p,y,p.inverse*y)
    @test old.value.valid
    @test maximum(abs.(direct.ions-old.value.ions))<.5
    # Quantisation only: no random electronics noise is injected.
    codes=digitise_adc(Xoshiro(0),forward.voltage_V,cal.adc).counts
    lsb=(Float64(cal.adc.maximum_V)-cal.adc.minimum_V)/2^cal.adc.bits
    drop=(mean(view(codes,p.groups[1])).-Float64.(codes)).*lsb
    quantised=SR.direct_ion_inverse(p,drop)
    @test quantised.valid
    # Existing electronics comparison target: 0.5 ENC, approximately 25 ions.
    @test maximum(abs.(quantised.ions-truth))<25.
    qy=[mean(view(drop,g))-mean(view(drop,p.groups[1])) for g in p.groups[2:end]]
    old_quantised=SR.fit_nonlinear(p,qy,p.inverse*qy)
    @test old_quantised.valid
    # Warm paths and time repeated direct inversions of the same observed vector.
    SR.direct_ion_inverse(p,x)
    times=Float64[];allocated=0
    for _ in 1:30
     timed=@timed SR.direct_ion_inverse(p,x);push!(times,timed.time);allocated=timed.bytes
    end
    push!(rows,(event_number=parse(Int,event),channel=parse(Int,split(name,'_')[2]),scale,
     buckets=length(q),samples=length(x),peak_voltage_drop_V=maximum(x),
     max_bucket_error_ions=maximum(abs.(direct.ions-truth)),total_error_ions=sum(direct.ions)-sum(truth),
     max_quantised_error_ions=maximum(abs.(quantised.ions-truth)),
     max_iterative_quantised_error_ions=maximum(abs.(old_quantised.ions-truth)),
     max_difference_from_fit_ions=maximum(abs.(direct.ions-old.value.ions)),
     direct_s=median(times),iterative_s=old.time,allocated_bytes=allocated,
     base_preparation_s=plan_time.time,direct_preparation_s=direct_setup.time))
    if scale==1
     examples[event]=Dict("input"=>truth,"direct"=>direct.ions,"iterative"=>old.value.ions,"quantised"=>quantised.ions,
       "time_ms"=>(collect(eachindex(q)).-.5).*it.bucket_ms)
    end
   end
  end
 end
 end
 budget=4*1024^2;GuardedFiles.space_check(outdir,budget,1024^3);tmp=target*".tmp"
 try
  h5open(tmp,"w") do f
   f["schema"]="itaca.direct_ion_inverse_validation.v1"
   f["source_sha256"]=bytes2hex(open(sha256,source))
   f["calibration_sha256"]=bytes2hex(open(sha256,joinpath(ROOT,card["nonlinear_calibration"])))
   f["added_noise"]=false;columns(create_group(f,"cases"),rows)
   write_control_tree(create_group(f,"examples"),examples)
   f["summary/cases"]=length(rows)
   f["summary/median_direct_s"]=median(r.direct_s for r in rows)
   f["summary/median_iterative_s"]=median(r.iterative_s for r in rows[2:end])
   f["summary/max_bucket_error_ions"]=maximum(r.max_bucket_error_ions for r in rows)
   f["summary/max_total_error_ions"]=maximum(abs(r.total_error_ions) for r in rows)
   f["summary/max_quantised_error_ions"]=maximum(r.max_quantised_error_ions for r in rows)
   f["summary/max_iterative_quantised_error_ions"]=maximum(r.max_iterative_quantised_error_ions for r in rows)
   f["summary/max_difference_from_fit_ions"]=maximum(r.max_difference_from_fit_ions for r in rows)
   files=[@__FILE__,joinpath(ROOT,"src/reconstruction/ions/direct_inverse.jl"),joinpath(ROOT,"src/reconstruction/ions/reconstruction.jl")]
   f["code_files"]=relpath.(files,ROOT);f["code_sha256"]=[bytes2hex(open(sha256,p)) for p in files]
   flush(f);GuardedFiles.check_size(tmp,budget)
  end
  mv(tmp,target)
 catch;isfile(tmp)&&rm(tmp);rethrow();end
 println(target)
end
validate_direct()
