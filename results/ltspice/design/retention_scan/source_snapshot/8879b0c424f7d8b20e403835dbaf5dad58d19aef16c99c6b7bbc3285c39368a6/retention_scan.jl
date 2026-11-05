module RetentionScan
using HDF5,TOML,SHA,Statistics,Test
include("pad_history_readout.jl")
const P=PadHistoryReadout
const D=P.D
const B=D.B
const A=D.A
sha(p)=open(io->bytes2hex(sha256(io)),p)
settings(root,card)=TOML.parsefile(joinpath(root,card))
attrs(g)=Dict(k=>read(attributes(g)[k]) for k in keys(attributes(g)))

"""Write charge-preserving uniform-arrival currents and independent long-tail inputs."""
function inputs(root,card)
 c=settings(root,card);dir=joinpath(root,c["output"]);mkpath(dir)
 target=joinpath(dir,"inputs.h5")
 if isfile(target)
  @assert h5open(h->read(attributes(h)["config_sha256"]),target,"r")==sha(joinpath(root,card))
  return target
 end
 d=TOML.parsefile(joinpath(root,c["detector"]));a=d["geometry"]["pad_width_m"]/2;g=d["geometry"]["gap_m"]
 v=d["gas"]["bulk_speed_m_s"]/d["bias"]["bulk_field_V_m"]*(d["bias"]["ic_voltage_V"]-d["bias"]["pad_voltage_V"])/g
 h5open(target,"w") do h
  at=attributes(h)
  for (k,value) in ("config_sha256"=>sha(joinpath(root,card)),"detector_sha256"=>sha(joinpath(root,c["detector"])),"generator_sha256"=>sha(@__FILE__),"waveform_sha256"=>sha(joinpath(@__DIR__,"bucket_study.jl")),"single_ion_charge_C"=>B.E*(1-B.phi(g,a)),"pad_halfwidth_m"=>a,"gap_m"=>g,"gap_speed_m_s"=>v)
   at[k]=value
  end
  at["assumptions"]="Central 0.5 mm square pad in an infinite grounded plane; pad weighting potential 1, plane and infinity 0; homogeneous gas; central constant-speed 3 m/s trajectory through the finite 5 mm gap. Uniform ion collection within each frame; finite-path charge retained without renormalization or a collection impulse."
  for frame in c["frames_s"]
   key="$(round(Int,1000frame))ms"
   for (r,step) in enumerate(c["max_steps_s"])
    stop=4frame+c["averaging_s"]
    t,i=B.pwl(frame,a,g,v,c["prehistory_s"],stop,step)
    for tail in (false,true)
     # The current is already zero at the end of its compact support.
     tt=tail ? vcat(t,frame+c["tail_duration_s"]+c["prehistory_s"]) : t
     ii=tail ? vcat(i,0.) : i
     group=create_group(h,"$key/$(tail ? "tail" : "fast")_r$r")
     group["time_s"]=tt;group["single_A_per_ion"]=ii
     @assert abs(B.integral(tt,ii)/read(at["single_ion_charge_C"])-1)<1e-9
    end
   end
  end
 end
 target
end

"""Select a rounded M4 bias that restores the frozen circuit's quiescent output."""
function tune(root,card,python,simulator)
 c=settings(root,card);dir=joinpath(root,c["output"])
 target=h5open(h->only(read(h["op/V(out)"])),joinpath(root,c["reference_results"]),"r")
 points=joinpath(dir,"operating_points.h5")
 function acquire(name,gate)
  script=joinpath(root,"scripts/electronics/feedback/run_retention_scan.py")
  args=[python,script,"--config",card,"--simulator",simulator,"--stage","op","--variant",name,"--m4-gate",string(gate)]
  result=strip(read(Cmd(Cmd(args);dir=root),String))
  h5open(h->read(attributes(h)["baseline_V"]),joinpath(root,result),"r")
 end
 for variant in c["variants"]
  name=variant["name"]
  isfile(points) && h5open(h->haskey(h,name),points,"r") && continue
  if name=="length_400um"
   gate=.6164;value=acquire(name,gate)
  else
   lo=.45;hi=.9;vlo=acquire(name,lo);vhi=acquire(name,hi)
   @assert (vlo-target)*(vhi-target)<0 (name,vlo,vhi,target)
   gate=NaN;value=NaN
   for _ in 1:22
    gate=(lo+hi)/2;value=acquire(name,gate)
    abs(value-target)<1e-5 && break
    if (value-target)*(vlo-target)>0;lo=gate;vlo=value;else;hi=gate;vhi=value;end
   end
   center=round(gate;digits=c["bias_decimal_places"]);step=10.0^(-c["bias_decimal_places"])
   gates=[center-step,center,center+step];values=[acquire(name,x) for x in gates]
   j=argmin(abs.(values.-target));gate=gates[j];value=values[j]
  end
  @assert abs(value-target)<.02 (name,gate,value,target)
  h5open(points,"cw") do h
   group=create_group(h,name);a=attributes(group)
   a["length_m"]=variant["length_m"];a["m4_gate_V"]=gate;a["baseline_V"]=value
   a["target_baseline_V"]=target;a["tuner_sha256"]=sha(@__FILE__)
  end
  println(name," M4gate=",gate," baseline=",value);flush(stdout)
 end
end

"""Reduce matched transients and stationary SPICE spectra using the frame estimator."""
function analyze(root,card)
 c=settings(root,card);dir=joinpath(root,c["output"]);target=joinpath(dir,"analysis.h5")
 @testset "Long-retention charge, fit and stationary-noise checks" begin
 h5open(target,"w") do out
  attributes(out)["config_sha256"]=sha(joinpath(root,card));attributes(out)["analysis_sha256"]=sha(@__FILE__)
  attributes(out)["averaging_s"]=c["averaging_s"]
  attributes(out)["assumptions"]="Same 180 nm transistor card, 1.25 pF external input capacitance, 3 fF feedback and 1.8 V supply as frozen reference. Only feedback length and rounded M4 bias are varied; gate 0 V, width 0.36 um. Small-signal stationary noise at the quiescent point; 5 ms centred averages and decay-corrected frame differences. All spectral values are one-sided."
  for variant in c["variants"]
   name=variant["name"];path=joinpath(dir,"variants",name*".h5")
   h5open(path,"r") do raw
    vg=create_group(out,name)
    for (k,v) in attrs(raw);attributes(vg)[k]=v;end
    attributes(vg)["raw_sha256"]=sha(path)
    f=read(raw["noise_400/frequency"]);asd=read(raw["noise_400/V(onoise)"])
    fc=read(raw["noise_200/frequency"]);ac=read(raw["noise_200/V(onoise)"])
    devices=sort(filter(n->occursin(r"^V\(m\w+\)$",n),keys(raw["noise_400"])))
    vg["frequency_Hz"]=f;vg["output_asd_V_sqrtHz"]=asd
    for frame in c["frames_s"]
     key="$(round(Int,1000frame))ms";src=raw[key];pre=c["prehistory_s"];width=c["averaging_s"]
     group=create_group(vg,key);at=attributes(group)
     tail(n,r)=begin;t,v=D.waveform(src["tail"],n,r);t.-pre,v;end
     t,v=tail(100,2);tc,vc=tail(100,1);tl,vl=tail(10,2)
     start=frame+.02;stop=frame+c["tail_duration_s"]-.01
     tau,fiterr,_,_,_=D.fit_tail(t,v,start,stop,.002)
     taucoarse,_,_,_,_=D.fit_tail(tc,vc,start,stop,.002)
     taulow,_,_,_,_=D.fit_tail(tl,vl,start,stop,.002)
     @test abs(taucoarse/tau-1)<1e-3
     @test abs(taulow/tau-1)<.01
     @test fiterr<.001
     windows=[(frame+.02,frame+1.),(frame+1.,frame+3.),(frame+3.,stop)]
     localtau=[D.fit_tail(t,v,a,b,.002)[1] for (a,b) in windows]
     @test maximum(abs.(localtau./tau.-1))<.01
     alpha=exp(-frame/tau)
     measure(t,v)=A.waveform_gain(t.-pre,v,width,frame,alpha)
     wave(n,r)=D.waveform(src,n,r)
     t100,v100=wave(100,2);tco,vco=wave(100,1)
     gain=measure(t100,v100)/100
     sigma=A.noise_rms(f,asd,width,frame,alpha)
     spectral=abs(A.noise_rms(fc,ac,width,frame,alpha)/sigma-1)
     quadrature=abs(A.noise_rms(f,asd,width,frame,alpha;order=16)/sigma-1)
     transient=abs(measure(tco,vco)/(100gain)-1)
     @test gain>0
     @test spectral<1e-3
     @test quadrature<1e-5
     @test transient<1e-3
     spectra=[A.noise_rms(f[f.>=low],asd[f.>=low],width,frame,alpha) for low in (.0001,.001,.1)]
     charge_errors=[];gain_ratios=[]
     q=h5open(h->read(attributes(h)["single_ion_charge_C"]),joinpath(dir,"inputs.h5"),"r")
     for n in filter(>(0),c["counts"]),r in (1,2)
      tt,vv=wave(n,r);push!(charge_errors,abs(B.integral(tt,read(src["single_N$(n)_r$r/I(Is1)"]))/(n*q)-1))
      push!(gain_ratios,measure(tt,vv)/(n*gain))
     end
     @test maximum(charge_errors)<2e-5
     @test maximum(abs.(gain_ratios.-1))<.01
     vars=[A.noise_rms(f,max.(read(raw["noise_400"][device]),1e-150),width,frame,alpha)^2 for device in devices]
     @test abs(sum(vars)/sigma^2-1)<1e-3
     for (k,x) in ("frame_s"=>frame,"tau_s"=>tau,"alpha"=>alpha,"best_index"=>1,"tail_fit_error"=>fiterr,"tau_refinement_error"=>abs(taucoarse/tau-1),"tau_amplitude_error"=>abs(taulow/tau-1),"transient_error"=>transient,"spectral_error"=>spectral,"quadrature_error"=>quadrature,"charge_error"=>maximum(charge_errors),"gain_change"=>maximum(abs.(gain_ratios.-1)),"ENC_e_rms"=>sigma/gain)
      at[k]=x
     end
     group["width_s"]=[width];group["gain_V_per_ion"]=[gain];group["noise_rms_V"]=[sigma];group["ENC_e_rms"]=[sigma/gain]
     group["device_names"]=devices;group["device_enc_e"]=sqrt.(vars)./gain;group["device_variance_fraction"]=vars./sigma^2
     group["local_tau_s"]=localtau;group["noise_low_Hz"]=[.0001,.001,.1];group["ENC_by_low_cut_e"]=spectra./gain
     group["time_s"]=t100.-pre;group["voltage_V"]=v100
     j=1:10:length(t);group["tail_time_s"]=t[j];group["tail_voltage_V"]=v[j]
     println(name," ",key," tau=",tau," ENC=",sigma/gain," oldbandENC=",spectra[3]/gain," gain=",gain," local_tau=",localtau);flush(stdout)
    end
   end
  end
 end
 end
 target
end

"""Write history and compression cards for the calibrated feedback variants."""
function history_cards(root,card)
 c=settings(root,card);dir=joinpath(root,c["output"]);mkpath(joinpath(dir,"cards"))
 original=TOML.parsefile(joinpath(root,c["history_card"]))
 viewer=TOML.parsefile(joinpath(root,c["compression_card"]))
 for variant in c["variants"]
  name=variant["name"];r=deepcopy(original)
  r["output"]=joinpath(c["output"],"screens",name)
  r["report"]=joinpath(c["report"],"histories",name)
  r["circuit_results"]=joinpath(c["output"],"variants",name*".h5")
  r["circuit_analysis"]=joinpath(c["output"],"analysis.h5");r["circuit_name"]=name
  r["source_decks"]=joinpath(c["output"],"decks",name)
  r["averaging_s"]=c["averaging_s"];r["prehistory_s"]=c["prehistory_s"]
  path=joinpath(c["output"],"cards",name*"_readout.toml")
  open(io->TOML.print(io,r),joinpath(root,path),"w")
  v=deepcopy(viewer);v["readout_card"]=path
  v["output"]=joinpath(c["output"],"compression",name)
  v["report"]=joinpath(c["report"],"compression",name)
  v["input_description"]="Direct collection, longer-retention feedback scan: "*name
  open(io->TOML.print(io,v),joinpath(dir,"cards",name*"_viewer.toml"),"w")
 end
end

"""Screen every stored pad, preserving empty buckets and source charge accounting."""
function screens(root,card)
 c=settings(root,card);history_cards(root,card)
 for variant in c["variants"]
  cp=joinpath(root,c["output"],"cards",variant["name"]*"_readout.toml")
  r=TOML.parsefile(cp)
  all(isfile(joinpath(root,r["output"],s["name"]*".h5")) for s in r["samples"]) && continue
  withenv("NAUSICAA_PAD_HISTORY_CARD"=>cp) do;P.screen(root);end
 end
end

"""Collect calibrated noise and completed nonlinear event statistics for presentation."""
function summary(root,card)
 c=settings(root,card);dir=joinpath(root,c["output"])
 readout=TOML.parsefile(joinpath(root,c["history_card"]))
 @testset "Retention comparison coverage" begin
 h5open(joinpath(dir,"analysis.h5"),"r") do analysis
 h5open(joinpath(dir,"presentation.h5"),"w") do out
  attributes(out)["analysis_sha256"]=sha(joinpath(dir,"analysis.h5"))
  attributes(out)["config_sha256"]=sha(joinpath(root,card))
  attributes(out)["summary_sha256"]=sha(@__FILE__)
  out["variant_names"]=[v["name"] for v in c["variants"]]
  for variant in c["variants"]
   name=variant["name"];HDF5.copy_object(analysis[name],out,name)
   vg=out[name]
   directory=joinpath(dir,"compression",name)
   statistics=[joinpath(directory,d,"statistics.h5") for d in readdir(directory) if isfile(joinpath(directory,d,"statistics.h5"))]
   @test length(statistics)==1
   h5open(only(statistics),"r") do stats
    attributes(vg)["statistics_source"]=relpath(only(statistics),root)
    attributes(vg)["statistics_sha256"]=sha(only(statistics))
    indices=findall(read(stats["criterion_noise_multiple"]).==1.)
    for sample in readout["samples"]
     sn=sample["name"];source=stats[sn];a=attributes(source)
     @test read(a["full_candidates_complete"])==1
     @test read(a["completed"])==read(a["planned"])
     g=create_group(vg,"histories/"*sn)
     for (k,value) in attrs(source);attributes(g)[k]=value;end
     g["loss_fraction"]=read(stats["criterion_loss_fraction"])[indices]
     g["event_count"]=[read(attributes(source["criterion_$j"])["event_count"]) for j in indices]
     g["event_fraction"]=[read(attributes(source["criterion_$j"])["event_fraction"]) for j in indices]
     screen=joinpath(dir,"screens",name,sn*".h5")
     h5open(screen,"r") do h
      for k in ("pads","minimum_linear_output_V","maximum_bucket_ions","pads_below_guard","events_with_guard_risk")
       attributes(g)[k]=read(attributes(h)[k])
      end
     end
    end
   end
  end
 end
 end
 end
end

end
