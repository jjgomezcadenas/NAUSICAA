module PadSequenceExample
using HDF5,TOML,Test
include("pad_history_readout.jl")
const P=PadHistoryReadout
const B=P.B
config(root)=TOML.parsefile(joinpath(root,"config/pad_sequence_example.toml"))

"""Save charge-preserving finite-gap currents for the example and paired controls."""
function build(root)
 c=config(root);dir=joinpath(root,c["output"]);mkpath(dir)
 destination=joinpath(dir,"inputs.h5");isfile(destination) && error("Preserving existing inputs")
 n=c["counts"];frame=c["frame_s"];pre=c["prehistory_s"]
 h5open(destination,"w") do h
  attributes(h)["config_sha256"]=P.sha(joinpath(root,"config/pad_sequence_example.toml"))
  attributes(h)["generator_sha256"]=P.sha(@__FILE__)
  @testset "Example input charge" begin
   for name in ("requested","without_probe","without_large","without_large_or_probe")
    counts=copy(n)
    if name in ("without_large","without_large_or_probe");counts[8]=0;end
    if name in ("without_probe","without_large_or_probe");counts[9]=0;end
    g=create_group(h,name);g["counts"]=counts
    g["edges_s"]=collect(0:length(n)).*frame
    for (r,step) in enumerate(c["max_steps_s"])
     u=collect(range(-.005/3,0.,length=ceil(Int,.005/3/step)+1))
     unit_t=vcat(u,frame.+u);unit_i=B.current.(unit_t,frame,.00025,.005,3.)
     active=findall(>(0),counts)
     times=sort!(unique(round.(vcat(0.,[pre+x+(k-1)*frame for k in active for x in unit_t],pre+length(n)*frame+c["recovery_tail_s"]);digits=12)))
     current=[sum(counts[k]*((first(unit_t)<=t-pre-(k-1)*frame<=last(unit_t)) ? B.interp(unit_t,unit_i,t-pre-(k-1)*frame) : 0.) for k in active) for t in times]
     expected=sum(counts)*B.E*(1-B.phi(.005,.00025))
     error=abs(B.integral(times,current)/expected-1)
     @test error<1e-8
     group=create_group(g,"r$r");group["time_s"]=times;group["current_A"]=current
     attributes(group)["charge_error"]=error
    end
   end
  end
 end
end

"""Average SPICE boundary readings and compare recovered counts and probe gain."""
function analyze(root)
 c=config(root);dir=joinpath(root,c["output"]);frame=c["frame_s"];pre=c["prehistory_s"]
 reference=P.settings(root);k=P.kernel(root,reference,frame);width=c["averaging_s"]
 gain=k.gain
 noise=h5open(joinpath(root,c["circuit_analysis"]),"r") do h
  group=h[c["circuit_name"]]["$(round(Int,frame*1000))ms"]
  j=read(attributes(group)["best_index"]);read(group["noise_rms_V"])[j]
 end
 h5open(joinpath(dir,"analysis.h5"),"w") do out
  for (name,value) in ("config_sha256"=>P.sha(joinpath(root,"config/pad_sequence_example.toml")),"analysis_sha256"=>P.sha(@__FILE__),"gain_V_per_ion"=>gain,"noise_rms_V"=>noise,"ENC_e"=>noise/gain,"tau_s"=>k.tau,"alpha"=>k.alpha,"baseline_V"=>k.baseline,"frame_s"=>frame,"averaging_s"=>width,"assumptions"=>"Uniform intra-bucket arrivals; on-axis 0.5 mm square pad in grounded plane; 5 mm gap; 3 m/s speed; finite-path charge retained; each run starts at quiescent baseline. No random noise added.")
   attributes(out)[name]=value
  end
  h5open(joinpath(dir,"inputs.h5"),"r") do inp
   @testset "Example transient validation" begin
    for name in keys(inp)
     source=joinpath(dir,name*".h5")
     h5open(source,"r") do raw
      counts=read(inp[name]["counts"]);edges=read(inp[name]["edges_s"]);times=edges.+pre
      function reduce(r)
       h=raw["history_r$r"];b=raw["baseline_r$r"]
       t=read(h["time"]);v=read(h["V(out)"]);tb=read(b["time"]);vb=read(b["V(out)"])
       delta=v.-[B.interp(tb,vb,x) for x in t]
       cumulative=vcat(0.,cumsum(diff(t).*(delta[1:end-1].+delta[2:end])./2))
       primitive(x)=begin
        j=clamp(searchsortedlast(t,x),1,length(t)-1);dx=x-t[j]
        cumulative[j]+delta[j]*dx+(delta[j+1]-delta[j])/(t[j+1]-t[j])*dx^2/2
       end
       avg=[(primitive(x+width/2)-primitive(x-width/2))/width for x in times]
       point=[B.interp(t,v,x) for x in times]
       signal=-(avg[2:end].-k.alpha.*avg[1:end-1])
       (;t,v,avg,point,signal)
      end
      fine=reduce(2);coarse=reduce(1)
      error=maximum(abs.(fine.signal.-coarse.signal))/maximum(abs.(fine.signal))
      charge=B.integral(fine.t,read(raw["history_r2/I(Is1)"]))
      expected=sum(counts)*B.E*(1-B.phi(.005,.00025))
      @test abs(charge/expected-1)<1e-8
      @test error<1e-4
      linear_v,linear_avg,linear_signal=P.response(counts,k)
      g=create_group(out,name)
      for (key,value) in ("counts"=>counts,"edges_s"=>edges,"output_at_edges_V"=>fine.point,"average_at_edges_V"=>fine.avg.+k.baseline,"frame_signal_V"=>fine.signal,"estimated_ions"=>fine.signal./gain,"linear_output_at_edges_V"=>linear_v[1:length(edges)],"linear_estimated_ions"=>linear_signal[1:length(counts)]./gain,"time_s"=>fine.t.-pre,"output_V"=>fine.v)
       g[key]=value
      end
      attributes(g)["raw_sha256"]=P.sha(source)
      attributes(g)["minimum_output_V"]=minimum(fine.v)
      attributes(g)["refinement_error"]=error
      attributes(g)["charge_error"]=abs(charge/expected-1)
      println(name," minV=",minimum(fine.v)," refinement=",error)
      for j in eachindex(counts)
       println(j," ",counts[j]," ",fine.point[j+1]," linear=",linear_signal[j]/gain," estimated=",fine.signal[j]/gain)
      end
     end
    end
   end
  end
  with_large=(read(out["requested/frame_signal_V"])[9]-read(out["without_probe/frame_signal_V"])[9])/gain
  without_large=(read(out["without_large/frame_signal_V"])[9]-read(out["without_large_or_probe/frame_signal_V"])[9])/gain
  attributes(out)["probe_estimated_with_large_e"]=with_large
  attributes(out)["probe_estimated_without_large_e"]=without_large
  attributes(out)["probe_relative_gain_change"]=with_large/without_large-1
  println("Marginal 200-ion probe: after5000=",with_large," without5000=",without_large," change=",with_large/without_large-1)
 end
end
end
