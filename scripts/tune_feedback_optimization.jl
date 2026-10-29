using HDF5,TOML,SHA
root=normpath(joinpath(@__DIR__,".."));python,sim=ARGS
cp=joinpath(root,"config/feedback_optimization.toml");c=TOML.parsefile(cp);dir=joinpath(root,c["output"])
target=h5open(h->read(attributes(h)["baseline_V"]),joinpath(dir,"reference.h5"),"r")
function acquire(name,gate;full=false)
 args=[python,joinpath(root,"scripts/run_feedback_retuned.py"),"--simulator",sim,"--variant",name,"--m4-gate",string(gate)]
 full || push!(args,"--op-only")
 result=read(Cmd(Cmd(args);dir=root),String)
 full && (print(result);return target)
 path=joinpath(root,strip(result))
 h5open(h->read(attributes(h)["baseline_V"]),path,"r")
end
pointsfile=joinpath(dir,"retuned_points.h5")
for v in c["variants"]
 name=v["name"];name=="reference" && continue
 done=isfile(pointsfile) && h5open(h->haskey(h,name),pointsfile,"r")
 if done
  gate=h5open(h->read(attributes(h[name])["m4_gate_V"]),pointsfile,"r")
 else
  lo=.45;hi=.7;vlo=acquire(name,lo);vhi=acquire(name,hi)
  @assert (vlo-target)*(vhi-target)<0 (name,vlo,vhi)
  gate=NaN;value=NaN
  for iteration in 1:32
   gate=(lo+hi)/2;value=acquire(name,gate)
   abs(value-target)<1e-4 && break
   if (value-target)*(vlo-target)>0;lo=gate;vlo=value;else;hi=gate;vhi=value;end
  end
  @assert abs(value-target)<1e-4 (name,gate,value,target)
  h5open(pointsfile,"cw") do h
   g=create_group(h,name);a=attributes(g)
   for (k,val) in v;a[k]=val;end
   a["m4_gate_V"]=gate;a["baseline_V"]=value;a["target_baseline_V"]=target
   a["tuner_sha256"]=bytes2hex(sha256(read(@__FILE__)))
  end
 end
 println("Retuned ",name," M4gate=",gate);flush(stdout)
 acquire(name,gate;full=true);flush(stdout)
end
