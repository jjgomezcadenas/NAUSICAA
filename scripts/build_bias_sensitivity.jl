using TOML,HDF5,SHA
root=normpath(joinpath(@__DIR__,".."));cp=joinpath(root,"config/bias_sensitivity.toml");c=TOML.parsefile(cp)
out=joinpath(root,c["output"]);mkpath(out)
points=NamedTuple[]
for family in ("paper","v2")
 p=c[family]
 for device in (family=="paper" ? ("m4","m3") : ("m4",))
  for (j,offset) in enumerate(c["offsets_V"])
   push!(points,(name="$(family)_$(device)_$(j)",family=family,kind=device,offset_V=offset,
     m3_gate_V=p["m3_gate_V"]+(device=="m3" ? offset : 0.),m4_gate_V=p["m4_gate_V"]+(device=="m4" ? offset : 0.)))
  end
 end
 for digits in c["rounding_digits"]
  push!(points,(name="$(family)_rounded_$(digits)",family=family,kind="rounded",offset_V=NaN,
     m3_gate_V=round(p["m3_gate_V"];digits=digits),m4_gate_V=round(p["m4_gate_V"];digits=digits)))
 end
end
h5open(joinpath(out,"points.h5"),"w") do h
 attributes(h)["config_sha256"]=bytes2hex(sha256(read(cp)))
 attributes(h)["generator_sha256"]=bytes2hex(sha256(read(@__FILE__)))
 for point in points
  g=create_group(h,point.name)
  for k in keys(point);attributes(g)[string(k)]=getproperty(point,k);end
 end
end
println(length(points)," bias points")
