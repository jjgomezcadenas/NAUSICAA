using HDF5,TOML,SHA
root=normpath(joinpath(@__DIR__,".."))
python=ARGS[1];simulator=ARGS[2]
dir=joinpath(root,"results/ltspice/rogelio_v2_retuned");mkpath(dir)
target=0.878931543
gates=Float64[];outputs=Float64[]
function point(gate)
 run(`$python $(joinpath(root,"scripts/run_v2_gate_point.py")) --simulator $simulator --gate $gate`)
 value=h5open(joinpath(dir,"tuning_raw.h5"),"r") do h
  names=collect(keys(h));j=only(filter(n->isapprox(read(attributes(h[n])["gate_V"]),gate;atol=1e-14,rtol=0),names))
  only(read(h[j]["V(out)"]))
 end
 push!(gates,gate);push!(outputs,value);println("gate=",gate," output=",value)
 value
end
lo=.64;hi=.8
@assert point(lo)>target>point(hi)
for iteration in 1:32
 global lo,hi
 mid=(lo+hi)/2;value=point(mid)
 if value>target;lo=mid;else;hi=mid;end
 abs(value-target)<1e-6 && break
end
best=argmin(abs.(outputs.-target));gate=gates[best]
h5open(joinpath(dir,"tuning.h5"),"w") do h
 h["gate_V"]=gates;h["output_V"]=outputs
 attributes(h)["target_output_V"]=target;attributes(h)["selected_gate_V"]=gate
 attributes(h)["selected_output_V"]=outputs[best]
 attributes(h)["source_sha256"]=bytes2hex(sha256(read(joinpath(dir,"tuning_raw.h5"))))
end
c=TOML.parsefile(joinpath(root,"config/rogelio_v2.toml"))
c["output"]="results/ltspice/rogelio_v2_retuned";c["retuned"]=true
c["variants"]=[Dict("name"=>"v2_retuned","capacitance_F"=>1.25e-12,"temperature_C"=>20.,"m4_gate_V"=>gate)]
path=joinpath(root,"config/rogelio_v2_retuned.toml")
buffer=IOBuffer();TOML.print(buffer,c);text=String(take!(buffer))
if isfile(path);@assert read(path,String)==text;else;write(path,text);end
println("Selected gate=",gate," Vout=",outputs[best])
