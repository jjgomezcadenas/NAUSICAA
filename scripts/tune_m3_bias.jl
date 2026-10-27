using HDF5
root=normpath(joinpath(@__DIR__,".."))
python,simulator=ARGS
raw=joinpath(root,"results/ltspice/m3_scan/raw.h5")
getoutput(g)=only(read(g["op/V(out)"]))
target=h5open(h->getoutput(h["scale_1"]),raw,"r")
function acquire(scale,gate;full=false)
    args=[python,joinpath(root,"scripts/run_m3_scan.py"),"--simulator",simulator,"--scales",string(scale),"--gate",string(gate)]
    full || push!(args,"--op-only")
    run(Cmd(Cmd(args);dir=root))
    tag="scale_$(scale)_gate_"*string(round(gate;digits=10))
    h5open(raw,"r") do h
        matches=filter(k->startswith(k,"scale_$(scale)_gate_") && abs(read(attributes(h[k])["gate_V"])-gate)<1e-10,keys(h))
        getoutput(h[only(matches)])
    end
end
for scale in (2,4,8)
    lo=.35;hi=.9
    vlo=acquire(scale,lo);vhi=acquire(scale,hi)
    @assert (vlo-target)*(vhi-target)<0
    gate=NaN;value=NaN
    for iteration in 1:32
        gate=(lo+hi)/2;value=acquire(scale,gate)
        abs(value-target)<1e-5 && break
        if (value-target)*(vlo-target)>0
            lo=gate;vlo=value
        else
            hi=gate;vhi=value
        end
    end
    @assert abs(value-target)<1e-5
    println("Tuned scale $scale: gate=$gate, output=$value, target=$target")
    acquire(scale,gate;full=true)
    run(Cmd(Cmd([python,joinpath(root,"scripts/run_m3_scan.py"),"--simulator",simulator,"--scales",string(scale),"--gate",string(gate),"--injection-only","--injection-counts","0","10","100","--refinements","1","2"]);dir=root))
end
