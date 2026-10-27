using HDF5,SHA
include(joinpath(@__DIR__,"../src/ion_signal/m4_scan.jl"))
root=normpath(joinpath(@__DIR__,".."));dir=joinpath(root,"results/ltspice/m3_scan")
h5open(joinpath(dir,"raw.h5"),"r") do h
 h5open(joinpath(dir,"injection_diagnostics.h5"),"w") do o
  attributes(o)["raw_sha256"]=bytes2hex(sha256(read(joinpath(dir,"raw.h5"))))
  attributes(o)["script_sha256"]=bytes2hex(sha256(read(@__FILE__)))
  attributes(o)["purpose"]="One-electron injection numerical diagnostics. Final calibration uses 100-electron injection, checked at 10 electrons; injection_precise tolerance trials are excluded from final ENC."
  for tag in sort(filter(k->startswith(k,"scale_") && haskey(h[k],"noise_400"),keys(h)))
    g=h[tag];out=create_group(o,tag)
    for kind in ("injection","injection_precise")
      for r in (1,2,4,8)
        haskey(g,"$(kind)_N1_r$r") || continue
        t,v=M4Scan.wave(g,kind,1,r)
        b=g["$(kind)_N0_r$r"];tb=read(b["time"]);vb=read(b["V(out)"])
        gain,peak=M4Scan.injection_gain(t,v,300e-6)
        a=attributes(create_group(out,"$(kind)_r$r"))
        for (k,x) in ("gain_V_per_e"=>gain,"peak_time_s"=>peak,"baseline_range_after_4ms_V"=>maximum(vb[tb.>.004])-minimum(vb[tb.>.004]),"signal_min_V"=>minimum(v),"signal_max_V"=>maximum(v))
          a[k]=x
        end
        println(tag," ",kind," r=",r," gain=",gain," baseline range=",read(a["baseline_range_after_4ms_V"])," signal min/max=",extrema(v))
      end
    end
  end
 end
end
