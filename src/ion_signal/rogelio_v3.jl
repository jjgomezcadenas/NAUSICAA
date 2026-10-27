module RogelioV3
using HDF5,SHA,LinearAlgebra
include("discharge_calibration.jl")
const D=DischargeCalibration
const B=D.B
function analyze(root)
 dir=joinpath(root,"results/ltspice/rogelio_v3")
 h5open(joinpath(dir,"raw.h5"),"r") do raw
  h5open(joinpath(dir,"analysis.h5"),"w") do out
   a=attributes(out);a["raw_sha256"]=D.sha(joinpath(dir,"raw.h5"));a["analysis_sha256"]=D.sha(@__FILE__)
   a["baseline_V"]=only(read(raw["op/V(out)"]))
   a["M4_display_Vout_over_Id_Ohm"]=a["baseline_V"][]/only(read(raw["op/Id(M4)"]))
   names=String[];logtext=read(attributes(raw["op"])["log"])
   for line in split(logtext,'\n')
    startswith(line,"Name:") && (names=split(line)[2:end])
    if startswith(line,"Gds:") && "m4" in names
     gd=parse(Float64,split(line)[findfirst(==("m4"),names)+1]);a["M4_gds_S"]=gd;a["M4_incremental_ro_Ohm"]=1/gd
    end
   end
   for name in ("original","original_uncompressed","exp_r1","exp_r2","uniform_1000_r1","uniform_1000_r2")
    h=raw[name];t=read(h["time"]);i=read(h["I(Is1)"]);v=read(h["V(out)"])
    g=create_group(out,name);at=attributes(g)
    charge=B.integral(t,i);voltage=first(v).-v
    at["charge_C"]=charge;at["equivalent_electrons"]=charge/B.E
    at["resolved_charge_integration"]=!startswith(name,"original")
    at["peak_voltage_V"]=maximum(voltage);at["peak_voltage_time_s"]=t[argmax(voltage)]
    at["current_max_A"]=maximum(i);at["first_positive_current_s"]=t[findfirst(>(0),i)]
    qkey=only(filter(k->lowercase(k)=="v(qmip_c)",collect(keys(h))))
    at["behavioral_charge_C"]=last(read(h[qkey]))
    g["time_s"]=t;g["current_A"]=i;g["voltage_drop_V"]=voltage
    if startswith(name,"exp")
     inds=findall(j->.20001<t[j]<.22 && i[j]>0,eachindex(t))
     x=t[inds].-.2;coef=hcat(ones(length(x)),x)\log.(i[inds])
     at["fitted_current_tau_s"]=-1/coef[2];at["fitted_current_at_0p2s_A"]=exp(coef[1])
    end
    println(name," charge fC=",charge*1e15," electrons=",charge/B.E," peak mV=",maximum(voltage)*1e3," t_peak=",t[argmax(voltage)]," idt fC=",read(at["behavioral_charge_C"])*1e15)
   end
   for prefix in ("exp","uniform_1000")
    coarse=attributes(out[prefix*"_r1"]);fine=attributes(out[prefix*"_r2"])
    fine["charge_refinement_error"]=abs(read(coarse["charge_C"])/read(fine["charge_C"])-1)
    fine["voltage_refinement_error"]=abs(read(coarse["peak_voltage_V"])/read(fine["peak_voltage_V"])-1)
    @assert read(fine["charge_refinement_error"])<1e-5
    @assert read(fine["voltage_refinement_error"])<1e-4
    @assert abs(read(fine["behavioral_charge_C"])/read(fine["charge_C"])-1)<1e-5
   end
   expa=attributes(out["exp_r2"])
   a["nominal_exp_area_C"]=48e-15*1.7e-3
   a["exp_area_relative_error"]=read(expa["charge_C"])/read(a["nominal_exp_area_C"])-1
   @assert abs(read(a["exp_area_relative_error"]))<1e-5
   a["uniform_to_exp_charge_ratio"]=read(attributes(out["uniform_1000_r2"])["charge_C"])/read(expa["charge_C"])
   a["uniform_to_exp_voltage_ratio"]=read(attributes(out["uniform_1000_r2"])["peak_voltage_V"])/read(expa["peak_voltage_V"])
  end
 end
end
end
