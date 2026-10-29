#!/usr/bin/env julia
# Evaluate Petit's electron width prescription at the configured pressure.
# sigma_T,L [mm] = d_T,L [mm sqrt(bar/cm)] sqrt(L [cm]/P [bar]).
using TOML,HDF5,SHA
root=normpath(joinpath(@__DIR__,"../.."));card=joinpath(root,"config/electron_transport.toml")
c=TOML.parsefile(card);p=Float32(c["gas"]["pressure_bar"])
lengths=Float32.(c["evaluation"]["drift_cm"])
output=joinpath(root,"results/electron_diffusion");mkpath(output)
h5open(joinpath(output,"electron_diffusion.h5"),"w") do f
    f["pressure_bar"]=p;f["helium_fraction"]=Float32(c["gas"]["helium_fraction"])
    f["drift_cm"]=lengths
    f["provenance/card_toml"]=read(card,String)
    f["provenance/code_sha256"]=bytes2hex(sha256(read(@__FILE__)))
    f["provenance/prescription"]="Petit src/itaca_functions.jl:sigma_t_mm/sigma_l_mm; Xe-He coefficients identified with 15% helium by user design input"
    for gas in ("pure_xe","xe_he"),axis in ("T","L")
        base=c["diffusion"][gas]["d_"*axis*"_mm_sqrtbar_per_sqrtcm"]
        effective=Float32(base)/sqrt(p)
        sigma=effective.*sqrt.(lengths)
        reference=base.*sqrt.(Float64.(lengths)./Float64(p))
        f[gas*"/d_"*axis*"_mm_per_sqrtcm"]=effective
        f[gas*"/sigma_"*axis*"_mm"]=sigma
        f[gas*"/max_Float64_error_"*axis*"_mm"]=maximum(abs.(Float64.(sigma).-reference))
        @assert maximum(abs.(Float64.(sigma).-reference))<2e-6
        println(gas," d",axis,"(15 bar)=",effective," mm/sqrt(cm); sigma=",sigma," mm")
    end
end
