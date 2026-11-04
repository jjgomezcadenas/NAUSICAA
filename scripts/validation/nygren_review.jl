# Analytical charge-noise and drift-time budgets for the local SeF6 paper review.
using HDF5, TOML, SHA
root=normpath(joinpath(@__DIR__,"../.."))
card=joinpath(root,"config/nygren_review.toml"); c=TOML.parsefile(card)
N=c["Q_eV"]/c["W_eV"]; F=c["fano"]; g=2sqrt(2log(2))
enc=c["enc_e"]; channels=c["independent_channels"]
resolution=[g*sqrt(F*N+n*e^2)/N for n in channels,e in enc]
budgets=[sqrt((r*N/g)^2-F*N) for r in c["target_fwhm_fraction"]]
channel_limits=[b^2/e^2 for b in budgets,e in enc]
durations=[z/v for z in c["z_extent_m"],v in c["drift_speed_m_s"]]
mobility=c["reduced_mobility_m2_V_s"].*(c["reference_pressure_Pa"]/c["pressure_Pa"])*
    (c["temperature_K"]/c["reference_temperature_K"])
speeds=mobility.*c["field_V_m"]
electron_charge=1.602176634e-19
shotnoise=sqrt.(c["leakage_A"].*c["integration_s"]./electron_charge)
@assert all(isfinite,resolution) && all(resolution.>g*sqrt(F/N))
@assert durations[2,:]≈2durations[1,:]
@assert isapprox(resolution[1,1],.00506;atol=2e-5) # reproduces CDR's rounded 0.51%
outdir=joinpath(root,"results/nygren_review");mkpath(outdir)
write(joinpath(outdir,".gitignore"),"*.h5\n")
h5open(joinpath(outdir,"budgets.h5"),"w") do h
    attributes(h)["schema"]="nausicaa.nygren_review.v1"
    for (k,v) in c
        k=="sources" && continue
        h[k]=v
    end
    h["primary_pairs"]=N;h["intrinsic_fwhm_fraction"]=g*sqrt(F/N)
    h["resolution_fwhm_fraction"]=resolution
    h["target_noise_budget_e"]=budgets;h["maximum_independent_channels"]=channel_limits
    h["arrival_duration_s"]=durations;h["mobility_speed_m_s"]=speeds
    h["shot_noise_e"]=shotnoise
    h["electronic_1000e_fwhm_fraction"]=g*sqrt(F*N+1000^2)/N
    h["card_toml"]=read(card,String)
    h["calculation_code"]=read(@__FILE__,String)
    h["calculation_sha256"]=bytes2hex(sha256(read(@__FILE__)))
    h["julia_version"]=string(VERSION)
    for (j,p) in enumerate(c["sources"])
        d=create_group(h,"sources/$j");d["path"]=p
        d["sha256"]=bytes2hex(sha256(read(joinpath(root,p))))
    end
end
println("Primary pairs: ",N,"; intrinsic FWHM: ",100g*sqrt(F/N),"%")
for (i,n) in enumerate(channels);println(n," channels: ",100resolution[i,:]," % FWHM");end
println("<1% FWHM channel limits for ENC ",enc,": ",channel_limits[3,:])
println("Mobility-derived speed range: ",speeds," m/s")
