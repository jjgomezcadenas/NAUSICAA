using HDF5,TOML,SHA
const START=time()
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis
include("fast_trigger_emin_functions.jl")
card=TOML.parsefile(joinpath(@__DIR__,"fast_trigger_emin.toml"))
readcard(key)=TOML.parsefile(joinpath(ROOT,card[key]))
r=readcard("realisation_card");d=readcard("detector_card");sensor=readcard("sipm_card");asic=readcard("electronics_card");o=d["optics"];gas=d[card["gas"]]
Y=Float64(el_photon_yield(Float32(gas["yield_gain_photons_per_electron_kV"]),Float32(gas["yield_threshold_kV_cm_bar"]),Float32(d["reduced_field_kV_cm_bar"]),Float32(d["pressure_bar"]),Float32(d["hole_depth_mm"])))
a=Float64(elm_light_acceptance(Float32(d["hole_diameter_mm"]/2),Float32(d["hole_depth_mm"]);tpb_efficiency=Float32(o["tpb_conversion_efficiency"]),forward_fraction=Float32(o["forward_fraction"]),guide_index=Float32(o["guide_refractive_index"]),gap_index=Float32(o["gap_refractive_index"]),propagation_transmission=Float32(o["propagation_transmission"]),coupling_transmission=Float32(o["coupling_transmission"]),electrode_transmission=Float32(o["electrode_transmission"])))
@assert d["electron_collection_efficiency"]==1
@assert !sensor["crosstalk"] && !sensor["afterpulsing"] && !sensor["gain_spread"]
@assert isapprox(asic["bucket_width_s"],d["bucket_width_us"]*1e-6)
rows=NamedTuple[]
for E in card["energies_keV"], L in card["electron_drift_mm"], K in card["summed_sipm_buckets"]
    survival=exp(-L/d["drift_speed_mm_us"]/(1000d["lifetime_ms"]))
    v=emin_moments(;energy_keV=E,w_eV=r["w_eV"],fano=r["fano_factor"],survival,yield_per_e=Y,J=d["j_variance_ratio"],acceptance=a,pde=sensor["pde"],cells=K,dark_per_cell=sensor["dark_reference_rate_Hz"]*asic["bucket_width_s"],noise_pe=asic["noise_rms_pe"])
    clip=asic["current_limit_A"]*asic["bucket_width_s"]/(1.602176634e-19*sensor["gain"])
    # Simultaneous arrivals in one bucket are a conservative occupancy check.
    M=sensor["illuminated_microcell_count"]
    occupied=M*(-expm1(-v.pe/M));loss=1-occupied/v.pe
    push!(rows,merge((energy_keV=E,drift_mm=L,summed_cells=K,survival,clip_pe=clip,instantaneous_occupancy_loss_fraction=loss),v))
end
path=joinpath(@__DIR__,"results/fast_trigger_emin.h5")
h5open(path,"w") do f
    f["schema"]="itaca.fast_trigger_emin.v1";f["card"]=read(joinpath(@__DIR__,"fast_trigger_emin.toml"),String)
    for key in ("realisation_card","detector_card","sipm_card","electronics_card")
        f[key]=read(joinpath(ROOT,card[key]),String)
    end
    f["yield_per_electron"]=Y;f["optical_acceptance"]=a;f["pde"]=sensor["pde"];f["asic_noise_rms_pe"]=asic["noise_rms_pe"]
    for k in keys(first(rows));f[string(k)]=[getproperty(v,k) for v in rows];end
    f["implementation_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,"fast_trigger_emin_functions.jl"))))
    f["elapsed_seconds"]=time()-START
end
println("Saved ",path," in ",round(time()-START;digits=2)," s")
