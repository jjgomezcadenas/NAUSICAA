#!/usr/bin/env julia
#=
Compute the analytical worked example in docs/elm.tex and save HDF5 results.

Run from the repository with `julia --project=. scripts/documentation/calculate_elm_example.jl`.
Detector and response inputs are read from the adopted ELM cards. Existing event-analysis functions
evaluate the EL yields and axial hole acceptance; the guide uses the agreed
meridional, isotropic-in-core approximation. Probabilities and mean photon counts
are deterministic. No event production or SiPM response is run.
=#
using HDF5, TOML, SHA
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis

detector_path = joinpath(ROOT,"config/specs/detector/electrons/elm_4mm_250ns_ep2_drift10cm.toml")
sensor_path = joinpath(ROOT,"config/specs/response/sipm/elm_pitch4mm_3015ps_tau12p5ns.toml")
electronics_path = joinpath(ROOT,"config/specs/response/sipm_electronics/readout_s14160_1315ps.toml")
detector, sensor, electronics = TOML.parsefile.((detector_path,sensor_path,electronics_path))
optics = detector["optics"]
p = Dict("pressure_bar"=>detector["pressure_bar"], "ep_kV_cm_bar"=>detector["reduced_field_kV_cm_bar"],
    "hole_depth_mm"=>detector["hole_depth_mm"], "hole_diameter_mm"=>detector["hole_diameter_mm"], "pitch_mm"=>detector["hole_pitch_mm"],
    "gas_gap_radial_mm"=>0.5, "guide_length_cm"=>10.0,
    "guide_index"=>1.50, "gas_index"=>1.00, "tpb_efficiency"=>0.65,
    "propagation_transmission"=>1.0, "coupling_transmission"=>1.0,
    "electrode_transmission"=>1.0, "electron_collection"=>1.0,
    "sensor_coverage"=>1.0, "cathode_kV"=>1.2, "drift_V_cm"=>200.0,
    "drift_length_cm"=>150.0, "transfer_V_cm"=>300.0, "transfer_gap_cm"=>1.0,
    "sensor_plane_kV"=>0.0, "sipm_side_mm"=>detector["sipm_side_mm"], "bucket_us"=>detector["bucket_width_us"],
    "window_index"=>1.57, "adhesive_index"=>1.52)
for (name, key, scale) in [("gas_gap_radial_mm","gas_gap_radial_mm",1.0),
    ("guide_length_cm","guide_length_mm",0.1), ("guide_index","guide_refractive_index",1.0),
    ("gas_index","gap_refractive_index",1.0), ("tpb_efficiency","tpb_conversion_efficiency",1.0),
    ("propagation_transmission","propagation_transmission",1.0),
    ("coupling_transmission","coupling_transmission",1.0),
    ("electrode_transmission","electrode_transmission",1.0),
    ("cathode_kV","cathode_potential_V",0.001),("sensor_plane_kV","sensor_plane_potential_V",0.001)]
    p[name]=optics[key]*scale
end
p["electron_collection"]=detector["electron_collection_efficiency"]
fitpath = joinpath(ROOT, "config/electron_transport.toml")
fit = TOML.parsefile(fitpath)["el_yield"]
r, h = p["hole_diameter_mm"]/2, p["hole_depth_mm"]
A = Float64(axial_light_acceptance(Float32(r), Float32(h)))
capture = 1-p["gas_index"]/p["guide_index"]
forward = optics["forward_fraction"]
optical = A*p["tpb_efficiency"]*forward*capture*p["propagation_transmission"]*
    p["coupling_transmission"]*p["electrode_transmission"]*p["sensor_coverage"]
vd = p["cathode_kV"]+p["drift_V_cm"]*p["drift_length_cm"]/1000
vin = vd+p["transfer_V_cm"]*p["transfer_gap_cm"]/1000
dv = p["pressure_bar"]*p["ep_kV_cm_bar"]*h/10
vout = vin+dv
grading = (vout-p["sensor_plane_kV"])/p["guide_length_cm"]

# Independent midpoint integration of solid angle, and the angular capture cone.
n = 100_000
aq = sum(0.5*(1-z/hypot(z,r)) for z in ((i-0.5)*h/n for i in 1:n))/n
theta = acos(p["gas_index"]/p["guide_index"])
cq = sum(sin((i-0.5)*theta/n) for i in 1:n)*theta/n
@assert abs(A-aq)<1e-7
@assert abs(capture-cq)<1e-9
@assert 0<optical<A<0.5
@assert p["hole_diameter_mm"]<=p["sipm_side_mm"]
@assert isapprox(vin*1000,optics["entry_potential_V"])
@assert isapprox(vout*1000,optics["exit_potential_V"])
@assert sensor["illuminated_microcell_count"]==round(Int,sensor["microcell_count"]*pi*r^2/p["sipm_side_mm"]^2)

out = joinpath(ROOT,"results/elm/elm_example.h5")
mkpath(dirname(out))
h5open(out,"w") do f
    attributes(f)["schema"] = "elm_analytical_example_v2"
    attributes(f)["model"] = "uniform axial EL; absorbing walls; isotropic-in-core meridional capture; ideal propagation/coupling"
    attributes(f)["calculation_sha256"] = bytes2hex(sha256(read(@__FILE__)))
    attributes(f)["yield_card_sha256"] = bytes2hex(sha256(read(fitpath)))
    attributes(f)["grid_readout_sha256"] = bytes2hex(sha256(read(joinpath(ROOT,"src/event_analysis/grid_readout.jl"))))
    provenance=create_group(f,"provenance")
    for (key,path) in [("detector",detector_path),("sensor",sensor_path),("electronics",electronics_path)]
        provenance[key*"_path"]=relpath(path,ROOT)
        provenance[key*"_toml"]=read(path,String)
        provenance[key*"_sha256"]=bytes2hex(sha256(read(path)))
    end
    inputs=create_group(f,"inputs")
    for (k,v) in p; inputs[k]=v; end
    vals=Dict("Acceptance"=>A,"Capture"=>capture,"OpticalFraction"=>optical,
        "OpticalPercent"=>100optical,"ForwardFactor"=>forward,
        "TrappingAngleDeg"=>rad2deg(theta),"DriftEndKV"=>vd,"EntryKV"=>vin,
        "ELDropKV"=>dv,"ExitKV"=>vout,"ELFieldKVcm"=>p["pressure_bar"]*p["ep_kV_cm_bar"],
        "GradingKVcm"=>grading,"ReducedGrading"=>grading/p["pressure_bar"],
        "FieldRatio"=>p["pressure_bar"]*p["ep_kV_cm_bar"]*1000/p["transfer_V_cm"],
        "SquareAperturePercent"=>100pi*r^2/p["pitch_mm"]^2,
        "HexAperturePercent"=>100pi*r^2/(sqrt(3)/2*p["pitch_mm"]^2),
        "ChannelDiameterMM"=>2r+2p["gas_gap_radial_mm"],
        "WebMM"=>p["pitch_mm"]-2r-2p["gas_gap_radial_mm"],
        "VUVExitPerThousand"=>1000A,
        "BlueTotalPerThousand"=>1000A*p["tpb_efficiency"],
        "BlueForwardPerThousand"=>1000A*p["tpb_efficiency"]*forward,
        "BlueSiPMPerThousand"=>1000optical,
        "OptimumGlueIndex"=>sqrt(p["guide_index"]*p["window_index"]),
        "PDE"=>sensor["pde"],"Gain"=>sensor["gain"],"Microcells"=>sensor["microcell_count"],
        "IlluminatedMicrocells"=>sensor["illuminated_microcell_count"],
        "DarkOnlyMicrocells"=>sensor["microcell_count"]-sensor["illuminated_microcell_count"],
        "DarkRateKCPS"=>sensor["dark_reference_rate_Hz"]/1000,
        "RecoveryNS"=>sensor["recovery_s"]*1e9,"BucketNS"=>p["bucket_us"]*1000,
        "NoisePE"=>electronics["noise_rms_pe"],"CurrentLimitMA"=>electronics["current_limit_A"]*1000,
        "ClipThresholdPE"=>electronics["current_limit_A"]*electronics["bucket_width_s"]/(1.602176634e-19*sensor["gain"]))
    for (label,gas) in [("Xe","pure_xe"),("He","xe_he")]
        c=fit[gas]
        y=Float64(el_photon_yield(Float32(c["gain_photons_per_electron_kV"]),
            Float32(c["threshold_kV_cm_bar"]),Float32(p["ep_kV_cm_bar"]),
            Float32(p["pressure_bar"]),Float32(h)))
        vals[label*"Yield"]=y
        vals[label*"ExitPhotons"]=y*A
        vals[label*"SiPMPhotons"]=y*optical
        vals[label*"PrimaryPE"]=y*optical*sensor["pde"]
    end
    for (label,ep) in [("Two",2.0),("TwoHalf",2.5),("Three",3.0)]
        drop=p["pressure_bar"]*ep*h/10
        vals[label*"ELDrop"]=drop
        vals[label*"Exit"]=vin+drop
        vals[label*"Length"]=(vin+drop-p["sensor_plane_kV"])/grading
    end
    values=create_group(f,"values")
    for (k,v) in vals; values[k]=v; end
    validation=create_group(f,"validation")
    validation["acceptance_quadrature_absolute_error"]=abs(A-aq)
    validation["capture_quadrature_absolute_error"]=abs(capture-cq)
    validation["quadrature_steps"]=n
    # Copy saved Julia-reduced pilot scalars; presentation performs no statistics.
    pilotpath=joinpath(ROOT,"results/elm_sipm_pitch4mm_pilot/elm_sipm_pitch4mm_pilot.h5")
    oldpath=joinpath(ROOT,"results/elm_sipm_pilot/elm_sipm_pilot.h5")
    for (key,path) in [("pilot_four",pilotpath),("pilot_five",oldpath)]
        provenance[key*"_path"]=relpath(path,ROOT)
        provenance[key*"_sha256"]=bytes2hex(sha256(read(path)))
    end
    pilots=create_group(f,"pilots")
    h5open(pilotpath,"r") do source
        for (label,prefix) in [("four",""),("five","reference/")], tag in ("tau_12p5","tau_25")
            g=create_group(pilots,label*"_"*tag)
            src=source[prefix*tag*"/xe_he15"]
            for key in keys(src)
                ndims(src[key])==0 || continue
                value=read(src[key])
                value isa Number && (g[key]=value)
            end
        end
    end
end
println(out)
println("Optical acceptance = ",optical,"; analytical quadrature checks passed.")
