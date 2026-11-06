"""Validated calibration loading and response-card resolution for CSA production."""
module IonCSANonlinearIO
using HDF5, TOML, SHA
include("ion_csa_nonlinear.jl")
const N = IonCSANonlinear
sha(path) = open(io -> bytes2hex(sha256(io)), path)
resolve(root,path) = isabspath(path) ? path : joinpath(root,path)
const DEFAULT_CALIBRATION = "results/ion_csa_nonlinear/v1/calibration.h5"

"""
    load_calibration(path; expected_sha256=nothing)

Read an immutable, independently validated calibration into owned arrays.
Require a completed successful validation and the exact numerical implementation
used during calibration. Return the model, current bound, frame widths, circuit
identity, averaging window and guard. A checksum mismatch or failed validation
raises `ArgumentError` before production starts.
"""
function load_calibration(path;expected_sha256=nothing)
    digest=sha(path)
    expected_sha256===nothing || digest==expected_sha256 || throw(ArgumentError("Nonlinear calibration checksum changed"))
    h5open(path,"r") do h
        String(read(attributes(h)["schema"]))=="nausicaa.ion_csa_nonlinear.v1" || throw(ArgumentError("Nonlinear calibration schema"))
        read(h["complete"]) && read(h["passed"]) || throw(ArgumentError("Nonlinear calibration has not passed validation"))
        String(read(h["model_code_sha256"]))==sha(joinpath(@__DIR__,"ion_csa_nonlinear.jl")) || throw(ArgumentError("Revalidate modified nonlinear implementation"))
        model=N.NonlinearCSA(read(h["knots_V"]),read(h["gain_V_per_e"]),read(h["recovery_V_s"]);
            baseline_V=read(h["baseline_V"]),step_s=read(h["step_s"]))
        (;model,sha256=digest,maximum_current_A=read(h["maximum_current_A"]),
          frame_ms=read(h["validated_frame_ms"]),averaging_s=read(h["averaging_s"]),
          guard_V=read(h["guard_V"]),maximum_error_ENC=read(h["maximum_error_ENC"]),
          circuit_sha256=String(read(h["circuit_sha256"])))
    end
end

"""
    resolve_card(root, card; default_mode="hybrid_spice")

Copy a response card and resolve `response_model` to `hybrid_spice` or
`hybrid_nonlinear`. The nonlinear choice pins the validated calibration by SHA256
and checks the selected circuit, guard and averaging window. The current source
uses the calibrated finite-gap prescription. New circuit/input prescriptions
require a corresponding validation. The production entry point supplies its
default; historical low-level acquisition cards retain their SPICE interpretation.
"""
function resolve_card(root,card;default_mode="hybrid_spice")
    c=copy(card);mode=get(c,"response_model",default_mode)
    execution=get(c,"mode","validation")
    execution in ("production","validation") || throw(ArgumentError("CSA mode must be production or validation"))
    execution=="production" && mode!="hybrid_nonlinear" && throw(ArgumentError("Production mode requires the validated hybrid_nonlinear response"))
    mode in ("hybrid_spice","hybrid_nonlinear") || throw(ArgumentError("Unknown CSA response_model: $mode"))
    # Preserve the canonical identity of historical cards without a mode key.
    mode=="hybrid_spice" && return c
    c["response_model"]=mode
    c["nonlinear_calibration"]=get(c,"nonlinear_calibration",DEFAULT_CALIBRATION)
    k=load_calibration(resolve(root,c["nonlinear_calibration"]);expected_sha256=get(c,"nonlinear_calibration_sha256",nothing))
    c["nonlinear_calibration_sha256"]=k.sha256
    circuit=resolve(root,joinpath(c["calibration_directory"],"decks",c["variant"],"corrected.inc"))
    sha(circuit)==k.circuit_sha256 || throw(ArgumentError("Nonlinear calibration circuit differs"))
    c["averaging_s"]==k.averaging_s && c["initial_spice_guard_V"]==k.guard_V || throw(ArgumentError("Nonlinear calibration aperture/guard differs"))
    c["linear_error_tolerance_ENC"]>=k.maximum_error_ENC || throw(ArgumentError("Requested accuracy requires new validation"))
    # These are the input-current conditions represented by the saved 20/30/40 ms histories.
    for (name,value) in ("input_timing"=>"uniform_mean_arrivals_with_finite_gap_current",
                        "pad_side_m"=>0.0005,"gap_m"=>0.005,"gap_speed_m_s"=>3.,"current_step_s"=>5e-6)
        c[name]==value || throw(ArgumentError("Nonlinear calibration input prescription differs: $name"))
    end
    c
end
end
