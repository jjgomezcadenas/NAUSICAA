root=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(root,"src/ion_signal/ion_csa_response.jl"))
const C=IonCSAResponse
stage=ARGS[1];card=length(ARGS)>1 ? ARGS[2] : "config/specs/runs/ion_csa_v3_pilot.toml"
if stage=="prepare";C.prepare(root,card)
else
 dir=abspath(ARGS[3])
 if stage=="reduce";C.reduce_cases(root,card,dir)
 elseif stage=="nonlinear";C.nonlinear_cases(root,card,dir)
 elseif stage=="assess";C.assess(root,card,dir)
 elseif stage=="validate";C.validate_cases(root,card,dir);C.validate_nonlinear_controls(root,card,dir)
 elseif stage=="finish";C.finish(root,card,dir)
 elseif stage=="summary";include(joinpath(root,"src/ion_signal/ion_csa_summary.jl"));IonCSASummary.run(C,root,card,dir)
 else;error("Unknown ion CSA stage")
 end
end
