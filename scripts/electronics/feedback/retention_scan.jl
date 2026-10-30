root=normpath(joinpath(@__DIR__,"../../.."))
include(joinpath(root,"src/ion_signal/retention_scan.jl"))
const R=RetentionScan
stage,card=ARGS[1:2]
if stage=="inputs";R.inputs(root,card)
elseif stage=="tune";R.tune(root,card,ARGS[3],ARGS[4])
elseif stage=="analyze";R.analyze(root,card)
elseif stage=="screens";R.screens(root,card)
elseif stage=="controls"
 include(joinpath(root,"src/ion_signal/retention_validation.jl"))
 RetentionValidation.run(R,root,card)
elseif stage=="summary";R.summary(root,card)
else;error("Unknown retention scan stage $stage");end
