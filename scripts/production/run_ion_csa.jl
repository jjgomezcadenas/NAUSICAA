"""Run or publish the ion CSA stage using a resolved production run card."""
started=time_ns()
root=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(root,"src/ion_signal/ion_csa_pipeline.jl"))
card=isempty(ARGS) ? "config/specs/runs/ion_csa_v3_pilot_all.toml" : only(ARGS)
manifest=IonCSAProduction.run(root,card)
seconds=(time_ns()-started)/1e9
IonCSAProduction.write_toml(joinpath(dirname(manifest),"timing.toml"),Dict(
    "entrypoint_elapsed_s"=>seconds,"includes"=>"module loading, calculation, reconstruction, presentation and publication"))
println("Complete ion CSA workflow: ",round(seconds;digits=2)," s")
