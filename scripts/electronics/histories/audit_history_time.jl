include(joinpath(@__DIR__,"../../../src/ion_signal/history_time_audit.jl"))
println(HistoryTimeAudit.run(isempty(ARGS) ? 3790 : parse(Int,ARGS[1])))
