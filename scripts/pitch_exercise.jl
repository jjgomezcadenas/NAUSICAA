include(joinpath(@__DIR__,"../src/ion_signal/pitch_exercise.jl"))
root=normpath(joinpath(@__DIR__,".."))
length(ARGS)==1 || error("Choose inputs, analyze or checks")
if ARGS[1]=="inputs";PitchExercise.inputs(root)
elseif ARGS[1]=="analyze";PitchExercise.analyze(root)
elseif ARGS[1]=="checks";PitchExercise.checks(root)
else;error("Unknown stage")
end
