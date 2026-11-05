using Test,TOML
include("../scripts/production/reduce_elm_sipm_pilot.jl")
const R=ELMSiPMReduction
@testset "Shared and legacy electron analysis bounds" begin
 r=TOML.parsefile(joinpath(R.ROOT,"config/specs/realisations/elm_matched_full.toml"))
 e=TOML.parsefile(joinpath(R.ROOT,"config/specs/detector/electrons/elm_4mm_250ns_ep2_shared.toml"))
 @test R.electron_reference_bounds(e,r)==(5.,145.,150.)
 legacy=TOML.parsefile(joinpath(R.ROOT,"config/specs/detector/electrons/elm_4mm_250ns_ep2_drift10cm.toml"))
 @test R.electron_reference_bounds(legacy,r)==(10.,150.,150.)
 original=TOML.parsefile(joinpath(R.ROOT,"config/specs/generated/df1414afcdc2/realisation.toml"))
 @test R.electron_reference_bounds(e,original)==(0.,140.,150.)
end
