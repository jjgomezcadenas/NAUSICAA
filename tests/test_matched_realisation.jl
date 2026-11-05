using Test,TOML
include("../scripts/production/specs.jl")
const S=ProductionSpecs
const ROOT=normpath(joinpath(@__DIR__,".."))
@testset "Shared realization and config-independent primary stream" begin
    old=TOML.parsefile(joinpath(ROOT,"config/specs/generated/df1414afcdc2/realisation.toml"))
    new=TOML.parsefile(joinpath(ROOT,"config/specs/realisations/elm_matched_full.toml"))
    legacyvalues=Dict(k=>v for (k,v) in old if !(k in ("samples","max_events_per_sample")))
    expected=S.content_hash(join(("bb0nu","source-sha",S.canonical(legacyvalues)),"|"))
    a=S.realisation_identities(ROOT,"bb0nu","source-sha",old)
    b=S.realisation_identities(ROOT,"bb0nu","source-sha",new)
    @test a.stream_identity==b.stream_identity==expected
    @test a.realisation_id!=b.realisation_id
    shifted=deepcopy(new);shifted["min_ion_drift_cm"]=10.;shifted["max_ion_drift_cm"]=140.
    shifted["chamber_length_cm"]=160.
    c=S.realisation_identities(ROOT,"bb0nu","source-sha",shifted)
    @test c.stream_identity==b.stream_identity
    @test c.realisation_id!=b.realisation_id
    broken=deepcopy(new);broken["fano_factor"]=.2
    @test_throws ErrorException S.realisation_identities(ROOT,"bb0nu","source-sha",broken)
    spec=S.load_specs(ROOT,"config/specs/runs/elm_matched_detector_smoke.toml")
    @test spec.ions.data["pad_pitch_mm"]==4
    @test spec.ions.data["bucket_width_ms"]==40
    @test !haskey(spec.electrons.data,"placement")
end
