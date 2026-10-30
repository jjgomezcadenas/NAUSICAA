using Test
include(joinpath(@__DIR__,"../../scripts/production/specs.jl"))
using .ProductionSpecs

@testset "Production TOML resolver and identities" begin
    root=normpath(joinpath(@__DIR__,"../.."))
    spec=load_specs(root,"config/specs/runs/production_pilot.toml")
    @test spec.realisation.data["chamber_length_cm"]==150.0
    @test spec.electrons.data["reduced_field_kV_cm_bar"]==3.0
    @test spec.electrons.data["j_variance_ratio"]==0.01
    @test spec.electrons.data["j_status"]=="oliveira_fig6_user_reading"
    @test isfile(spec.sample_paths["bb0nu"])
    @test canonical(Dict("b"=>2,"a"=>1))==canonical(Dict("a"=>1,"b"=>2))
    @test content_hash("abc")=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    @test event_seed("identity",1,"pairs")==event_seed("identity",1,"pairs")
    @test event_seed("identity",1,"pairs")!=event_seed("identity",1,"ion_phase")
    @test pair_checksum(Int64[1,2])==pair_checksum(Int64[1,2])
    @test pair_checksum(Int64[1,2])!=pair_checksum(Int64[2,1])
end
