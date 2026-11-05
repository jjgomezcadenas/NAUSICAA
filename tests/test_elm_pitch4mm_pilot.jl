using Test
include("test_elm_sipm_pilot.jl")
@testset "Selected-gas replay and 4 mm ELM cards" begin
    @test Replay.selected_electron_gases(Dict())==["pure_xe","xe_he15"]
    @test Replay.selected_electron_gases(Dict("selected_gases"=>["xe_he15"]))==["xe_he15"]
    @test_throws ArgumentError Replay.selected_electron_gases(Dict("selected_gases"=>String[]))
    @test_throws ArgumentError Replay.selected_electron_gases(Dict("selected_gases"=>["argon"]))
    spec=Replay.load_specs(Prod.ROOT,"config/specs/runs/elm_photons_pitch4mm_ep2.toml")
    @test spec.electrons.data["hole_pitch_mm"]==4
    @test spec.electrons.data["hole_diameter_mm"]==2.5
    @test spec.electrons.data["optics"]["gas_gap_radial_mm"]==.5
    @test Replay.electron_optics(spec.electrons.data).acceptance < Replay.elm_light_acceptance(1.5f0,5f0;tpb_efficiency=.65f0,forward_fraction=.5f0,guide_index=1.5f0,gap_index=1f0)
    pilot=Prod.Pilot.settings(Prod.ROOT,"config/specs/runs/elm_pitch4mm_sipm_response_tau12p5ns.toml")
    @test pilot.sensor.dark_rate_Hz≈700000
    @test pilot.c["illuminated_microcell_count"]==round(Int,39984pi*1.25^2/9)==21808
    @test pilot.run["gases"]==["xe_he15"]
end
