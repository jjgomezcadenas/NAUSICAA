using Test
include("../../scripts/production/run_framed_detector_response.jl")
@testset "Event 604 main-track frame" begin
    c=TOML.parsefile(joinpath(ROOT,"config/readout/single_electrons_detector_response_3mm_mainframe20.toml"))
    h,info=read_mc_event(joinpath(ROOT,c["hits"]),Int64(604));h=center_hits(h,info)
    particles=read_mc_particles(expanduser(c["particles"]),Int64[604])[604]
    rng=Xoshiro(c["seed"])
    ids=h5open(joinpath(ROOT,c["hits"]),"r") do f;read(f["events/event_id"]);end
    row=findfirst(==(604),ids)
    placements=[place_event(rng,id,1500f0) for id in ids[1:row]];p=last(placements)
    ft=FastTriggerSettings(1500f0,50f0,1450f0,(20f0,20f0,20f0),3f0)
    o=mc_fast_trigger(h,info,particles,p,ft)
    @test o.decision==trigger_accepted
    @test o.frame.lower_mm[1]≈-59.855103f0
    @test all(i->contains(o.frame,(h.x_mm[i],h.y_mm[i],h.z_mm[i])),eachindex(h.x_mm))
    selected=hits_in_frame(h,o.frame)
    @test selected.energy_keV==h.energy_keV
    s=event_settings(c,p)
    pairs=sample_fano_pairs(Xoshiro(seedfor(c["seed"],604,"pairs")),selected,s.W_eV,s.fano_factor)
    result=sim_ion_detector_response(Xoshiro(seedfor(c["seed"],604,"ions")),selected,pairs,p,o.frame,s)
    @test result.losses.excluded_total==0
    @test sum(result.ions.counts)==sum(pairs)
end
