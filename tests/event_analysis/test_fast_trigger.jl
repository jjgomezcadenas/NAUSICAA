using Test,Random,HDF5
@testset "MC fast trigger" begin
    info=MCEventInfo(Int64(7),(0f0,0f0,0f0),Int32[1])
    particles=[MCParticle(Int32(1),Int32(0),"e-"),MCParticle(Int32(2),Int32(1),"gamma"),MCParticle(Int32(3),Int32(2),"e-"),MCParticle(Int32(4),Int32(3),"gamma"),MCParticle(Int32(5),Int32(1),"e-")]
    h=MCEventHits(Int64(7),mc_source_frame,Float32[0,1,30,31,2],zeros(Float32,5),zeros(Float32,5),Float32[100,100,1,2,5],Int32[1,1,3,4,5],Int32[0,1,0,0,0])
    hc=center_hits(h,info);s=FastTriggerSettings(1500f0,50f0,1450f0,(10f0,10f0,10f0),3f0);p=EventPlacement(Int64(7),75f0,75f0)
    roots=gamma_map(reverse(particles),info)
    @test roots==Dict(1=>0,2=>2,3=>2,4=>2,5=>0)
    gs=gamma_deposits(hc,roots)
    @test length(gs)==1 && gs[1].energy_keV==3f0
    @test mc_fast_trigger(hc,info,particles,p,s).decision==trigger_accepted
    hc.energy_keV[4]=nextfloat(2f0)
    hc.energy_keV[3]=nextfloat(1f0)
    @test mc_fast_trigger(hc,info,particles,p,s).decision==outside_frame_gamma
    @test h.energy_keV[4]==2f0
    hc.x_mm[3:4].=11f0
    @test mc_fast_trigger(hc,info,particles,p,s).decision==trigger_accepted
    @test contains(make_frame(primary_box(hc,info),s.margin_mm),(11f0,0f0,0f0))
    @test mc_fast_trigger(hc,info,particles,EventPlacement(Int64(7),5f0,145f0),s).decision==trigger_accepted
    @test mc_fast_trigger(hc,info,particles,EventPlacement(Int64(7),4f0,146f0),s).decision==outside_vertex_fiducial
    hc.z_mm[3]=1600f0
    @test mc_fast_trigger(hc,info,particles,EventPlacement(Int64(7),4f0,146f0),s).decision==outside_detector
    stats=FastTriggerStatistics();o=mc_fast_trigger(hc,info,particles,p,s)
    stats=update_statistics(stats,o)
    @test stats.outside_detector_events==stats.input_events==1
    @test merge_statistics(stats,stats).input_events==2
    @test place_event(Xoshiro(42),Int64(7),1500f0)==place_event(Xoshiro(42),Int64(7),1500f0)
    @test_throws ArgumentError gamma_map([MCParticle(Int32(1),Int32(0),"e-"),MCParticle(Int32(2),Int32(99),"e-")],info)
    @test_throws ArgumentError gamma_map(vcat(particles,particles[1:1]),info)
    @test_throws ArgumentError gamma_map([particles[1],MCParticle(Int32(2),Int32(3),"e-"),MCParticle(Int32(3),Int32(2),"e-")],info)
    @test_throws ArgumentError FastTriggerSettings(1500f0,1450f0,50f0,(10f0,10f0,10f0),3f0)
    @test_throws ArgumentError mc_fast_trigger(h,info,particles,p,s)
end
@testset "MC event reader" begin
    mktempdir() do dir
        path=joinpath(dir,"hits.h5")
        h5open(path,"w") do f
            attributes(f)["schema"]=EventAnalysis.EventIO.SCHEMA
            f["events/event_id"]=Int64[7];f["events/first_row"]=Int64[1];f["events/n_hits"]=Int64[2]
            f["events/first_primary_row"]=Int64[1];f["events/n_primaries"]=Int32[1];f["events/primary_particle_ids"]=Int32[9]
            for axis in ("x","y","z");f["events/vertex_$(axis)_mm"]=Float32[10];f["hits/$(axis)_mm"]=Float32[10,11];end
            f["hits/energy_MeV"]=Float32[.001,.002];f["hits/particle_id"]=Int32[9,9];f["hits/hit_id"]=Int32[0,1]
        end
        h,info=read_mc_event(path,Int64(7))
        @test h.energy_keV==Float32[1,2]
        @test info.primary_particle_ids==Int32[9]
        @test info.generation_vertex_mm==(10f0,10f0,10f0)
        @test center_hits(h,info).z_mm==Float32[0,1]
        @test_throws KeyError read_mc_event(path,Int64(8))
        @test_throws ArgumentError read_mc_event(path,Int64(8);event_row=1)
    end
end

@testset "Main-track frame includes delta electrons" begin
    info=MCEventInfo(7,(0f0,0f0,0f0),Int32[10,20])
    particles=[MCParticle(10,0,"e-"),MCParticle(20,0,"e-"),MCParticle(30,10,"e-"),
        MCParticle(40,30,"e-"),MCParticle(50,40,"gamma"),MCParticle(60,50,"e-"),MCParticle(70,60,"e-")]
    # Use exact Int32 fields in the existing concrete MCParticle constructor.
    h=MCEventHits(7,mc_vertex_frame,Float32[0,2,-35,-45,100,120],zeros(Float32,6),zeros(Float32,6),
        Float32[100,100,20,30,10,5],Int32[10,20,30,40,60,70],Int32[0,0,0,0,0,0])
    b=main_track_box(h,info,reverse(particles))
    @test b.lower_mm==(-45f0,0f0,0f0)
    @test b.upper_mm==(2f0,0f0,0f0)
    frame=make_frame(b,(20f0,20f0,20f0))
    @test all(i->contains(frame,(h.x_mm[i],h.y_mm[i],h.z_mm[i])),1:4)
    @test !contains(frame,(100f0,0f0,0f0))
    gs=gamma_deposits(h,gamma_map(particles,info))
    @test only(gs).energy_keV==15f0
    @test gamma_veto(gs,frame,3f0)
    # A gamma inside the new frame is accepted although outside the direct-primary frame.
    h.x_mm[5:6].=-40f0
    settings=FastTriggerSettings(1500f0,50f0,1450f0,(20f0,20f0,20f0),3f0)
    @test mc_fast_trigger(h,info,particles,EventPlacement(7,75f0,75f0),settings).decision==trigger_accepted
end
