using Test
@testset "Events generated as a gamma" begin
    # Generated gamma 1 from the wall. It produces: the photoelectron 2 (main electron) with its delta electron 3 and a
    # bremsstrahlung gamma 7 absorbed as electron 8; the K X-ray 4 absorbed as electron 5; a Compton electron 6.
    particles=[MCParticle(Int32(1),Int32(0),"gamma"),MCParticle(Int32(2),Int32(1),"e-"),MCParticle(Int32(3),Int32(2),"e-"),MCParticle(Int32(4),Int32(1),"gamma"),
               MCParticle(Int32(5),Int32(4),"e-"),MCParticle(Int32(6),Int32(1),"e-"),MCParticle(Int32(7),Int32(2),"gamma"),MCParticle(Int32(8),Int32(7),"e-")]
    wall=(1800f0,0f0,300f0);start=(200f0,-100f0,600f0)                       # generation point and creation point of electron 2
    generated=MCEventInfo(Int64(4),wall,Int32[1]);lead=LeadingElectron(Int32(2),start,2413f0);info=leading_electron_info(generated,lead)
    @test info.generation_vertex_mm==start && info.primary_particle_ids==Int32[2] && info.event_id==4

    # The main electron and its electron descendants are the main track; every other product of the gamma is its own deposit
    roots=gamma_map(particles,info)
    @test roots[2]==0 && roots[3]==0
    @test roots[1]==1 && roots[4]==4 && roots[5]==4 && roots[6]==6 && roots[7]==7 && roots[8]==7
    @test_throws ArgumentError gamma_map(particles,generated)                # the gamma itself cannot be the primary of the trigger

    # hits: main electron along X, its delta, the gamma's own deposit at the absorption point, the X-ray 10 mm away
    L=1300f0;s=FixedFrameTriggerSettings(L,1800f0,160f0,3f0,true,50f0)
    event(extra_x,extra_e;pid=Int32(6))=begin
        local x=Float32[0,40,80,20,0,10,extra_x];local y=Float32[0,0,0,15,0,0,0];local z=Float32[0,0,0,0,0,0,0]
        local e=Float32[800,800,800,13,5,30,extra_e];local ids=Int32[2,2,2,3,1,5,pid]
        center_hits(MCEventHits(Int64(4),mc_source_frame,x.+start[1],y.+start[2],z.+start[3],e,ids,Int32[0,1,2,0,0,0,0]),info)
    end
    h=event(30f0,0f0);p=placement_from_vertex(info,L)
    @test p.ion_drift_cm==60f0 && p.electron_drift_cm==70f0                  # drift from the track, not from the wall
    box=main_track_box(h,info,particles);@test box.lower_mm==(0f0,0f0,0f0) && box.upper_mm==(80f0,15f0,0f0)
    o=mc_fast_trigger(h,info,particles,p,s);@test o.decision==trigger_accepted && contains(o.frame,box)
    deposits=gamma_deposits(h,roots);@test [d.gamma_particle_id for d in deposits]==Int32[1,4] && [d.energy_keV for d in deposits]==Float32[5,30]

    # A Compton electron far from the track: vetoed above the energy cut, ignored below it
    @test mc_fast_trigger(event(500f0,400f0),info,particles,p,s).decision==outside_frame_gamma
    @test mc_fast_trigger(event(500f0,2f0),info,particles,p,s).decision==trigger_accepted
    # The same deposit inside the frame is kept for the component count
    @test mc_fast_trigger(event(60f0,400f0),info,particles,p,s).decision==trigger_accepted
    # A far deposit below the cut does not combine with the X-ray: each product of the gamma is its own deposit
    @test mc_fast_trigger(event(500f0,2f0),info,particles,p,s).decision==trigger_accepted

    # An event generated as electrons is mapped as before: no particle maps to itself unless it is a gamma
    electrons=[MCParticle(Int32(1),Int32(0),"e-"),MCParticle(Int32(2),Int32(1),"gamma"),MCParticle(Int32(3),Int32(2),"e-"),MCParticle(Int32(5),Int32(1),"e-")]
    r=gamma_map(electrons,MCEventInfo(Int64(9),(0f0,0f0,0f0),Int32[1]));@test r[1]==0 && r[5]==0 && r[2]==2 && r[3]==2
end
