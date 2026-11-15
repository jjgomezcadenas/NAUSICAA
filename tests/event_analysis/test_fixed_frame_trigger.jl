using Test
@testset "Fixed-frame fast trigger" begin
    L=1300f0;R=1800f0;s=FixedFrameTriggerSettings(L,R,200f0,3f0)
    @test_throws ArgumentError FixedFrameTriggerSettings(L,R,0f0,3f0)
    @test_throws ArgumentError FixedFrameTriggerSettings(L,50f0,200f0,3f0)
    # particles: primary 1; gamma 2 with electron 3; delta electron 5 of the primary
    particles=[MCParticle(Int32(1),Int32(0),"e-"),MCParticle(Int32(2),Int32(1),"gamma"),MCParticle(Int32(3),Int32(2),"e-"),MCParticle(Int32(5),Int32(1),"e-")]
    event(v;x=Float32[0,60,0,0],y=Float32[0,0,0,40],z=Float32[0,0,30,-50],e=Float32[100,100,5,20])=begin
        local vertex_info=MCEventInfo(Int64(9),v,Int32[1])
        local source_hits=MCEventHits(Int64(9),mc_source_frame,x.+v[1],y.+v[2],z.+v[3],e,Int32[1,1,3,5],Int32[0,1,0,0])
        center_hits(source_hits,vertex_info),vertex_info
    end
    trig(h,info)=mc_fast_trigger(h,info,particles,placement_from_vertex(info,L),s)

    # Placement: the vertex as generated; drifts add up to the chamber length
    h,info=event((100f0,-200f0,400f0));p=placement_from_vertex(info,L)
    @test p.ion_drift_cm==40f0 && p.electron_drift_cm==90f0 && p.event_id==9
    @test 10f0*(p.ion_drift_cm+p.electron_drift_cm)≈L
    @test_throws ArgumentError placement_from_vertex(MCEventInfo(Int64(9),(0f0,0f0,-1f0),Int32[1]),L)
    @test_throws ArgumentError placement_from_vertex(MCEventInfo(Int64(9),(0f0,0f0,1301f0),Int32[1]),L)
    @test_throws ArgumentError mc_fast_trigger(h,info,particles,EventPlacement(Int64(9),50f0,80f0),s)   # not the generated vertex

    # Accepted: the frame is the cube centred on the main-track box (primary and delta electron)
    o=trig(h,info)
    @test o.decision==trigger_accepted
    @test main_track_box(h,info,particles).lower_mm==(0f0,0f0,-50f0) && main_track_box(h,info,particles).upper_mm==(60f0,40f0,0f0)
    @test o.frame.lower_mm==(-70f0,-80f0,-125f0) && o.frame.upper_mm==(130f0,120f0,75f0)
    @test all(o.frame.upper_mm[i]-o.frame.lower_mm[i]==200f0 for i in 1:3)
    @test box_side(main_track_box(h,info,particles))==60f0

    # Track larger than the frame: exactly 200 mm fits, more does not
    h2,i2=event((0f0,0f0,650f0);x=Float32[0,200,0,0]);@test trig(h2,i2).decision==trigger_accepted
    h2,i2=event((0f0,0f0,650f0);x=Float32[0,200.5,0,0]);@test trig(h2,i2).decision==track_larger_than_frame
    h2,i2=event((0f0,0f0,650f0);z=Float32[0,0,30,-201]);@test trig(h2,i2).decision==track_larger_than_frame

    # Frame outside the detector: at both ends of Z and at the wall
    h2,i2=event((0f0,0f0,100f0));@test trig(h2,i2).decision==frame_outside_detector            # cube reaches z = -25 mm
    h2,i2=event((0f0,0f0,125f0));@test trig(h2,i2).decision==trigger_accepted                  # cube starts at z = 0
    h2,i2=event((0f0,0f0,1230f0));@test trig(h2,i2).decision==frame_outside_detector           # cube ends at z = 1305 mm
    h2,i2=event((0f0,0f0,1225f0));@test trig(h2,i2).decision==trigger_accepted                 # cube ends at z = 1300 mm
    h2,i2=event((1700f0,0f0,650f0));@test trig(h2,i2).decision==frame_outside_detector         # corner at x = 1830 mm, beyond R
    h2,i2=event((1500f0,0f0,650f0));@test trig(h2,i2).decision==trigger_accepted

    # Gamma veto: a deposit above the cut must be fully inside the cube
    h2,i2=event((0f0,0f0,650f0);z=Float32[0,0,80,-50]);@test trig(h2,i2).decision==outside_frame_gamma   # gamma at z = 80, cube ends at 75
    h2,i2=event((0f0,0f0,650f0);z=Float32[0,0,80,-50],e=Float32[100,100,3,20]);@test trig(h2,i2).decision==trigger_accepted   # 3 keV is not above the cut
    h2,i2=event((0f0,0f0,650f0);z=Float32[0,0,75,-50]);@test trig(h2,i2).decision==trigger_accepted      # on the edge

    # Outside the detector: main-track hits, and gamma deposits above the cut, beyond the radius or the ends
    h2,i2=event((1790f0,0f0,650f0);x=Float32[0,20,0,0]);@test trig(h2,i2).decision==outside_detector     # primary hit at r = 1810 mm
    h2,i2=event((0f0,0f0,20f0));@test trig(h2,i2).decision==outside_detector                              # delta electron at z = -30 mm
    h2,i2=event((1700f0,0f0,650f0);x=Float32[0,20,150,0]);@test trig(h2,i2).decision==outside_detector   # 5 keV gamma deposit at r = 1850 mm
    h2,i2=event((1700f0,0f0,650f0);x=Float32[0,20,150,0],e=Float32[100,100,2,20]);@test trig(h2,i2).decision==frame_outside_detector   # 2 keV deposit ignored; the cube still crosses the wall

    # Statistics: disjoint counts, merged by addition; legacy decisions are refused
    st=FixedFrameStatistics()
    for (hh,ii) in (event((0f0,0f0,650f0)),event((0f0,0f0,100f0)),event((0f0,0f0,650f0);x=Float32[0,300,0,0]),event((0f0,0f0,650f0);z=Float32[0,0,80,-50]),event((0f0,0f0,20f0)))
        st=update_statistics(st,trig(hh,ii))
    end
    @test (st.input_events,st.accepted_events,st.frame_outside_detector_events,st.track_larger_than_frame_events,st.outside_frame_gamma_events,st.outside_detector_events)==(5,1,1,1,1,1)
    @test merge_statistics(st,st).input_events==10
    @test_throws ArgumentError update_statistics(st,FastTriggerOutput(outside_vertex_fiducial,o.frame))

    # Simulation box: main track and contained gamma deposits plus the margin, limited to the cube
    @test simulation_box(h,info,particles,o.frame,nothing)==o.frame
    b=simulation_box(h,info,particles,o.frame,20f0)
    @test b.lower_mm==(-20f0,-20f0,-70f0) && b.upper_mm==(80f0,60f0,50f0)          # gamma deposit at z = 30 included
    b=simulation_box(h,info,particles,o.frame,500f0);@test b.lower_mm==o.frame.lower_mm && b.upper_mm==o.frame.upper_mm

    # Inputs are preserved
    @test h.z_mm==Float32[0,0,30,-50] && info.generation_vertex_mm==(100f0,-200f0,400f0)
end

@testset "Fixed frame moved inside the detector" begin
    L=1300f0;R=1800f0;centred=FixedFrameTriggerSettings(L,R,160f0,3f0);inside=FixedFrameTriggerSettings(L,R,160f0,3f0,true)
    @test !centred.inside_detector && inside.inside_detector
    @test_throws ArgumentError FixedFrameTriggerSettings(L,R,2600f0,3f0,true)
    particles=[MCParticle(Int32(1),Int32(0),"e-"),MCParticle(Int32(2),Int32(1),"gamma"),MCParticle(Int32(3),Int32(2),"e-"),MCParticle(Int32(5),Int32(1),"e-")]
    event(v;x=Float32[0,60,0,0],y=Float32[0,0,0,40],z=Float32[0,0,30,-50],e=Float32[100,100,5,20])=begin
        local vertex_info=MCEventInfo(Int64(9),v,Int32[1])
        local source_hits=MCEventHits(Int64(9),mc_source_frame,x.+v[1],y.+v[2],z.+v[3],e,Int32[1,1,3,5],Int32[0,1,0,0])
        center_hits(source_hits,vertex_info),vertex_info
    end
    trig(h,info,s)=mc_fast_trigger(h,info,particles,placement_from_vertex(info,L),s)
    detector(f,i)=(ntuple(k->f.lower_mm[k]+i.generation_vertex_mm[k],3),ntuple(k->f.upper_mm[k]+i.generation_vertex_mm[k],3))

    # Away from the boundary the cube is the centred one
    h,i=event((0f0,0f0,650f0));@test trig(h,i,inside).frame==trig(h,i,centred).frame && trig(h,i,inside).decision==trigger_accepted

    # Near the TopMetal end: track at detector Z 10..60 mm; the centred cube reaches Z = -45 mm
    h,i=event((0f0,0f0,60f0))
    @test trig(h,i,centred).decision==frame_outside_detector
    o=trig(h,i,inside);lo,hi=detector(o.frame,i)
    @test o.decision==trigger_accepted && 0f0<=lo[3]<0.01f0 && isapprox(hi[3],160f0;atol=0.01)
    @test all(hi[k]-lo[k]≈160f0 for k in 1:3) && contains(o.frame,main_track_box(h,i,particles))
    # Near the EL end: track at detector Z 1240..1290 mm (the 2 keV gamma deposit at 1320 mm is below the cut)
    h,i=event((0f0,0f0,1290f0);e=Float32[100,100,2,20]);o=trig(h,i,inside);lo,hi=detector(o.frame,i)
    @test trig(h,i,centred).decision==frame_outside_detector && o.decision==trigger_accepted && 1299.99f0<hi[3]<=1300f0

    # Near the wall: track at detector X 1700..1760; the centred cube has a corner at radius 1813 mm
    h,i=event((1700f0,0f0,650f0))
    @test trig(h,i,centred).decision==frame_outside_detector
    o=trig(h,i,inside);lo,hi=detector(o.frame,i)
    @test o.decision==trigger_accepted
    @test maximum(hypot(x,y) for x in (lo[1],hi[1]), y in (lo[2],hi[2])) <= R
    @test isapprox(maximum(hypot(x,y) for x in (lo[1],hi[1]), y in (lo[2],hi[2])),R;atol=0.01)
    # Track touching the wall: no cube that contains it fits inside the detector
    h,i=event((1720f0,0f0,650f0);x=Float32[0,79,0,0])
    @test trig(h,i,inside).decision==frame_outside_detector

    # The gamma veto uses the moved cube: the gamma deposit at Z + 30 mm is inside it
    h,i=event((0f0,0f0,60f0);z=Float32[0,0,95,-50])            # gamma at detector Z = 155 mm
    @test trig(h,i,inside).decision==trigger_accepted
    h,i=event((0f0,0f0,60f0);z=Float32[0,0,105,-50])           # gamma at detector Z = 165 mm, beyond the cube
    @test trig(h,i,inside).decision==outside_frame_gamma
end

@testset "Main track away from the EL" begin
    L=1300f0;R=1800f0;s=FixedFrameTriggerSettings(L,R,160f0,3f0,true,50f0)
    @test s.el_exclusion_mm==50f0 && FixedFrameTriggerSettings(L,R,160f0,3f0,true).el_exclusion_mm==0f0
    @test_throws ArgumentError FixedFrameTriggerSettings(L,R,160f0,3f0,true,-1f0)
    particles=[MCParticle(Int32(1),Int32(0),"e-"),MCParticle(Int32(2),Int32(1),"gamma"),MCParticle(Int32(3),Int32(2),"e-"),MCParticle(Int32(5),Int32(1),"e-")]
    event(v;z=Float32[0,0,-30,-50],e=Float32[100,100,5,20])=begin
        local vertex_info=MCEventInfo(Int64(9),v,Int32[1])
        local source_hits=MCEventHits(Int64(9),mc_source_frame,Float32[0,60,0,0].+v[1],Float32[0,0,0,40].+v[2],z.+v[3],e,Int32[1,1,3,5],Int32[0,1,0,0])
        center_hits(source_hits,vertex_info),vertex_info
    end
    trig(h,i,set=s)=mc_fast_trigger(h,i,particles,placement_from_vertex(i,L),set)
    # The main track spans detector Z = vertex - 50 .. vertex
    h,i=event((0f0,0f0,1250f0));@test trig(h,i).decision==trigger_accepted            # top of the track exactly 50 mm from the EL
    h,i=event((0f0,0f0,1251f0));@test trig(h,i).decision==main_track_near_el
    @test trig(h,i,FixedFrameTriggerSettings(L,R,160f0,3f0,true)).decision==trigger_accepted   # no exclusion: accepted
    # The cut looks at the whole main track, not at the vertex: vertex at 1200 mm, a primary hit at 1260 mm
    h,i=event((0f0,0f0,1200f0);z=Float32[0,60,-30,-50]);@test trig(h,i).decision==main_track_near_el
    # A gamma deposit near the EL does not trigger it
    h,i=event((0f0,0f0,1240f0);z=Float32[0,0,40,-50]);@test trig(h,i).decision==trigger_accepted
    st=update_statistics(FixedFrameStatistics(),trig(event((0f0,0f0,1251f0))...))
    @test st.main_track_near_el_events==1 && st.input_events==1 && st.accepted_events==0
end
