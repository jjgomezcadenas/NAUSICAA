@testset "MC selection observables" begin
    # The low-energy component is longer than the leading-energy component.
    h=EventHits(Int64(5),source_frame,Float32[0,1,2,3,100,101],zeros(Float32,6),
        zeros(Float32,6),Float32[.01,.01,.01,.01,.5,.5],zeros(Float32,6),
        ones(Int32,6),Int32.(0:5))
    r=analyze_mc_event(h,MCGraphSettings())
    @test r.event.n_components==2 && !r.event.accepted
    @test r.event.track_length_mm==3f0
    @test sort(r.component_length_mm)==Float32[1,3]
    @test r.event.voxel_energy_MeV≈sum(h.energy_MeV)
    @test r.event.voxel_reduction_error_MeV<1e-6
    r.voxel_energy_MeV[1]=0
    @test h.energy_MeV[1]==.01f0
    r=analyze_mc_event(h,MCGraphSettings(energy_threshold_keV=100))
    @test r.event.accepted && r.event.n_components==1
    @test r.event.track_length_mm==1f0 && r.event.n_retained_voxels==2
    @test r.event.retained_energy_MeV==1f0
    r=analyze_mc_event(h,MCGraphSettings(energy_threshold_keV=1000))
    @test !r.event.accepted && r.event.n_components==0 && r.event.track_length_mm==0f0
end
