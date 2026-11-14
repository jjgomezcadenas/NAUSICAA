using Test,Statistics
include("../../src/event_analysis/EventAnalysis.jl")
include("../../src/detector_response/DetectorResponse.jl")
include("../../src/electronics/ElectronicsResponse.jl")
include("../../src/reconstruction/SignalReconstruction.jl")
using .EventAnalysis,.DetectorResponse,.ElectronicsResponse,.SignalReconstruction
@testset "Signed event charge, positive hits and shared Z" begin
 ch=ReadoutChannels(Int32[0],Int32[0],Float32[2],Float32[2])
 checks=IonReconstructionChecks(Float32[.9],falses(1),trues(1),falses(1),Float32[0])
 ions=ReconstructedIons(1,ch,IonReadoutTime(100f0,20f0,2),reshape(Float32[10,-3],2,1),trues(2,1),checks)
 pe=ReconstructedPE(1,ch,ElectronReadoutTime(2465f0,5f0,2),2475f0,reshape(Float32[-2,20],2,1),trues(2,1),falses(2,1))
 p=EventPlacement(1,1.65f0,148.35f0)
 ds=(W_eV=21.9f0,electron_transport=(electron_drift_velocity_mm_us=.6f0,),ion_transport=(ion_drift_velocity_mm_ms=.15f0,),electron_lifetime_ms=30f0,el_photons_per_electron=100f0,optical_acceptance=.1f0,chamber_length_mm=1500f0)
 q=event_charge_estimates(ions,pe,p,ds,.32)
 @test q.ions_total==7
 @test q.pe_total≈18*exp(1483.5/.6/1000/30) rtol=1e-6
 ih=reconstruct_hits(ions,p,ds;keV_per_count=q.ion_keV_per_count)
 eh=reconstruct_hits(pe,p,ds;charge_correction=q.attachment_correction,keV_per_count=q.pe_keV_per_count)
 @test ih.charge==Float32[10]
 @test length(eh.charge)==1
 @test ih.z_mm[1]≈0 atol=1e-4
 @test eh.z_mm[1]≈0 atol=1e-4
 @test ih.energy_keV[1]≈.219 rtol=1e-6
 pe.valid[1]=false
 @test isnan(event_charge_estimates(ions,pe,p,ds,.32).pe_total)
 @test event_charge_estimates(ions,pe,p,ds,.32).invalid_pe_buckets==1
end
include("../../scripts/production/track_controls.jl")
@testset "Hit cut and rebuilt kNN" begin
 values=reshape(Float32[4,1,4],3,1);valid=trues(3,1);input=ones(Int64,3,1)
 channels=ReadoutChannels(Int32[0],Int32[0],Float32[0],Float32[0])
 hits=ReconstructedHits(zeros(Float32,3),zeros(Float32,3),Float32[0,10,20],Float32[4,1,4],Float32[4,1,4])
 stats=hit_cut_statistics(values,valid,input,1.,3.)
 @test stats.kept_hits==2
 @test stats.input_energy_before_keV==3
 @test stats.input_energy_after_keV==2
 topology=hit_cut_topology(values,valid,input,channels,hits,1.,3.,4.,1,6.)
 @test topology["components_before"]==1
 @test topology["components_after"]==1
 @test topology["max_distance_mm"]==24
 counts_only=hit_cut_topology(values,valid,input,channels,hits,1.,3.,4.,1,6.;compute_endpoints=false)
 @test counts_only["components_before"]==topology["components_before"]
 @test counts_only["signal_components_after"]==topology["signal_components_after"]
 @test topology["endpoint_max_shift_mm"]==0
 # Distance follows diffusion, even though the hit coordinates are unchanged.
 short=hit_cut_topology(values,valid,input,channels,hits,1.,3.,2.,1,6.)
 @test short["components_after"]==2
end
