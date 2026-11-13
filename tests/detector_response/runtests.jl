using Test, Random, Statistics
include("../../src/event_analysis/EventAnalysis.jl")
using .EventAnalysis
include("../../src/detector_response/DetectorResponse.jl")
using .DetectorResponse
function settings(;sigma=0f0,acceptance=1f0,yield=10f0,J=0f0,lifetime=1f20)
 DetectorResponseSettings(1500f0,20f0,.15f0,Xe_85_He_15,ReadoutGrid(4f0,(0f0,0f0)),ReadoutGrid(4f0,(0f0,0f0)),
  1f0,1f0,IonTransportSettings(1f0,sigma,sigma),ElectronTransportSettings(1f0,sigma,sigma),
  lifetime,yield,J,acceptance,:blue_tpb)
end
hits(x,y,z,e=fill(1f0,length(x)))=MCEventHits(1,mc_vertex_frame,Float32.(x),Float32.(y),Float32.(z),Float32.(e),fill(Int32(1),length(x)),Int32.(1:length(x)))
p=EventPlacement(1,75f0,75f0)
f=EventFrame(mc_vertex_frame,(-1f0,-1f0,-2f0),(5f0,5f0,2f0))
@testset "Frame, whole channels and times" begin
 ch=readout_channels(f,settings().ion_grid)
 @test length(ch.ix)==9
 @test Set(ch.ix)==Set(Int32[-1,0,1])
 @test length(readout_channels(EventFrame(mc_vertex_frame,(0f0,0f0,0f0),(4f0,4f0,4f0)),settings().ion_grid).ix)==1
 edge=EventFrame(mc_vertex_frame,(-1f0,-1f0,-1000f0),(5f0,5f0,1000f0))
 clipped=acquisition_frame(edge,p,1500f0)
 @test clipped.lower_mm[3]==-750f0
 @test clipped.upper_mm[3]==750f0
 it=ion_readout_time(f,p,1500f0,1f0,1f0)
 et=electron_readout_time(f,p,1500f0,1f0,1f0)
 @test it.start_ms==et.start_us==748f0
 @test it.bucket_count==et.bucket_count==4
 h=hits([0,0],[0,0],[-1,1]);pairs=Int64[10,20]
 i=sim_ion_detector_response(Xoshiro(1),h,pairs,p,f,settings())
 e=sim_electron_detector_response(Xoshiro(2),h,pairs,p,f,settings())
 @test vec(sum(i.ions.counts;dims=2))==[0,10,0,20]
 @test vec(sum(e.photons.counts;dims=2))==[0,200,0,100]
 @test count(==(0),vec(sum(i.ions.counts;dims=1)))==8
 # Independent known channel counts: a position beyond x=5 still falls in acquired cell 1.
 cells=Dict((Int64(1),Int64(0),Int64(1))=>Int64(7),(Int64(2),Int64(0),Int64(1))=>Int64(3),
  (Int64(2),Int64(0),Int64(0))=>Int64(2))
 v,l,o=DetectorResponse._vectors(cells,ch,it)
 @test sum(v.counts)==7
 @test (l.excluded_xy,l.excluded_time,l.excluded_total)==(5,2,5)
 @test (l.channels_with_signal,l.outside_channels_with_signal)==(2,1)
 @test sum(o.count)==5
 @test_throws ArgumentError DetectorResponse._vectors(cells,ch,it;max_vector_bytes=1)
 @test electron_readout_time(f,p,1500f0,1f0,125f0).start_us/1000 ≈ ion_readout_time(f,p,1500f0,1000f0,.125f0).start_ms atol=1e-6
end
@testset "Frame selection and charge" begin
 h=hits([0,20],[0,0],[0,0],[5,2.9])
 inside=hits_in_frame(h,f)
 @test inside.energy_keV==[5f0]
 @test h.energy_keV==[5f0,2.9f0]
 pairs=sample_fano_pairs(Xoshiro(3),inside,20f0,.15f0)
 @test length(pairs)==1
 r=sim_ion_detector_response(Xoshiro(4),inside,pairs,p,f,settings(sigma=5f0))
 @test sum(r.ions.counts)+r.losses.excluded_total==sum(pairs)
 @test r.losses.excluded_total>0
 r2=sim_ion_detector_response(Xoshiro(4),inside,pairs,p,f,settings(sigma=5f0))
 @test r2.ions.counts==r.ions.counts
 e=sim_electron_detector_response(Xoshiro(5),inside,pairs,p,f,settings(sigma=5f0,acceptance=.3f0,lifetime=30f0))
 @test sum(e.photons.counts)+e.losses.excluded_total+e.uncollected_photons==e.EL_photons
 @test e.e_lost_attachment<=sum(pairs)
 @test vector_statistics(r.ions).total==sum(r.ions.counts)
 @test_throws ArgumentError sim_ion_detector_response(Xoshiro(1),h,Int64[1,1],p,f,settings())
end
@testset "Sampling moments" begin
 rng=Xoshiro(42);n=20000
 a=[sample_integer_moments(rng,100.,.15) for _ in 1:n]
 @test abs(mean(a)-100)<.12
 @test abs(var(a)-15)<.6
 b=[sample_el_photons(rng,Int64(10),100f0,.03f0) for _ in 1:n]
 @test abs(mean(b)-1000)<.2
 @test abs(var(b)-30)<1.2
 c=[sample_incident_photons(rng,Int64(1000),.1f0) for _ in 1:n]
 @test abs(mean(c)-100)<.3
 @test abs(var(c)-90)<3
 a=sample_attached_electrons(rng,fill(Int64(100),n),3f0,30f0)
 prob=exp(-.1)
 @test abs(mean(a)-100prob)<.1
 @test abs(var(a)-100prob*(1-prob))<.4
end
@testset "Production collection comparison" begin
 h=hits([1,2],[1,2],[-1,1]);pairs=Int64[1000,2000]
 wide=EventFrame(mc_vertex_frame,(-100f0,-100f0,-100f0),(100f0,100f0,100f0))
 oldh=EventHits(1,source_frame,h.x_mm,h.y_mm,h.z_mm,h.energy_keV./1000f0,zeros(Float32,2),h.particle_id,h.hit_id)
 oldgrid=GridSpec(4f0,1f0,(0f0,0f0,0f0),(0f0,0f0),Int8(1))
 old=collect_ion_pad_history(Xoshiro(12),oldh,pairs,p,oldgrid;sigma_t_mm=1f0,sigma_l_mm=1f0)
 new=sim_ion_detector_response(Xoshiro(12),h,pairs,p,wide,settings(sigma=1f0))
 lookup=Dict((Int64(x),Int64(y))=>j for (j,(x,y)) in enumerate(zip(new.ions.channels.ix,new.ions.channels.iy)))
 @test sum(new.ions.counts)==sum(old.cells.count)==3000
 @test new.losses.excluded_total==0
 for j in eachindex(old.cells.count)
   @test new.ions.counts[old.cells.bucket[j]+101,lookup[(old.cells.pixel_ix[j],old.cells.pixel_iy[j])]]==old.cells.count[j]
 end
end

@testset "Detector edge and zero signal" begin
 pe=EventPlacement(1,0f0,150f0)
 fe=EventFrame(mc_vertex_frame,(-4f0,-4f0,0f0),(4f0,4f0,4f0))
 he=hits([0],[0],[0])
 r=sim_ion_detector_response(Xoshiro(5),he,Int64[10000],pe,fe,settings(sigma=1f0))
 @test 4500 < r.losses.excluded_time < 5500
 @test sum(r.ions.counts)+r.losses.excluded_total==10000
 zero=sim_ion_detector_response(Xoshiro(5),he,Int64[0],pe,fe,settings())
 @test all(iszero,zero.ions.counts)
 @test vector_statistics(zero.ions).mean==0
end
@testset "Independent ion settings" begin
 h=hits([0,1],[0,1],[0,1]);pairs=Int64[1000,2000];s=settings(sigma=1f0,acceptance=.3f0,lifetime=30f0)
 other=DetectorResponseSettings((name==:ion_grid ? ReadoutGrid(3f0,(1f0,1f0)) : getfield(s,name) for name in fieldnames(DetectorResponseSettings))...)
 a=sim_electron_detector_response(Xoshiro(19),h,pairs,p,f,s)
 b=sim_electron_detector_response(Xoshiro(19),h,pairs,p,f,other)
 @test a.photons.counts==b.photons.counts
 @test a.e_lost_attachment==b.e_lost_attachment
 @test a.EL_photons==b.EL_photons
end
