using Test,Random,Statistics
include("../../scripts/production/frame_noise_hits.jl")
@testset "Noise hits of the unsimulated frame cells" begin
    pitch=3f0;grid=ReadoutGrid(pitch,(0f0,0f0));v=0.165f0;dt=20f0;L=1300f0
    detector=(ion_grid=grid,chamber_length_mm=L,ion_transport=(ion_drift_velocity_mm_ms=v,))
    p=EventPlacement(Int64(1),65f0,65f0)                                   # vertex at z = 650 mm
    cube=EventFrame(mc_vertex_frame,(-100f0,-100f0,-100f0),(100f0,100f0,100f0))
    box=EventFrame(mc_vertex_frame,(-30f0,-30f0,-30f0),(30f0,30f0,30f0))
    ri=(channels=readout_channels(box,grid),time=ion_readout_time(box,p,L,v,dt))
    pads=readout_channels(cube,grid);time=ion_readout_time(cube,p,L,v,dt)
    sim_pads=Set(zip(ri.channels.ix,ri.channels.iy));free=length(pads.ix)*time.bucket_count-length(ri.channels.ix)*ri.time.bucket_count
    @test length(pads.ix)==67^2 || length(pads.ix)==68^2
    sigma=32.6;scale=21.9/1000;cut=1.5                                      # cut at 2.1 sigma, to have statistics
    tail=ccdf(Normal(),cut/scale/sigma);expected=free*tail
    counts=Int[];all_hits=nothing
    for seed in 1:40
        h=frame_noise_hits(Xoshiro(seed),cube,ri,p,detector,sigma,scale,cut);push!(counts,length(h.x_mm));seed==1&&(all_hits=h)
    end
    @test abs(mean(counts)-expected)<4*sqrt(expected/40)                    # Poisson mean
    @test abs(var(counts)/expected-1)<0.6
    h=all_hits
    @test all(>(Float32(cut)),h.energy_keV)
    @test all(cube.lower_mm[1]-pitch.<=h.x_mm.<=cube.upper_mm[1]+pitch) && all(cube.lower_mm[3]-dt*v.<=h.z_mm.<=cube.upper_mm[3]+dt*v)
    # tail shape: mean excess over the cut of a Gaussian tail
    t0=cut/scale/sigma;mean_tail=sigma*scale*pdf(Normal(),t0)/tail
    @test abs(mean(h.energy_keV)-mean_tail)<5*std(h.energy_keV)/sqrt(length(h.energy_keV))
    # no hit in a simulated cell: pads of the box within the simulated buckets
    k(z)=floor(Int,(Float64(z)+650)/(dt*v));ksim=round(Int,ri.time.start_ms/dt):(round(Int,ri.time.start_ms/dt)+ri.time.bucket_count-1)
    ix(x)=floor(Int32,x/pitch)
    @test !any((ix(h.x_mm[j]),ix(h.y_mm[j])) in sim_pads && k(h.z_mm[j]) in ksim for j in eachindex(h.x_mm))
    @test any((ix(h.x_mm[j]),ix(h.y_mm[j])) in sim_pads for j in eachindex(h.x_mm))   # simulated pads still get hits outside the simulated buckets
    # hits lie on pad centres and bucket centres
    @test all(isapprox(mod(x/pitch,1),0.5;atol=1e-3) for x in h.x_mm)
    @test all(isapprox(mod((Float64(zv)+650)/(dt*v),1),0.5;atol=1e-3) for zv in h.z_mm)
    # reproducible, and no hits when the whole frame is simulated or the cut is far in the tail
    @test frame_noise_hits(Xoshiro(3),cube,ri,p,detector,sigma,scale,cut).energy_keV==frame_noise_hits(Xoshiro(3),cube,ri,p,detector,sigma,scale,cut).energy_keV
    whole=(channels=pads,time=time);@test isempty(frame_noise_hits(Xoshiro(1),cube,whole,p,detector,sigma,scale,cut).x_mm)
    @test isempty(frame_noise_hits(Xoshiro(1),cube,ri,p,detector,sigma,scale,30.0).x_mm)
    a=join_hits(h,h);@test length(a.x_mm)==2length(h.x_mm) && a.energy_keV[end]==h.energy_keV[end]
    @test_throws ArgumentError frame_noise_hits(Xoshiro(1),cube,ri,p,detector,0.0,scale,cut)
end
