using Random, Statistics, DataFrames, HDF5

"""One-particle test hits with arbitrary physical coordinates and energies."""
function diffusion_hits(x=Float32[0],y=zero(x),z=zero(x);energy=ones(Float32,length(x)))
    return EventHits(42,source_frame,x,y,z,energy,zeros(Float32,length(x)),
                     fill(Int32(7),length(x)),Int32.(0:length(x)-1))
end


mutable struct ReplayNormals <: AbstractRNG
    draws::Vector{Float32}
    cursor::Int
end
function Random.randn(r::ReplayNormals,::Type{Float32})
    r.cursor += 1
    return r.draws[r.cursor]
end

@testset "Ion production and thermal scaling" begin
    h = diffusion_hits(Float32[0,1,2,3];energy=Float32[0,0.25,0.75,1.25])
    p = ion_counts(h,500000f0)
    @test p.expected == Float32[0,0.5,1.5,2.5]
    @test p.sampled == Int64[0,0,2,2] # round nearest, ties to even
    @test p.total_expected === 4.5f0 && p.rounding_difference === -0.5f0
    @test p.mean_per_hit === 1.125f0 && p.event_energy_MeV === 2.25f0
    @test p.expected !== h.energy_MeV
    @test_throws ArgumentError ion_counts(h,0f0)
    @test_throws ArgumentError ion_counts(h,NaN32)
    @test_throws ArgumentError ion_counts(h,floatmin(Float32))
    s = sigma_t_ion_mm(293.15f0,100f0,200f0)
    si = sqrt(2*1.380649e-23/1.602176634e-19*293.15*1/20000)*1000
    @test s ≈ si rtol=3e-6 # Petit rounded kB/e versus exact SI constants
    @test sigma_t_ion_mm(293.15f0,400f0,200f0) ≈ 2f0*s
    @test sigma_t_ion_mm(293.15f0,100f0,800f0) ≈ s/2f0
    @test sigma_t_ion_mm(1172.6f0,100f0,200f0) ≈ 2f0*s
    @test sigma_t_ion_mm(293.15f0,0f0,200f0) === 0f0
    @test_throws ArgumentError sigma_t_ion_mm(293.15f0,-1f0,200f0)
    @test_throws ArgumentError sigma_t_ion_mm(293.15f0,1f0,0f0)
end

@testset "Explicit uniform RNG boundaries" begin
    a,b = Xoshiro(52),Xoshiro(52)
    @test sample_drift_cm(a,10f0,10f0) === 10f0
    @test rand(a,UInt64) == rand(b,UInt64)
    a,b = Xoshiro(123),Xoshiro(123)
    @test [sample_drift_cm(a,10f0,150f0) for _ in 1:20] == [sample_drift_cm(b,10f0,150f0) for _ in 1:20]
    drift = [sample_drift_cm(a,10f0,150f0) for _ in 1:20000]
    @test all(10f0 .<= drift .< 150f0)
    @test abs(mean(drift)-80f0) < 1.5f0
    @test var(drift) ≈ (140f0^2/12f0) rtol=0.04
    offsets = [sample_grid_offset_mm(a,2f0) for _ in 1:20000]
    x,y = first.(offsets),last.(offsets)
    @test all(0f0 .<= x .< 2f0) && all(0f0 .<= y .< 2f0)
    @test abs(mean(x)-1f0) < 0.03f0 && abs(mean(y)-1f0) < 0.03f0
    @test abs(cor(x,y)) < 0.03f0
    @test_throws ArgumentError sample_drift_cm(a,20f0,10f0)
    @test_throws ArgumentError sample_grid_offset_mm(a,0f0)
end

@testset "Petit diffusion replay and accounting" begin
    h5open(joinpath(@__DIR__,"fixtures/petit_diffusion_reference.h5"),"r") do f
        h = diffusion_hits(Float32.(read(f["input/x"])),Float32.(read(f["input/y"])),Float32.(read(f["input/z"])))
        rng = ReplayNormals(Float32.(read(f["draws"])),0)
        r = diffuse_xyz_image_mc(rng,h,read(f["input/electrons"]);sigma_t_mm=0.5f0,sigma_l_mm=0.5f0,
            nbins=8,nsigma=2f0,ionisation_energy_eV=20f0)
        @test rng.cursor == 24
        for name in (:x,:y,:z,:energy)
            @test r.bins[!,name] ≈ read(f["output/"*string(name)]) rtol=2e-7
            @test eltype(r.bins[!,name]) == Float32
        end
        @test r.bins.ions == read(f["output/electrons"])
        @test r.sampled_ions == 8 && r.retained_ions == 6 && r.outside_ions == 2
        @test r.retained_ions == sum(r.bins.ions)
        @test h.x_mm == Float32[-1,1]
        @test length(r.x_edges_mm) == 9 && eltype(r.x_edges_mm) == Float32
    end
    edges = Float32[-1,0,1]
    @test [EventAnalysis._histogram_index(x,edges) for x in
        Float32[prevfloat(-1f0),-1,prevfloat(0f0),0,1,nextfloat(1f0)]] == [0,1,1,2,2,0]
    h = diffusion_hits()
    kw = (;sigma_t_mm=1f0,sigma_l_mm=1f0,nbins=20,nsigma=3f0,ionisation_energy_eV=21.9f0)
    r = diffuse_xyz_image_mc(ReplayNormals(Float32[3,0,0, nextfloat(3f0),0,0],0),h,Int64[2];kw...)
    @test r.retained_ions == 1 && r.outside_ions == 1 # final edge included
    a,b = Xoshiro(17),Xoshiro(17)
    empty = diffuse_xyz_image_mc(a,h,Int64[0];kw...)
    @test isempty(empty.bins) && empty.sampled_ions == empty.retained_ions == empty.outside_ions == 0
    @test rand(a,UInt64) == rand(b,UInt64)
    @test eltype(empty.bins.ions) == Int64
    @test diffuse_xyz_image_mc(Xoshiro(17),h,Int64[1000];kw...).bins ==
          diffuse_xyz_image_mc(Xoshiro(17),h,Int64[1000];kw...).bins
    @test_throws DimensionMismatch diffuse_xyz_image_mc(a,h,Int64[];kw...)
    @test_throws ArgumentError diffuse_xyz_image_mc(a,h,Int64[-1];kw...)
    @test_throws ArgumentError diffuse_xyz_image_mc(a,h,Int64[1];kw...,sigma_t_mm=-1f0)
    @test_throws ArgumentError diffuse_xyz_image_mc(a,h,Int64[1];kw...,nbins=0)
    @test_throws ArgumentError diffuse_xyz_image_mc(a,h,Int64[1];kw...,nsigma=Inf32)
end

@testset "Gaussian statistics and finite-window losses" begin
    h = diffusion_hits()
    n = Int64(200000)
    kw = (;sigma_t_mm=1f0,sigma_l_mm=2f0,nbins=100,ionisation_energy_eV=21.9f0)
    r = diffuse_xyz_image_mc(Xoshiro(59),h,[n];kw...,nsigma=5f0)
    d = r.bins
    weights = Float64.(d.ions)/r.retained_ions
    xyz = Float64.(Matrix(d[:,[:x,:y,:z]]))
    mu = vec(sum(xyz.*weights;dims=1))
    centered = xyz .- mu'
    covariance = centered'*(centered.*weights)
    @test all(abs.(mu) .< [0.02,0.02,0.04])
    @test diag(covariance) ≈ [1,1,4] rtol=0.015
    @test abs(covariance[1,2]) < 0.02 && abs(covariance[1,3]) < 0.04 && abs(covariance[2,3]) < 0.04
    @test r.outside_ions < 10
    # Same seed gives the same physical draws for every padding and bin count.
    small = diffuse_xyz_image_mc(Xoshiro(59),h,[n];kw...,nsigma=3f0)
    wide = diffuse_xyz_image_mc(Xoshiro(59),h,[n];kw...,nsigma=4f0)
    @test small.outside_ions >= wide.outside_ions >= r.outside_ions
    expected_loss = 1-0.9973002039367398^3 # exact independent Gaussian ±3sigma box
    @test abs(small.outside_ions/n-expected_loss) < 6sqrt(expected_loss*(1-expected_loss)/n)
    coarse = diffuse_xyz_image_mc(Xoshiro(59),h,[n];kw...,nsigma=3f0,nbins=20)
    @test coarse.retained_ions == small.retained_ions
    @test sum(Float64,small.bins.energy) ≈ small.retained_ions*Float64(21.9f0)/1e6 rtol=2e-7
    @test all(small.bins.ions .> 0)
end

@testset "Zero diffusion, degenerate axes and input validation" begin
    h = diffusion_hits(Float32[1,1],Float32[0,0],Float32[-32,-32])
    r = ReplayNormals(zeros(Float32,15),0)
    kw = (;sigma_t_mm=0f0,sigma_l_mm=0f0,nbins=100,nsigma=3f0,ionisation_energy_eV=20f0)
    out = diffuse_xyz_image_mc(r,h,Int64[2,3];kw...)
    @test r.cursor == 15
    @test out.retained_ions == 5 && out.outside_ions == 0
    @test out.bins.x == Float32[1] && out.bins.y == Float32[0] && out.bins.z == Float32[-32]
    @test out.bins.ions == Int64[5]
    @test out.x_edges_mm == Float32[prevfloat(1f0),nextfloat(1f0)]
    @test out.y_edges_mm == Float32[prevfloat(0f0),nextfloat(0f0)]
    longitudinal = diffuse_xyz_image_mc(Xoshiro(9),h,Int64[100,100];kw...,sigma_t_mm=1f0)
    @test all(longitudinal.bins.z .== -32f0)
    @test length(longitudinal.z_edges_mm) == 2
    a,b = Xoshiro(12),Xoshiro(12)
    @test_throws OverflowError diffuse_xyz_image_mc(a,h,Int64[typemax(Int64),1];kw...)
    @test rand(a,UInt64) == rand(b,UInt64)
    far = diffusion_hits(Float32[1f8])
    @test_throws ArgumentError diffuse_xyz_image_mc(a,far,Int64[1];kw...,sigma_t_mm=0.1f0)
end
