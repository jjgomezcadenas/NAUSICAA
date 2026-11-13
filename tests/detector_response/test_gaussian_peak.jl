using Test, Random, Statistics
include("../../src/analysis/gaussian_peak.jl")
using .GaussianPeak
rng=Xoshiro(710)
x=2458 .+ 3 .* randn(rng,100000)
@testset "Truncated Gaussian likelihood" begin
 for w in (2.,2.5,3.)
  r=fit_peak(x;window_sigma=w)
  @test abs(r.mu-2458)<.04
  @test abs(r.sigma-3)<.04
  @test r.lower<r.mu<r.upper
 end
 tail=vcat(x,fill(2200.,2000))
 r=fit_peak(tail)
 @test abs(r.sigma-3)<.05
 @test abs(r.mu-2458)<.04
 @test r.entries<length(tail)
end
