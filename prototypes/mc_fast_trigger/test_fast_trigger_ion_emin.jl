using Test,LinearAlgebra
include("fast_trigger_ion_emin_functions.jl")
@testset "Linear ion S/N and temporal covariance" begin
    cv=Matrix{Float64}(I,10,10)*100
    a=ion_emin_moments(1.,20.,.15,cv,4;arrangement="independent_pads")
    b=ion_emin_moments(1.,20.,.15,cv,4;arrangement="consecutive_buckets")
    @test a.ions==50
    @test a.noise_rms_ions==20
    @test a.sn_noise==2.5
    @test a.fano_variance==7.5
    @test a==b
    cv+=Matrix(Tridiagonal(fill(-25.,9),zeros(10),fill(-25.,9)))
    c=ion_emin_moments(1.,20.,.15,cv,4;arrangement="consecutive_buckets")
    @test c.noise_variance==4*100-2*3*25
    @test ion_emin_moments(1.,20.,.15,cv,1;arrangement="consecutive_buckets").noise_rms_ions==10
    @test_throws ArgumentError ion_emin_moments(1.,20.,.15,cv,11;arrangement="consecutive_buckets")
end
