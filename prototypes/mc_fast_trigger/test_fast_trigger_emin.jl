using Test
include("fast_trigger_emin_functions.jl")
base=(energy_keV=10.,w_eV=20.,fano=.15,survival=1.,yield_per_e=100.,J=.03,acceptance=.1,pde=.3,cells=1,dark_per_cell=.2,noise_pe=.5)
@testset "Emin moments" begin
    x=emin_moments(;base...)
    @test x.primary_electrons==500
    @test x.pe==1500
    @test x.attachment_variance==0
    @test x.optical_variance+x.pde_variance≈500*100*.03*(1-.03)
    y=emin_moments(;merge(base,(cells=4,))...)
    @test y.pe==x.pe
    @test y.noise_rms_pe≈2x.noise_rms_pe
    @test y.sn_noise≈x.sn_noise/2
    z=emin_moments(;merge(base,(survival=.5,))...)
    @test z.pe==x.pe/2
    @test z.attachment_variance≈500*.5*.5*3^2
    ideal=emin_moments(;merge(base,(fano=0.,J=0.,acceptance=1.,pde=1.,dark_per_cell=0.,))...)
    @test ideal.total_rms_pe==.5
end
