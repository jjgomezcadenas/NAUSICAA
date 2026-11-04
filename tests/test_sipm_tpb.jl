using Test,Random,Distributions,Statistics
include("../scripts/production/run_sipm_tpb_pilot.jl")
const T=TPBPilot
const EA=T.EA
@testset "Square aperture and TPB optical fluctuations" begin
    # Independent polar-angle integral of the depth-averaged solid angle:
    # at azimuth phi the square boundary is r=a/max(|cos phi|,|sin phi|).
    a=.65;h=5.;n=100000
    polar=sum(begin
        phi=2pi*(i-.5)/n;r=a/max(abs(cos(phi)),abs(sin(phi)))
        (h+r-hypot(h,r))/(2h)
    end for i in 1:n)/n
    square=EA.axial_square_acceptance(1.3,5.,2.5)
    @test square≈polar rtol=1e-6
    @test square≈EA.axial_square_acceptance(13.,50.,25.)
    @test EA.axial_light_acceptance(.65f0,5f0)<square<EA.axial_light_acceptance(Float32(.65sqrt(2)),5f0)
    @test_throws ArgumentError EA.axial_square_acceptance(5.,5.,2.5)
    c=Dict("shape"=>"square","side_m"=>.0013,"pde"=>.32,
        "tpb"=>Dict("conversion_efficiency"=>.62,"forward_fraction"=>.5))
    @test T.PIPE.sensor_area_m2(c)≈1.69e-6
    @test T.PIPE.effective_pde(c)≈.0992
    c["tpb"]["conversion_efficiency"]=1.2
    @test_throws ArgumentError T.PIPE.effective_pde(c)
    # Nested aperture + TPB/PDE must give the full direct binomial variance.
    rng=Xoshiro(17);trials=10000;ac=Float64(EA.axial_light_acceptance(2.5f0,5f0));p=.0992
    values=[rand(rng,Binomial(rand(rng,Binomial(rand(rng,Binomial(trials,ac)),square/ac)),p)) for _ in 1:50000]
    target=trials*square*p;variance=target*(1-square*p)
    @test abs(mean(values)-target)<6sqrt(variance/length(values))
    @test var(values)≈variance rtol=.04
    grid=EA.GridSpec(10f0,.25f0,(0f0,0f0,0f0),(0f0,0f0),Int8(-1))
    cells=EA.SparseCells(Int64[0,1],Int64[0,1],Int64[-1,2],Int64[300,700])
    hist=EA.PhotonSiPMHistory(EA.EventPlacement(Int64(1),50f0,100f0),:pure_xe,grid,cells,Int64(100),Int64(99),Int64(10000))
    thin=T.thin_history(Xoshiro(12),hist,.3)
    @test sum(thin.cells.count)<=sum(cells.count)
    @test thin.emitted_photons==hist.emitted_photons
    @test thin.initial_pairs==hist.initial_pairs
    @test T.thin_history(Xoshiro(12),hist,0.).cells.count==Int64[]
    @test T.thin_history(Xoshiro(12),hist,1.).cells.count==cells.count
    @test cells.count==[300,700]
end
