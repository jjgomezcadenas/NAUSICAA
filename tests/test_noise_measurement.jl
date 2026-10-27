using HDF5, Test
include("../src/ion_signal/averaging_scan.jl")
@testset "Stationary noise measurement" begin
    # White-noise RMS has a known rectangular-band integral.
    f=[.1,1e6]; asd=[2e-6,2e-6]
    @test AveragingScan.noise_rms(f,asd,0.,.02,0.) ≈ 2e-6*sqrt(f[2]-f[1]) rtol=1e-12
    # Independent cosine quadrature checks the covariance sign and normalization.
    hi=1000.;lo=.1;alpha=.99;dt=.02
    expected=2e-6*sqrt((1+alpha^2)*(hi-lo)-2alpha*(sin(2pi*hi*dt)-sin(2pi*lo*dt))/(2pi*dt))
    @test AveragingScan.noise_rms([lo,hi],asd,0.,dt,alpha) ≈ expected rtol=1e-12
    h5open(joinpath(@__DIR__,"../results/ltspice/noise_measurement/analysis.h5"),"r") do h
        for name in keys(h)
            a=attributes(h[name*"/400"]);b=attributes(h[name*"/200"])
            for key in ("raw_rms_V","average_rms_V","difference_rms_V")
                @test abs(read(a[key])/read(b[key])-1)<2e-6
            end
            @test 0 < read(a["averaged_noise_correlation"]) < 1
        end
        a=attributes(h["v3/400"])
        @test read(a["ltspice_measured_rms_V"]) ≈ read(a["raw_rms_V"]) rtol=1e-5
    end
end
