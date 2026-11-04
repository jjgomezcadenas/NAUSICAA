using Test
include("../scripts/production/reduce_electron_operating_points.jl")
@testset "Energy correction and distribution" begin
    sensor=Dict("dark_reference_rate_Hz"=>0.,"diameter_m"=>.005,"dark_reference_area_m2"=>9e-6,"recovery_s"=>73.4e-9,"microcell_count"=>7854)
    electron=Dict("drift_speed_mm_us"=>1.,"lifetime_ms"=>15.)
    @test corrected_total(100.,10,250e-9,sensor,electron,0.) == (100.,0.)
    @test corrected_total(100.,10,250e-9,sensor,electron,150.)[1] ≈ 100exp(.1)
    sensor["dark_reference_rate_Hz"]=1e6
    corrected,dark=corrected_total(100.,10,250e-9,sensor,electron,0.)
    @test corrected+dark≈100
    @test 0<dark<10*1e6*(pi*.005^2/4)/9e-6*250e-9
    mktemp() do path,io
        close(io)
        h5open(path,"w") do f
            distribution(f,"test",[1.,2.,3.])
            @test sum(read(f["test/counts"]))==3
            @test read(f["test/rms"])≈sqrt(2/3)
        end
    end
end
