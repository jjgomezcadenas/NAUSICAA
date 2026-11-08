using Test
include("../scripts/production/reduce_response_overview.jl")
@testset "Paired pulse support" begin
    @test response_metrics([0,1,100], [999,-5,90]) == (2.0,85.0,42.5,90.0)
    @test response_metrics([1,1], [-3,-7]) == (2.0,-10.0,-5.0,-3.0)
    @test response_metrics([1,3], [1,3]) == (2.0,4.0,2.0,3.0)
    @test isnan(response_metrics([0],[2])[3])
    @test_throws DimensionMismatch response_metrics([1],[1,2])
end

@testset "Retired overview preserves its saved reduction" begin
    mktempdir(joinpath(ROOT,"work")) do dir
        output=joinpath(dir,"response_overview.h5")
        write(output,"saved reduction fixture")
        card=joinpath(dir,"run.toml")
        open(card,"w") do io
            TOML.print(io,Dict("historical_sources_retired"=>true,"output"=>dir,
                "sample"=>"bb0nu","sipm_manifest"=>joinpath(dir,"missing.toml"),
                "ions"=>[Dict("path"=>joinpath(dir,"missing.h5"),"pitch_mm"=>3)]))
        end
        @test_throws ArgumentError main(card)
        @test read(output,String)=="saved reduction fixture"
    end
end
