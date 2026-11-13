using Test
include("../../scripts/production/run_framed_detector_response.jl")
@testset "Electron display rebin" begin
    ch=ReadoutChannels(Int32[0,1],Int32[0,0],Float32[2,6],Float32[2,2])
    counts=reshape(Int64.(1:46),23,2)
    v=DetectorVectors(ch,ElectronReadoutTime(1000f0,.25f0,Int32(23)),counts)
    r=(placement=EventPlacement(1,75f0,75f0),excluded=SparseCells(Int64[],Int64[],Int64[],Int64[]))
    g=ReadoutGrid(4f0,(0f0,0f0))
    raw=projections(r,v,g,.6f0,false)
    grouped=projections(r,v,g,.6f0,false;rebin=20)
    @test grouped["xy_count"]==raw["xy_count"]
    @test sum(grouped["xz_count"])==sum(raw["xz_count"])==sum(counts)
    @test length(grouped["xz_count"])==4
    @test grouped["xz_count"][1]==sum(counts[1:20,1])
    @test grouped["xz_count"][2]==sum(counts[21:23,1])
    @test grouped["display_z_bin_mm"]≈3f0
    @test grouped["xz_y"][1]>grouped["xz_y"][2]
    @test v.counts==reshape(Int64.(1:46),23,2)
    @test v.time.bucket_us==.25f0
end
