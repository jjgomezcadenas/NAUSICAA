using Test
include(joinpath(@__DIR__,"../src/ion_signal/pixel_area.jl"))
using .PixelArea
root=normpath(joinpath(@__DIR__,".."))
@testset "Area inputs and physical geometry" begin
    devices=read_devices(joinpath(root,"LTSpice/variants/m4_x4_m3_x4/noise.cir"))
    @test length(devices)==9
    @test count(d->d.instance=="M2",devices)==3
    @test only(filter(d->d.name=="M5",devices)).W_m≈220e-6
    @test only(filter(d->d.name=="M5",devices)).L_m≈.5e-6
    @test sum(d.W_m*d.L_m for d in devices)≈(110+3*13.5+400+4.1472+110+6+36)*1e-12
    for margin in (5e-6,10e-6,20e-6)
        boxes,height=pack_devices(devices;width=200e-6,margin=margin,gap=margin,outer=2margin,passive_side=50e-6)
        @test validate_boxes(boxes,200e-6,height)
        @test height<1e-3
        @test 200e-6*height>=sum(b.w*b.h for b in boxes)
        @test sum(b.w*b.h for b in boxes)>sum(d.W_m*d.L_m for d in devices)
        @test Set(b.name for b in boxes)==union(Set(d.name for d in devices),Set(["Cf reserve"]))
    end
    @test_throws ErrorException pack_devices(devices;width=10e-6,margin=20e-6,gap=10e-6,outer=0.,passive_side=50e-6)
    collision=[(name="a",x=0.,y=0.,w=1.,h=1.),(name="b",x=.5,y=.5,w=1.,h=1.)]
    @test_throws AssertionError validate_boxes(collision,2.,2.)
    @test_throws AssertionError validate_boxes(collision[1:1],.5,.5)
end
