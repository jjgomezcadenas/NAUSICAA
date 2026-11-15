# Main-track box and zoom cube of the event display (scripts/production/track_controls.jl).
using Test
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"));using .EventAnalysis
include(joinpath(ROOT,"scripts/production/track_controls.jl"))
voxels(x,y,z)=(voxels=EventAnalysis.DataFrame(x=Float32.(x),y=Float32.(y),z=Float32.(z),energy=ones(Float32,length(x))),)
@testset "main_track_box" begin
    big=voxels([0,3,6,9],[0,0,3,3],[10,10,13,16]);small=voxels([40,43],[40,40],[40,40])
    box=main_track_box([small,big],3)            # the component with most voxels, not the first
    @test box==Float32[-1.5 10.5;-1.5 4.5;8.5 17.5]
    @test main_track_box([small],3)==Float32[38.5 44.5;38.5 41.5;38.5 41.5]
    @test main_track_box([],3)===nothing
end
@testset "zoom_cube" begin
    frame=Float32[-80 80;-80 80;-80 80]
    a=Float32[-1.5 10.5;-1.5 4.5;8.5 17.5];b=Float32[0 20;0 2;0 2]
    z=zoom_cube([a,b],frame)
    side=z[:,2].-z[:,1]
    @test all(side.≈side[1])                       # a cube
    @test side[1]≈21.5+20                          # largest side of the union (x: -1.5..20) plus 10 mm each side
    for i in 1:3;@test z[i,1]<=min(a[i,1],b[i,1])&&z[i,2]>=max(a[i,2],b[i,2]);end
    @test all(z[:,1].>=frame[:,1]) && all(z[:,2].<=frame[:,2])
    # margin of 10 % when the track is large: 100 mm plus 10 mm each side; the cube is shifted to stay inside the frame
    c=Float32[-75 25;-10 10;-10 10];zc=zoom_cube([c],frame)
    @test all((zc[:,2].-zc[:,1]).≈120)
    @test zc[1,1]≈-80 && zc[1,2]≈40            # centred at -25 it would start at -85; shifted by 5 mm
    @test zc[2,:]≈Float32[-60,60] && zc[3,:]≈Float32[-60,60]
    # larger than the frame on an axis: the frame itself, still a cube
    d=Float32[-79 79;0 1;0 1];zd=zoom_cube([d],frame)
    @test zd==frame
    # larger than the frame and a box beyond it: a cube of the widened extent, centred on the frame, containing the box
    g=Float32[-83 79;0 1;0 1];zg=zoom_cube([g],frame)
    @test all((zg[:,2].-zg[:,1]).≈163) && zg[1,1]≈-83 && zg[1,2]≈80 && zg[2,:]≈Float32[-81.5,81.5]
    # a box reaching beyond the frame widens it: the box stays contained
    e=Float32[-83 -23;0 10;0 10];ze=zoom_cube([e],frame)
    @test ze[1,1]<=-83 && ze[1,2]>=-23 && all((ze[:,2].-ze[:,1]).≈80)
    # nothing entries are skipped; no boxes gives the frame; no frame gives no clipping
    @test zoom_cube([nothing,a],frame)==zoom_cube([a],frame)
    @test zoom_cube([nothing],frame)==frame
    znf=zoom_cube([a],nothing);@test all((znf[:,2].-znf[:,1]).≈(12+20))   # largest side 12 mm (x), margins 10 mm
    @test znf[1,:]≈Float32[-1.5-10,10.5+10]
end
