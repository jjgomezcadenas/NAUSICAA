using Test
include("../../scripts/network/network_events_builder.jl")
@testset "Readout lattice inside the frame" begin
    @test lattice_cells(160,3)==55 && lattice_cells(160,4)==42 && lattice_cells(160,3.3)==50
    # Pads every 3 mm with the first centre 1 mm inside the frame, given out of order
    x=Float32[-79+3k for k in (5,0,53,2)]
    @test lattice_indices(x,-80,3.0,55)==UInt8[5,0,53,2]
    # A pad centred just below the lower face still belongs to cell 0
    @test lattice_indices(Float32[-81.4,-78.4],-80,3.0,55)==UInt8[0,1]
    # The upper edge: a centre half a pitch above the frame is the last cell
    @test lattice_indices(Float32[-81.4+3k for k in 0:54],-80,3.0,55)==UInt8.(0:54)
    @test_throws ErrorException lattice_indices(Float32[-79,-77.5],-80,3.0,55)      # off the lattice
    @test_throws ErrorException lattice_indices(Float32[-79,-79+3*60],-80,3.0,55)   # outside the grid
    # Non-integer pitch, as the ion frame of 3.3 mm
    z=[10.2+3.3k for k in 0:40];@test lattice_indices(z,10.0,3.3,50)==UInt8.(0:40)
    @test isempty(lattice_indices(Float32[],0,3.0,55))
    # A shorter last time group: one slice displaced by 0.4 pitch towards the others; the lattice is that of the majority
    zs=vcat(fill(10.2+0.4*2.842,3),[10.2+2.842k for k in 1:20]);@test lattice_indices(zs,10.0,2.842,58;tolerance=0.5)==vcat(fill(0x00,3),UInt8.(1:20))
    @test_throws ErrorException lattice_indices(zs,10.0,2.842,58)
end
@testset "Pitches from the resolved card" begin
    card=Dict("detector"=>Dict("ions"=>Dict("pad_pitch_mm"=>3.0,"ion_bucket_ms"=>20.0),"drift"=>Dict("ion_drift_velocity_cm_s"=>16.5,"electron_drift_velocity_mm_us"=>0.5684),
        "el"=>Dict("hole_pitch_mm"=>4.0),"sipm"=>Dict("electron_bucket_us"=>0.25)),"run"=>Dict("electrons"=>Dict("rebin_factor"=>20)))
    p=lattice_pitches(card);@test p.itrack[1]==3 && p.itrack[3]≈3.3 && p.etrack[1]==4 && p.etrack[3]≈2.842
end
@testset "ROI cut" begin
    # Strict bounds, compared in double precision as in the selection tables
    @test in_roi(Float32[2451.9,2452.0,2452.1,2463.5,2463.6,2463.7],(2452.0,2463.6))==[false,false,true,true,false,false]
    lower=2452.0387;e=Float32[prevfloat(Float32(lower)),nextfloat(Float32(lower))];@test in_roi(e,(lower,2463.6))==(Float64.(e).>lower)
    mktempdir() do dir
        path=joinpath(dir,"comparison.h5");h5open(f->(f["roi/lower_keV"]=2452.0;f["roi/upper_keV"]=2463.6),path,"w");@test roi_of(path)==(2452.0,2463.6)
    end
end
@testset "One hit per cell" begin
    ix=UInt8[1,2,1,3];iy=UInt8[0,0,0,0];iz=UInt8[5,5,5,5]
    @test one_hit_per_cell(ix,iy,iz,Float32[3.0,7.0,3.1,4.0])==[2,3,4]
    @test one_hit_per_cell(ix[1:2],iy[1:2],iz[1:2],Float32[1,2])==[1,2]
end
