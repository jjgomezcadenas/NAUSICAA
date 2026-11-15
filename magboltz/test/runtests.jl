using Test,HDF5
include("../MagboltzDriver.jl")
using .MagboltzDriver

@testset "Magboltz input and density" begin
    p=read_point(joinpath(@__DIR__,"..","pilot.toml"))
    @test p.gas==:Xe_85_He_15
    @test effective_pressure_bar(p)==17.5
    @test isapprox(reduced_field_Td(p),0.06864;rtol=1e-3)
    lines=split(chomp(input_card(p)),'\n')
    @test lines[1]=="2 2 0 1 0.0"
    @test lines[2]=="7 3 0 0 0 0"
    third=parse.(Float64,split(lines[3]))
    @test third[1:6]==[85,15,0,0,0,0]
    @test isapprox(third[7],16.85;atol=1e-12)
    @test isapprox(third[8],17.5*750.061683;rtol=1e-12)
    @test lines[4]=="300.0 0 0"
    @test lines[5]=="0 0 0 0 0"
    corrected=GasPoint(density_mode=:partial_number)
    @test isapprox(density_factor(corrected),1.0957786611209301;rtol=1e-12)
    @test isapprox(effective_pressure_bar(corrected),19.176126569616276;rtol=1e-12)
    @test isapprox(reduced_field_Td(corrected),0.06264;rtol=1e-3)
    geant4=GasPoint(density_mode=:geant4_material)
    @test isapprox(density_factor(geant4),1.11209763;rtol=1e-6)
    @test density_factor(geant4)>density_factor(corrected)
    pure=GasPoint(gas=:Xe_100)
    pure_lines=split(chomp(input_card(pure)),'\n')
    @test pure_lines[1:2]==["1 2 0 1 0.0","7 0 0 0 0 0"]
    @test parse.(Float64,split(pure_lines[3]))[1:6]==[100,0,0,0,0,0]
end

@testset "Magboltz output and diffusion units" begin
    output="""
      Z DRIFT VELOCITY = 0.6120E+00 MICRONS/NANOSECOND  +-    0.10%
      TRANSVERSE DIFFUSION   = 0.2500D+03 +-    0.80%
          =  0.01234 EV. +- 0.83%
          =  250.000 MICRONS/CENTIMETER**0.5  +-    0.60%
      LONGITUDINAL DIFFUSION = 0.6400D+03 +-    1.10%
          =  0.01345 EV. +- 1.12%
          =  400.000 MICRONS/CENTIMETER**0.5  +-    0.75%
      NUMBER OF COLLISIONS IN FINAL ENERGY BIN =    25.0
    """
    x=parse_output(output)
    @test x.velocity_mm_us==0.612
    @test x.transverse_diffusion_cm2_s==250
    @test x.longitudinal_um_sqrtcm==400
    @test x.final_energy_bin_collisions==25
    t=diffusion_units(x.transverse_um_sqrtcm,17.5)
    @test t.mm_sqrtm==2.5
    @test isapprox(t.card_mm_sqrtbar_sqrtcm,0.25sqrt(17.5))
    @test_throws ArgumentError parse_output("missing results")
    mktempdir() do dir
        input=joinpath(dir,"input.card");raw=joinpath(dir,"magboltz.out")
        err=joinpath(dir,"magboltz.err");source=joinpath(dir,"source.f")
        binary=joinpath(dir,"binary");table=joinpath(dir,"table.h5")
        write(input,"test input");write(raw,output);write(err,"")
        write(source,"source");write(binary,"binary")
        row=(point=GasPoint(),result=x,elapsed_s=1.0,input=input,output=raw,error=err)
        save_table(table,[row];source_path=source,binary_path=binary)
        h5open(table,"r") do f
            @test read(f["schema"],String)=="itaca.magboltz.transport.v1"
            @test read(f["results/velocity_mm_us"])[1]==0.612
            @test read(f["results/transverse_mm_sqrtm"])[1]==2.5
            @test read(f["collision_target"])[1]==80_000_000
            @test read(f["raw_outputs"])[1]==output
        end
    end
    mktempdir() do dir
        fake=joinpath(dir,"fake_magboltz.sh")
        write(fake,"#!/bin/sh\ncat >/dev/null\ncat <<'RESULT'\n"*output*"RESULT\n")
        chmod(fake,0o755)
        row=run_point(GasPoint(),fake,joinpath(dir,"run"))
        @test row.result.velocity_mm_us==0.612
        @test isfile(row.input)
        @test isfile(row.output)
    end
end
