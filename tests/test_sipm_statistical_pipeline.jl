using Test,HDF5,Random
include("../scripts/production/sipm_production.jl")
const Prod=SiPMProduction
const Pilot=Prod.Pilot
const EA=Prod.EA
const ROOT=normpath(joinpath(@__DIR__,".."))
@testset "Statistical card and PE-bucket stage integration" begin
    pilot=Pilot.settings(ROOT,"config/specs/runs/sipm_response_statistical_pilot.toml")
    @test pilot.c["response_model"]=="statistical"
    @test pilot.c["statistical_substeps"]==128
    reference=Pilot.settings(ROOT,"config/specs/runs/sipm_response_ep3.toml")
    @test pilot.sensor_model!=reference.sensor_model
    @test pilot.electronics_model==reference.electronics_model
    counts=Int64[0,1000,50000,0,100]
    dispatched=Pilot.simulate_sensor_channel(Xoshiro(33),counts,250e-9,pilot)
    direct=Prod.P.simulate_channel_statistical(Xoshiro(33),counts,250e-9,pilot.sensor;
        options=Prod.P.StatisticalOptions(128,16))
    @test dispatched==direct
    grid=EA.GridSpec(10f0,.25f0,(0f0,0f0,0f0),(0f0,0f0),Int8(-1))
    cells=EA.SparseCells(Int64[-1,-1,2],Int64[3,3,1],Int64[-2,2,0],Int64[1000,50000,100])
    history=EA.PhotonSiPMHistory(EA.EventPlacement(Int64(0),50f0,100f0),:pure_xe,grid,cells,
        Int64(100),Int64(90),Int64(100000))
    mkpath(joinpath(ROOT,"work/sipm_tests"))
    mktempdir(joinpath(ROOT,"work/sipm_tests")) do dir
        input=joinpath(dir,"photons.h5")
        h5open(input,"w") do f
            g=create_group(f,"events/0");EA.write_photon_history(g,history)
            g["pair_checksum"]="fixture-pairs"
        end
        src=(;path=input,sourcehash=Prod.sha(input),sample="fixture",gas="pure_xe",dt=250e-9,stream="fixture-stream")
        branch=(;src,sensor_dir=joinpath(dir,"sensor"),asic_dir=joinpath(dir,"asic"),sensorid="sensor-fixture",asicid="asic-fixture")
        settings=(;pilot,codehash="fixture-code")
        sensor=Prod.sensor_shard(settings,branch,Int64[0],1)
        asic=Prod.asic_shard(settings,branch,sensor,Int64[0],1)
        savedhash=Prod.sha(sensor)
        @test Prod.sensor_shard(settings,branch,Int64[0],1)==sensor
        @test Prod.sha(sensor)==savedhash
        h5open(sensor,"r") do s
            h5open(asic,"r") do a
                sg=s["events/0"];ag=a["events/0"]
                @test read(s["metadata/complete"]) && read(a["metadata/complete"])
                @test read(s["metadata/signal_unit"])=="PE"
                @test read(sg["response_model"])=="statistical"
                @test read(sg["incident_photons_total"])==sum(cells.count)
                @test read(sg["pair_checksum"])==read(ag["pair_checksum"])=="fixture-pairs"
                @test read(sg["bucket"])==read(ag["bucket"])
                @test minimum(read(sg["bucket"]))==-6
                @test read(sg["pixel_ix"])==read(ag["pixel_ix"])
                @test read(sg["pixel_iy"])==read(ag["pixel_iy"])
                @test read(a["metadata/source_sha256"])==savedhash
                @test all(isfinite,read(ag["readout_pe_equiv"]))
                @test sum(read(sg["sensor_pe_equiv"]))≈read(sg["sensor_pe_total"])
                @test read(sg["statistical_state_repairs"])==0
                @test read(sg["statistical_variance_limits"])==0
            end
        end
        # New ASIC noise card reuses the complete sensor state and charge history.
        altpilot=merge(pilot,(;electronics=Prod.E.Electronics(250e-9,1.,1e-3),electronics_model="alternate-noise"))
        altsettings=merge(settings,(;pilot=altpilot))
        altbranch=merge(branch,(;asic_dir=joinpath(dir,"alternate_asic"),asicid="alternate-asic"))
        Prod.asic_shard(altsettings,altbranch,sensor,Int64[0],1)
        @test Prod.sha(sensor)==savedhash
    end
end
