using Test,TOML,HDF5
include("../scripts/production/reduce_elm_sipm_pilot.jl")
const R=ELMSiPMReduction
@testset "Streaming reduction across split shard entries" begin
    original=joinpath(R.ROOT,"results/elm_sipm_pitch4mm_pilot")
    mkpath(joinpath(R.ROOT,"work/elm_reduction_tests"))
    mktempdir(joinpath(R.ROOT,"work/elm_reduction_tests")) do dir
        files=TOML.parsefile(joinpath(original,"files.toml"))
        for (tag,response) in files["responses"]
            manifest=TOML.parsefile(response["manifest"])
            shards=Dict{String,Any}[]
            for entry in manifest["shards"]
                for ids in Iterators.partition(entry["event_ids"],50)
                    partial=copy(entry);partial["event_ids"]=collect(ids);push!(shards,partial)
                end
            end
            manifest["shards"]=shards
            target=joinpath(dir,"manifest_$tag.toml")
            open(io->TOML.print(io,manifest),target,"w")
            response["manifest"]=target
            # Report regeneration must succeed after the original cards disappear.
            response["card"]=joinpath(dir,"missing-working-tree-card.toml")
        end
        open(io->TOML.print(io,files),joinpath(dir,"files.toml"),"w")
        R.main(output_dir=dir,report_name="streamed",store_inclusive=false)
        h5open(joinpath(dir,"streamed.h5"),"r") do fresh
          h5open(joinpath(original,"elm_sipm_pitch4mm_pilot.h5"),"r") do prior
            for tag in ("12p5","25")
                a=fresh["tau_$tag/xe_he15"];b=prior["tau_$tag/xe_he15"]
                @test length(read(a["sensor_paths"]))==2
                @test read(a["ledger_checks"])==100
                @test TOML.parse(String(read(a["realisation_card"]))) isa Dict
                @test !haskey(a,"inclusive_detected_pe")
                for key in ("detected_total","optical_charge_total","readout_total","event_mean_pe",
                    "event_max_pe","clip_loss","retained_fraction","corrected_readout_rms_percent",
                    "retained_readout_rms_percent","largest_detected_pulse","largest_sensor_pulse")
                    @test read(a[key])==read(b[key])
                end
                @test read(a["quiet_noise_mean_pe"])≈read(b["quiet_noise_mean_pe"])
                @test read(a["quiet_noise_rms_pe"])≈read(b["quiet_noise_rms_pe"])
            end
          end
        end
    end
end
@testset "Disjoint four-worker shard assignment" begin
    assigned=[[n for n in 1:80 if mod(n-1,4)==w-1] for w in 1:4]
    @test all(length(x)==20 for x in assigned)
    @test sort(vcat(assigned...))==collect(1:80)
end
