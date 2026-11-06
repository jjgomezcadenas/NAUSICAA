using Test,HDF5
module IonReplay
include("../scripts/production/run_ion_readout.jl")
end
module ElectronReplay
include("../scripts/production/run_electron_250ns.jl")
end
@testset "Both arms require the shared primary implementation" begin
    mktempdir() do dir
        path=joinpath(dir,"ledger.h5")
        h5open(path,"w") do f
            attributes(f)["schema"]="nausicaa.event_realisation.v2"
            m=create_group(f,"metadata")
            for (k,v) in pairs((;complete=true,event_count=2,source_sha256="source",card_hash="card",
                stream_identity="stream",seed=42,primary_code_hash=IonReplay.primary_code_hash(IonReplay.ROOT)))
                m[string(k)]=v
            end
        end
        paths=(;source_sha="source",realisation_id="card",stream_identity="stream")
        checks=(IonReplay.check_realisation,ElectronReplay.check_common_realisation)
        for check in checks
            @test h5open(f->check(f,paths,2,42),path,"r")===nothing
        end
        h5open(path,"r+") do f
            delete_object(f["metadata"],"primary_code_hash");f["metadata/primary_code_hash"]="wrong"
        end
        for check in checks
            @test_throws ErrorException h5open(f->check(f,paths,2,42),path,"r")
        end
        h5open(path,"r+") do f;delete_object(f["metadata"],"primary_code_hash");end
        for check in checks
            @test_throws ErrorException h5open(f->check(f,paths,2,42),path,"r")
        end
    end
end
