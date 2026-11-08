#!/usr/bin/env julia
using HDF5,SHA
module PhotonPilot
include("run_electron_250ns.jl")
end
module ResponsePilot
include("benchmark_sipm_production.jl")
end
include(joinpath(@__DIR__,"../../src/ion_signal/guarded_files.jl"))
const ROOT=PhotonPilot.ROOT

"""Archive pilot photon files beside their parent realisation and save the response input index.

Photon contents and identities are preserved. The copied index points to the
archived files and records the original pilot index as provenance.
"""
function archive_photon_pilot(source_index::AbstractString)
    target=joinpath(ROOT,"results/simulations/production_pilot/electron_250ns_ep2_pilot100.h5")
    GuardedFiles.checked_copy(source_index,target*".tmp";force=true)
    h5open(target*".tmp","r+") do f
        g=f["bb0nu"]
        base=dirname(String(read(g["realisation_path"])))
        dir=joinpath(base,"detector","photons",String(read(g["electron_id"])));mkpath(dir)
        for gas in ("pure_xe","xe_he15")
            key="$(gas)_path";source=String(read(g[key]))
            dest=joinpath(dir,basename(source))
            if isfile(dest)
                sha256(read(dest))==sha256(read(source)) || error("Archived photon contents differ")
            else
                GuardedFiles.checked_copy(source,dest)
            end
            delete_object(g,key);g[key]=dest
        end
        f["metadata/pilot_cache_index"]=String(source_index)
    end
    mv(target*".tmp",target;force=true)
    target
end

"""Run 100 bb0nu events at E/P=2, J=0.030 through photons, statistical SiPM and ASIC.

Both gas prescriptions use the existing common primary realisation and 250 ns
buckets. Output cards bind the response to the archived 100-event photon index.
"""
function run_ep2_pilot()
    index=PhotonPilot.main_electrons("config/specs/runs/electron_250ns_ep2.toml",100)
    archive_photon_pilot(index)
    ResponsePilot.benchmark_response(tags=("ep2",),
        output_dir=joinpath(ROOT,"results/sipm_statistical_production_pilot/ep2"))
end
if abspath(PROGRAM_FILE)==@__FILE__;run_ep2_pilot();end
