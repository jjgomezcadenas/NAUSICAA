include(joinpath(@__DIR__,"../../../src/ion_signal/bucket_study.jl"))
using HDF5,TOML,SHA
root=normpath(joinpath(@__DIR__,"../../.."))
configpath=joinpath(root,"config/bucket_study.toml")
cfg=TOML.parsefile(configpath)
for frame in cfg["bucket_sizes_s"]
    dir=joinpath(root,cfg["output"],"$(round(Int,frame*1000))ms")
    h5open(joinpath(dir,"raw.h5"),"r") do h
        @assert read(attributes(h)["config_sha256"])==bytes2hex(sha256(read(configpath))) "Configuration differs from saved acquisition"
        @assert read(attributes(h)["frame_s"])==frame
        @assert read(attributes(h)["input_sha256"])==bytes2hex(sha256(read(joinpath(dir,"input.h5"))))
    end
end
BucketStudy.analyze(root)
