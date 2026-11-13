using HDF5,TOML,SHA,LinearAlgebra
const START=time()
include("fast_trigger_ion_emin_functions.jl")
root=normpath(joinpath(@__DIR__,"../.."));cardpath=joinpath(@__DIR__,"fast_trigger_ion_emin.toml")
c=TOML.parsefile(cardpath);r=TOML.parsefile(joinpath(root,c["realisation_card"]))
m=TOML.parsefile(joinpath(root,c["production_manifest"]));source=m["ion_response"]
key="covariance_by_frames/$(c["covariance_frames"])/covariance_e2"
cv,metadata=h5open(source,"r") do f
    @assert read(f["metadata/complete"])==true
    @assert read(f["metadata/pitch_mm"])==c["expected_pitch_mm"]
    @assert read(f["metadata/frame_s"])==c["expected_frame_s"]
    (read(f[key]),Dict(k=>read(f["metadata/"*k]) for k in ["model_id","response_card_toml","reconstruction_card_toml"]))
end
@assert issymmetric(cv) && isposdef(Symmetric(cv))
rows=NamedTuple[]
for E in c["energies_keV"], K in c["summed_pad_buckets"], arrangement in ["independent_pads","consecutive_buckets"]
    v=ion_emin_moments(E,r["w_eV"],r["fano_factor"],cv,K;arrangement)
    push!(rows,merge((energy_keV=E,summed_cells=K,arrangement=arrangement),v))
end
out=joinpath(@__DIR__,"results/fast_trigger_ion_emin.h5")
h5open(out,"w") do f
    f["schema"]="itaca.fast_trigger_ion_emin.v1";f["card"]=read(cardpath,String)
    f["source_response"]=source;f["source_covariance_dataset"]=key;f["covariance_e2"]=cv
    f["covariance_sha256"]=bytes2hex(sha256(reinterpret(UInt8,vec(cv))))
    for (k,v) in metadata;f[k]=v;end
    f["realisation_card"]=read(joinpath(root,c["realisation_card"]),String)
    for k in keys(first(rows));f[string(k)]=[getproperty(v,k) for v in rows];end
    f["frame_s"]=c["expected_frame_s"];f["pitch_mm"]=c["expected_pitch_mm"]
    f["single_bucket_enc_ions"]=sqrt(cv[cld(size(cv,1),2),cld(size(cv,1),2)])
    f["elapsed_seconds"]=time()-START
end
println("Saved ",out," in ",round(time()-START;digits=3)," s")
