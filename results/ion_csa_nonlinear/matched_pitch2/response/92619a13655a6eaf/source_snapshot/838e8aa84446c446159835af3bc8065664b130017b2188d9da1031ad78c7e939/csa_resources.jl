"""Content checks for the physical resources used by CSA response calculations."""
module CSAResources
using HDF5, SHA
sha(path)=open(io->bytes2hex(sha256(io)),path)
resolve(root,path)=normpath(isabspath(path) ? path : joinpath(root,path))
require(ok,message)=ok || throw(ArgumentError(message))

"""
    verify(root, card)::Dict{String,String}

Verify the selected linear waveform, production deck and transistor model against
hashes recorded when the circuit calibration was acquired. The analysis pins the
waveform file; that file pins its acquisition deck and model. Return their actual
SHA-256 digests for response identities and saved provenance. Missing resources or
any content mismatch raise `ArgumentError` before cache reuse or simulation.
"""
function verify(root,card)
    base=resolve(root,card["calibration_directory"]);variant=card["variant"]
    paths=Dict("analysis"=>joinpath(base,"analysis.h5"),
        "kernel"=>joinpath(base,"variants",variant*".h5"),
        "deck"=>joinpath(base,"decks",variant,"30ms_single_N100_r2.cir"),
        "model"=>joinpath(base,"decks",variant,"corrected.inc"))
    for p in values(paths);require(isfile(p),"CSA resource is absent: $p");end
    hashes=Dict(k*"_sha256"=>sha(p) for (k,p) in paths)
    h5open(paths["analysis"],"r") do h
        a=attributes(h[variant])
        for (key,expected) in (("raw_sha256",hashes["kernel_sha256"]),("model_sha256",hashes["model_sha256"]))
            require(haskey(a,key) && String(read(a[key]))==expected,
                "CSA $key differs from acquired calibration; regenerate and validate calibration")
        end
    end
    h5open(paths["kernel"],"r") do h
        a=attributes(h["30ms/single_N100_r2"])
        require(haskey(a,"deck_sha256") && String(read(a["deck_sha256"]))==hashes["deck_sha256"],
            "CSA production deck differs from acquired calibration; regenerate and validate calibration")
        require(String(read(attributes(h)["model_sha256"]))==hashes["model_sha256"],
            "CSA waveform transistor model differs")
    end
    hashes
end
end
