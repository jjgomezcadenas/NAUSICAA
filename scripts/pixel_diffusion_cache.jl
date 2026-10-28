# Reusable transport cache, shared by every collection-grid configuration.

"""Open a diffusion-only cache identified by source and physical transport inputs.

The fingerprint includes source content/path (the RNG identity), W, drift range,
thermal widths' temperature/field inputs, histogram resolution/padding, base seed,
Julia version and transport implementation hashes. Pixel pitch, frame size,
collection direction, histogram selection and requested event count are excluded.
The cache is appendable by event ID. Each event stores losslessly compressed
Float32 occupied diffusion-bin centres/energies, Int64 ions, explicit histogram
edges, evaluated drift/widths, seed and exact produced/retained/outside counts.
A completion marker is written last; incomplete event groups are rebuilt.
"""
function open_pixel_diffusion_cache(c,input,sample,det,root)
    tc=c["transport"]
    code=["scripts/pixel_diffusion_cache.jl","src/event_analysis/production.jl","src/event_analysis/diffusion.jl"]
    hashes=[bytes2hex(sha256(read(joinpath(root,p)))) for p in code]
    source_hash=open(io->bytes2hex(sha256(io)),input)
    physics=join(string.((source_hash,realpath(input),tc["ionisation_energy_eV"],tc["min_cm"],tc["max_cm"],
        tc["nbins"],tc["nsigma"],det["gas"]["temperature_K"],det["bias"]["bulk_field_V_m"],
        c["run"]["seed"],VERSION,"thermal-isotropic-rng-streams-v1",hashes...)),"|")
    fingerprint=bytes2hex(sha256(physics))
    directory=joinpath(root,get(c["run"],"diffusion_cache_directory","results/diffused_tracks"));mkpath(directory)
    path=joinpath(directory,sample*"_"*fingerprint[1:16]*".h5")
    cache=h5open(path,isfile(path) ? "r+" : "w")
    if haskey(cache,"metadata/fingerprint")
        read(cache["metadata/fingerprint"])==fingerprint || error("Diffusion cache fingerprint mismatch")
    else
        attributes(cache)["schema"]="nausicaa.diffused_tracks.v1"
        cache["metadata/fingerprint"]=fingerprint;cache["metadata/source_path"]=input
        cache["metadata/source_sha256"]=source_hash;cache["metadata/transport_inputs"]=physics
        cache["metadata/code_paths"]=code;cache["metadata/code_sha256"]=hashes
        cache["metadata/julia_version"]=string(VERSION)
        create_group(cache,"events")
    end
    return cache,path
end

"""Read and validate a completed cached cloud; return nothing for an uncached event.

Return the same owned diffusion NamedTuple needed by collection, plus saved drift,
width, seed and timing. Validate exact ion totals and preserve Float32 positions
and Int64 counts without revoxelisation or Gaussian draws.
"""
function read_cached_diffusion(cache,id)
    key="events/"*string(id)
    haskey(cache,key*"/complete") || return nothing
    g=cache[key]
    bins=DataFrame(:event_id=>fill(Int64(id),length(g["ions"])),
        (Symbol(k)=>read(g[k]) for k in ("x","y","z","energy","ions"))...)
    sampled=read(g["sampled_ions"]);retained=read(g["retained_ions"]);outside=read(g["outside_ions"])
    @assert sum(bins.ions)==retained && sampled==retained+outside
    d=(;bins,sampled_ions=sampled,retained_ions=retained,outside_ions=outside,
        x_edges_mm=read(g["x_edges_mm"]),y_edges_mm=read(g["y_edges_mm"]),z_edges_mm=read(g["z_edges_mm"]))
    return (;d,drift=read(g["drift_cm"]),sigma=read(g["sigma_mm"]),seed=parse(UInt64,read(g["seed"])))
end

"""Save a complete occupied diffusion cloud with lossless compression and a final marker."""
function save_cached_diffusion(cache,id,d,drift,sigma,seed)
    key="events/"*string(id)
    haskey(cache,key) && HDF5.delete_object(cache,key)
    g=create_group(cache,key)
    for k in ("x","y","z","energy","ions")
        values=d.bins[!,Symbol(k)]
        dataset=create_dataset(g,k,eltype(values),(length(values),);chunk=(max(1,min(length(values),4096)),),shuffle=true,deflate=3)
        write(dataset,values)
    end
    for k in (:sampled_ions,:retained_ions,:outside_ions,:x_edges_mm,:y_edges_mm,:z_edges_mm)
        g[string(k)]=getproperty(d,k)
    end
    g["drift_cm"]=drift;g["sigma_mm"]=sigma;g["seed"]=string(seed)
    g["complete"]=true
end
