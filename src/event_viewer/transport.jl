# Application adapter around the tested event-analysis transport functions.
using Random

"""Transport one event with deterministic file/event/realisation random streams.

Drift is uniform within the configured cm limits, or the supplied fixed cm value.
A SHA-256-derived seed identifies each canonical input path, event ID, base seed
and realisation. Separate streams keep the random drift fixed when transport
settings are edited. Reuse ion_counts and the Petit Gaussian ion sampler with
thermal longitudinal/transverse widths. Return the occupied ion table, explicit
axes, exact charge ledger, evaluated parameters and seed strings. Inputs retain
source-frame mm/MeV values. Detector fields are read in SI and converted explicitly.
"""
function transport_event(hits,source,mode,fixed_cm,realisation)
    c=CONFIG["transport"];det=TOML.parsefile(joinpath(ROOT,c["detector_card"]))
    mode in ("random","fixed") || error("Drift mode must be random or fixed")
    isfinite(fixed_cm) && c["min_cm"]<=fixed_cm<=c["max_cm"] || error("Fixed drift must be within the configured range")
    realisation>=0 || error("Realisation must be nonnegative")
    token=bytes2hex(sha256(join((realpath(source),string(hits.event_id),string(c["seed"]),string(realisation)),"|")))
    seed=parse(UInt64,token[1:16];base=16)
    drift=mode=="fixed" ? Float32(fixed_cm) : sample_drift_cm(Xoshiro(seed),Float32(c["min_cm"]),Float32(c["max_cm"]))
    sigma=sigma_t_ion_mm(Float32(det["gas"]["temperature_K"]),drift,Float32(det["bias"]["bulk_field_V_m"])/100f0)
    w=Float32(c["ionisation_energy_eV"]);n=ion_counts(hits,w)
    d=diffuse_xyz_image_mc(Xoshiro(seed⊻0x9e3779b97f4a7c15),hits,n.sampled;
        sigma_t_mm=sigma,sigma_l_mm=sigma,nbins=Int(c["nbins"]),nsigma=Float32(c["nsigma"]),ionisation_energy_eV=w)
    return (;d,drift_cm=drift,sigma_t_mm=sigma,sigma_l_mm=sigma,mode,fixed_cm=Float32(fixed_cm),realisation,
        seed=string(seed),ionisation_energy_eV=w,temperature_K=Float32(det["gas"]["temperature_K"]),
        field_V_cm=Float32(det["bias"]["bulk_field_V_m"])/100f0,
        detector_toml=read(joinpath(ROOT,c["detector_card"]),String))
end

"""Save transported ions and exact Int64 display projections without reconstruction cuts.

The displayed cells are the occupied intermediate diffusion bins. Their centres,
widths and charge counts are saved explicitly; each ion contributes to every
projection once. HDF5 projection axes follow the MC display convention.
"""
function save_transport(f,t)
    base="transport/";d=t.d
    for key in (:drift_cm,:sigma_t_mm,:sigma_l_mm,:mode,:fixed_cm,:realisation,:seed,:ionisation_energy_eV,:temperature_K,:field_V_cm,:detector_toml)
        f[base*string(key)]=getproperty(t,key)
    end
    for key in (:sampled_ions,:retained_ions,:outside_ions)
        f[base*string(key)]=getproperty(d,key)
    end
    save_columns(create_group(f,base*"bins"),d.bins)
    xe,ye,ze=d.x_edges_mm,d.y_edges_mm,d.z_edges_mm
    for (a,edges) in zip(("x","y","z"),(xe,ye,ze));f[base*"display/axes/"*a*"_edges_mm"]=edges;end
    xy=zeros(Int64,length(xe)-1,length(ye)-1);xz=zeros(Int64,length(xe)-1,length(ze)-1);yz=zeros(Int64,length(ye)-1,length(ze)-1)
    occupied=EventAnalysis.DataFrame(ix=Int64[],iy=Int64[],iz=Int64[],ions=Int64[],hit_count=Int64[])
    for r in eachrow(d.bins)
        i=EventAnalysis._histogram_index(r.x,xe);j=EventAnalysis._histogram_index(r.y,ye);k=EventAnalysis._histogram_index(r.z,ze)
        xy[i,j]+=r.ions;xz[i,k]+=r.ions;yz[j,k]+=r.ions
        push!(occupied,(i,j,k,r.ions,r.ions))
    end
    for (name,values) in (("xy",xy),("xz",xz),("yz",yz))
        f[base*"display/projections/"*name*"_ions"]=values
        f[base*"display/projections/"*name*"_hit_count"]=values
    end
    save_columns(create_group(f,base*"display/occupied"),occupied)
end
