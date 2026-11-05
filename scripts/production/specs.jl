"""Script-side TOML resolution and deterministic production identities."""
module ProductionSpecs
using TOML, SHA

export load_specs, canonical, content_hash, event_seed, pair_checksum

"""Canonical representation of parsed TOML values, independent of key order/comments."""
function canonical(value)
    if value isa AbstractDict
        return "{"*join((repr(k)*":"*canonical(value[k]) for k in sort!(collect(keys(value)))),",")*"}"
    elseif value isa AbstractVector
        return "["*join(canonical.(value),",")*"]"
    elseif value isa AbstractString || value isa Number || value isa Bool
        return repr(value)
    end
    throw(ArgumentError("Unsupported TOML value $(typeof(value))"))
end

"""SHA-256 hex digest of a string or file, used for source and configuration IDs."""
content_hash(value::AbstractString;file::Bool=false)=
    bytes2hex(sha256(file ? read(value) : codeunits(value)))

"""Derive one UInt64 RNG seed from an identity, event ID and stage label."""
function event_seed(identity::AbstractString,event_id::Int64,label::AbstractString)::UInt64
    digest=content_hash(join((identity,string(event_id),label),"|"))
    return parse(UInt64,digest[1:16];base=16)
end

"""Hash ordered Int64 pair counts for branch-by-branch identity checks."""
pair_checksum(counts::Vector{Int64})=
    bytes2hex(sha256(Vector{UInt8}(reinterpret(UInt8,counts))))

"""Read the five pilot cards, check their identities and return resolved paths and values."""
function load_specs(root::AbstractString,run_card::AbstractString)
    run_path=realpath(isabspath(run_card) ? run_card : joinpath(root,run_card))
    run_spec=TOML.parsefile(run_path)
    get(run_spec,"schema",nothing)=="nausicaa.detector_run.v1" ||
        error("Run card schema mismatch")
    get(run_spec,"pilot_only",false)==true || error("Pilot runner requires pilot_only=true")
    selected=get(run_spec,"selected_samples",["bb0nu","single_electrons"])
    selected isa AbstractVector && !isempty(selected) &&
        all(name->name in ("bb0nu","single_electrons"),selected) && allunique(selected) ||
        error("Select one or both known samples without duplicates")
    config_path=realpath(joinpath(root,run_spec["detector_configuration"]))
    configuration=TOML.parsefile(config_path)
    get(configuration,"schema",nothing)=="nausicaa.detector_configuration.v1" ||
        error("Detector configuration schema mismatch")
    load(key,expected)=begin
        path=realpath(joinpath(root,configuration[key]))
        data=TOML.parsefile(path)
        get(data,"schema",nothing)==expected || error("$key card schema mismatch")
        (;path,data)
    end
    realisation=load("realisation","nausicaa.realisation_spec.v1")
    ions=load("ions","nausicaa.ion_detector_spec.v1")
    electrons=load("electrons","nausicaa.electron_detector_spec.v1")
    r=realisation.data;i=ions.data;e=electrons.data
    0<=r["min_ion_drift_cm"]<r["chamber_length_cm"] || error("Invalid drift range")
    0<=r["fano_factor"]<1 && r["w_eV"]>0 && r["max_events_per_sample"]>0 ||
        error("Invalid realisation card")
    all(>(0),(i["temperature_K"],i["drift_field_V_cm"],i["drift_speed_m_s"],
               i["pad_pitch_mm"],i["bucket_width_ms"])) || error("Invalid ion card")
    all(>(0),(e["pressure_bar"],e["drift_speed_mm_us"],e["lifetime_ms"],
               e["hole_pitch_mm"],e["hole_diameter_mm"],e["hole_depth_mm"],
               e["bucket_width_us"])) || error("Invalid electron card")
    0<=e["helium_fraction"]<1 && 1.5<=e["reduced_field_kV_cm_bar"]<=3 &&
        0<=e["j_variance_ratio"]<1 || error("Invalid gas/EL parameters")
    sensor_size = get(e,"sipm_shape","circle")=="square" ? e["sipm_side_mm"] : e["sipm_diameter_mm"]
    sensor_size>=e["hole_diameter_mm"] || error("SiPM must cover hole")
    j_anchors=Dict(1.5=>0.1,2.0=>0.05,3.0=>0.01)
    j_status=get(e,"j_status",nothing)
    j_status in ("oliveira_fig6_user_reading","user_specified") ||
        error("Pilot J must identify its source")
    if j_status=="oliveira_fig6_user_reading"
        expected=get(j_anchors,e["reduced_field_kV_cm_bar"],nothing)
        expected!==nothing && isapprox(e["j_variance_ratio"],expected;atol=1e-9) ||
            error("Figure 6 J reading applies only at the recorded E/P anchor")
    end
    Set(keys(r["samples"]))==Set(("bb0nu","single_electrons")) ||
        error("Pilot requires both sample paths")
    paths=Dict(name=>realpath(joinpath(root,path)) for (name,path) in r["samples"])
    return (;run_spec,run_path,configuration,config_path,realisation,ions,electrons,
            sample_paths=paths,results_dir=joinpath(root,run_spec["results_dir"]),
            report_dir=joinpath(root,run_spec["report_dir"]))
end
end
