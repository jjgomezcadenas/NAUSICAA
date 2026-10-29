# Reproducible numerical checkpoint: Julia calculations -> HDF5 -> Python report.
using HDF5, DataFrames, TOML, SHA, Dates, Random, Statistics, LinearAlgebra
include(joinpath(@__DIR__,"../../src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis
const ROOT = normpath(joinpath(@__DIR__,"../.."))
const CARD = isempty(ARGS) ? joinpath(ROOT,"config/event_analysis.toml") : abspath(ARGS[1])
const CONFIG = TOML.parsefile(CARD)
const DETECTOR = TOML.parsefile(joinpath(ROOT,CONFIG["inputs"]["detector_card"]))
const OUT = joinpath(ROOT,"results/event_analysis/step2_validation.h5")

"""Save a Julia dataframe as typed HDF5 columns with the original column order."""
function write_table(parent,name,df)
    g = create_group(parent,name)
    attributes(g)["columns"] = names(df)
    for column in names(df); g[column] = df[!,column]; end
    return g
end

"""Independent Float64 moments of retained histogram centres for validation."""
function moments(d)
    xyz = Float64.(Matrix(d[:,[:x,:y,:z]]))
    weights = Float64.(d.ions)/sum(d.ions)
    mean_xyz = vec(sum(xyz.*weights;dims=1))
    centered = xyz.-mean_xyz'
    covariance = centered'*(centered.*weights)
    return mean_xyz,covariance
end

"""One point used to test the diffusion kernel against an analytic Gaussian."""
point_hits() = EventHits(0,source_frame,Float32[0],Float32[0],Float32[0],Float32[1],
                        Float32[0],Int32[7],Int32[0])

function main()
    tests = TOML.parsefile(joinpath(ROOT,"work/event_analysis/step2_validation.toml"))
    tests["status"] == "passed" || error("Run the step-two tests first")
    for (path,hash) in tests["files"]
        bytes2hex(sha256(read(joinpath(ROOT,path)))) == hash || error("Changed after testing: $path")
    end
    w = Float32(CONFIG["ionization"]["ionisation_energy_eV"])
    temperature = Float32(DETECTOR["gas"]["temperature_K"])
    field = Float32(DETECTOR["bias"]["bulk_field_V_m"])/100f0
    nbins = Int(CONFIG["diffusion"]["nbins"])
    padding = Float32(CONFIG["diffusion"]["nsigma"])
    longitudinal = CONFIG["diffusion"]["longitudinal"]
    longitudinal in ("thermal","zero") || error("Unknown longitudinal model")
    seed = UInt64(CONFIG["random"]["seed"])
    n = Int64(CONFIG["validation"]["synthetic_ions"])
    n > 0 || error("synthetic_ions must be positive")
    h = point_hits()
    kw = (;sigma_t_mm=1f0,sigma_l_mm=1f0,ionisation_energy_eV=w)
    # Compile before timing; the histogram still allocates on each measured call.
    diffuse_xyz_image_mc(Xoshiro(seed),h,Int64[100];kw...,nbins=nbins,nsigma=padding)
    mkpath(dirname(OUT))
    h5open(OUT,"w") do f
        attributes(f)["schema"] = "nausicaa.event.step2.validation.v1"
        attributes(f)["created_utc"] = string(now(UTC))
        attributes(f)["julia_version"] = string(VERSION)
        attributes(f)["scope"] = "Ion production, thermal diffusion and histogram convergence; source-frame hits; common L per event; per-hit rounding."
        f["configuration/card"] = read(CARD,String)
        f["configuration/detector"] = read(joinpath(ROOT,CONFIG["inputs"]["detector_card"]),String)
        f["configuration/ionisation_energy_eV"] = w
        f["configuration/temperature_K"] = temperature
        f["configuration/field_V_cm"] = field
        f["configuration/nbins"] = nbins
        f["configuration/nsigma"] = padding
        f["tests/passed"] = tests["counts"]["passed"]
        paths = sort(collect(keys(tests["files"])))
        append!(paths,[relpath(CARD,ROOT),"scripts/validation/validate_event_step2.jl"])
        f["provenance/paths"] = paths
        f["provenance/sha256"] = [bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in paths]

        # Independent analytic box retention from the standard normal CDF.
        rows = NamedTuple[]
        for bins in sort(unique(vcat(Int.(CONFIG["validation"]["bin_counts"]),nbins))), ns in sort(unique(vcat(Float32.(CONFIG["validation"]["padding_sigma"]),padding)))
            measured = @timed diffuse_xyz_image_mc(Xoshiro(seed),h,[n];kw...,nbins=bins,nsigma=ns)
            r = measured.value
            mu,covariance = moments(r.bins)
            erf_value = ccall((:erf,Base.Math.libm),Float64,(Float64,),Float64(ns)/sqrt(2.0))
            push!(rows,(nbins=bins,nsigma=ns,seed=seed,sampled=n,retained=r.retained_ions,
                outside=r.outside_ions,expected_outside_fraction=1-erf_value^3,
                outside_fraction=r.outside_ions/n,
                expected_variance_mm2=1-2Float64(ns)*exp(-Float64(ns)^2/2)/(sqrt(2pi)*erf_value),
                mean_x_mm=mu[1],variance_x_mm2=covariance[1,1],cov_xy_mm2=covariance[1,2],
                seconds=measured.time,allocated_bytes=measured.bytes))
            @assert r.retained_ions+r.outside_ions == n
            if bins == nbins && ns == padding
                g = create_group(f,"gaussian_example")
                g["x_edges_mm"] = r.x_edges_mm
                centres = r.x_edges_mm[1:end-1]./2f0 .+ r.x_edges_mm[2:end]./2f0
                marginal = zeros(Int64,length(centres))
                for row in eachrow(r.bins)
                    marginal[EventAnalysis._histogram_index(row.x,r.x_edges_mm)] += row.ions
                end
                g["x_centres_mm"] = centres
                g["observed_density_per_mm"] = Float64.(marginal)./(r.retained_ions.*Float64.(diff(r.x_edges_mm)))
                g["expected_density_per_mm"] = exp.(-Float64.(centres).^2 ./ 2)./(sqrt(2pi)*erf_value)
            end
        end
        write_table(f,"synthetic_scan",DataFrame(rows))

        for (sample_index,relative_input) in enumerate(CONFIG["inputs"]["samples"])
            input = joinpath(ROOT,relative_input)
            label = replace(basename(input),"_hitdf.h5"=>"")
            h5open(input,"r") do src
                g = create_group(f,"samples/"*label)
                g["input_sha256"] = bytes2hex(open(sha256,input))
                md = EventAnalysis.EventIO.read_event_metadata(src)
                allhits = EventAnalysis.EventIO.read_hitdf(src)
                totals = NamedTuple[]
                for row in eachrow(md)
                    range = row.first_row:row.first_row+row.n_hits-1
                    event = EventHits(row.event_id,source_frame,allhits.x_mm[range],allhits.y_mm[range],
                        allhits.z_mm[range],allhits.energy_MeV[range],allhits.time_ns[range],
                        allhits.particle_id[range],allhits.hit_id[range])
                    counts = ion_counts(event,w)
                    # Independent reference quantifies rounding versus Float32 arithmetic.
                    reference = sum(Float64,event.energy_MeV)*1e6/Float64(w)
                    push!(totals,(event_id=row.event_id,energy_MeV=counts.event_energy_MeV,
                        expected=counts.total_expected,sampled=sum(counts.sampled),
                        rounding_difference=counts.rounding_difference,
                        relative_rounding_difference=counts.total_expected > 0 ? counts.rounding_difference/counts.total_expected : 0f0,
                        relative_expectation_error=reference > 0 ? (counts.total_expected-reference)/reference : 0.0))
                end
                write_table(g,"production",DataFrame(totals))
                nread = min(Int(CONFIG["validation"]["events_per_sample"]),nrow(md))
                nread > 0 || error("events_per_sample must be positive")
                indices = unique(round.(Int,range(1,nrow(md);length=nread)))
                drift_rng = Xoshiro(seed+UInt64(sample_index))
                cases = NamedTuple[]
                for (case_index,index) in enumerate(indices)
                    hits,info = read_event(src,md.event_id[index])
                    counts = ion_counts(hits,w)
                    sampled_L = sample_drift_cm(drift_rng,Float32(CONFIG["drift"]["min_cm"]),Float32(CONFIG["drift"]["max_cm"]))
                    lengths = vcat(Float32.(CONFIG["validation"]["drift_cm"]),sampled_L)
                    event_seed = seed+UInt64(1000sample_index+case_index)
                    for L in lengths
                        sigma = sigma_t_ion_mm(temperature,L,field)
                        sigma_z = longitudinal == "thermal" ? sigma : 0f0
                        timed = @timed diffuse_xyz_image_mc(Xoshiro(event_seed),hits,counts.sampled;
                            sigma_t_mm=sigma,sigma_l_mm=sigma_z,nbins,nsigma=padding,ionisation_energy_eV=w)
                        r = timed.value
                        @assert r.sampled_ions == r.retained_ions+r.outside_ions == sum(counts.sampled)
                        @assert sum(r.bins.ions) == r.retained_ions
                        push!(cases,(event_id=hits.event_id,L_cm=L,sigma_mm=sigma,seed=event_seed,
                            sampled_ions=r.sampled_ions,retained_ions=r.retained_ions,outside_ions=r.outside_ions,
                            occupied_bins=nrow(r.bins),seconds=timed.time,allocated_bytes=timed.bytes))
                        if case_index == 1 && L == first(lengths)
                            ex = write_table(g,"example_diffusion",r.bins)
                            ex["x_edges_mm"] = r.x_edges_mm;ex["y_edges_mm"] = r.y_edges_mm;ex["z_edges_mm"] = r.z_edges_mm
                            attributes(ex)["L_cm"] = L
                        end
                    end
                end
                write_table(g,"diffusion_cases",DataFrame(cases))
                t = DataFrame(totals); c = DataFrame(cases)
                g["summary/events"] = nrow(t)
                g["summary/mean_rounding_difference"] = mean(t.rounding_difference)
                g["summary/max_absolute_relative_rounding"] = maximum(abs,t.relative_rounding_difference)
                g["summary/max_absolute_relative_expectation_error"] = maximum(abs,t.relative_expectation_error)
                g["summary/mean_seconds"] = mean(c.seconds)
                g["summary/max_outside_fraction"] = maximum(c.outside_ions./c.sampled_ions)
                println(label,": production ",nrow(t)," events; diffusion ",nrow(c)," cases; mean ",round(mean(c.seconds)*1000;digits=1)," ms")
            end
        end
    end
    println(OUT)
end
main()
