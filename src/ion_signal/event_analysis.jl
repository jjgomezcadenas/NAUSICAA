module EventAnalysis

using HDF5, SHA, TOML
include("event_io.jl")
using .EventIO

export event_statistics, common_grid, bin_index, project_event, analyze_sample, select_hits,
       write_analysis, run_analysis

const ROOT = normpath(joinpath(@__DIR__, "../.."))
const SCHEMA = "nausicaa.ion_tracks.v1"
sha(path) = open(io -> bytes2hex(sha256(io)), path)

"""Contiguous event slices, preserving event IDs and the original row order."""
function event_ranges(ids)
    isempty(ids) && throw(ArgumentError("An empty sample has no ensemble mean"))
    starts = vcat(1, findall(diff(ids) .!= 0) .+ 1)
    stops = vcat(starts[2:end] .- 1, length(ids))
    event_ids = ids[starts]
    allunique(event_ids) || error("Input event rows must be contiguous")
    return event_ids, [a:b for (a,b) in zip(starts, stops)]
end

"""Expectations in ions, not integer realizations; input coordinates stay in mm."""
function event_statistics(df, vertices; w_eV=21.9)
    isfinite(w_eV) && w_eV > 0 || throw(ArgumentError("w_eV must be positive and finite"))
    all(isfinite, df.energy_MeV) && all(>=(0), df.energy_MeV) || error("Invalid hit energy")
    for name in (:x_mm, :y_mm, :z_mm)
        all(isfinite, df[!,name]) || error("Invalid hit coordinate")
    end
    ids, ranges = event_ranges(df.event_id)
    allunique(vertices.event_id) || error("Repeated generation vertices")
    vertex_rows = Dict(id=>j for (j,id) in enumerate(vertices.event_id))
    all(id -> haskey(vertex_rows,id),ids) || error("Missing generation vertex")
    n = length(ids)
    energy = zeros(n); cx = zeros(n); cy = zeros(n); cz = zeros(n)
    extent = zeros(3)
    nhits = length.(ranges)
    for (j, rows) in enumerate(ranges)
        energy[j] = sum(i -> Float64(df.energy_MeV[i]), rows)
        energy[j] > 0 || error("Event $(ids[j]) has zero energy")
        for (axis, column, origin_column, center) in zip(1:3, (:x_mm,:y_mm,:z_mm),
                (:vertex_x_mm,:vertex_y_mm,:vertex_z_mm), (cx,cy,cz))
            center[j] = Float64(vertices[vertex_rows[ids[j]],origin_column])
            isfinite(center[j]) || error("Invalid generation vertex")
            extent[axis] = max(extent[axis], maximum(i -> abs(Float64(df[i,column])-center[j]), rows))
        end
    end
    expected_ions = energy .* (1e6/w_eV)
    return (; ids, ranges, nhits, energy, expected_ions,
            ions_per_hit=expected_ions./nhits, cx, cy, cz, extent)
end

"""Square XY grid common to all samples; symmetric Z grid; retain all hits."""
function common_grid(stats; pitch_mm=2.0, z_bin_mm=2.0, margin=0.1, extra_bins=1)
    all(x -> isfinite(x) && x > 0, (pitch_mm,z_bin_mm)) || throw(ArgumentError("Invalid bin width"))
    isfinite(margin) && margin >= 0 && extra_bins >= 1 || throw(ArgumentError("Invalid grid margin"))
    isempty(stats) && throw(ArgumentError("At least one sample is required"))
    xy_extent = maximum(max(s.extent[1],s.extent[2]) for s in stats)
    z_extent = maximum(s.extent[3] for s in stats)
    hxy = ceil(Int, (1+margin)*xy_extent/pitch_mm)+extra_bins
    hz = ceil(Int, (1+margin)*z_extent/z_bin_mm)+extra_bins
    x = collect(-hxy:hxy) .* pitch_mm
    z = collect(-hz:hz) .* z_bin_mm
    return (; x_edges_mm=x, y_edges_mm=copy(x), z_edges_mm=z)
end

"""Left-closed/right-open bins, except the final bin includes its upper edge."""
function bin_index(value, edges)
    isfinite(value) || throw(DomainError(value, "Non-finite coordinate"))
    first(edges) <= value <= last(edges) || throw(DomainError(value, "Hit outside grid"))
    return min(searchsortedlast(edges, value), length(edges)-1)
end

function centered_event(df, stats, j, w_eV)
    r = stats.ranges[j]
    return (x_mm=Float64.(df.x_mm[r]).-stats.cx[j],
            y_mm=Float64.(df.y_mm[r]).-stats.cy[j],
            z_mm=Float64.(df.z_mm[r]).-stats.cz[j],
            ions=Float64.(df.energy_MeV[r]).*(1e6/w_eV))
end

"""Select the trace after shifting to its generation vertex; never recenter it."""
function select_hits(df,stats; xy_halfwidth_mm=Inf)
    xy_halfwidth_mm > 0 || throw(ArgumentError("Cut halfwidth must be positive"))
    selected = falses(size(df,1)); extent = zeros(3)
    for j in eachindex(stats.ids), i in stats.ranges[j]
        x = Float64(df.x_mm[i])-stats.cx[j]
        y = Float64(df.y_mm[i])-stats.cy[j]
        z = Float64(df.z_mm[i])-stats.cz[j]
        selected[i] = abs(x)<=xy_halfwidth_mm && abs(y)<=xy_halfwidth_mm
        if selected[i]
            extent .= max.(extent,(abs(x),abs(y),abs(z)))
        end
    end
    any(selected) || error("Spatial cut removed the entire sample")
    return (; selected,extent)
end

function add_event!(xy, z, event, grid)
    for i in eachindex(event.ions)
        ix = bin_index(event.x_mm[i], grid.x_edges_mm)
        iy = bin_index(event.y_mm[i], grid.y_edges_mm)
        iz = bin_index(event.z_mm[i], grid.z_edges_mm)
        xy[ix,iy] += event.ions[i]
        z[iz] += event.ions[i]
    end
end

function project_event(event, grid)
    xy = zeros(length(grid.x_edges_mm)-1, length(grid.y_edges_mm)-1)
    z = zeros(length(grid.z_edges_mm)-1)
    add_event!(xy,z,event,grid)
    return (; xy, z)
end

function summary(values)
    avg = sum(values)/length(values)
    # Population dispersion of the saved sample, not uncertainty on its mean.
    return [avg, sqrt(sum(x -> (x-avg)^2, values)/length(values)), minimum(values), maximum(values)]
end

"""Two marginal projections averaged over ALL events, including empty bins."""
function analyze_sample(df, stats, grid; w_eV=21.9, example_event_id=0, tolerance=1e-11,
                        xy_halfwidth_mm=Inf, diagnostic_grid=grid)
    example_index = findfirst(==(example_event_id), stats.ids)
    example_index === nothing && throw(KeyError(example_event_id))
    xy = zeros(length(grid.x_edges_mm)-1, length(grid.y_edges_mm)-1)
    z = zeros(length(grid.z_edges_mm)-1)
    errors = Float64[]
    selection = select_hits(df,stats; xy_halfwidth_mm)
    selected_energy = zeros(length(stats.ids)); selected_hits = zeros(Int,length(stats.ids))
    hx = zeros(length(diagnostic_grid.x_edges_mm)-1); hy = zeros(length(diagnostic_grid.y_edges_mm)-1)
    cx = zero(hx); cy = zero(hy)
    # Check every event by independent grouping in each projection, without
    # allocating a full dense matrix for each small track.
    for j in eachindex(stats.ids)
        event = centered_event(df, stats, j, w_eV)
        keep = selection.selected[stats.ranges[j]]
        selected_hits[j] = count(keep)
        selected_energy[j] = sum(i -> Float64(df.energy_MeV[i]),stats.ranges[j][keep]; init=0.0)
        xy_sparse = Dict{Tuple{Int,Int},Float64}()
        z_sparse = Dict{Int,Float64}()
        for i in eachindex(event.ions)
            ix_raw = bin_index(event.x_mm[i],diagnostic_grid.x_edges_mm)
            iy_raw = bin_index(event.y_mm[i],diagnostic_grid.y_edges_mm)
            hx[ix_raw] += event.ions[i]; hy[iy_raw] += event.ions[i]
            cx[ix_raw] += 1; cy[iy_raw] += 1
            keep[i] || continue
            ix = bin_index(event.x_mm[i],grid.x_edges_mm)
            iy = bin_index(event.y_mm[i],grid.y_edges_mm)
            iz = bin_index(event.z_mm[i],grid.z_edges_mm)
            xy_sparse[(ix,iy)] = get(xy_sparse,(ix,iy),0.0)+event.ions[i]
            z_sparse[iz] = get(z_sparse,iz,0.0)+event.ions[i]
        end
        for ((ix,iy), value) in xy_sparse
            xy[ix,iy] += value
        end
        for (iz, value) in z_sparse
            z[iz] += value
        end
        total = selected_energy[j]*(1e6/w_eV)
        actual = [sum(values(xy_sparse);init=0.0),sum(values(z_sparse);init=0.0),sum(event.ions[keep])]
        push!(errors, total>0 ? maximum(abs.(actual./total .- 1)) : maximum(abs,actual))
    end
    n = length(stats.ids)
    xy ./= n; z ./= n
    selected_ions = selected_energy.*(1e6/w_eV)
    mean_ions = sum(selected_ions)/n
    xy_error = abs(sum(xy)/mean_ions-1)
    z_error = abs(sum(z)/mean_ions-1)
    maxerr = max(maximum(errors),xy_error,z_error)
    maxerr <= tolerance || error("Ion conservation failed: relative error $maxerr")
    example = centered_event(df,stats,example_index,w_eV)
    keep = selection.selected[stats.ranges[example_index]]
    example = NamedTuple{keys(example)}(Tuple(values(example)[i][keep] for i in 1:4))
    example_hist = project_event(example,grid)
    return (; xy, z, example, example_hist, example_index, errors, xy_error, z_error,
            per_hit_ions=Float64.(df.energy_MeV).*(1e6/w_eV), selected=selection.selected,
            selected_energy,selected_ions,selected_hits,
            diagnostic_grid,diagnostic_x_ions=hx./n,diagnostic_y_ions=hy./n,
            diagnostic_x_hits=cx./n,diagnostic_y_hits=cy./n)
end

function put(group,name,value,unit="1"; axes="")
    if value isa AbstractArray && !isempty(value)
        chunk = ntuple(i -> min(size(value,i), ndims(value)==1 ? 4096 : 64), ndims(value))
        group[name, chunk=chunk, shuffle=true, deflate=3] = value
    else
        group[name] = value
    end
    attributes(group[name])["units"] = unit
    isempty(axes) || (attributes(group[name])["axis_order_on_disk"] = axes)
end

"""HDF5 matrices read as [x,y] in Julia and [y,x] in h5py; both are labelled."""
function write_analysis(path, input, sample, stats, grid, result, params, config_path, detector_path)
    mkpath(dirname(path))
    mkpath(joinpath(ROOT,"work"))
    mktempdir(joinpath(ROOT,"work"); prefix="ion-analysis-") do scratch
        temporary = joinpath(scratch,"analysis.h5")
        h5open(temporary,"w") do f
            attributes(f)["schema"] = SCHEMA
            meta = create_group(f,"metadata")
            for (name,value) in ("sample_id"=>sample["id"], "sample_label"=>sample["label"],
                "input_file"=>relpath(input,ROOT), "input_sha256"=>sha(input),
                "config_sha256"=>sha(config_path), "detector_config_sha256"=>sha(detector_path),
                "analysis_code_sha256"=>sha(@__FILE__), "manifest_sha256"=>sha(joinpath(ROOT,"Manifest.toml")),
                "julia_version"=>string(VERSION), "centering"=>"Common primary generation vertex; translations only; no rotations or barycenter calculation",
                "normalization"=>"Ions per bin per saved event; includes zero contributions; no per-track normalization",
                "assumptions"=>"Mean ionization production only. No Fano fluctuations, recombination, diffusion, focusing or collection losses. Raw totals and spatially selected projections are recorded separately.",
                "time_interpretation"=>"Longitudinal coordinate relative to the generation vertex divided by uniform speed; not an absolute arrival time or a transport simulation",
                "boundary_rule"=>"Left-closed/right-open; last bin includes upper edge; out-of-range hits are errors")
                meta[name] = value
            end
            h5open(input,"r") do src
                meta["original_file"] = read(src["metadata/source_file"])
                meta["original_sha256"] = read(src["metadata/source_sha256"])
                if haskey(src,"metadata/configuration")
                    cfg = create_group(meta,"source_configuration")
                    cfg["param_key"] = read(src["metadata/configuration/param_key"])
                    cfg["param_value"] = read(src["metadata/configuration/param_value"])
                end
            end
            p = create_group(f,"parameters")
            for (key,value) in pairs(params)
                p[string(key)] = value
            end
            ev = create_group(f,"events")
            for (name,value,unit) in (("event_id",stats.ids,"1"), ("first_row",first.(stats.ranges),"1"),
                ("n_hits",stats.nhits,"1"), ("energy_MeV",stats.energy,"MeV"),
                ("expected_ions",stats.expected_ions,"ions"), ("mean_ions_per_hit",stats.ions_per_hit,"ions/hit"),
                ("vertex_x_mm",stats.cx,"mm"), ("vertex_y_mm",stats.cy,"mm"), ("vertex_z_mm",stats.cz,"mm"))
                put(ev,name,value,unit)
            end
            attributes(ev)["index_base"] = 1
            put(ev,"selected_energy_MeV",result.selected_energy,"MeV")
            put(ev,"selected_expected_ions",result.selected_ions,"ions")
            put(ev,"selected_n_hits",result.selected_hits,"hits")
            put(ev,"discarded_energy_MeV",stats.energy.-result.selected_energy,"MeV")
            put(ev,"discarded_expected_ions",stats.expected_ions.-result.selected_ions,"ions")
            hits = create_group(f,"hits")
            put(hits,"expected_ions",result.per_hit_ions,"ions")
            put(hits,"selected",UInt8.(result.selected),"boolean")
            attributes(hits)["row_alignment"] = "Exactly the input hitdf row order; identifiers and coordinates stay in input"
            bins = create_group(f,"bins")
            for (name,edges) in pairs(grid)
                put(bins,string(name),edges,"mm")
            end
            put(bins,"relative_time_edges_ms",grid.z_edges_mm./params.drift_speed_mm_ms,"ms")
            means = create_group(f,"mean")
            put(means,"xy",result.xy,"ions/bin/event"; axes="y,x")
            put(means,"z",result.z,"ions/bin/event"; axes="z")
            ex = create_group(f,"example")
            put(ex,"event_id",stats.ids[result.example_index])
            for (name,value) in pairs(result.example)
                put(ex,string(name),value,name==:ions ? "ions" : "mm")
            end
            put(ex,"xy",result.example_hist.xy,"ions/bin"; axes="y,x")
            put(ex,"z",result.example_hist.z,"ions/bin"; axes="z")
            totals = create_group(f,"summary")
            totals["statistic_order"] = ["mean","population_std","min","max"]
            put(totals,"n_events",length(stats.ids)); put(totals,"n_hits",sum(stats.nhits))
            put(totals,"energy_MeV",summary(stats.energy),"MeV")
            put(totals,"expected_ions",summary(stats.expected_ions),"ions/event")
            put(totals,"mean_ions_per_hit",summary(stats.ions_per_hit),"ions/hit")
            put(totals,"n_hits_per_event",summary(stats.nhits),"hits/event")
            put(totals,"total_expected_ions",sum(stats.expected_ions),"ions")
            put(totals,"total_energy_MeV",sum(stats.energy),"MeV")
            put(totals,"pooled_mean_ions_per_hit",sum(stats.expected_ions)/sum(stats.nhits),"ions/hit")
            sel = create_group(f,"selection")
            sel["description"] = "Keep |x - vertex_x| <= XY halfwidth and |y - vertex_y| <= XY halfwidth; no Z cut or recentering"
            put(sel,"xy_halfwidth_mm",get(params,:xy_halfwidth_mm,Inf),"mm")
            put(sel,"n_selected_hits",sum(result.selected_hits),"hits")
            put(sel,"n_discarded_hits",sum(stats.nhits)-sum(result.selected_hits),"hits")
            put(sel,"empty_events",count(==(0),result.selected_hits),"events")
            put(sel,"mean_expected_ions",sum(result.selected_ions)/length(stats.ids),"ions/event")
            put(sel,"total_selected_ions",sum(result.selected_ions),"ions")
            put(sel,"total_discarded_ions",sum(stats.expected_ions)-sum(result.selected_ions),"ions")
            put(sel,"discarded_ion_fraction",1-sum(result.selected_ions)/sum(stats.expected_ions))
            put(sel,"discarded_hit_fraction",1-sum(result.selected_hits)/sum(stats.nhits))
            put(sel,"maximum_event_discarded_fraction",maximum(1 .- result.selected_ions./stats.expected_ions))
            put(sel,"energy_MeV",summary(result.selected_energy),"MeV")
            put(sel,"expected_ions",summary(result.selected_ions),"ions/event")
            diag = create_group(f,"diagnostics")
            for (name,value,unit) in (("x_edges_mm",result.diagnostic_grid.x_edges_mm,"mm"),
                ("y_edges_mm",result.diagnostic_grid.y_edges_mm,"mm"),("uncut_x_ions",result.diagnostic_x_ions,"ions/bin/event"),
                ("uncut_y_ions",result.diagnostic_y_ions,"ions/bin/event"),("uncut_x_hits",result.diagnostic_x_hits,"hits/bin/event"),
                ("uncut_y_hits",result.diagnostic_y_hits,"hits/bin/event"))
                put(diag,name,value,unit)
            end
            valid = create_group(f,"validation")
            put(valid,"event_max_relative_error",maximum(result.errors))
            put(valid,"xy_mean_relative_error",result.xy_error)
            put(valid,"z_mean_relative_error",result.z_error)
            put(valid,"out_of_range_hits",0)
            put(valid,"max_centered_extent_mm",stats.extent,"mm")
        end
        # Check serialized shapes and normalization independently of the writer.
        h5open(temporary,"r") do f
            a = read(f["mean/xy"]); z = read(f["mean/z"])
            size(a) == (length(grid.x_edges_mm)-1,length(grid.y_edges_mm)-1) || error("XY shape mismatch")
            for total in (sum(a),sum(z))
                isapprox(total, sum(result.selected_ions)/length(stats.ids); rtol=params.tolerance) || error("Saved charge mismatch")
            end
        end
        Base.Filesystem.rename(temporary,path)
    end
    return path
end

function run_analysis(config_path=joinpath(ROOT,"config/ion_tracks.toml"))
    config = TOML.parsefile(config_path)
    detector_path = joinpath(ROOT,"config/detector.toml")
    detector = TOML.parsefile(detector_path)
    cfg = config["binning"]
    cfg["centering"] == "generation_vertex" || error("Unsupported centering")
    pitch = 1000*detector["geometry"]["pitch_m"]
    speed = detector["gas"]["bulk_speed_m_s"] # 1 m/s = 1 mm/ms, interface conversion.
    params = (w_eV=Float64(config["ionization"]["w_eV"]), pitch_mm=pitch,
              frame_ms=Float64(cfg["frame_ms"]), drift_speed_mm_ms=speed,
              z_bin_mm=speed*cfg["frame_ms"], margin_fraction=Float64(cfg["margin_fraction"]),
              xy_halfwidth_mm=Float64(config["selection"]["xy_halfwidth_mm"]),
              extra_bins=Int(cfg["extra_bins"]), tolerance=Float64(config["numerics"]["charge_relative_tolerance"]),
              w_status=config["ionization"]["w_status"], reference_key=config["ionization"]["reference_key"])
    isfinite(params.tolerance) && params.tolerance>0 || error("Invalid tolerance")
    samples = config["samples"]
    allunique(s["id"] for s in samples) || error("Duplicate sample names")
    all(s -> occursin(r"^[A-Za-z0-9_-]+$",s["id"]),samples) || error("Invalid sample name")
    stats = Any[]; hashes = String[]; retained_extents = Any[]
    for sample in samples
        input = joinpath(ROOT,sample["input"])
        df = read_hitdf(input)
        s = event_statistics(df,read_event_metadata(input); w_eV=params.w_eV)
        sample["example_event_id"] in s.ids || throw(KeyError(sample["example_event_id"]))
        push!(stats,s); push!(hashes,sha(input))
        push!(retained_extents,select_hits(df,s; xy_halfwidth_mm=params.xy_halfwidth_mm))
        println(sample["id"], ": ",length(s.ids)," events; centered extents (mm) ",s.extent)
        flush(stdout)
    end
    diagnostic_grid = common_grid(stats; pitch_mm=params.pitch_mm,z_bin_mm=params.z_bin_mm,
                       margin=params.margin_fraction,extra_bins=params.extra_bins)
    grid = common_grid(retained_extents; pitch_mm=params.pitch_mm,z_bin_mm=params.z_bin_mm,
                       margin=params.margin_fraction,extra_bins=params.extra_bins)
    println("Common grid: ",length(grid.x_edges_mm)-1," x ",length(grid.y_edges_mm)-1,
            "; Z bins: ",length(grid.z_edges_mm)-1)
    for (sample,s,hash) in zip(samples,stats,hashes)
        input = joinpath(ROOT,sample["input"])
        sha(input) == hash || error("Input changed between analysis passes")
        result = analyze_sample(read_hitdf(input),s,grid; w_eV=params.w_eV,
                                example_event_id=sample["example_event_id"],tolerance=params.tolerance,
                                xy_halfwidth_mm=params.xy_halfwidth_mm,diagnostic_grid)
        output = joinpath(ROOT,"results/ion_tracks",sample["id"],"analysis.h5")
        write_analysis(output,input,sample,s,grid,result,params,config_path,detector_path)
        println("Saved ",relpath(output,ROOT),"; mean retained ions/event = ",sum(result.selected_ions)/length(s.ids),
                "; max conservation error = ",max(maximum(result.errors),result.xy_error,result.z_error))
        flush(stdout)
    end
end

end
