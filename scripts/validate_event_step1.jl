# Step-one checkpoint: real-file reader, arithmetic and deterministic-core checks.
# Run after make test-events test-event-core test-ion-analysis.
using HDF5, DataFrames, TOML, SHA, Dates
include(joinpath(@__DIR__,"../src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis
const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT,"results/event_analysis/step1_validation.h5")
mkpath(dirname(OUT))
core = TOML.parsefile(joinpath(ROOT,"work/event_analysis/step1_validation.toml"))
core["status"] == "passed" || error("Run the core tests successfully first")

"""Read test-suite pass counts from the captured Julia test summaries."""
function test_passes(file)
    lines = readlines(file)
    any(line -> occursin("ERROR:",line) || occursin("Test Failed",line),lines) && error("Failed test log: $file")
    counts = [parse(Int,m.captures[1]) for line in lines
              for m in [match(r"^.+\|\s+(\d+)\s+\d+\s+[0-9.]+s$",line)] if m !== nothing]
    isempty(counts) && error("No test summary: $file")
    return sum(counts)
end

h5open(OUT,"w") do f
    attributes(f)["schema"] = "nausicaa.event.step1.validation.v1"
    attributes(f)["created_utc"] = string(now(UTC))
    attributes(f)["julia_version"] = string(VERSION)
    attributes(f)["scope"] = "Reader and deterministic reconstruction checkpoint. Raw-event smoke checks use 2 mm bins, kNN k=10, 4 mm cutoff, zero threshold; these are software checks."
    provenance = create_group(f,"provenance")
    for (path,hash) in core["files"]
        bytes2hex(sha256(read(joinpath(ROOT,path)))) == hash || error("Changed after test: $path")
    end
    provenance["paths"] = collect(keys(core["files"]))
    provenance["sha256"] = String[v for v in values(core["files"])]
    tests = create_group(f,"tests")
    tests["suite"] = ["Event core","Event IO","Legacy ion analysis"]
    tests["passed"] = [core["counts"]["passed"],
        test_passes(joinpath(ROOT,"work/event_analysis/step1_io.log")),
        test_passes(joinpath(ROOT,"work/event_analysis/step1_legacy.log"))]
    for sample in ("0nubb","electron_0nubb_energy")
        input = joinpath(ROOT,"events",sample*"_hitdf.h5")
        h5open(input,"r") do src
            g = create_group(f,"samples/"*sample)
            md = EventAnalysis.EventIO.read_event_metadata(src)
            allhits = EventAnalysis.EventIO.read_hitdf(src)
            g["events"] = nrow(md); g["hits"] = nrow(allhits)
            g["bytes"] = filesize(input)
            g["source_sha256"] = read(src["metadata/source_sha256"])
            g["reduced_sha256"] = bytes2hex(open(sha256,input))
            g["primary_counts"] = sort(unique(md.n_primaries))
            relative_errors = Float64[]; absolute_errors = Float64[]
            for (start,count,total) in zip(md.first_row,md.n_hits,md.energy_MeV)
                reference = sum(Float64,view(allhits.energy_MeV,start:start+count-1))
                err = abs(Float64(total)-reference)
                push!(absolute_errors,err)
                push!(relative_errors,reference == 0 ? 0.0 : err/reference)
            end
            # Float64 is used solely as the independent validation reference.
            g["max_energy_absolute_error_MeV"] = maximum(absolute_errors)
            g["max_energy_relative_error"] = maximum(relative_errors)
            g["mean_energy_relative_error"] = sum(relative_errors)/length(relative_errors)
            ids = md.event_id[unique(round.(Int,range(1,nrow(md);length=5)))]
            smoke = create_group(g,"raw_event_smoke")
            smoke["event_id"] = ids
            ncomponents = Int64[]; ntruth = Int64[]; energies = Float32[]
            # Warm the indexed reader before measuring five independent reads.
            read_event(src,first(ids))
            elapsed = @elapsed for id in ids
                read_event(src,id)
            end
            g["mean_indexed_read_ms"] = 1000*elapsed/length(ids)
            for id in ids
                hits,info = read_event(src,id)
                bins = voxelize_event(to_petit_dataframe(hits),2f0)
                tracks = make_tracks(bins;max_distance_mm=4f0,k=10)
                truth = compute_mc_path(hits,info,0.5f0)
                for t in tracks
                    path = walk_track_from_extremes(t)
                    isfinite(path.total_length) && isfinite(path.confidence) || error("Nonfinite backbone")
                    raw = get_raw_path(t,path.path_indices)
                    isempty(raw.s) || isapprox(last(raw.s),path.total_length) || error("Path length mismatch")
                end
                push!(ncomponents,length(tracks));push!(ntruth,nrow(truth));push!(energies,sum(bins.energy))
            end
            smoke["components"] = ncomponents; smoke["truth_rows"] = ntruth
            smoke["binned_energy_MeV"] = energies
            println(sample,": ",nrow(md)," events, max relative Float32 energy error ",maximum(relative_errors))
        end
    end
end
println(OUT)
