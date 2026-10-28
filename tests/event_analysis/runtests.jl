using Test, TOML, SHA, Dates

const ROOT = normpath(joinpath(@__DIR__, "../.."))
const STEP4 = "--step4" in ARGS
const STEP3 = STEP4 || "--step3" in ARGS
const STEP2 = STEP3 || "--step2" in ARGS
const STEP = STEP4 ? 4 : STEP3 ? 3 : STEP2 ? 2 : 1
const REPORT = joinpath(ROOT, "work/event_analysis/step$(STEP)_validation.toml")
const INPUT_FILES = sort(vcat(["Project.toml", "Manifest.toml", "src/ion_signal/event_io.jl"],
    [relpath(joinpath(dir,file),ROOT) for base in ("src/event_analysis", "tests/event_analysis")
     for (dir,_,files) in walkdir(joinpath(ROOT,base)) for file in files]))
report = Dict{String,Any}(
    "schema" => "nausicaa.event_core.validation.v1", "status" => "running",
    "julia_version" => string(VERSION), "started_utc" => string(now(UTC)),
    "files" => Dict(path => bytes2hex(sha256(read(joinpath(ROOT, path)))) for path in INPUT_FILES))
mkpath(dirname(REPORT))
open(io -> TOML.print(io, report), REPORT, "w")

include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis
include("fixtures.jl")

suite = @testset "Event analysis through step $STEP" begin
    include("test_types.jl")
    include("test_coordinates.jl")
    include("test_docstrings.jl")
    include("test_reconstruction.jl")
    STEP2 && include("test_diffusion.jl")
    STEP3 && include("test_collection.jl")
    STEP4 && include("test_ensemble.jl")
    include("test_display.jl")
    include("test_mc_inspection.jl")
    include("test_mc_selection.jl")
    include("test_pixel_charge.jl")
    include("test_diffusion_cache.jl")
end

function counts(testset)
    c = Test.get_test_counts(testset)
    return Dict("passed" => c.passes + c.cumulative_passes,
                "failed" => c.fails + c.cumulative_fails,
                "errors" => c.errors + c.cumulative_errors,
                "broken" => c.broken + c.cumulative_broken)
end

report["status"] = "passed"
report["counts"] = counts(suite)
report["groups"] = [merge(Dict("name" => t.description), counts(t))
                    for t in suite.results if t isa Test.DefaultTestSet]
report["completed_utc"] = string(now(UTC))
open(io -> TOML.print(io, report), REPORT, "w")
println("Validation record: ", relpath(REPORT, ROOT))
