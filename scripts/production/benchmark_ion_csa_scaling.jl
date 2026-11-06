# Measure fresh production outputs after warm-up; preserve exact cards and timings.
using TOML,HDF5,Statistics
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/ion_signal/ion_csa_pipeline.jl"))
const P=IonCSAProduction
base=joinpath(ROOT,"results/ion_csa_nonlinear/scaling",string(time_ns()))
mkpath(base);write(joinpath(base,".gitignore"),"execution/\n*.h5\n*.tmp\n")
template=TOML.parsefile(joinpath(ROOT,"config/specs/runs/ion_csa_nonlinear_production_100.toml"))
response=TOML.parsefile(joinpath(ROOT,template["response_card"]))
records=Dict{String,Any}[]
for (ordinal,n) in enumerate((10,100,300))
    c=copy(response)
    # A recorded fresh noise seed gives each benchmark an uncached response identity.
    c["noise_seed"]=Int(time_ns()%UInt64(1_000_000_000))+ordinal
    cp=P.write_toml(joinpath(base,"response_$n.toml"),c)
    run=copy(template);run["response_card"]=relpath(cp,ROOT);run["events_per_sample"]=n
    run["results_directory"]=relpath(joinpath(base,"events_$n"),ROOT)
    card=P.write_toml(joinpath(base,"run_$n.toml"),run)
    s=P.settings(ROOT,card)
    # Check the actual filtered response identity, not just the outer run directory.
    execution=joinpath(s.dir,"execution/pitch4")
    raw=P.write_toml(joinpath(execution,"filtered_run.toml"),Dict("schema"=>"nausicaa.ion_csa_run.v1",
        "response_card"=>relpath(P.write_toml(joinpath(execution,"resolved_response.toml"),s.c),ROOT),
        "ion_index"=>run["ion_index"],"samples"=>run["samples"],"pitches_mm"=>[4],
        "events_per_sample"=>n,"selection"=>run["selection"],
        "results_directory"=>relpath(joinpath(execution,"filtered"),ROOT),
        "report_directory"=>relpath(joinpath(execution,"filtered_report"),ROOT)))
    probe=P.R.settings(ROOT,raw)
    @assert all(x->!isfile(x.response),probe.sources) "Benchmark would reuse a filtered response"
    maxframes=64
    for source in probe.sources
        h5open(source.path,"r") do h
            for id in source.selected
                buckets=read(h["events/$id/bucket"])
                maxframes=max(maxframes,maximum(buckets)-minimum(buckets)+probe.c["quiet_before_frames"]+ceil(Int,probe.c["recovery_tail_s"]/source.frame)+1)
            end
        end
    end
    calibration_start=time_ns();P.R.calibration(probe,maxframes)
    calibration_s=(time_ns()-calibration_start)/1e9
    started=time_ns();manifest=P.run(ROOT,card);elapsed=(time_ns()-started)/1e9
    timing=TOML.parsefile(joinpath(dirname(manifest),"execution/pitch4/timing.toml"))
    row=Dict{String,Any}("events"=>n,"warmup"=>ordinal==1,"elapsed_s"=>elapsed,
        "work_without_presentation_s"=>elapsed-timing["presentation_s"],"manifest"=>manifest,
        "fresh_filtered_outputs"=>true,"noise_seed"=>c["noise_seed"],
        "calibration_setup_s"=>calibration_s,"maximum_frames"=>maxframes)
    merge!(row,timing);push!(records,row)
    P.write_toml(joinpath(base,"measurements.toml"),Dict("measurements"=>records))
    println("BENCHMARK ",n," events: ",row["work_without_presentation_s"]," s excluding presentation; warmup=",ordinal==1);flush(stdout)
end
a,b=records[2:3];slope=(b["work_without_presentation_s"]-a["work_without_presentation_s"])/(b["events"]-a["events"])
fixed=a["work_without_presentation_s"]-slope*a["events"]
summary=Dict("seconds_per_event"=>slope,"fixed_warm_cost_s"=>fixed,
    "estimated_8000_event_work_s"=>fixed+8000slope,"measurements"=>records)
P.write_toml(joinpath(base,"scaling.toml"),summary)
P.write_toml(joinpath(dirname(base),"latest.toml"),Dict("summary"=>joinpath(base,"scaling.toml")))
println("SCALING ",slope," s/event; fixed ",fixed," s; 8000 events ",fixed+8000slope," s")
