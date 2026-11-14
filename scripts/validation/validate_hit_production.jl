#!/usr/bin/env julia
include("../production/run_hit_production.jl")
using Test
# Fail immediately if production enters any diagnostic constructor/fill routine.
control_histogram(x,edges)=error("Control histogram called in production")
histogram_edges(values,bins)=error("Histogram edges called in production")
channel_histograms(values,valid,input,bins)=error("Noise histograms called in production")
hit_debug(c;kwargs...)=error("Debug entry called in production")
function event_digest(event)
    io=IOBuffer()
    write(io,event.info.event_number,event.info.ion_drift_cm,event.info.electron_drift_cm,
      event.info.ion_sigma_t_mm,event.info.electron_sigma_t_mm)
    for bounds in (event.info.frame.lower_mm,event.info.frame.upper_mm),v in bounds;write(io,v);end
    for key in fieldnames(EventEnergy);write(io,getfield(event.energy,key));end
    for arm in (:itrack,:etrack),key in fieldnames(SavedHits);write(io,getfield(getfield(event,arm),key));end
    bytes2hex(sha256(take!(io)))
end
card=abspath(ARGS[1]);TOML.parsefile(card)["controls"]["mode"]=="none"||error("Validation requires controls=none")
expected=Dict{Int64,String}()
path=run_hit_production(card;on_event=e->expected[e.info.event_number]=event_digest(e))
@testset "Production-only output and exact read-back" begin
    h5open(path,"r") do f
        @test Set(keys(f))==Set(["schema","events","itrack","etrack","configuration"])
        @test length(read(f["events/event_number"]))==length(expected)
        for arm in ("itrack","etrack")
            @test Set(keys(f[arm]))==Set(["offsets","x_mm","y_mm","z_mm","energy_keV"])
        end
    end
    reader=open_hit_input(path)
    seconds=@elapsed for (id,digest) in expected
        @test event_digest(read_hit_event(reader,id))==digest
    end
    close_hit_input(reader)
    @test !isfile(path*".tmp")
    @test !isdir(joinpath(dirname(path),"debug"))
    println("READ_BACK seconds=$seconds events=$(length(expected))")
end
