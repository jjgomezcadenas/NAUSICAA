#!/usr/bin/env julia
#=
Compare two detector-response control files event by event.

  julia --project=. scripts/validation/compare_detector_controls.jl reference.h5 candidate.h5

Every dataset under /events of the candidate is compared with the reference on
the common event numbers (the candidate may hold the first N events of the
reference). Exact equality is expected when the two runs share the seeds and
the physics; the script prints, per dataset, the number of differing events and
the largest absolute difference, and exits with status 1 if any differ.
=#
using HDF5
length(ARGS) == 2 || error("Usage: compare_detector_controls.jl reference.h5 candidate.h5")
ref, cand = ARGS
status = 0
h5open(ref, "r") do a
    h5open(cand, "r") do b
        ida = read(a["events/event_number"]); idb = read(b["events/event_number"])
        pos = Dict(id => i for (i, id) in enumerate(ida))
        all(haskey(pos, id) for id in idb) || error("Candidate holds events absent from the reference")
        rows = [pos[id] for id in idb]
        println("events compared: ", length(idb), " of ", length(ida), " in the reference")
        for name in sort(collect(keys(b["events"])))
            name == "response_s" && continue           # timing, not physics
            haskey(a["events"], name) || (println(rpad(name, 40), " absent from the reference"); continue)
            x = read(a["events"][name])[rows]; y = read(b["events"][name])
            differ = count(x .!= y)
            maxdiff = eltype(x) <: Number ? maximum(abs.(Float64.(x) .- Float64.(y)); init = 0.0) : NaN
            println(rpad(name, 40), " differing events ", lpad(differ, 4), "   max |difference| ", maxdiff)
            differ > 0 && (global status = 1)
        end
        for name in ("source_sha256", "particles_sha256")
            haskey(a, name) && haskey(b, name) && println(rpad(name, 40), read(a[name]) == read(b[name]) ? " equal" : " DIFFERENT")
        end
    end
end
println(status == 0 ? "IDENTICAL" : "DIFFERENCES FOUND")
exit(status)
