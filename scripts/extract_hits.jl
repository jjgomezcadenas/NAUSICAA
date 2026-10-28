#!/usr/bin/env julia

include(joinpath(@__DIR__, "../src/ion_signal/event_io.jl"))
using .EventIO

function main(args)
    if isempty(args) || "--help" in args || "-h" in args
        println("""
        Usage: julia --project=. scripts/extract_hits.jl [--overwrite] INPUT.h5 [INPUT.h5 ...]

        Convert NEXUS /MC/hits to NAME_hitdf.h5 beside each input.
        Keep event_id, particle_id, hit_id, x_mm, y_mm, z_mm, energy_MeV, time_ns.
        Preserve Float32, mm, MeV, ns, all particles and source row order.
        Store primary IDs, the common vertex, an event index and provenance; verify before publishing each output.
        Existing output files require --overwrite. Original inputs are never replaced.
        """)
        return
    end
    overwrite = "--overwrite" in args
    inputs = filter(!=("--overwrite"), args)
    isempty(inputs) && throw(ArgumentError("At least one input is required"))
    any(startswith("-"), inputs) && throw(ArgumentError("Unknown option; see --help"))
    for input in inputs
        output = default_output(input)
        println("Converting $input -> $output (mm, MeV, ns; Float32)")
        flush(stdout)
        t = time()
        result = convert_hits(input; overwrite)
        println("Verified $(result.n_events) events, $(result.n_hits) hits; ",
                round(result.output_bytes/1e6; digits=2), " MB; ",
                round(time()-t; digits=1), " s")
        flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main(ARGS)
    catch e
        showerror(stderr, e)
        println(stderr)
        exit(1)
    end
end
