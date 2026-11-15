#!/usr/bin/env julia
include("MagboltzDriver.jl")
using .MagboltzDriver

function main(args)
    length(args) in (2,5) || error("Usage: julia --project=. magboltz/run.jl prepare CARD | run CARD BINARY RUN_DIR RESULT_H5")
    mode=args[1];point=read_point(args[2])
    if mode=="prepare" && length(args)==2
        print(input_card(point))
    elseif mode=="run" && length(args)==5
        source=joinpath(@__DIR__,"source","magboltz-11.21.f")
        row=run_point(point,args[3],args[4])
        save_table(args[5],[row];source_path=source,binary_path=args[3])
        println("Wrote ",args[5]," in ",round(row.elapsed_s;digits=2)," s")
    else
        error("Unknown command or argument count")
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main(ARGS)
