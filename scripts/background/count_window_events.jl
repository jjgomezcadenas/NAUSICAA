#!/usr/bin/env julia
#=
Count the events of NEXUS files whose energy deposited in the gas lies in a window.

  julia --project=. scripts/background/count_window_events.jl LO_keV HI_keV OUTPUT.toml FILE.h5 ...

For each file the energy of an event is the sum of MC/hits (MeV to keV), read
in blocks. Writes OUTPUT.toml with the number of events read and in the
window, in total and per file, and the saved window recorded in the
configuration of the first file. Used for the eITACA 208Tl production, saved
on a wider window than the 2420-2490 keV of the copper-shield runs.
=#
using HDF5, TOML

function count_window(lo::Float64, hi::Float64, files::Vector{String})
    started = time(); read_events = 0; in_window = 0; per_file = Dict{String,Any}(); saved = Dict{String,Any}()
    for (i, path) in enumerate(files)
        h5open(path, "r") do f
            if i == 1
                for r in read(f["MC/configuration"])
                    k = String(strip(String(r.param_key)))
                    occursin("DefaultEventAction", k) && (saved[k] = String(strip(String(r.param_value))))
                end
            end
            ds = f["MC/hits"]; n = length(ds); energy = Dict{Int64,Float64}(); block = 2_000_000
            for first in 1:block:n, h in ds[first:min(first + block - 1, n)]
                e = Int64(h.event_id); energy[e] = get(energy, e, 0.0) + Float64(h.energy) * 1000
            end
            ne = length(energy); nw = count(v -> lo <= v <= hi, values(energy))
            read_events += ne; in_window += nw; per_file[basename(path)] = [ne, nw]
        end
        i % 50 == 0 && (println("$i files: $read_events events read, $in_window in the window, $(round(time() - started; digits = 0)) s"); flush(stdout))
    end
    Dict("window_keV" => [lo, hi], "files" => length(files), "read_events" => read_events, "in_window" => in_window,
         "saved_window_of_files" => saved, "per_file" => per_file, "time_s" => time() - started)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 4 || error("Usage: count_window_events.jl LO_keV HI_keV OUTPUT.toml FILE.h5 ...")
    result = count_window(parse(Float64, ARGS[1]), parse(Float64, ARGS[2]), String.(ARGS[4:end]))
    open(ARGS[3], "w") do io; TOML.print(io, result); end
    println("COUNT $(ARGS[3]): $(result["in_window"]) of $(result["read_events"]) events in $(result["window_keV"]) keV, $(result["files"]) files")
end
