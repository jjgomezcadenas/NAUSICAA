module ProjectPaths
using TOML
export resolve_recorded_path

"""Resolve a recorded repository path through the migration maps without changing cards."""
function resolve_recorded_path(root::AbstractString, path::AbstractString)
    root = abspath(root)
    value = isabspath(path) ? relpath(path, root) : String(path)
    startswith(value, "../") && return abspath(path)
    moves = Pair{String,String}[]
    mapping = joinpath(root, "docs/path_migrations.toml")
    if isfile(mapping)
        for entry in TOML.parsefile(mapping)["moves"]
            push!(moves, entry["old"] => entry["new"])
        end
    end
    sort!(moves; by=x -> length(first(x)), rev=true)
    for _ in 1:8
        found = findfirst(x -> value == first(x) || startswith(value, first(x)*"/"), moves)
        found === nothing && return joinpath(root, value)
        old, new = moves[found]
        replacement = new * value[length(old)+1:end]
        replacement == value && return joinpath(root, value)
        value = replacement
    end
    error("Cyclic relocation mapping for $path")
end
end
