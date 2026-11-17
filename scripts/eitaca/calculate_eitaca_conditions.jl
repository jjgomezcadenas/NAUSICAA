#!/usr/bin/env julia
#=
eITACA detector: resolve the default set (variant IE_L1500_D400), compute the
derived quantities and save both for the detector page.

  julia --project=. scripts/eitaca/calculate_eitaca_conditions.jl

Output: results/eitaca/detector/conditions.h5 (schema eitaca.conditions.v1):
the card values as dotted keys, the derived quantities with unit and note,
and the SHA-256 of the cards. scripts/eitaca/present_eitaca_detector.py writes
notes/html/eitaca/index.html from it. An existing file is preserved.
=#
using HDF5, TOML
const ROOT = normpath(joinpath(@__DIR__, "../.."))
include(joinpath(ROOT, "src/event_analysis/EventAnalysis.jl"))
include(joinpath(ROOT, "src/itaca/ItacaConditions.jl"))
include(joinpath(ROOT, "src/eitaca/EItaca.jl")); using .EItaca

function main()
    r = resolve_detector(joinpath(ROOT, "config/eitaca/eitaca.toml"); variant = joinpath(ROOT, "config/eitaca/variants/IE_L1500_D400.toml"), root = ROOT)
    out = joinpath(ROOT, "results/eitaca/detector"); mkpath(out); target = joinpath(out, "conditions.h5")
    isfile(target) && error("Preserve existing conditions: $target")
    h5open(target * ".tmp", "w") do f
        f["schema"] = "eitaca.conditions.v1"; f["variant"] = r.variant; f["date"] = r.set.date
        for (k, v) in ItacaConditions.flatten(r.detector)
            f["inputs/$k"] = v isa Real && !(v isa Bool) ? Float64(v) : v isa AbstractDict ? sprint(io -> TOML.print(io, v)) : string(v)
        end
        for (i, (g, n, v, u, note)) in enumerate(derived_conditions(r.detector; root = ROOT))
            d = create_group(f, "derived/$g.$n"); d["value"] = v; d["unit"] = u; d["note"] = note; d["order"] = i
        end
        for (k, v) in r.set.sha256; f["sha256/$k"] = v; end
    end
    Base.rename(target * ".tmp", target); println("OUTPUT $target")
end
main()
