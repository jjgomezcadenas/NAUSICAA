#!/usr/bin/env julia
# Process the input batches of a run card in sequence: detector response, then
# hits, for each batch; then chain the hit files.
#   julia --project=. scripts/production/run_batches.jl run_card.toml FIRST LAST [--no-merge]
# Batch b reads the inputs with {batch} = b0000b and writes under
# results/readout/<run_label>/batches/b0000b/. A batch whose hit file exists is
# skipped, so an interrupted run can be resumed. The chained file is
# results/readout/<run_label>/<output_subdir>/hits.h5.
include("run_hit_production.jl")
include("merge_hit_batches.jl")

function run_batches(card::AbstractString,first::Int,last::Int;merge::Bool=true)
    0<=first<=last||error("Give FIRST <= LAST, both from 0")
    run=TOML.parsefile(card);label=run["run_label"];subdir=run["hit_production"]["output_subdir"]
    base=joinpath(ROOT,"results/readout",label)
    hit_path(b)=joinpath(base,"batches",batch_tag(b),subdir,"hits.h5")
    control_path(b)=joinpath(base,"batches",batch_tag(b),"detector_response/controls.h5")
    for b in first:last
        if isfile(hit_path(b));println("BATCH $(batch_tag(b)): hits exist, skipped");continue;end
        t=time()
        isfile(control_path(b))||main(card,0;batch=b)
        run_hit_production(card;batch=b)
        println("BATCH $(batch_tag(b)) done in $(round(time()-t;digits=1)) s");flush(stdout)
    end
    merge||return nothing
    files=[hit_path(b) for b in first:last];target=joinpath(base,subdir,"hits.h5")
    merge_hit_batches(files,target)
    counts=merge_selection_counts([control_path(b) for b in first:last],joinpath(base,subdir,"selection_counts.toml"))
    println("CHAINED $(length(files)) batches: $target ($(filesize(target)) bytes); selection $(counts)")
    target
end

if abspath(PROGRAM_FILE)==@__FILE__
    args=filter(!=("--no-merge"),ARGS)
    length(args)==3||error("Usage: run_batches.jl run_card.toml FIRST LAST [--no-merge]")
    run_batches(abspath(args[1]),parse(Int,args[2]),parse(Int,args[3]);merge=!("--no-merge" in ARGS))
end
