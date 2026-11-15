#!/usr/bin/env julia
# Classify the detached graph components of a run's batches against the MC
# truth and write the counts as TOML.
#   julia --project=. scripts/production/classify_run_batches.jl OUTPUT.toml "DESCRIPTION" [--ion-links=7,8] RUN_CARD [RUN_CARD ...]
# Needs the batch hit files under results/readout/<run_label>/batches/ and the
# MC inputs of the run cards (NAUSICAA_DATA). The itrack is classified at each
# listed connection radius (default: the card's); the etrack at the card's.
include("classify_ion_components.jl")

function classify_run_batches(output,description,cards;ion_links=nothing)
    isfile(output)&&error("Preserve existing classification: $output")
    rows=Dict{String,Any}[];scratch=mktempdir()
    for card in cards
        run=TOML.parsefile(card);label=run["run_label"];sample=run["sample"];tc=run["track_analysis"];subdir=run["hit_production"]["output_subdir"]
        base=joinpath(ROOT,"results/readout",label,"batches");tags=sort(filter(t->occursin(r"^b\d{5}$",t),readdir(base)))
        links=ion_links===nothing ? [Float64(tc["ion_connection_mm"])] : ion_links
        for (arm,radii) in (("itrack",links),("etrack",[Float64(tc["electron_connection_mm"])])), radius in radii
            graph=joinpath(scratch,"graph.toml");settings=copy(tc);arm=="itrack"&&(settings["ion_connection_mm"]=radius)
            open(io->TOML.print(io,Dict("track_analysis"=>settings)),graph,"w")
            events=0;multi=0;gammas=0
            for tag in tags
                hits=joinpath(base,tag,subdir,"hits.h5");events+=h5open(f->length(read(f["events/event_number"])),hits,"r")
                inputs=Dict(k=>input_path(ROOT,replace(run[k],"{batch}"=>tag)) for k in ("hits","particles"));out=joinpath(scratch,"$(label)_$(arm)_$(radius)_$tag.h5")
                try
                    classify(hits,graph,inputs["hits"],inputs["particles"],out;arm)
                    h5open(out,"r") do f;multi+=length(read(f["event_number"]));gammas+=count(read(f["all_secondary_gamma"]));end
                catch e
                    e isa BoundsError||e isa ArgumentError||rethrow()   # a batch with no multi-component event
                end
            end
            push!(rows,Dict{String,Any}("sample"=>sample,"run_label"=>label,"track"=>arm,"link_mm"=>radius,"events"=>events,
                "one_component"=>events-multi,"floating_gammas"=>gammas,"broken"=>multi-gammas))
        end
    end
    result=Dict{String,Any}("schema"=>"itaca.component_classification_summary.v1","sample_description"=>description,
        "rule"=>"A detached component is a floating gamma when at least 80% of its energy lies within the connection radius of Monte Carlo gamma hits and nearer to them than to the main track; an event with a detached component that fails this is a broken track.",
        "rows"=>rows)
    mkpath(dirname(output));open(io->TOML.print(io,result),output,"w")
    println("CLASSIFICATION $output");output
end

function classify_from_arguments(arguments)
    links=nothing;args=String[]
    for a in arguments;startswith(a,"--ion-links=") ? (links=parse.(Float64,split(split(a,"=")[2],","))) : push!(args,a);end
    length(args)>=3||error("Usage: classify_run_batches.jl OUTPUT.toml DESCRIPTION [--ion-links=7,8] RUN_CARD [RUN_CARD ...]")
    classify_run_batches(abspath(args[1]),args[2],abspath.(args[3:end]);ion_links=links)
end

if abspath(PROGRAM_FILE)==@__FILE__
    classify_from_arguments(ARGS)
end

